"""Command line interface: `skinpipe <command>`."""
from __future__ import annotations

import logging
from enum import Enum
from pathlib import Path
from typing import Optional

import typer
from rich.console import Console
from rich.markup import escape as rich_escape
from rich.table import Table

from . import report as report_mod
from . import triage
from .filters import safety as safety_f
from .filters import stability as stability_f
from .filters import topical as topical_f
from .filters.admet_exposure import (APPLIED_PERCENT_MAX, APPLIED_PERCENT_MIN,
                                     PRODUCT_USE, validate_applied_percent)
from .registry import Registry
from .types import SEVERITY_ORDER, Grade

app = typer.Typer(add_completion=False, help="Skin-targeted small-molecule discovery pipeline.")

def _note_faster_platform() -> None:
    """Say once that MM/GBSA could run on a GPU, without choosing one.

    The default stays CPU on purpose: MM/GBSA energies move in the last decimals
    between platforms, so a table assembled from two of them is not one table,
    and picking whatever the host happens to have would make the platform depend
    on the machine rather than on the run. But the measured cost is 2.8 minutes
    per compound on CPU against a GPU sitting at 0% -- silence there is not
    neutral either, it is an hour per twenty compounds nobody chose to spend.
    So: name the flag and let the caller decide.
    """
    try:
        import openmm
        names = {openmm.Platform.getPlatform(i).getName()
                 for i in range(openmm.Platform.getNumPlatforms())}
    except Exception:                            # openmm absent or misbuilt
        return
    faster = [n for n in ("CUDA", "OpenCL") if n in names]
    if faster:
        console.print(
            f"[yellow]MM/GBSA will run on CPU. This host also offers "
            f"{' and '.join(faster)}; `--platform {faster[0]}` is roughly ten "
            "times faster. Do not mix platforms within one table.[/]")

console = Console()

_GRADE_STYLE = {Grade.A: "green", Grade.B: "yellow", Grade.C: "dark_orange", Grade.D: "red"}


# ---------------------------------------------------------------- vocabularies
#
# A closed vocabulary passed as a bare string is not an input contract, it is a
# guess. `--admet-gate hihg` used to be accepted and run the whole funnel with
# the safety gate silently off: the comparison is `severity <= gate` over a
# fixed severity list, an unknown gate is not in that list, so nothing was ever
# excluded and a critical hazard reached the table captioned "worth ordering".
# The same shape of mistake in `--product` made every compound come back "ADMET
# not assessed", which also drops the gate.
#
# So both are click.Choice now: an unknown value exits non-zero and names the
# values that exist. `none` stays a real, distinct value -- it means "run no
# ADMET gate", which is a decision someone may want to make out loud, and it
# must not be confusable with a typo.

class DockEngine(str, Enum):
    'Which docking program runs.'
    vina = "vina"
    hybrid = "hybrid"
    unidock = "unidock"


class ProteinAction(str, Enum):
    any = "any"
    inhibit = "inhibit"
    activate = "activate"


class AdmetGate(str, Enum):
    """Severity at which the ADMET gate drops a compound.

    The members are the severities the hazard comparison actually understands
    (`types.SEVERITY_ORDER`), plus `none`, the explicit "run no gate" value.
    """
    critical = "critical"
    high = "high"
    moderate = "moderate"
    low = "low"
    context = "context"
    none = "none"

    @property
    def gate(self) -> Optional[str]:
        """The value to hand the screen: None when the gate is switched off."""
        return None if self is AdmetGate.none else self.value


#: Built from the producer's own table, so the CLI cannot drift from the product
#: formats `admet_exposure` can actually price an exposure for.
ProductFormat = Enum(  # type: ignore[misc]
    "ProductFormat", {k: k for k in sorted(PRODUCT_USE)}, type=str, module=__name__)

_DEFAULT_PRODUCT = ProductFormat("body_lotion")

# If a severity is ever added upstream without being offered here, this CLI
# would start refusing a legitimate value. Say so loudly at import. Not an
# `assert`: `python -O` strips those, and a guard that cannot fire in the
# configuration it ships under is the failure this file is being fixed for.
_missing = sorted(set(SEVERITY_ORDER) - {m.value for m in AdmetGate})
if _missing:
    raise RuntimeError(
        "types.SEVERITY_ORDER has gained a severity the CLI does not offer: "
        f"{_missing}. Add it to AdmetGate or --admet-gate will refuse it.")


def _applied_percent(value: float) -> float:
    """Validate `--applied-percent` with the exposure model's own rule.

    `--applied-percent -1` used to start a whole screen and report every
    compound at SED = -1.3 mg/kg bw/day in the `negligible` band; `101` and
    `nan` were accepted the same way. The check belongs here, at the option,
    rather than after `run_screen` returns: a use level that is not a
    concentration should cost nothing and write nothing.

    Deliberately not `typer.Option(..., min=0, max=100)`. Click implements a
    range as `value < min or value > max`, both of which are False for NaN, so
    a range option would refuse -1 and 101 and then hand `nan` straight to the
    screen -- a guard that cannot fire on the input that most needs it. Calling
    the producer's validator instead tests finiteness first and keeps one rule
    for the CLI and the Python API.
    """
    try:
        return validate_applied_percent(value)
    except ValueError as e:
        raise typer.BadParameter(str(e)) from e


def _topical_cutoff(value: float) -> float:
    """Click's range accepts NaN; use the same finite gate as run_screen."""
    from .pipeline import validate_topical_cutoff
    try:
        return validate_topical_cutoff(value)
    except ValueError as exc:
        raise typer.BadParameter(str(exc)) from exc


def _setup(verbose: bool) -> None:
    logging.basicConfig(level=logging.INFO if verbose else logging.WARNING,
                        format="%(levelname)s %(name)s: %(message)s")


#: `--core auto`, the value that asks the CLI to look the core up rather than
#: be handed one. Kept as a constant because it is compared in two places and
#: is also a legal, if perverse, SMARTS-shaped string otherwise.
_CORE_AUTO = "auto"


def _resolve_core(value: Optional[str], target: Optional[str]):
    """The pharmacophore core for a `--core` value, or None for no anchoring.

    Why `skinpipe catalog` needs this at all: a 2D fingerprint measures the
    whole molecule, and on a small one the bulk wins. Screening PubChem by
    similarity to 4-butylresorcinol returns simple alkylphenols -- 2-butyl-,
    4-pentyl-, 4-hexylphenol -- because the butylbenzene half of the query is
    most of its fingerprint and dropping one hydroxyl barely moves it. That
    hydroxyl is the entire mechanism: the resorcinol 1,3-diol chelates the
    binuclear copper and a phenol has nothing to chelate with. So the search
    retrieves the scaffold, discards the pharmacophore, and every hit is
    worthless in a way no downstream filter catches, because each one is a
    perfectly reasonable small molecule.

    `skinpipe.pharmacophore` was written to fix precisely that and was reachable
    from nothing in the package, so the shipped command kept running the search
    the module's own docstring calls worthless. This function is the wire.

    A bad value raises `typer.BadParameter`, which exits 2 before any network
    call -- the same contract `--admet-gate` and `--product` hold, and for the
    same reason: a funnel that runs and reports an empty list is
    indistinguishable from a genuine zero-hit screen.
    """
    from . import pharmacophore as ph

    if not value:
        return None
    if value.strip().lower() == _CORE_AUTO:
        gene = (target or "").strip().upper()
        core = ph.core_for(gene=gene) if gene else None
        if core is None:
            named = f"--target {target}" if target else "no --target was given"
            raise typer.BadParameter(
                f"--core auto needs a target this package declares a binding "
                f"core for, and {named}. Declared for: "
                f"{', '.join(sorted(ph.KNOWN_CORES))}. Otherwise pass the core "
                f"as SMARTS, e.g. --core '{ph.KNOWN_CORES['TYR']['smarts']}'.",
                param_hint="--core")
        return core
    core = ph.core_for(smarts=value)
    if core is None or not core.valid:
        raise typer.BadParameter(
            f"{value!r} is not a SMARTS RDKit can compile, so it cannot anchor "
            "a search. Pass a valid SMARTS, or 'auto' to use the declared core "
            "for --target.", param_hint="--core")
    return core


def _spec_or_exit(reg: Registry, marker: str):
    """Resolve a marker, or exit 1 with the message the registry already wrote.

    `Registry._raw` raises a `KeyError` that names every marker in the config
    file -- which is the right message and was reaching the user as a rich
    traceback with that message quoted inside a `KeyError(...)` repr, because
    nothing on the path between it and the terminal was catching it. A mistyped
    marker is ordinary user error, not a crash.
    """
    try:
        return reg.spec(marker)
    except KeyError as e:
        # KeyError stringifies to the repr of its argument, quotes included.
        console.print(f"[red]{rich_escape(str(e.args[0] if e.args else e))}[/]")
        raise typer.Exit(1) from e


@app.command("targets")
def targets_cmd(category: Optional[str] = typer.Option(None, "--category", "-c"),
                resolve: bool = typer.Option(True, help="query UniProt/ChEMBL/PDB"),
                budget: float = typer.Option(
                    300.0, "--budget", min=0.0,
                    help="seconds to spend on live lookups before reporting the "
                         "rest as unresolved; 0 removes the ceiling"),
                verbose: bool = typer.Option(False, "--verbose", "-v")) -> None:
    """Triage every marker in the registry."""
    _setup(verbose)
    # The ceiling is on lookups, not on the command: markers already resolved
    # are still printed, and the ones past the line are printed as "?" -- the
    # same cell an API failure produces, because they mean the same thing to a
    # reader (no counts were retrieved) and differ only in what to do about it,
    # which the footer says.
    reg = Registry(budget_s=budget or None)
    specs = reg.specs(category=category, resolve=resolve, offline=not resolve)
    tbl = Table(title="Skin marker triage", header_style="bold")
    for col in ("marker", "category", "dir", "mode", "node", "PDB", "ligands", "grade"):
        tbl.add_column(col, justify="right" if col in ("PDB", "ligands") else "left")
    n_blind = 0
    for s in specs:
        t = triage.assess(s)
        node = next((n for n in s.nodes if n.gene == t.screening_node), None)
        g = t.grade
        if t.counts_unresolved:
            # A grade computed from counts nobody fetched is not a grade. "D"
            # in this column reads as "not a small-molecule target", which is a
            # conclusion, and there is no evidence for it here.
            n_blind += 1
            cell, pdb, lig = "[dim]?[/]", "—", "—"
        else:
            cell = f"[{_GRADE_STYLE[g]}]{g.value}{'*' if t.provisional else ''}[/]"
            pdb = str(node.n_pdb) if node else "0"
            lig = str(node.n_ligands_p6) if node else "0"
        tbl.add_row(s.marker, s.category, s.direction.value, s.mode.value,
                    t.screening_node or "—", pdb, lig, cell)
    console.print(tbl)
    console.print("[dim]* provisional: graded on raw ChEMBL counts. "
                  "Run `skinpipe screen` to re-grade on topically deliverable ligands.[/dim]")
    if n_blind:
        console.print(f"[yellow]? {n_blind} marker(s) could not be resolved; their "
                      f"counts were never retrieved. Not a verdict -- re-run "
                      f"(or drop --no-resolve) once the APIs answer.[/]")
    starved = reg.skipped_for_budget
    if starved:
        # Said separately because the two silences call for different actions.
        # An API failure means wait for the endpoint; a spent budget means the
        # run stopped early and a longer one will get further, immediately.
        console.print(
            f"[yellow]{len(starved)} gene(s) were never looked up at all: the "
            f"{budget:.0f}s budget ran out first ({', '.join(starved[:6])}"
            f"{', ...' if len(starved) > 6 else ''}). That is this run giving "
            f"up, not an answer about them. What did resolve is cached, so "
            f"`--budget 0` or a plain re-run picks up where this stopped.[/]")


@app.command("explain")
def explain_cmd(marker: str, verbose: bool = typer.Option(False, "--verbose", "-v")) -> None:
    """Why a marker is or is not a small-molecule target."""
    _setup(verbose)
    reg = Registry()
    spec = _spec_or_exit(reg, marker)
    t = triage.assess(spec)
    console.print(f"[bold]{spec.marker}[/] — {spec.category}, wants to {spec.direction.value}")
    console.print(Registry.explain(spec))
    console.print()
    for n in spec.nodes:
        console.print(f"  [cyan]{n.gene}[/] {n.uniprot or ''}  "
                      f"PDB={n.n_pdb}  ligands={n.n_ligands_p6}\n     {n.rationale}")
    console.print(f"\n[bold]Recommendation[/] {t.recommendation}")
    for w in t.warnings:
        console.print(f"  [yellow]![/] {w}")


@app.command("route")
def route_cmd(marker: str, verbose: bool = typer.Option(False, "--verbose", "-v")) -> None:
    """Which discovery route this marker deserves, and why."""
    _setup(verbose)
    from . import strategy
    reg = Registry()
    spec = _spec_or_exit(reg, marker)
    d = strategy.choose_track(spec, triage.assess(spec))
    console.print(strategy.summarise(d))


@app.command("assess")
def assess_cmd(smiles: str,
               target: Optional[str] = typer.Option(None, "--target", "-t",
                                                    help="gene, for mechanism-aware alerts"),
               name: str = typer.Option("query", "--name")) -> None:
    """Skin verdict on one molecule: permeation, sensitisation, stability."""
    from .types import LigandRecord
    from rdkit import Chem
    from rdkit.Chem import Descriptors

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        console.print("[red]could not parse that SMILES[/]")
        raise typer.Exit(1)
    if mol.GetNumHeavyAtoms() == 0:
        # `mol is None` is not the whole test: RDKit parses "" into a perfectly
        # valid molecule that simply has no atoms, and every descriptor below
        # then returns 0.0. `skinpipe assess ''` exited 0 and reported MW 0.0,
        # "permeates with a suitable vehicle", "no electrophilic alert" and a
        # topical score of 9.0/10 -- a confident skin verdict about nothing,
        # and the most persuasive one the tool can print. `library.standardize`
        # already refuses a zero-heavy-atom structure; this is the same rule at
        # the other end of the package.
        console.print(f"[red]{rich_escape(repr(smiles))} parses but has no heavy "
                      "atoms, so there "
                      "is nothing to assess. An empty molecule scores 0 on every "
                      "descriptor, which reads as a clean, permeating compound.[/]")
        raise typer.Exit(1)
    props = {"mw": round(Descriptors.MolWt(mol), 2), "logp": round(Descriptors.MolLogP(mol), 2),
             "tpsa": round(Descriptors.TPSA(mol), 2), "hbd": Descriptors.NumHDonors(mol),
             "hba": Descriptors.NumHAcceptors(mol)}
    lig = LigandRecord(ligand_id=name, smiles=smiles, name=name, props=props)
    skin = topical_f.assess(lig)
    sens = safety_f.summarise(smiles, target_gene=target)
    stab = stability_f.summarise(smiles)

    console.print(f"[bold]{rich_escape(name)}[/]  MW {props['mw']}  logP {props['logp']}  TPSA {props['tpsa']}")
    console.print(f"  permeation   log Kp {skin.log_kp}  [{skin.permeation_class}]  "
                  f"topical score {skin.topical_score}/10")
    console.print(f"               {skin.verdict}")
    console.print(f"  sensitisation burden {sens['burden']}/10 — {sens['verdict']}")
    for a in sens["alerts"]:
        console.print(f"     [yellow]{a['severity']}[/] {a['name']} ({a['domain']}) {a['note']}")
    console.print(f"  stability    {stab['verdict']}")
    for a in stab["alerts"]:
        console.print(f"     {a['mode']}: {a['name']} → {a['mitigation']}")


@app.command("audit")
def audit_cmd(gene: str, pchembl_min: float = 6.0, max_records: int = 2000,
              verbose: bool = typer.Option(False, "--verbose", "-v")) -> None:
    """Assay provenance audit: how much of a target's data is really human?"""
    _setup(verbose)
    from .fetch import chembl, provenance
    u = chembl.uniprot_for_gene(gene)
    if not u:
        console.print(f"[red]no reviewed UniProt entry for {rich_escape(gene)}[/]")
        raise typer.Exit(1)
    tid = chembl.chembl_target_for_accession(u["accession"])
    if not tid:
        console.print(f"[red]no human ChEMBL target for {rich_escape(gene)} ({u['accession']})[/]")
        raise typer.Exit(1)
    acts, meta = chembl.fetch_activities(tid, pchembl_min=pchembl_min,
                                         max_records=max_records, target_gene=gene,
                                         with_meta=True)
    a = provenance.audit(acts, gene)
    console.print(f"[bold]{rich_escape(gene)}[/] {u['accession']} → {tid}   {a['n']} activities "
                  f"at pChEMBL ≥ {pchembl_min}")
    if not meta["complete"]:
        # Percentages computed from a truncated pull look exactly like
        # percentages computed from a complete one.
        console.print(f"[red]incomplete pull: ChEMBL reports {meta['total_count']} "
                      f"activities and only {meta['retrieved']} were retrieved. "
                      f"The percentages below are unreliable; re-run.[/]")
    elif meta["capped_at_max_records"]:
        console.print(f"[yellow]capped at --max-records {max_records} of "
                      f"{meta['total_count']} available[/]")
    tbl = Table(header_style="bold")
    tbl.add_column("provenance")
    tbl.add_column("n", justify="right")
    tbl.add_column("%", justify="right")
    for k, v in sorted(a["counts"].items(), key=lambda t: -t[1]):
        tbl.add_row(k, str(v), f"{100*v/max(a['n'],1):.1f}")
    console.print(tbl)
    console.print(f"human-target activity eligible: [bold]{a['trustworthy_fraction']*100:.1f}%[/]")
    console.print(f"human target in nonhuman cellular systems: {a.get('n_human_target_nonhuman_system', 0)}; "
                  f"annotated human skin functional context: {a.get('n_human_skin_function', 0)}. "
                  "These context classifications do not verify binding or skin efficacy.")
    if a["warning"]:
        console.print(f"[yellow]{a['warning']}[/]")


@app.command("catalog")
def catalog_cmd(query: str,
                threshold: int = typer.Option(80, "--threshold", "-t", min=0, max=100,
                                              help="PubChem 2D similarity cutoff (0-100)"),
                max_records: int = typer.Option(500, "--max-records", min=1),
                min_topical: float = typer.Option(6.0, "--min-topical", min=0.0, max=10.0,
                                                  callback=_topical_cutoff,
                                                  help="topical score cutoff (0-10)"),
                min_vendors: int = typer.Option(1, "--min-vendors", min=0),
                shape_top: int = typer.Option(60, "--shape-top", min=1,
                                              help="how many survivors get 3D treatment"),
                n_confs: int = typer.Option(10, "--confs", min=1),
                admet_gate: AdmetGate = typer.Option(
                    AdmetGate.critical, "--admet-gate", case_sensitive=False,
                    help="drop compounds at this ADMET severity or worse; "
                         "'none' runs no ADMET gate at all"),
                target: Optional[str] = typer.Option(
                    None, "--target",
                    help="gene symbol the search is for, e.g. TYR. A hazard the "
                         "target itself activates cannot be judged without it: at "
                         "TYR/TYRP1 a 4-substituted phenol is escalated from a "
                         "background alert to the mechanism of harm."),
                core_opt: Optional[str] = typer.Option(
                    None, "--core",
                    help="anchor the search on a pharmacophore core: 'auto' for "
                         "the declared core of --target, or a SMARTS. Without "
                         "it the search is a bare 2D similarity query, in which "
                         "the scaffold outvotes the pharmacophore -- a query on "
                         "4-butylresorcinol returns alkylphenols that cannot "
                         "chelate the copper and that no later filter rejects."),
                out: Optional[Path] = typer.Option(None, "--out", "-o"),
                top_n: int = typer.Option(25, "--top", min=1),
                verbose: bool = typer.Option(False, "--verbose", "-v")) -> None:
    """Find purchasable compounds that match a known binder's shape.

    For targets a structure-based screen cannot serve: no usable experimental
    structure, or seeds too large to hop by fragment generation. Give it a
    validated binder and it returns compounds a vendor sells.

    With `--core` the first stage is a substructure query on the pharmacophore
    instead of a similarity query on the whole molecule, so nothing lacking the
    binding motif can enter the funnel at all.
    """
    _setup(verbose)
    from . import catalog as cat
    from . import pharmacophore as ph

    # Settled before anything is fetched: a core that does not compile should
    # cost nothing, like every other closed-vocabulary value in this CLI.
    core = _resolve_core(core_opt, target)
    try:
        if core is not None:
            hits = ph.screen(core, query_smiles=query, max_records=max_records,
                             min_topical=min_topical, min_vendors=min_vendors,
                             shape_top_n=shape_top, n_confs=n_confs,
                             admet_gate=admet_gate.gate, target_gene=target)
        else:
            hits = cat.screen(query, threshold=threshold, max_records=max_records,
                              min_topical=min_topical, min_vendors=min_vendors,
                              shape_top_n=shape_top, n_confs=n_confs,
                              admet_gate=admet_gate.gate, target_gene=target)
    except cat.PubChemBusy as e:
        console.print(f"[yellow]{rich_escape(str(e))}[/]")
        raise typer.Exit(2)
    except ValueError as e:
        # e.g. a query SMILES that will not parse, or one with no heavy atoms.
        # Refusing it is the point: running the funnel and printing "nothing
        # survived" would be a transport-shaped answer to an input error.
        console.print(f"[red]{rich_escape(str(e))}[/]")
        raise typer.Exit(2)



    _ctx = (f"ADMET gate: {admet_gate.value}"
            + (f"  ·  target: {target}" if target else "  ·  target: none given"))








    _core_ctx = (f"anchored on core {rich_escape(core.smarts)} ({core.source})"
                 if core else "unanchored 2D similarity search")
    if not hits:
        console.print("[yellow]no purchasable compound survived the funnel[/]")
        console.print(f"[dim]{_ctx}[/dim]")
        console.print(f"[dim]{_core_ctx}[/dim]")
        raise typer.Exit(1)

    has_2d = any(h.get("tanimoto_2d") is not None for h in hits)
    tbl = Table(title=f"purchasable shape matches ({len(hits)})", header_style="bold")
    for col, just in (("CID", "left"), ("shape", "right"), ("2D", "right"),
                      ("MW", "right"), ("log Kp", "right"), ("vendors", "right"),
                      ("hazard", "left")):
        tbl.add_column(col, justify=just)
    for h in hits[:top_n]:
        # An anchored run has no 2D Tanimoto to the query, and an overlay that
        # failed has no shape. Neither absence is a zero, and printing one as
        # 0.000 would rank a compound on a number nobody measured.
        t2d = h.get("tanimoto_2d")
        shape_val = h.get("shape_tanimoto")
        novel = t2d is not None and (t2d or 1.0) < 0.6
        shape = f"{shape_val:.3f}" if shape_val is not None else "—"
        highlight = novel and shape_val is not None and shape_val >= 0.7
        tbl.add_row(str(h["cid"]),
                    f"[green]{shape}[/]" if highlight else shape,
                    f"{t2d:.2f}" if t2d is not None else "—",
                    f"{h['mw']:.0f}" if h.get("mw") is not None else "—",
                    f"{h['log_kp']:.2f}" if h.get("log_kp") is not None else "—",
                    str(h["n_vendors"]), h["admet_hazard"])
    console.print(tbl)
    if has_2d:
        console.print("[dim]green = shape-similar but 2D-distant, which is the "
                      "combination worth ordering[/dim]")
    console.print(f"[dim]{_ctx}[/dim]")
    console.print(f"[dim]{_core_ctx}[/dim]")
    # The N-13 state, surfaced. A dash in the shape column is a measurement
    # that failed, not a compound with no shape, and those rows are pinned to
    # the bottom of the table rather than scored on anything else.
    n_unmeasured = sum(1 for h in hits if h.get("shape_status") == ph.SHAPE_FAILED)
    if n_unmeasured:
        console.print(f"[yellow]{n_unmeasured} compound(s) show '—' for shape: the "
                      f"3D overlay failed. They are ranked last as unmeasured, "
                      f"not scored on their skin fitness instead.[/]")
    if core is None and target and target.strip().upper() in ph.KNOWN_CORES:
        # The C-05 defect, said out loud on the run that has it. A similarity
        # search on a small molecule is outvoted by its own scaffold, and this
        # package knows the motif that does the binding for this target.
        console.print(f"[yellow]this was an unanchored similarity search, but a "
                      f"binding core is declared for {target.strip().upper()}: "
                      f"{ph.KNOWN_CORES[target.strip().upper()]['why']}. Re-run "
                      f"with --core auto to exclude compounds that lack it.[/]")
    if not target:
        console.print("[yellow]no --target given: hazards were judged out of "
                      "context. A 4-substituted phenol reads as a background "
                      "alert here and is the mechanism of harm on a TYR/TYRP1 "
                      "programme.[/]")
    if out:
        import csv as _csv
        out = Path(out)
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w", newline="") as fh:
            w = _csv.DictWriter(fh, fieldnames=list(hits[0].keys()))
            w.writeheader()
            for h in hits:
                w.writerow(h | {"vendors": "; ".join(h["vendors"])})
        console.print(f"[green]written[/] {out} ({len(hits)} rows)")


@app.command("natural")
def natural_cmd(
    smiles: Optional[str] = typer.Argument(
        None, help="a structure to ask about; omit to report index status only"),
    search: bool = typer.Option(
        False, "--search",
        help="find natural products SIMILAR to this structure instead of asking "
             "whether this exact one is known"),
    threshold: float = typer.Option(0.35, "--threshold", min=0.0, max=1.0),
    limit: int = typer.Option(20, "--limit", min=1),
    pathway: Optional[str] = typer.Option(None, "--pathway",
                                          help="NPClassifier pathway filter"),
    max_mass: Optional[float] = typer.Option(None, "--max-mass"),
    build: bool = typer.Option(False, "--build",
                               help="build the index from the LOTUS export"),
    force: bool = typer.Option(False, "--force", help="rebuild even if current"),
) -> None:
    'Ask LOTUS whether a structure occurs in nature, and in what.'
    _setup(False)
    from skinpipe.fetch import lotus_local as L

    if build:
        try:
            info = L.build_index(force=force)
        except Exception as e:                       # noqa: BLE001
            console.print(f"[red]index build failed: {type(e).__name__}: {e}[/]")
            raise typer.Exit(1) from e
        console.print(f"[green]index built[/] {info}")
        return

    status = L.verify_index()
    console.print(f"index: [bold]{status.get('status')}[/] "
                  f"rows={status.get('rows')} {L.index_path()}")
    if not status.get("ok"):
        console.print(f"[yellow]{status.get('reason', '')}[/]")
        console.print("[yellow]`skinpipe natural --build` builds it from the "
                      f"LOTUS export (${L.ENV_VAR} points at the archive).[/]")
        if smiles:
            raise typer.Exit(1)
        return
    if not smiles:
        return

    if search:
        hits = L.similar_to(smiles, threshold=threshold, limit=limit,
                            pathway=pathway, max_mass=max_mass)
        if not hits:
            console.print(f"no natural product within Tanimoto {threshold:g}. "
                          "That is a statement about this index at this "
                          "threshold, not about nature.")
            return
        t = Table(title=f"{len(hits)} natural product(s) at Tanimoto >= {threshold:g}")
        for col in ("name", "similarity", "pathway", "organisms"):
            t.add_column(col)
        for h in hits[:limit]:
            orgs = h.get("organisms") or []
            t.add_row(str(h.get("name") or "")[:38],
                      f"{h.get('similarity', 0):.3f}",
                      str(h.get("pathway") or "")[:24],
                      ", ".join(orgs[:3])[:52])
        console.print(t)
        return

    key = _inchikey(smiles)
    if key is None:
        console.print("[red]could not read that SMILES[/]")
        raise typer.Exit(2)
    found = L.lookup([key])
    hit = found.get(key)
    if not hit:
        # An InChIKey carries stereochemistry in its second block, so a query
        # built from a SMILES typed without stereocentres misses a compound the
        # index holds -- and reads exactly like a real absence. Ask the skeleton
        # before reporting nothing.
        near = L.lookup_skeleton([key]).get(key) or []
        if near:
            console.print(f"[green]this skeleton occurs in nature[/] "
                          f"({key.split('-')[0]}), in {len(near)} recorded "
                          "form(s). The query's stereochemistry matches none of "
                          "them exactly, and stereoisomers are different "
                          "compounds:")
            for h in near[:6]:
                orgs = ", ".join((h.get("organisms") or [])[:3])
                console.print(f"  - {h.get('name') or '(unnamed)'}  "
                              f"[{h['inchikey']}]  {orgs}")
            return
        console.print(f"[yellow]not found[/] in LOTUS ({key}), and no compound "
                      "shares its skeleton. An occurrence database has finite "
                      "coverage, so this is 'no record here', not 'not a "
                      "natural product'.")
        return
    orgs = hit.get("organisms") or []
    console.print(f"[green]natural product[/]: {hit.get('name') or '(unnamed)'}")
    console.print(f"  pathway    : {hit.get('pathway') or '-'}")
    console.print(f"  superclass : {hit.get('superclass') or '-'}")
    console.print(f"  organisms  : {hit.get('n_organisms')} reported"
                  + (f" -- {', '.join(orgs[:6])}" if orgs else ""))


def _inchikey(smiles: str) -> Optional[str]:
    'InChIKey for a SMILES, or None when RDKit cannot read it.'
    from rdkit import Chem, RDLogger
    from rdkit.Chem import inchi
    RDLogger.DisableLog("rdApp.*")
    mol = Chem.MolFromSmiles(smiles)
    return inchi.MolToInchiKey(mol) if mol is not None else None


@app.command("screen")
def screen_cmd(marker: str,
               protein: bool = typer.Option(False, "--protein", help="resolve a human gene, primary UniProt accession or ChEMBL target ID directly"),
               action: ProteinAction = typer.Option(ProteinAction.any, "--action", help="desired action for --protein; a request, not validated skin biology"),
               out: Path = typer.Option(Path("runs"), "--out", "-o"),
               dock: bool = typer.Option(True, help="run AutoDock Vina"),
               engine: DockEngine = typer.Option(
                   DockEngine.vina, "--engine", case_sensitive=False,
                   help="Docking engine: vina uses CPU AutoDock Vina; hybrid "
                        "uses Uni-Dock for ligands with up to 8 rotatable bonds "
                        "and Vina for the remainder. GPU scores pass the configured "
                        "plausibility checks. unidock uses GPU for benchmarking."),
               unidock_bin: Optional[str] = typer.Option(
                   None, "--unidock-bin",
                   help="path to the `unidock` binary for --engine hybrid or "
                        "unidock. Defaults to $UNIDOCK_BIN, then PATH. It ships "
                        "as a separate conda package with a CUDA build, so it "
                        "is usually in its own environment."),
               gbsa: bool = typer.Option(False, help="MM/GBSA rescore the top poses"),
               gbsa_top: int = typer.Option(20, "--gbsa-top", min=1),
               pdb: Optional[str] = typer.Option(None, help="force a specific PDB ID"),
               ligand_code: Optional[str] = typer.Option(
                   None, "--ligand-code",
                   help="Cocrystal ligand that defines the docking box. "
                        "Choose explicitly when more than one ligand is present."),
               pchembl_min: Optional[float] = typer.Option(
                   None, "--pchembl-min", min=0.0, max=14.0,
                   help="pChEMBL floor. Default: the screening node's own "
                        "potency_floor from the registry, else 5.0. A kinase "
                        "and a copper oxidase do not share a threshold."),
               max_chembl: int = typer.Option(
                   0, "--max-chembl", min=0,
                   help="activities to pull. 0 (the default) requests every live "
                        "page, including when the advisory count is unavailable. "
                        "The acquisition receipt records coverage. A positive "
                        "cap smaller than the query is an activity-ID-ordered "
                        "slice; it does not select the most potent compounds."),
               dock_top: int = typer.Option(
                   900, "--dock-top", min=0,
                   help="candidates to dock, selected by --dock-selection. "
                        "The default balanced policy reserves exploration slots. "
                        "Docking is the wall clock -- ~5 ligands/min at "
                        "exhaustiveness 16 on 32 cores -- so 900 is about three "
                        "hours. 0 docks every candidate that passed the gates; "
                        "any cap is disclosed on the report."),
               dock_selection: str = typer.Option("balanced", "--dock-selection",
                   help="balanced reserves unmeasured/diverse candidates; potency reproduces legacy selection"),
               final_top: int = typer.Option(
                   300, "--final-top", min=0,
                   help="rows in candidate_summary.csv. 0 includes every ranked candidate. "
                        "candidates.csv always preserves the full acquired discovery population."),
               min_topical: float = typer.Option(4.0, "--min-topical", min=0.0, max=10.0,
                                                 callback=_topical_cutoff,
                                                 help="topical score cutoff (0-10)"),
               exhaustiveness: int = typer.Option(16, min=1),
               admet_gate: AdmetGate = typer.Option(
                   AdmetGate.critical, "--admet-gate", case_sensitive=False,
                   help="drop compounds at this ADMET severity or worse; "
                        "'none' runs no ADMET gate at all"),
               applied_percent: float = typer.Option(
                   1.0, "--applied-percent", callback=_applied_percent,
                   help="use level for the exposure estimate, as a percent of "
                        f"the formulation ({APPLIED_PERCENT_MIN:g}-"
                        f"{APPLIED_PERCENT_MAX:g}); SED is linear in it"),
               product: ProductFormat = typer.Option(
                   _DEFAULT_PRODUCT, "--product", case_sensitive=False,
                   help="product format the exposure estimate is priced for"),
               purchasable_only: bool = typer.Option(
                   False, help="keep only compounds a vendor actually sells"),
               exclude_patented: bool = typer.Option(
                   False, help="drop compounds found in the patent literature"),
               no_sourcing: bool = typer.Option(
                   False, help="skip the vendor and patent lookup"),
               allow_surrogate: bool = typer.Option(
                   False, help="keep activities that failed provenance QC"),
               library_file: Optional[Path] = typer.Option(None, "--library", "-l"),
               platform: Optional[str] = typer.Option(
                   None, "--platform",
                   help="OpenMM platform for MM/GBSA: CUDA, OpenCL, CPU, "
                        "Reference. Default CPU. Platform changes can materially "
                        "change energies; do not mix platforms in one table."),
               gbsa_cpu_profile: str = typer.Option(
                   "auto", "--gbsa-cpu-profile",
                   help="MM/GBSA CPU profile: auto (one thread unless explicitly "
                        "set), serial (one thread), parallel (split cores). "
                        "Multi-thread numerical repeats may differ."),
               gbsa_threads: Optional[int] = typer.Option(
                   None, "--gbsa-threads", min=1,
                   help="Explicit CPU threads per MM/GBSA pose."),
               gbsa_charge_method: str = typer.Option(
                   "am1bcc-fixed", "--gbsa-charge-method",
                   help="MM/GBSA charges. Default am1bcc-fixed uses a recorded "
                        "SQM protocol and fails without substituting another method. "
                        "Use auto explicitly for the legacy adaptive ladder."),
               gbsa_ligand_forcefield: str = typer.Option(
                   "openff-2.2.0", "--gbsa-ligand-forcefield",
                   help="Versioned ligand force field for MM/GBSA; default openff-2.2.0."),
               top_n: int = typer.Option(25, "--top", min=1),
               verbose: bool = typer.Option(False, "--verbose", "-v")) -> None:
    """Run the full screen for one marker."""
    _setup(verbose)
    from .mmgbsa.engine import cpu_threads
    from .mmgbsa.ligand import validate_screen_protocol
    try:
        cpu_threads(gbsa_threads, gbsa_cpu_profile)
        gbsa_charge_method, gbsa_ligand_forcefield = validate_screen_protocol(
            gbsa_charge_method, gbsa_ligand_forcefield)
    except (TypeError, ValueError) as exc:
        raise typer.BadParameter(str(exc), param_hint="MM/GBSA settings") from exc
    from .pipeline import run_screen
    registry_args = {}
    if protein:
        from .fetch.target_identity import resolve, TargetIdentityError
        from .protein_registry import ProteinRegistry
        try:
            identity = resolve(marker)
        except TargetIdentityError as exc:
            console.print(f"[red]{rich_escape(str(exc))}[/]")
            for choice in exc.choices:
                console.print(rich_escape(str(choice)))
            raise typer.Exit(2) from exc
        marker = identity["query"]
        registry_args["registry"] = ProteinRegistry(identity, action=action.value)
    elif action is not ProteinAction.any:
        console.print("[red]--action requires --protein; registered marker goals come from the registry.[/]")
        raise typer.Exit(2)
    workdir = Path(out) / marker.replace("/", "_")
    if platform is None and gbsa:
        _note_faster_platform()
    res = run_screen(marker, workdir, dock=dock, gbsa=gbsa, gbsa_top_n=gbsa_top,
                     dock_top=dock_top or None, final_top=final_top or None,
                     dock_selection=dock_selection,
                     engine=engine.value, unidock_bin=unidock_bin,
                     ligand_code=ligand_code,
                     pdb_id=pdb, pchembl_min=pchembl_min, max_chembl=max_chembl,
                     min_topical=min_topical, exhaustiveness=exhaustiveness,
                     trustworthy_only=not allow_surrogate,
                     admet_gate=admet_gate.gate,
                     applied_percent=applied_percent, product=product.value,
                     check_sourcing=not no_sourcing,
                     purchasable_only=purchasable_only,
                     exclude_patented=exclude_patented,
                     platform=platform,
                     gbsa_cpu_profile=gbsa_cpu_profile, gbsa_threads=gbsa_threads,
                     gbsa_charge_method=gbsa_charge_method,
                     gbsa_ligand_forcefield=gbsa_ligand_forcefield,
                     extra_files=[library_file] if library_file else None, **registry_args)
    md = report_mod.to_markdown(res, workdir / "report.md", top_n=top_n)
    report_mod.candidate_exports(res, workdir, summary_top=final_top or None)
    _, n_order = report_mod.orderable_to_csv(res, workdir / "orderable.csv")
    report_mod.research_lists_to_csv(res, workdir)
    report_mod.to_json(res, workdir / "result.json")




    report_mod.to_manifest(
        res, workdir / "run_manifest.json",
        parameters={"pdb_id": pdb, "dock": dock, "gbsa": gbsa,
                    "target_kind": "protein" if protein else "marker", "requested_action": action.value,
                    "gbsa_top_n": gbsa_top, "exhaustiveness": exhaustiveness,
                    "pchembl_min": pchembl_min, "max_chembl": max_chembl,
                    "dock_top": dock_top, "final_top": final_top, "dock_selection": dock_selection,
                    "engine": engine.value, "unidock_bin": unidock_bin,
                    "gbsa_cpu_profile": gbsa_cpu_profile, "gbsa_threads": gbsa_threads,
                    "gbsa_charge_method": gbsa_charge_method,
                    "gbsa_ligand_forcefield": gbsa_ligand_forcefield,
                    "ligand_code": ligand_code,
                    "min_topical": min_topical, "admet_gate": admet_gate,
                    "applied_percent": applied_percent, "product": product,
                    "purchasable_only": purchasable_only,
                    "exclude_patented": exclude_patented,
                    "trustworthy_only": not allow_surrogate,
                    "extra_files": [str(library_file)] if library_file else None,
                    "check_sourcing": not no_sourcing, "platform": platform})
    console.print(md)
    # All four, because `orderable.csv` is the one a formulator opens and it was
    # the one the message left out -- a file nobody is told about is a file nobody
    # reads. Its row count travels with it: an empty orderable list is a finding
    # only if the gates ran, and `report.md` is where that is said.
    console.print(
        f"\n[green]written[/] {workdir}/report.md, candidates.csv, "
        "ranked_candidates.csv, candidate_summary.csv, all_candidates.csv, "
        f"orderable.csv ({n_order} row{'' if n_order == 1 else 's'}), "
        f"result.json, run_manifest.json")


def main() -> None:
    app()


@app.command("resolve-protein")
def resolve_protein(query: str) -> None:
    """Print a verified human protein identity or explicit ambiguity as JSON."""
    import json
    from .fetch.target_identity import resolve, TargetIdentityError
    try:
        result = resolve(query)
    except TargetIdentityError as exc:
        print(json.dumps(exc.to_dict(), ensure_ascii=False))
        raise typer.Exit(2) from exc
    print(json.dumps(result, ensure_ascii=False))


@app.command("compare-hops")
def compare_hops(parent: str = typer.Option(..., "--parent"),
                 candidate: str = typer.Option(..., "--candidate")) -> None:
    """Describe and colour the exact common structure of a parent and hop."""
    from .hop_compare import compare
    import json
    try:
        comparison = compare(parent, candidate)
    except ValueError as exc:
        raise typer.BadParameter(str(exc)) from exc
    typer.echo(json.dumps(comparison, ensure_ascii=False, allow_nan=False))


@app.command("hop-screen")
def hop_screen(source: Path,
               parent_id: str = typer.Option(..., "--parent-id"),
               out: Path = typer.Option(..., "--out"),
               brics_products: int = typer.Option(200, "--brics-products"),
               max_3d: int = typer.Option(50, "--max-3d"),
               conformers: int = typer.Option(10, "--conformers"),
               max_sa: float = typer.Option(6.0, "--max-sa")) -> None:
    """Generate saved-parent hypotheses and preserve the complete proposal pool.

    BRICS has a raw product budget; bioisostere edits are separately counted.
    max-sa limits optional 3D comparisons, never the exported candidate pool.
    Output must be a new directory. No activity or patent claim is inherited.
    """
    from .hop_export import export_hops
    import json
    try:
        result = export_hops(source, parent_id, out, brics_products=brics_products,
                             max_3d=max_3d, conformers=conformers, max_sa=max_sa)
    except (ValueError, OSError) as exc:
        raise typer.BadParameter(str(exc)) from exc
    typer.echo(json.dumps(result, ensure_ascii=False, allow_nan=False))


if __name__ == "__main__":
    main()
