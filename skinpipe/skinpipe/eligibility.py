"""Can this go in a cosmetic at all? One verdict, with its basis attached.

The judgement currently arrives from three unrelated places -- an element in
the formula, a name on a statutory list, a structural alert read across from a
related substance -- and they are not the same kind of answer. Mixing them
produces advice that is either falsely confident or uselessly cautious.

So every verdict here carries its **basis**, and there are only two that matter:

* **Definitive.** The substance is named on a prohibition list, or it contains
  an element that cosmetics may not contain. Nothing further is needed and no
  data package changes the answer. Mercury is not going into a cream.
* **Read-across.** A structural feature resembles something restricted and
  prompts compound-specific review. An organoboron B-C motif does not inherit
  the classification of boric acid or borates. The screening layer retains its
  annotation without a blanket critical CMR flag; actual compound identity,
  toxicity and cosmetic eligibility still require their own evidence.

A third case sits outside both: the molecule is fine and the **claim** is the
problem. Minoxidil is barred as a drug substance; describing a perfectly legal
botanical as blocking DHT is barred as a claim. Same outcome, different fix.

The statutory list is loaded from ``config/kr_prohibited_ingredients.yaml``
when present. The module works without it and says so, because a missing list
must never read as a clean bill of health -- and "says so" has to be more than
a caveat stapled to a pass. When the list is unusable (``ListStatus`` separates
*absent*, *unreadable* and *parsed but empty*, which have different fixes), a
verdict that rests on the absence of a list hit comes back ``DATA_REQUIRED`` on
a ``NOT_ASSESSED`` basis: no lookup happened, so there is nothing to clear.
Verdicts resting on their own positive evidence -- a prohibited element, a
notified 고시원료, a drug substance -- stand, and carry the caveat. A list that
is *present but
unmatchable in kind* gets the same treatment: the 별표 1 export carries a
Korean name for every row, a CAS for 617 of 624, and a structure for almost
none, so a compound identified by SMILES or an English name is being checked
against keys the list mostly does not have. When nothing matches, the verdict
says which keys were available rather than implying the substance was cleared.
"""
from __future__ import annotations

import functools
import logging
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Optional, Sequence, Union

import yaml
from rdkit import Chem, RDLogger
from rdkit.Chem import inchi

from .regulatory import index_names, name_variants
from .types import concerns_of

RDLogger.DisableLog("rdApp.*")
log = logging.getLogger(__name__)

_PROHIBITED = Path(__file__).resolve().parent.parent / "config" / "kr_prohibited_ingredients.yaml"


class Verdict(str, Enum):
    PROHIBITED = "prohibited"        # barred outright
    RESTRICTED = "restricted"        # permitted under a stated limit
    DRUG_ONLY = "drug_only"          # a drug substance, not a cosmetic ingredient
    DATA_REQUIRED = "data_required"  # structural concern; commission data
    ELIGIBLE = "eligible"
    UNKNOWN = "unknown"


class Basis(str, Enum):
    STATUTE = "statute"              # named on a legal list
    ELEMENT = "element"              # contains a prohibited element
    READ_ACROSS = "read_across"      # structural resemblance only
    NOT_ASSESSED = "not_assessed"    # the list never loaded; nothing was checked
    NONE = "none"

    @property
    def definitive(self) -> bool:
        return self in (Basis.STATUTE, Basis.ELEMENT)


class ListStatus(str, Enum):
    """Why the statutory list is, or is not, available to match against.

    These were one boolean, and the boolean lost the distinction that matters:
    a path that holds no file, a file that cannot be parsed, and a file that
    parses to no entries are three different operational failures with three
    different fixes, and only one of them is even visible in a log. Worse, the
    third used to set ``list_loaded=True`` -- an empty list matches nothing, so
    every compound came back "not on the list", which is exactly how a clean
    pass looks.
    """
    LOADED = "loaded"            # entries are present and searchable
    MISSING = "missing"          # nothing at the path
    UNREADABLE = "unreadable"    # present, but could not be read or parsed
    EMPTY = "empty"              # parsed, but carries no prohibited/restricted rows

    @property
    def usable(self) -> bool:
        """True only when a lookup against this list means anything."""
        return self is ListStatus.LOADED


# Element tables below are checked against 「화장품 안전기준 등에 관한 규정」
# 식약처고시 제2026-19호 (시행 2026.3.18.) [별표 1], read as the law.go.kr text.
#
# An element-atom screen can only express one shape of rule: "this element and
# all of its compounds". Several 별표 1 entries are not that shape. `크롬 ;
# 크로믹애씨드 및 그 염류` names an oxidation state. `셀렌 및 그 화합물(셀레늄
# 아스파테이트는 제외)` carries an exemption. Elsewhere the list writes
# `옥시염화비스머스 이외의 비스머스화합물` and exempts 바륨설페이트 -- set
# differences over substances, not over atoms. No atom test can represent any
# of those, so substance-level matching (InChIKey / CAS / exact name) is the
# only faithful representation, and the atom screen is kept only where the
# regulation really does say "및 그 화합물".

#: Elements whose 별표 1 entry is a blanket ban on the element and all its
#: compounds. Verified one by one against the text; an atom is sufficient.
PROHIBITED_ELEMENTS: dict[str, str] = {
    "Hg": "수은 및 그 화합물", "Pb": "납 및 그 화합물", "As": "비소 및 그 화합물",
    "Cd": "카드뮴 및 그 화합물", "Sb": "안티몬 및 그 화합물", "Tl": "탈륨 및 그 화합물",
    "Be": "베릴륨 및 그 화합물", "Te": "텔루륨 및 그 화합물",
}

#: Elements whose 별표 1 entry is NARROWER than "all compounds of this element",
#: so the presence of the atom is a reason to look the substance up, never a
#: verdict. Reported as DATA_REQUIRED / READ_ACROSS with the statutory text
#: quoted, because the alternative -- dropping them -- would let an unflagged
#: nickel compound through, and that is worse than an over-flagged one.
#: symbol -> (the 별표 1 entry as written, what it does and does not reach)
WATCH_ELEMENTS: dict[str, tuple[str, str]] = {
    "Cr": ("크롬 ; 크로믹애씨드 및 그 염류",
           "금지 대상은 크롬 금속과 크로믹애씨드(6가 크롬) 계열입니다. 3가 크롬 안료는 "
           "오히려 명시적으로 허용되어 있습니다 — 「화장품의 색소 종류와 기준 및 시험방법」 "
           "별표 1 제75호 크로뮴옥사이드그린(CI 77288), 제76호 크로뮴하이드로사이드그린"
           "(CI 77289), 사용제한 '-'."),
    "Ni": ("니켈, 그리고 개별 등재된 니켈 화합물 10종(디하이드록사이드·디옥사이드·"
           "모노옥사이드·설파이드·설페이트·카보네이트·니켈(Ⅱ)트리플루오로아세테이트·"
           "디니켈트리옥사이드·테트라카보닐니켈·트리니켈디설파이드)",
           "별표 1에 '니켈 화합물' 전체를 금지하는 항목은 없습니다. 비의도적 잔류는 제6조"
           "제2항이 한도로 관리합니다 — 눈 화장용 35 μg/g, 색조 30 μg/g, 기타 10 μg/g."),
    "Se": ("셀렌 및 그 화합물(셀레늄아스파테이트는 제외)",
           "조문에 예외가 명시되어 있으므로 셀레늄아스파테이트는 금지 대상이 아닙니다."),
}

#: Removed from the atom screen entirely after checking the current 별표 1:
#: Pt is absent from 별표 1 and 별표 2 alike (and from EU Annex II), so there is
#: no regulatory source for a platinum ban in either jurisdiction. U, Po and Ra
#: are absent by name and are reachable only through the conditional entry
#: 방사성물질(다만, 제품에 포함된 방사능의 농도 등이 「생활주변방사선 안전관리법」
#: 제15조의 규정에 적합한 경우 제외). That is a condition, not an element ban, so
#: it is reported as a caveat. This is not an exhaustive radionuclide screen.
RADIOACTIVE_ELEMENTS: frozenset[str] = frozenset({"U", "Po", "Ra"})


@dataclass
class EligibilityVerdict:
    verdict: Verdict = Verdict.UNKNOWN
    basis: Basis = Basis.NONE
    reason: str = ""
    matched: Optional[dict[str, Any]] = None
    limit_text: Optional[str] = None
    caveats: list[str] = field(default_factory=list)
    #: The 식약처 화장품 원료성분정보 row this substance matches, if any.
    #:
    #: Deliberately not a `verdict` and not a `basis`. Appearing in that
    #: inventory means the regulator recognises the substance as a cosmetic raw
    #: material and has given it a standard name and a CAS; it does not permit
    #: anything. 별표 1 still prohibits what it prohibits -- and the inventory
    #: and the prohibition list are not disjoint, since a substance can be a
    #: recognised raw material whose use is barred. Folding this into the
    #: verdict would turn "this chemotype is already used in cosmetics", which
    #: is a discovery signal, into "this may be used", which is a legal claim
    #: nobody here is entitled to make.
    inventory: Optional[dict[str, Any]] = None

    @property
    def definitive(self) -> bool:
        return self.basis.definitive

    @property
    def usable(self) -> bool:
        """True only when nothing bars it. DATA_REQUIRED is not a yes."""
        return self.verdict in (Verdict.ELIGIBLE, Verdict.RESTRICTED)


def _cas_key(value: Optional[str]) -> str:
    """CAS registry numbers as a comparable key: digits and hyphens only."""
    if not value:
        return ""
    return "".join(ch for ch in str(value) if ch.isdigit() or ch == "-").strip("-")


#: The 식약처 inventory, built by `scripts/build_mfds_inventory.py`. Optional:
#: the file needs a per-account data.go.kr key to produce, so a checkout that
#: has never run that script simply has no inventory signal. Absent is not
#: "nothing matched" -- `inventory_coverage["rows"] == 0` is how a caller tells
#: those apart, and `assess` only reports a hit, never an absence.
_INVENTORY = Path(__file__).resolve().parent.parent / "config" / "kr_cosmetic_ingredients.yaml"


# Both loaders below are cached on the path. They parse committed YAML that
# does not change while the process runs, and `CosmeticEligibility()` is
# constructed per screen and per test -- uncached, the inventory alone cost
# 9.8 s per construction and turned a 57 s test run into a ten-minute one.
# The returned dicts are treated as read-only by every caller here.










_MFDS_RESTRICTIONS = (Path(__file__).resolve().parent.parent
                      / "config" / "mfds_restrictions.yaml")

#: The jurisdiction this module's verdicts are about. The file covers eleven;
#: only Korea decides `Verdict.PROHIBITED` here, and the rest travel on the
#: verdict as context, because a substance banned in the EU and legal in Korea
#: is a fact a formulator needs and is not a Korean prohibition.
_HOME_COUNTRY = "한국"


def _other_jurisdictions(hit: dict) -> str:
    """The same substance's status outside Korea, when the row carries it.

    A formulator shipping to more than one market needs this and cannot get it
    from 별표 1, which is silent about everywhere else. It is appended to the
    reason rather than folded into the verdict: the verdict is a statement
    about Korean law, and a ban in Brazil does not change it.
    """
    proh = [c for c in (hit.get("prohibited_in") or []) if c != _HOME_COUNTRY]
    lim = [c for c in (hit.get("limited_in") or []) if c != _HOME_COUNTRY]
    parts = []
    if proh:
        parts.append(f"금지: {', '.join(proh)}")
    if lim:
        parts.append(f"한도: {', '.join(lim)}")
    return f" 국외 규제 — {' / '.join(parts)}." if parts else ""


@functools.lru_cache(maxsize=4)
def _load_mfds_restrictions(path: Optional[Path] = None
                            ) -> tuple[dict[str, dict], dict[str, dict], dict[str, int]]:
    """``(by_cas, by_name, coverage)`` for the official restriction list."""
    path = Path(path or _MFDS_RESTRICTIONS)
    if not path.is_file():
        log.warning("MFDS restriction list at %s is missing; official API coverage is unavailable", path)
        return {}, {}, {"rows": 0, "cas": 0, "present": 0, "kr_prohibited": 0}
    try:
        data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except (OSError, yaml.YAMLError) as exc:
        log.warning("MFDS restriction list at %s is unreadable: %s", path, exc)
        return {}, {}, {"rows": 0, "cas": 0, "present": 0, "kr_prohibited": 0}
    if not isinstance(data, dict) or not isinstance(data.get("substances"), list):
        log.warning("MFDS restriction list at %s has no valid substances array", path)
        return {}, {}, {"rows": 0, "cas": 0, "present": 0, "kr_prohibited": 0}
    rows: list[dict] = []
    for r in data.get("substances") or []:
        if not isinstance(r, dict):
            continue
        proh = r.get("prohibited_in") or []
        lim = r.get("limited_in") or []
        kr_type = r.get("kr_type") or ""
        kr_conditions = r.get("kr_conditions") or []
        # `_section` is the shape `_list_hit`'s callers already read, and it is
        # what turns a row into a definitive verdict, so it is set only where
        # the row supports one.
        #
        # A 한도 row with no Korean condition text does not. 58 of the 655 such
        # rows are substances on the 색소 table whose 사용제한 column reads "-":
        # 크로뮴옥사이드그린 is listed, and 「화장품의 색소 종류와 기준 및
        # 시험방법」 별표 1 제75호 permits it outright. Calling that RESTRICTED
        # on a STATUTE basis would have the pipeline refuse a pigment the MFDS
        # allows, on the strength of a type code that means "appears in the
        # table", not "is limited". Those rows still travel -- the country
        # lists and the text are on the record -- they just decide nothing.
        if _HOME_COUNTRY in proh or kr_type == "금지":
            section = "prohibited"
        elif (_HOME_COUNTRY in lim or kr_type) and kr_conditions:
            section = "restricted"
        else:
            section = None
        rows.append(dict(r) | {"_section": section, "_source": "mfds_api"})
    by_cas: dict[str, dict] = {}
    for r in rows:
        key = _cas_key(r.get("cas"))
        if key and r["_section"]:
            by_cas.setdefault(key, r)
    by_name = index_names([r for r in rows if r["_section"]], ("en", "ko"))
    if not rows:
        log.warning("MFDS restriction list at %s contains no valid substance rows; official API coverage is unavailable", path)
    return by_cas, by_name, {
        "rows": len(rows), "cas": len(by_cas), "present": int(bool(rows)),
        "kr_prohibited": sum(1 for r in rows if r["_section"] == "prohibited"),
    }


@functools.lru_cache(maxsize=4)
def _load_inventory(path: Optional[Path] = None
                    ) -> tuple[dict[str, dict], dict[str, dict], dict[str, int]]:
    """``(by_cas, by_name, coverage)`` for the registered-raw-material list."""
    path = Path(path or _INVENTORY)
    if not path.is_file():
        return {}, {}, {"rows": 0, "cas": 0, "present": 0}
    try:
        data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except (OSError, yaml.YAMLError) as exc:
        log.warning("cosmetic-ingredient inventory at %s is unreadable: %s",
                    path, exc)
        return {}, {}, {"rows": 0, "cas": 0, "present": 0}
    rows = [r for r in (data.get("ingredients") or []) if isinstance(r, dict)]
    by_cas: dict[str, dict] = {}
    for r in rows:
        key = _cas_key(r.get("cas"))
        if key:
            by_cas.setdefault(key, r)
    by_name = index_names(rows, ("en", "ko"))
    return by_cas, by_name, {"rows": len(rows), "cas": len(by_cas), "present": 1}


class CosmeticEligibility:
    def __init__(self, prohibited_path: Optional[Path] = None,
                 regulations: Optional[Any] = None):
        self.path = Path(prohibited_path or _PROHIBITED)
        self._d: dict[str, Any] = {}
        #: Why the list is or is not usable. Never collapse this to a boolean
        #: at the call site: see `ListStatus`.
        self.list_status: ListStatus = ListStatus.MISSING
        #: The operational detail behind `list_status`, in the caveat's words.
        self.list_detail: str = f"'{self.path}' 경로에 파일이 없습니다"
        if self.path.is_file():
            try:
                self._d = yaml.safe_load(self.path.read_text()) or {}
                self.list_status = ListStatus.LOADED
                self.list_detail = ""
            except (yaml.YAMLError, OSError, UnicodeDecodeError) as exc:
                self.list_status = ListStatus.UNREADABLE
                self.list_detail = f"{type(exc).__name__}: {exc}"
                log.error("could not parse %s: %s", self.path.name, exc)
        # A YAML file is not obliged to hold a mapping. `safe_load` returns a
        # list or a bare string just as happily, and every `_d.get(...)` below
        # would then raise from inside a module whose whole job is to keep a
        # failed lookup from reading as a pass.
        if not isinstance(self._d, dict):
            self.list_status = ListStatus.UNREADABLE
            self.list_detail = (f"최상위 구조가 mapping이 아닙니다 "
                                f"({type(self._d).__name__})")
            self._d = {}
        # Three lookup keys, in descending order of reliability: InChIKey, CAS,
        # then an exact match on an expanded name variant. The coverage of each
        # is measured here rather than assumed, because it is the honest answer
        # to "is this substance on the list?" when the answer comes back no.
        self._by_key: dict[str, dict] = {}
        self._by_cas: dict[str, dict] = {}
        rows: list[dict] = []
        for section in ("prohibited", "restricted"):
            for item in self._d.get(section) or []:
                item = dict(item) | {"_section": section}
                rows.append(item)
                smi = item.get("smiles")
                if smi:
                    mol = Chem.MolFromSmiles(smi)
                    if mol is not None:
                        self._by_key[inchi.MolToInchiKey(mol)] = item
                    else:
                        log.warning("unparseable SMILES on %s", item.get("ko"))
                cas = _cas_key(item.get("cas"))
                if cas:
                    self._by_cas.setdefault(cas, item)
        self._by_name = index_names(rows, ("en", "ko"))
        if self.list_status is ListStatus.LOADED and not rows:
            self.list_status = ListStatus.EMPTY
            self.list_detail = ("파일은 읽혔으나 prohibited·restricted 항목이 0건입니다"
                                f" ('{self.path}')")
        self.coverage = {"rows": len(rows), "structural": len(self._by_key),
                         "cas": len(self._by_cas),
                         "en": sum(1 for r in rows if r.get("en")),
                         "name_corrupted": sum(1 for r in rows
                                               if r.get("name_corrupted"))}
        self._regs = regulations
        # The path is passed rather than defaulted, because the loaders are
        # cached on their argument: called as `_load_inventory()` the key is
        # always `None` while the file they actually read comes from a module
        # global, so the first load would answer for every later one and a
        # test that points the global somewhere else would silently get the
        # earlier file.
        self._inv_by_cas, self._inv_by_name, self.inventory_coverage = \
            _load_inventory(_INVENTORY)
        self._mfds_by_cas, self._mfds_by_name, self.mfds_coverage = \
            _load_mfds_restrictions(_MFDS_RESTRICTIONS)

    # -------------------------------------------------------- list state
    @property
    def list_usable(self) -> bool:
        """True only when a "not on the list" answer carries any information."""
        return self.list_status.usable

    @property
    def list_loaded(self) -> bool:
        """The old boolean, kept for callers that only ask yes/no.

        Deliberately False for `EMPTY` as well: a list with no entries clears
        every compound it is shown, which is the failure this property exists
        to prevent, and answering True for it would be the old bug with a new
        name.
        """
        return self.list_usable

    def list_caveat(self) -> Optional[str]:
        """The sentence a verdict must carry when the list could not be used.

        None when the list loaded. The state and the operational detail are
        both in the text, because this string is what survives into the CSV,
        the JSON and the Markdown, and "the list is missing" without saying
        *how* sends the reader back to the pipeline to find out.
        """
        if self.list_usable:
            return None
        return ("별표 1 사용금지 원료 목록이 로드되지 않았습니다 "
                f"(상태 {self.list_status.value}: {self.list_detail}). 원소 검사와 "
                "구조 경보만 적용됐으므로 '해당 없음'을 '사용 가능'으로 읽으면 안 됩니다. "
                "이 화면의 규제 판정은 미검증 상태이며, 목록을 복구해 재실행해야 합니다.")

    def _note_list_state(self, v: "EligibilityVerdict") -> "EligibilityVerdict":
        """Attach the list's state to a verdict that did not need the list.

        An element prohibition or a drug-substance finding stands on its own
        evidence and the verdict does not change. The reader still has to know
        that the statutory list behind the rest of this screen never loaded,
        and these paths used to return before any caveat was attached.
        """
        c = self.list_caveat()
        if c:
            v.caveats.append(c)
        if not self.mfds_coverage.get("present"):
            v.caveats.append("식약처 화장품 사용제한 원료정보 API 목록을 불러오지 못했습니다. "
                             "큐레이션 목록의 판정은 유지하지만 공식 API에만 있는 원료의 제한은 "
                             "확인되지 않았습니다.")
        v.caveats = list(dict.fromkeys(v.caveats))
        return v

    # ------------------------------------------------------------- helpers
    def _element_hit(self, mol: Chem.Mol) -> Optional[tuple[str, str]]:
        for atom in mol.GetAtoms():
            sym = atom.GetSymbol()
            if sym in PROHIBITED_ELEMENTS:
                return sym, PROHIBITED_ELEMENTS[sym]
        return None

    def _watch_hit(self, mol: Chem.Mol) -> Optional[tuple[str, str, str]]:
        for atom in mol.GetAtoms():
            sym = atom.GetSymbol()
            if sym in WATCH_ELEMENTS:
                return (sym,) + WATCH_ELEMENTS[sym]
        return None

    def _list_hit(self, mol: Chem.Mol, name: Optional[str],
                  cas: Optional[Union[str, Sequence[str]]] = None) -> Optional[dict]:
        """Structure, then CAS, then an exact match on an expanded name.

        CAS is the route that actually works on this data: 619 of the 626
        별표 1 rows carry a registry number, while only 15 carry a SMILES and
        those 15 collapse to 12 distinct InChIKeys. Structure is still tried
        first, because an InChIKey match is the only exact one available, but
        it can never hit more than those 12.

        The name route requires exact equality on a complete name, for the
        reasons set out in ``regulatory.name_variants``.
        """
        hit = self._by_key.get(inchi.MolToInchiKey(mol))
        if hit:
            return hit
        # Every registry number, not the first. PubChem lists several synonyms
        # for one structure -- resveratrol has 501-36-0 and 31100-06-8,
        # niacinamide has three -- and which of them 별표 1 filed the substance
        # under is not knowable in advance. Checking one and stopping would make
        # the answer depend on synonym ordering.
        for one in ([cas] if isinstance(cas, str) or cas is None else list(cas)):
            hit = self._by_cas.get(_cas_key(one))
            if hit:
                return hit
        for variant in name_variants(name, split_groups=False):
            hit = self._by_name.get(variant)
            if hit:
                return hit
        # The official list, after the curated one. Order matters: where both
        # hold a substance the curated row is the richer record -- it carries
        # the structure and the hand-checked reading of the clause -- and this
        # one is what answers for the 1,554 CAS numbers it does not hold at all.
        for one in ([cas] if isinstance(cas, str) or cas is None else list(cas)):
            hit = self._mfds_by_cas.get(_cas_key(one))
            if hit:
                return hit
        for variant in name_variants(name, split_groups=False):
            hit = self._mfds_by_name.get(variant)
            if hit:
                return hit
        return None

    # -------------------------------------------------------------- assess
    def assess(self, smiles: str, name: Optional[str] = None,
               admet_profile: Optional[Any] = None, *,
               cas: Optional[Union[str, Sequence[str]]] = None
               ) -> EligibilityVerdict:
        """The verdict, plus the inventory note, on every path.

        `_assess` below returns early on five of them -- a prohibited element,
        a 별표 1 hit, a 고시원료, a drug substance, an unparseable structure --
        and the inventory lookup used to sit after the last of them, so
        `verdict.inventory` was filled in only when nothing decisive happened.
        A field that is None both for "not a registered raw material" and for
        "we never got that far" is the missing-versus-absent confusion this
        module exists to avoid, so the annotation happens here, once, where
        every return passes through.
        """
        v = self._assess(smiles, name, admet_profile, cas=cas)
        inv = self._inventory_hit(name, cas)
        if inv and v.inventory is None:
            v.inventory = inv
            label = inv.get("ko") or inv.get("en") or ""
            v.caveats.append(
                f"식약처 화장품 원료성분정보에 등재된 원료입니다({label}). 등재는 "
                "표준명·CAS가 부여되어 화장품 원료로 인지된다는 사실이며, 사용 "
                "허가가 아닙니다. 사용 가부는 위 판정을 따릅니다.")
        return self._note_list_state(v)

    def _assess(self, smiles: str, name: Optional[str] = None,
                admet_profile: Optional[Any] = None, *,
                cas: Optional[Union[str, Sequence[str]]] = None
                ) -> EligibilityVerdict:
        """Definitive checks first, then read-across, then eligible.

        ``admet_profile`` may be an :class:`~skinpipe.types.ADMETProfile` or
        the ``to_dict()`` form of one -- the pipeline stores the dict on the
        candidate, and a gate that reads only the object shape is a gate that
        is off for every real screen.
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return EligibilityVerdict(Verdict.UNKNOWN, Basis.NONE,
                                      "구조를 해석할 수 없어 판정 불가")

        v = EligibilityVerdict()

        # 1. Element -- structural, exact, needs no list. Only the eight
        #    elements whose 별표 1 entry really is "및 그 화합물" reach this path.
        el = self._element_hit(mol)
        if el:
            sym, ko = el
            return self._note_list_state(EligibilityVerdict(
                Verdict.PROHIBITED, Basis.ELEMENT,
                f"{sym}를 포함합니다 ({ko}). 화장품 사용금지 원소이므로 자료로 해결되지 "
                "않습니다."))

        # 2. Statutory list.
        hit = self._list_hit(mol, name, cas)
        if hit:
            section = hit.get("_section")
            # Where the hit came from is part of the finding. A row from the
            # official service has clean names, a CAS and the other
            # jurisdictions; a row from the curated parse has the structure and
            # the read of the clause. Saying which one answered lets a reader
            # go to the right place to check it.
            other = _other_jurisdictions(hit)
            # `_note_list_state` even on a hit. A hit used to imply the curated
            # list had loaded, because it was the only list; the official
            # service is a second one and can answer while the curated file is
            # missing. The verdict stands on its own evidence, but the reader
            # still has to be told that the element tables and the exemption
            # readings -- which live only in the curated file -- did not run.
            if section == "prohibited":
                return self._note_list_state(EligibilityVerdict(
                    Verdict.PROHIBITED, Basis.STATUTE,
                    f"{hit.get('ko') or hit.get('en')}은(는) 「화장품 안전기준 등에 관한 "
                    "규정」 별표 1 사용금지 원료입니다."
                    + (" (출처: 식약처 화장품 사용제한 원료정보 API)" if hit.get("_source") == "mfds_api" else "")
                    + other, matched=hit))
            v = self._note_list_state(EligibilityVerdict(
                Verdict.RESTRICTED, Basis.STATUTE,
                f"{hit.get('ko') or hit.get('en')}은(는) 별표 2 사용제한 원료입니다."
                + (" (출처: 식약처 화장품 사용제한 원료정보 API)" if hit.get("_source") == "mfds_api" else "")
                + other,
                matched=hit, limit_text=hit.get("limit_text_ko") or "; ".join(hit.get("kr_conditions") or [])))

        # 3. Explicit regulatory permission outranks a structural alert.
        #    A 고시원료 has been examined by the regulator and permitted at a
        #    stated concentration. Retinol and retinyl palmitate are notified
        #    anti-wrinkle ingredients; the retinoid structural alert fires on
        #    them exactly as it does on tretinoin, and letting it win would have
        #    the pipeline reject an ingredient the MFDS has already approved.
        #    Retinoic acid is excluded by a different mechanism -- it is a drug
        #    substance, not a prohibited cosmetic raw material.
        if self._regs is not None and v.verdict is Verdict.UNKNOWN:
            notified = self._regs.find_notified(smiles=smiles, name=name)
            if notified:
                # Survives an unusable 별표 1 on purpose: a 고시원료 is a
                # positive permission the regulator issued at a stated
                # concentration, not an absence of evidence, and a substance
                # cannot be both notified and on the prohibition list. The
                # caveat still rides along, and the screen-level gate is still
                # marked unverified by the pipeline.
                return self._note_list_state(EligibilityVerdict(
                    Verdict.RESTRICTED, Basis.STATUTE,
                    f"{notified.get('ko')}은(는) {notified.get('category')} 고시원료로 "
                    f"{notified.get('concentration')} 범위에서 허용됩니다. 구조 경보가 "
                    "있더라도 규제기관이 명시적으로 허가한 원료입니다.",
                    matched=notified, limit_text=notified.get("concentration")))

        # 4. Drug substance, from the functional-cosmetics data.
        if self._regs is not None and v.verdict is Verdict.UNKNOWN:
            drug = self._regs.find_prohibited(name, smiles=smiles)
            if drug:
                return self._note_list_state(EligibilityVerdict(
                    Verdict.DRUG_ONLY, Basis.STATUTE,
                    f"{drug.get('ko')}은(는) 의약품 성분이라 화장품에 배합할 수 없습니다.",
                    matched=drug))

        # 5. An element whose statutory entry is narrower than the atom. Not a
        #    verdict on its own, so it lands on the read-across path with the
        #    entry quoted, and the user learns what the exception actually is.
        watch = self._watch_hit(mol)
        if v.verdict is Verdict.UNKNOWN and watch:
            sym, entry_ko, scope_ko = watch
            v = EligibilityVerdict(
                Verdict.DATA_REQUIRED, Basis.READ_ACROSS,
                f"{sym}을(를) 포함합니다. 별표 1의 해당 항목은 '{entry_ko}'이며, 원소 "
                f"전체를 금지하지 않습니다. {scope_ko} 이 물질이 그 항목에 해당하는지는 "
                "물질 단위로 확인해야 합니다.")

        # 6. Read-across from a structural alert. Never definitive.
        #    `concerns_of` because the caller may hand over either an
        #    `ADMETProfile` or the `to_dict()` form the pipeline actually
        #    stores on a candidate. The old `getattr(admet_profile,
        #    "concerns", [])` returned [] for the dict -- silently, with no
        #    error and no caveat -- so tavaborole, the benzoxaborole this
        #    module's own docstring is written about, came back
        #    eligible/none instead of data_required/read_across. CMR feature
        #    warnings remain data requirements even when their specificity
        #    review removes a blanket critical flag. Displaying an exploratory
        #    candidate is independent of clearing it for cosmetic use.
        review = [c for c in concerns_of(admet_profile)
                  if c.severity == "critical" or
                  (c.domain == "cmr" and c.severity in {"high", "moderate", "context"})]
        if review:
            c = review[0]
            alert = (f"{c.name} 구조 경보이며 법적 분류나 독성 측정 결과가 아닙니다. "
                     f"물질별 자료와 적용 조건을 확인해야 합니다. {c.mechanism[:140]}")
            if v.verdict is Verdict.UNKNOWN:
                v = EligibilityVerdict(Verdict.DATA_REQUIRED, Basis.READ_ACROSS, alert)
            elif v.basis is Basis.READ_ACROSS:
                v.caveats.append(alert)

        if v.verdict is Verdict.UNKNOWN:
            v = EligibilityVerdict(Verdict.ELIGIBLE, Basis.NONE,
                                   "사용금지 목록과 금지 원소에 해당하지 않습니다.")

        if not self.list_usable:
            # The whole point of the module, and for a long time it was only a
            # caveat: "not on the list" when there is no list is not a finding
            # about the substance, it is a finding about the pipeline. ELIGIBLE
            # here rests entirely on a lookup that never happened, so it is
            # demoted to the state that means "no judgement was possible".
            # DATA_REQUIRED is not usable, so `ScreenResult.orderable` holds it
            # back from the buy list instead of shipping it as cleared.
            if v.verdict is Verdict.ELIGIBLE:
                v = EligibilityVerdict(
                    Verdict.DATA_REQUIRED, Basis.NOT_ASSESSED,
                    "대조할 사용금지 원료 목록이 없어 목록 대조 자체가 이뤄지지 않았습니다. "
                    "원소 검사와 구조 경보만 통과했을 뿐이므로 '사용 가능'이 아니라 "
                    "'미검증'이며, 목록을 복구해 재평가해야 합니다.",
                    caveats=list(v.caveats))
            self._note_list_state(v)
        else:
            if hit is None:
                v.caveats.append(self._coverage_caveat())
            if self._d.get("completeness"):
                v.caveats.append(f"목록 완전성: {self._d['completeness']}")
        if self._d.get("retrieved"):
            v.caveats.append(f"목록 수집일 {self._d['retrieved']}. 고시는 개정됩니다.")
        rad = sorted({a.GetSymbol() for a in mol.GetAtoms()} & RADIOACTIVE_ELEMENTS)
        if rad:
            v.caveats.append(
                f"{'/'.join(rad)}: 별표 1에 원소명으로 등재된 사용금지 원료가 아닙니다. "
                "다만 조건부 항목 '방사성물질(다만, 제품에 포함된 방사능의 농도 등이 「생활주변"
                "방사선 안전관리법」 제15조의 규정에 적합한 경우 제외)'의 적용 여부를 확인해야 "
                "합니다. 이는 조건부 규정이지 원소 금지가 아닙니다.")
        # Last, and never touching `verdict` or `basis`. A hit here is the fact
        # that the regulator already lists this substance as a cosmetic raw
        # material, which is what a discovery screen wants to know and is not a
        # permission -- the prohibition verdict above stands whatever this says,
        # and a substance can be both an inventory row and a barred one. Only a
        # hit is reported: with no inventory file there is nothing to say, and
        # "not in the inventory" would read as a finding when it is silence.
        return v

    def _inventory_hit(self, name: Optional[str],
                       cas: Optional[Union[str, Sequence[str]]]) -> Optional[dict]:
        """The inventory row this substance matches, by CAS then by exact name."""
        for one in ([cas] if isinstance(cas, str) or cas is None else list(cas)):
            row = self._inv_by_cas.get(_cas_key(one))
            if row:
                return row
        for variant in name_variants(name, split_groups=False):
            row = self._inv_by_name.get(variant)
            if row:
                return row
        return None

    def _coverage_caveat(self) -> str:
        """Why "not on the list" is a statement about the lookup, not the risk.

        The module docstring promises that a missing list never reads as a
        clean bill of health. A list that is *present but unmatchable in kind*
        needs the same promise kept: the 별표 1 export carries a structure for
        only 15 of its 626 rows, and those 15 collapse to the 12 distinct
        InChIKeys this method reports below, so a compound identified by SMILES
        or an English name is being checked against a mostly Korean-name-only
        index and can miss a real hit. The number in the caveat is measured
        from the loaded file, which is why it cannot drift away from this
        paragraph unnoticed.
        """
        c = self.coverage
        return (f"조회 범위 한계: 별표 1·2 {c['rows']}건 중 구조식 보유 {c['structural']}건, "
                f"CAS 보유 {c['cas']}건, 영문명 보유 {c['en']}건(이름이 변환 과정에서 손상된 "
                f"행 {c['name_corrupted']}건은 이름 조회에서 제외). 이 판정이 실제로 대조한 "
                "키는 위에 적힌 것뿐입니다. 쓰이지 않은 키가 있다면 '목록에 없음'은 안전하다는 "
                "뜻이 아니라 그 키로 대조한 적이 없다는 뜻입니다.")


def summarise(verdict: EligibilityVerdict) -> str:
    tag = "확정" if verdict.definitive else "판단"
    return f"[{verdict.verdict.value} / {tag}] {verdict.reason}"
