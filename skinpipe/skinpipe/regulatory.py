"""Korean functional-cosmetic (기능성화장품) context for a discovery hit.

A binder is not a product. In Korea a cosmetic that claims a function falls
under 화장품법, and the route to market forks sharply:

* **보고** -- the formulation uses a notified ingredient (고시원료) at the
  notified concentration, in a permitted dosage form. No efficacy review.
* **심사** -- anything else, including every novel ingredient. A full dossier
  with human efficacy data.

The gap between those two is years and a great deal of money, so it is the
single most consequential thing a discovery programme can get wrong. A pipeline
that proposes a novel chemical entity without saying what that costs is giving
half an answer.

Three traps this module exists to catch.

**The wrong category.** Adenosine is a notified anti-wrinkle ingredient. It is
not a whitening ingredient, and proposing it for pigmentation forfeits the 보고
route entirely.

**A claim that is actually a drug claim.** Hair-loss work drifts here
constantly: the moment the mechanism is described as blocking DHT or androgen,
the product reads as 의약품 regardless of what is in it. Minoxidil and
finasteride are drugs and cannot be formulated into a cosmetic at all.

**A category that does not exist.** Several markers in this project's table --
sensitivity and damage repair among them -- have no statutory functional
category. Work on them can still make an excellent product; it simply cannot
carry a 기능성 claim, and that has to be known before the programme starts.

Data lives in ``config/kr_functional_cosmetics.yaml``, which carries a
``source_url`` per item and a ``verified`` flag where the research could not
confirm a value. Unverified items are surfaced, never silently used.
"""
from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Optional

import yaml
from rdkit import Chem, RDLogger
from rdkit.Chem import inchi

RDLogger.DisableLog("rdApp.*")
log = logging.getLogger(__name__)

_DATA = Path(__file__).resolve().parent.parent / "config" / "kr_functional_cosmetics.yaml"


# ----------------------------------------------------------- name matching
#
# Why this is exact matching over expanded variants, and not substring matching
# -----------------------------------------------------------------------------
# The previous implementation accepted a hit when either string contained the
# other. Against a statutory list that legitimately contains 헥산 (n-hexane),
# 페놀, 벤젠, 니켈, 붕산 and 브롬, that is not an occasional near miss; it is a
# machine for manufacturing false *definitive* verdicts. 1,2-헥산디올, a legal
# humectant, read as n-hexane. 폴리페놀 read as phenol. 소듐도데실벤젠설포네이트
# read as benzene. And it failed in the permissive direction too, which is
# worse: 알파-알부틴 matched the 별표 4 notified entry 알부틴, handing back a
# regulatory permission for an ingredient that has its own specification.
# kr_functional_cosmetics_NOTES.md section 5 item 6 records that exact
# conflation as a known cause of a 부적합 판정.
#
# The replacement expands each side into a set of *complete* alternative
# spellings and requires exact equality between them. Every generated variant
# is a whole name, never a fragment, so a hit always means "these two strings
# name the same substance" and a near miss is simply not a hit. Korean chemical
# names carry no spaces, so word-boundary logic would buy nothing here; exact
# equality after normalisation is the only safe default. There is deliberately
# no fuzzy fallback: this function feeds definitive verdicts.

_PAREN = re.compile(r"\(([^()]*)\)")
_INDEX = re.compile(r"\[\s*\d+\s*\]")          # the [1] / [2] group markers
_SPLIT = re.compile(r"[;·]|\sand\s", re.IGNORECASE)
_DROP = re.compile(r"[\s()\[\]{}_.,;:/\\'\"`~!?*+·‐-―‘’“”-]+")

#: Korean transliteration pairs that denote the same syllable. Applied to both
#: sides, so folding can only ever merge two spellings of one word, never two
#: different words. Deliberately minimal: a pair is added only when a primary
#: source shows both spellings used for the same substance. 식약처 고시 text
#: writes 히드로퀴논 and -에텔; trade and INCI-Korean usage writes 하이드로퀴논
#: and -에테르, and a name index that knows only one of them silently misses.
KO_SPELLING_FOLDS: tuple[tuple[str, str], ...] = (("하이드로", "히드로"),
                                                  ("에테르", "에텔"))

#: Collective suffixes. 별표 1 writes "X 및 그 염류"; a formulator writes "X".
#: Stripping these yields another complete name, so it stays exact matching.
_GROUP_SUFFIXES = ("및그염류및유도체", "그염류및유도체", "및그유도체", "및그염류",
                   "및그화합물", "및그염", "anditssalts", "andtheirsalts",
                   "anditssalt", "andits")

#: A parenthetical is either a synonym worth indexing -- 카테콜(피로카테콜) -- or
#: a clause that is not a name at all -- (2023년 추가), (셀레늄아스파테이트는 제외).
#: Indexing a clause as if it were a name is how a condition becomes a verdict.
_CLAUSE_MARKERS = ("제외", "추가", "시행", "경우", "이하", "이상", "까지", "포함",
                   "한함", "기준", "다만", "로서", "제품", "%", "ppm")


def normalise_name(value: Optional[str]) -> str:
    """Casefold, drop separators and punctuation, fold spelling variants."""
    if not value:
        return ""
    text = _INDEX.sub("", str(value)).strip().casefold()
    text = _DROP.sub("", text)
    for variant, canonical in KO_SPELLING_FOLDS:
        text = text.replace(variant, canonical)
    return text


def name_variants(value: Optional[str], *, split_groups: bool = True) -> set[str]:
    """Every complete name that ``value`` could reasonably be written as.

    ``split_groups`` splits grouped list entries -- ``과붕산나트륨·과붕산나트륨
    일수화물``, ``SODIUM PERBORATE and SODIUM PERBORATE MONOHYDRATE`` -- into
    their members. That is a convention of the statutory text, so it is applied
    when indexing the list and *not* to a query: taking one component out of a
    blend name and declaring the blend prohibited is the false-positive this
    whole function exists to prevent.
    """
    if not value:
        return set()
    raw = _INDEX.sub("", str(value)).strip()
    pieces = [p for p in _SPLIT.split(raw) if p.strip()] if split_groups else [raw]
    if len(pieces) > 1:
        pieces.append(raw)

    out: set[str] = set()
    for piece in pieces:
        out.add(normalise_name(piece))
        out.add(normalise_name(_PAREN.sub(" ", piece)))
        for inner in _PAREN.findall(piece):
            if len(inner.strip()) > 1 and not any(m in inner for m in _CLAUSE_MARKERS):
                out.add(normalise_name(inner))
                if split_groups:
                    # List synonyms inside parentheses, without splitting
                    # chemical locants such as the bare name 1,2-hexanediol.
                    out.update(normalise_name(alias) for alias in re.split(r",\s*", inner)
                               if alias.strip() and not re.match(r"^\d+[,\-]", inner))
    for name in list(out):
        for suffix in _GROUP_SUFFIXES:
            if name.endswith(suffix) and len(name) > len(suffix) + 1:
                out.add(name[: -len(suffix)])
    out.discard("")
    return out


def index_names(items: Iterable[dict], fields: tuple[str, ...]) -> dict[str, dict]:
    """Build variant -> item. Rows flagged ``name_corrupted`` are left out.

    47 rows of the 별표 1 export carry a truncated HWP->PDF fragment instead of
    a name ('단)', '콜)', '로서0.1%이하인경우는제외)'). Their CAS is intact and
    still keys them, but a fragment must never be able to produce a verdict by
    name, so it never enters the name index.
    """
    index: dict[str, dict] = {}
    for item in items:
        if item.get("name_corrupted"):
            continue
        for field_name in fields + ("aliases",):
            value = item.get(field_name)
            for one in (value if isinstance(value, list) else [value]):
                for variant in name_variants(one):
                    index.setdefault(variant, item)
    return index

#: Pipeline marker category -> statutory 기능성화장품 category.
#: ``None`` means there is no functional category for this kind of claim.
CATEGORY_MAP: dict[str, Optional[str]] = {
    "pigmentation": "whitening",
    "photoaging": "anti_wrinkle",
    "hair": "hair_loss",
    "acne_rosacea": "acne",
    "barrier": "skin_barrier",
    "sensitive": None,
    "repair": None,
}


@dataclass
class RegulatoryAssessment:
    marker: str
    pipeline_category: str
    kr_category: Optional[str] = None
    kr_category_ko: Optional[str] = None
    has_functional_category: bool = False
    notified_ingredients: list[dict] = field(default_factory=list)
    is_notified: bool = False
    notified_match: Optional[dict] = None
    prohibited: Optional[dict] = None
    route: str = "심사"
    dosage_form_restriction: Optional[str] = None
    claim_traps: list[dict] = field(default_factory=list)
    caveats: list[str] = field(default_factory=list)
    verdict: str = ""
    # A miss in this finite index is not a finding that no restriction exists.
    prohibited_lookup_status: str = "unassessed"
    prohibited_match_method: Optional[str] = None


class Regulations:
    def __init__(self, path: Optional[Path] = None):
        self.path = Path(path or _DATA)
        self._d: dict[str, Any] = yaml.safe_load(self.path.read_text())
        self.retrieved = self._d.get("retrieved")
        self._by_key = {c["key"]: c for c in self._d.get("categories", [])}
        self._no_notified = {
            c.get("key") if isinstance(c, dict) else c
            for c in self._d.get("categories_without_notified_ingredients", [])
        }
        self._inchikeys: dict[str, dict] = {}
        for n in self._d.get("notified_ingredients", []):
            smi = n.get("smiles")
            if not smi:
                continue
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                log.warning("notified ingredient %s has an unparseable SMILES", n.get("ko"))
                continue
            self._inchikeys[inchi.MolToInchiKey(mol)] = n
        self._notified_names = index_names(self._d.get("notified_ingredients", []),
                                           ("inci", "ko"))
        prohibited = self._d.get("prohibited_as_cosmetic", [])
        self._drug_names = index_names(prohibited, ("en", "ko"))
        self._drug_inchikeys: dict[str, dict] = {}
        for entry in prohibited:
            key = self._structure_key(entry.get("smiles"))
            if key:
                if entry.get("inchikey") and entry["inchikey"] != key:
                    raise ValueError("prohibited structure and recorded InChIKey disagree")
                self._drug_inchikeys[key] = entry

    # ------------------------------------------------------------- lookups
    def category_for(self, pipeline_category: str) -> Optional[str]:
        return CATEGORY_MAP.get(pipeline_category)

    def notified_for(self, kr_category: str) -> list[dict]:
        return [n for n in self._d.get("notified_ingredients", [])
                if n.get("category") == kr_category]

    def find_notified(self, smiles: Optional[str] = None,
                      name: Optional[str] = None) -> Optional[dict]:
        """Match a compound against the notified list by structure, then by name.

        Structure first: a name can be spelled a dozen ways and an InChIKey
        cannot. Only seven notified entries carry a SMILES, since most are
        extracts or salts, so the name route still does real work.

        The name route is exact over expanded variants. A notified match is a
        regulatory *permission* that outranks structural alerts downstream, so
        an approximate hit here is the most expensive kind of error in the
        module: 알파-알부틴 must not come back as the 별표 4 entry 알부틴.
        """
        structure_key = self._structure_key(smiles)
        if structure_key:
            hit = self._inchikeys.get(structure_key)
            if hit:
                return hit
        hit = self._by_name(self._notified_names, name)
        if hit and structure_key and hit.get("smiles"):
            reference = Chem.MolFromSmiles(hit["smiles"])
            if reference is None or inchi.MolToInchiKey(reference) != structure_key:
                log.warning("notified ingredient name %r conflicts with supplied structure", name)
                return None
        return hit

    @staticmethod
    def _structure_key(smiles: Optional[str]) -> Optional[str]:
        if not isinstance(smiles, str) or not smiles.strip():
            return None
        mol = Chem.MolFromSmiles(smiles)
        if mol is None or mol.GetNumHeavyAtoms() == 0:
            return None
        return inchi.MolToInchiKey(mol) or None

    def _prohibited_lookup(self, smiles: Optional[str], name: Optional[str]):
        key = self._structure_key(smiles)
        hit = self._drug_inchikeys.get(key) if key else None
        if hit:
            return hit, "matched", "exact_structure"
        by_name = self._by_name(self._drug_names, name)
        if by_name:
            # Do not prohibit a different molecule merely because its label
            # conflicts with an exact supplied structure. Hold the identity.
            reference = self._structure_key(by_name.get("smiles"))
            if key and reference and key != reference:
                return None, "identity_conflict", None
            if smiles is not None and not key:
                return None, "unassessed", None
            return by_name, "matched", "exact_name"
        if smiles is not None and not key:
            return None, "unassessed", None
        if not key and not name_variants(name, split_groups=False):
            return None, "unassessed", None
        return None, "not_matched", None

    def find_prohibited(self, name: Optional[str] = None, *,
                        smiles: Optional[str] = None) -> Optional[dict]:
        """Exact structure, then exact name. None is a lookup miss, not permission.

        The first positional argument remains the name for existing callers.
        Structures retain full stereochemistry; no family/salt/fuzzy expansion
        is implied by a hit on this small curated index.
        """
        return self._prohibited_lookup(smiles, name)[0]

    @staticmethod
    def _by_name(index: dict[str, dict], name: Optional[str]) -> Optional[dict]:
        for variant in name_variants(name, split_groups=False):
            hit = index.get(variant)
            if hit:
                return hit
        return None

    def traps_for(self, kr_category: Optional[str]) -> list[dict]:
        return [t for t in self._d.get("claim_traps", [])
                if kr_category is None or t.get("category") in (kr_category, None)]

    def route_detail(self, key: str) -> Optional[dict]:
        return next((r for r in self._d.get("routes", []) if r.get("key") == key), None)

    # ------------------------------------------------------------- assess
    def assess(self, marker: str, pipeline_category: str,
               smiles: Optional[str] = None,
               compound_name: Optional[str] = None) -> RegulatoryAssessment:
        kr = self.category_for(pipeline_category)
        cat = self._by_key.get(kr) if kr else None
        a = RegulatoryAssessment(marker=marker, pipeline_category=pipeline_category,
                                 kr_category=kr,
                                 kr_category_ko=cat.get("ko") if cat else None,
                                 has_functional_category=kr is not None)

        # Whether an ingredient may be used at all is a property of the
        # ingredient. Whether a functional claim can be made about it is a
        # property of the category. This lookup used to sit *inside* the
        # branch that runs only when a category exists, so the same molecule
        # got opposite answers depending on what was being claimed:
        #
        #     Minoxidil / hair          -> 불가        (prohibited found)
        #     Minoxidil / pigmentation  -> 불가        (prohibited found)
        #     Minoxidil / sensitive     -> 해당없음    (never looked)
        #     Minoxidil / repair        -> 해당없음    (never looked)
        #
        # and the 해당없음 verdict ends "design it as a general cosmetic",
        # which reads as permission for a 의약품 성분 that may not go into a
        # cosmetic at all. Dropping the functional claim is what the reader
        # would do on being told the category does not exist.
        prohibited, lookup_status, match_method = self._prohibited_lookup(smiles, compound_name)
        a.prohibited_lookup_status = lookup_status
        a.prohibited_match_method = match_method

        if kr is not None:
            a.dosage_form_restriction = (cat.get("product_restriction_ko")
                                         if cat else None)
            a.notified_ingredients = self.notified_for(kr)
            a.claim_traps = self.traps_for(kr)
        else:
            a.claim_traps = self.traps_for(None)

        if prohibited:
            a.prohibited = prohibited
            a.route = "불가"
            reason = prohibited.get("reason_ko") or "화장품 사용금지 원료 목록의 해당 근거를 확인하세요."
            a.verdict = (f"{prohibited.get('ko')}은(는) 화장품 배합 불가로 분류됩니다. "
                         f"근거: {reason} 기전 벤치마크로만 사용하세요.")
            if kr is None:
                # Both facts, and the order matters: the restriction binds
                # whatever the claim, so it is not softened by the missing
                # category -- it is the reason the missing category is moot.
                a.verdict += (f" ({pipeline_category}에는 대응하는 기능성화장품 "
                              "카테고리도 없으나, 일반 화장품으로 돌려도 배합은 "
                              "불가합니다.)")
            return a

        match = self.find_notified(smiles=smiles, name=compound_name)
        # Positive notified evidence may establish a category-specific route.
        # Otherwise an unmatched name/structure establishes no permission.
        if lookup_status in ("unassessed", "identity_conflict") or match is None:
            a.route = "미확인"
            a.verdict = ("성분의 이름·정확 구조 조회에서 규제 지위를 확인하지 못했습니다. "
                         "조회 실패·불일치는 제한 없음이나 화장품 배합 허가를 뜻하지 않습니다. "
                         "물질 식별과 현행 사용금지·사용제한 원문 확인이 필요합니다.")
            a.caveats.append("이 조회는 소수의 고시·의약품 항목만 포함하며, 전체 금지 목록이나 "
                             "염·유도체의 포괄 판정을 대신하지 않습니다.")
            if lookup_status == "identity_conflict":
                a.caveats.append("입력 이름과 구조가 다른 물질을 가리킵니다.")
            return a

        if kr is None:
            a.route = "해당없음"
            a.verdict = (f"{pipeline_category}에는 대응하는 기능성화장품 카테고리가 없습니다. "
                         "좋은 제품은 만들 수 있으나 기능성 표방은 불가하므로, 일반 화장품 "
                         "소구로 설계하거나 다른 카테고리로 프레이밍해야 합니다.")
            return a

        if match and match.get("category") == kr:
            a.is_notified = True
            a.notified_match = match
            a.route = "보고"
            a.verdict = (f"{match.get('ko')}은(는) 이 카테고리의 고시원료입니다 "
                         f"({match.get('concentration')}). 고시 함량과 제형을 지키면 "
                         "심사 없이 보고로 갑니다.")
        elif match:
            a.notified_match = match
            a.route = "심사"
            a.verdict = (f"{match.get('ko')}은(는) 고시원료이지만 카테고리가 다릅니다"
                         f"({match.get('category')}). 이 표방으로 쓰려면 심사가 필요합니다.")
        elif kr in self._no_notified:
            a.route = "심사"
            a.verdict = (f"{cat.get('ko') if cat else kr} 카테고리에는 고시원료가 없습니다. "
                         "모든 제품이 심사 대상입니다.")
        else:
            a.route = "심사"
            n_notified = len(a.notified_ingredients)
            a.verdict = (f"고시원료가 아니므로 심사 대상입니다. 이 카테고리에는 고시원료가 "
                         f"{n_notified}종 있으니, 속도가 중요하면 그중 하나를 쓰는 처방과 "
                         "비교해 판단하세요.")

        if match and match.get("verified") is False:
            a.caveats.append(
                f"{match.get('ko')}의 고시 내용은 원문 확인이 되지 않았습니다"
                + (f": {match.get('note_ko')}" if match.get("note_ko") else "")
                + ". 의약품안전나라에서 현행 고시를 확인하세요.")
        if self.retrieved:
            a.caveats.append(
                f"고시 데이터 수집일 {self.retrieved}. 고시는 개정되므로 실제 신청 전에는 "
                "law.go.kr에서 현행본을 재확인하세요.")
        return a


def novel_ingredient_cost_note() -> str:
    """What choosing a novel ingredient actually commits you to."""
    return ("신규 원료로 기능성화장품을 상업화하려면 인체적용시험을 포함한 전체 안전성·"
            "유효성 자료를 갖춘 심사를 거쳐야 하고, 현실적으로 수 년과 상당한 비용이 "
            "듭니다. 고시원료를 쓰면 보고로 끝납니다. 신규성과 특허를 얻는 대신 무엇을 "
            "내주는지 계산하고 결정하세요.")
