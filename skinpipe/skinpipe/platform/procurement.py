"""Recorded purchase availability; supplier checkout and stock are separate."""
from datetime import date, datetime
from urllib.parse import urlparse


def _dated_offer(offer):
    if not isinstance(offer, dict) or offer.get("available") is not True:
        return False
    source, checked = offer.get("source"), offer.get("checked_at")
    if not isinstance(source, str) or not source.strip() or not isinstance(checked, str):
        return False
    try:
        # ISO dates or timestamps, including an explicit UTC Z suffix.
        if len(checked) == 10:
            date.fromisoformat(checked)
        else:
            datetime.fromisoformat(checked.replace("Z", "+00:00"))
    except ValueError:
        return False
    return True


def has_purchase_record(row):
    """Accept affirmative flags or dated, attributed available supplier offers.

    Missing flags and truthy strings are not availability evidence. An imported
    affirmative flag remains a catalog annotation, not a live-stock guarantee.
    """
    from .scout import offer_evidence_status
    if offer_evidence_status(row.get("catalog_evidence") or row) != "complete":
        return False
    offers = row.get("offers")
    offers = offers if isinstance(offers, list) else []
    return (row.get("purchasable") is True or row.get("purchase_available") is True
            or any(_dated_offer(offer) for offer in offers))


def documented_offers(row, *, today=None):
    """Dated exact-form evidence for a final purchase list; stock is separate.

    Catalog flags remain useful for retrieval but cannot prove a specific
    product matches stereochemistry, isotope, charge and salt form.
    """
    from .identity import chemical_identity
    identity = chemical_identity(row.get("smiles"))
    exact = identity.get("exact_structure_id")
    if exact is None:
        return []
    result = []
    for offer in row.get("offers") or []:
        if not isinstance(offer, dict):
            continue
        supplier = offer.get("supplier") or offer.get("vendor")
        link = offer.get("product_url") or offer.get("url")
        checked = offer.get("checked_at")
        if not all(isinstance(v, str) and v.strip() for v in (supplier, link, checked)):
            continue
        parsed = urlparse(link)
        if parsed.scheme not in {"http", "https"} or not parsed.hostname or parsed.username or parsed.password:
            continue
        try:
            checked_date = date.fromisoformat(checked) if len(checked) == 10 else datetime.fromisoformat(checked.replace("Z", "+00:00")).date()
        except ValueError:
            continue
        current = today or date.today()
        if checked_date > current:
            continue
        offered = chemical_identity(offer.get("smiles") or offer.get("exact_smiles"))
        if offered.get("exact_structure_id") != exact:
            continue
        # Explicit unavailability is retained in the original offer, but cannot
        # supply a currently orderable purchase shortlist.
        if offer.get("available") is False:
            continue
        result.append({**offer, "supplier": supplier, "product_url": link,
                       "exact_structure_id": exact, "age_days": (current-checked_date).days,
                       "age_warning": (current-checked_date).days > 180,
                       "stock_status": offer.get("stock_status", "unconfirmed"),
                       "lead_time_status": offer.get("lead_time_status", "unconfirmed"),
                       "order_status": offer.get("order_status", "unconfirmed")})
    return result


def procurement_state(row):
    from .scout import offer_evidence_status
    supplier_status = offer_evidence_status(row.get("catalog_evidence") or row)
    if supplier_status != "complete":
        return {"catalog_record": None, "documented_exact_form": False,
                "status": supplier_status, "offers": [], "stock_confirmed": False,
                "reason": "supplier_evidence_incomplete"}
    offers = documented_offers(row)
    return {"catalog_record": has_purchase_record(row), "documented_exact_form": bool(offers),
            "status": "documented" if offers else "procurement_unconfirmed", "offers": offers,
            "stock_confirmed": any(o.get("stock_status") == "confirmed" for o in offers)}
