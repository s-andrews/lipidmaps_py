import re
from typing import Any, Dict, List, Optional, Set

from ..models.reaction import ReactionData, CompoundComponent
from .lipid_reaction_rules import lipid_reaction_rules


class ReactionEvaluator:
    """Evaluate ReactionData objects against rule metadata and optional dataset evidence.

    Minimal, self-contained evaluator that annotates ReactionData with
    `possible` (bool) and `possible_explanation` (str).
    """

    def __init__(self, rules: Optional[Dict[str, Any]] = None):
        self.rules = rules or lipid_reaction_rules

    def annotate_reactions(self, reactions: List[ReactionData], dataset: Optional[Any] = None) -> None:
        for r in reactions:
            try:
                res = self.evaluate_reaction(r, dataset=dataset)
            except Exception as e:
                # ensure reactions are always annotated even if evaluator errors
                setattr(r, "possible", False)
                setattr(r, "possible_explanation", f"error_evaluating:{e}")
                setattr(r, "evaluation", {"possible": False, "explanation": f"error_evaluating:{e}"})
            else:
                setattr(r, "possible", res.get("possible", False))
                setattr(r, "possible_explanation", res.get("explanation", ""))
                # also attach the full evaluation dict for UI selection / detailed view
                try:
                    setattr(r, "evaluation", res)
                except Exception:
                    # best-effort: evaluation may be non-serializable; fall back to basic fields
                    setattr(r, "evaluation", {"possible": res.get("possible", False), "explanation": res.get("explanation", "")})

    def evaluate_reaction(self, reaction: ReactionData, dataset: Optional[Any] = None) -> Dict[str, Any]:
        reasons: List[str] = []
        possible = True

        dataset_lm_ids: Set[str] = set()
        lipids = []
        if dataset:
            lipids = getattr(dataset, "lipids", [])
            for lip in lipids:
                if getattr(lip, "lm_id", None):
                    dataset_lm_ids.add(lip.lm_id)
                if getattr(lip, "generic_lm_id", None):
                    dataset_lm_ids.add(lip.generic_lm_id)

        def _get_val(comp, field: str):
            if comp is None:
                return None
            if isinstance(comp, dict):
                return comp.get(field)
            return getattr(comp, field, None)

        def _get_headgroup_from_comp(comp) -> Optional[str]:
            if not comp:
                return None
            # Prefer explicit compound headgroup field
            hg = _get_val(comp, "compound_headgroup")
            if hg:
                return hg
            # Accept QuantifiedLipid-style fields (headgroup, generic_lm_id, standardized_name, main_class)
            q_hg = _get_val(comp, "headgroup") or _get_val(comp, "generic_lm_id") or _get_val(comp, "standardized_name") or _get_val(comp, "main_class")
            if q_hg:
                return q_hg

            # Generic LM IDs may encode the class (legacy compound field)
            gl = _get_val(comp, "compound_generic_lm_id")
            if gl:
                return gl

            # Try abbrev / full_struct / name fields and extract leading headgroup token.
            for fld in ("compound_abbrev", "compound_full_struct", "compound_name"):
                txt = _get_val(comp, fld)
                if not txt:
                    continue
                t = str(txt).strip()
                # match headgroup optionally followed by space + P- or O- (e.g., 'LPC P-')
                m = re.match(r"^([A-Za-z0-9]+(?:\s(?:P|O)-)?)", t)
                if m:
                    return m.group(1)
            return None

        def _get_linkage_from_comp(comp) -> Optional[str]:
            if not comp:
                return None
            txt = _get_val(comp, "compound_abbrev") or _get_val(comp, "compound_headgroup") or _get_val(comp, "compound_name") or ""
            if not txt:
                return None
            t = str(txt).strip()
            if t.startswith("P-") or " P-" in t:
                return "ether_vinyl"
            if t.startswith("O-") or " O-" in t:
                return "ether_alkyl"
            if t.startswith("Cer") or t.startswith("SM") or str(t).lower().startswith("dhcer"):
                return "amide"
            return "ester"

        def _chain_count_from_comp(comp) -> Optional[int]:
            if not comp:
                return None
            cab = _get_val(comp, "compound_abbrev_chains")
            if cab:
                try:
                    return int(cab)
                except Exception:
                    pass
            text = _get_val(comp, "compound_abbrev") or _get_val(comp, "compound_name") or _get_val(comp, "compound_full_struct") or ""
            if text and "/" in str(text):
                return str(text).count("/") + 1
            if text and re.search(r"\d{1,2}:\d", str(text)):
                return len(re.findall(r"\d{1,2}:\d", str(text)))
            # fallback to dataset lipids
            for lid in dataset_lm_ids:
                for l in lipids:
                    if getattr(l, "lm_id", None) == lid or getattr(l, "generic_lm_id", None) == lid:
                        s = getattr(l, "standardized_name", None) or ""
                        if s and "/" in s:
                            return s.count("/") + 1
                        if s and re.search(r"\d{1,2}:\d", s):
                            return len(re.findall(r"\d{1,2}:\d", s))
            return None

        def _chain_carbon_list_from_comp(comp) -> List[int]:
            """Return list of carbon counts found in component (e.g., [16,18] for 16:0/18:1)."""
            if not comp:
                return []
            text = " ".join([str(_get_val(comp, f) or "") for f in ("compound_abbrev", "compound_full_struct", "compound_name")])
            if not text:
                return []
            found = re.findall(r"(\d{1,2}):\d", text)
            return [int(x) for x in found]

        def _total_carbons(comps: List[CompoundComponent]) -> Optional[int]:
            totals = []
            for c in comps:
                lst = _chain_carbon_list_from_comp(c)
                if lst:
                    totals.extend(lst)
            if not totals:
                return None
            return sum(totals)

        def _chain_db_list_from_comp(comp) -> List[int]:
            """Return list of double-bond counts found in component (e.g., [0,1] for 16:0/18:1)."""
            if not comp:
                return []
            text = " ".join([str(_get_val(comp, f) or "") for f in ("compound_abbrev", "compound_full_struct", "compound_name")])
            if not text:
                return []
            found = re.findall(r"\d{1,2}:(\d)", text)
            return [int(x) for x in found]

        def _total_dbs(comps: List[CompoundComponent]) -> Optional[int]:
            totals = []
            for c in comps:
                lst = _chain_db_list_from_comp(c)
                if lst:
                    totals.extend(lst)
            if not totals:
                return None
            return sum(totals)

        def _has_sphingoid_evidence(comp) -> bool:
            txt = " ".join([str(_get_val(comp, f) or "") for f in ("compound_abbrev", "compound_full_struct", "compound_name")])
            if re.search(r"\b[dts]\d{1,2}:\d\b", txt):
                return True
            for lid in dataset_lm_ids:
                for l in lipids:
                    if getattr(l, "lm_id", None) == lid or getattr(l, "generic_lm_id", None) == lid:
                        s = getattr(l, "standardized_name", None) or ""
                        if re.search(r"\b[dts]\d{1,2}:\d\b", s):
                            return True
            return False

        pairs_info: List[Dict[str, Any]] = []

        for react in reaction.reactants:
            for prod in reaction.products:
                hr = _get_headgroup_from_comp(react)
                hp = _get_headgroup_from_comp(prod)
                if not hr or not hp:
                    continue
                rules = self.rules.get("headgroup_reactions", {}).get(hr, {})
                conv = rules.get("conversion_rules", {}) or {}
                rule = conv.get(hp)
                if not rule:
                    continue

                # linkage
                if rule.get("require_same_linkage", True):
                    lr = _get_linkage_from_comp(react) or rules.get("linkage_type")
                    lp = _get_linkage_from_comp(prod) or self.rules.get("headgroup_reactions", {}).get(hp, {}).get("linkage_type")
                    if lr and lp and lr != lp:
                        possible = False
                        reasons.append(f"linkage_mismatch:{hr}({lr})->{hp}({lp})")

                # acyl chains
                req_ch = rule.get("required_acyl_chains")
                if req_ch is not None:
                    rc = _chain_count_from_comp(react)
                    pc = _chain_count_from_comp(prod)
                    if rc is None and pc is None:
                        possible = False
                        reasons.append(f"missing_acyl_info:{hr}->{hp} needs {req_ch}")
                    else:
                        if not ((rc is not None and rc == req_ch) or (pc is not None and pc == req_ch)):
                            possible = False
                            reasons.append(f"acyl_count_mismatch:{hr}->{hp} need {req_ch} reactant={rc} product={pc}")

                # sphingoid
                src_rule = self.rules.get("headgroup_reactions", {}).get(hr, {})
                tgt_rule = self.rules.get("headgroup_reactions", {}).get(hp, {})
                if src_rule.get("has_sphingoid") or tgt_rule.get("has_sphingoid"):
                    has_src = _has_sphingoid_evidence(react)
                    has_tgt = _has_sphingoid_evidence(prod)
                    if not (has_src or has_tgt):
                        possible = False
                        reasons.append(f"missing_sphingoid:{hr}->{hp}")

                # total-carbon conservation check (apply when a conversion rule defines required_acyl_chains)
                req_ch = rule.get("required_acyl_chains")
                if req_ch is not None:
                    # compute total carbons across all reactants/products in this pair
                    reactants_all = [react]  # pairwise; in future consider multi-component sums
                    products_all = [prod]
                    r_total = _total_carbons(reactants_all)
                    p_total = _total_carbons(products_all)
                    # include any free fatty acid components present in reaction's reactants/products
                    # scan full reaction for FA components
                    fa_total = 0
                    for c in (reaction.reactants + reaction.products):
                        hg = _get_headgroup_from_comp(c) or ""
                        if hg == "FA" or (_get_val(c, "compound_type") and _get_val(c, "compound_type") == "fa"):
                            lst = _chain_carbon_list_from_comp(c)
                            if lst:
                                fa_total += sum(lst)
                    # adjust p_total by adding any FA components found among products (they represent released tails)
                    if p_total is None and fa_total:
                        # if product lacks main-chain info but FA present, use FA total as product contribution
                        p_total = fa_total

                    if r_total is None or p_total is None:
                        possible = False
                        reasons.append(f"missing_chain_info:{hr}->{hp} reactant_total={r_total} product_total={p_total} fa_total={fa_total}")
                    else:
                        if r_total != p_total:
                            possible = False
                            reasons.append(f"carbon_mismatch:{hr}->{hp} reactant_total={r_total} product_total={p_total} delta={r_total - p_total}")

                # gather detailed per-pair info for UI/display
                try:
                    react_carbons = _chain_carbon_list_from_comp(react)
                    prod_carbons = _chain_carbon_list_from_comp(prod)
                    react_dbs = _chain_db_list_from_comp(react)
                    prod_dbs = _chain_db_list_from_comp(prod)
                    react_chain_count = _chain_count_from_comp(react)
                    prod_chain_count = _chain_count_from_comp(prod)
                    lr = _get_linkage_from_comp(react) or rules.get("linkage_type")
                    lp = _get_linkage_from_comp(prod) or self.rules.get("headgroup_reactions", {}).get(hp, {}).get("linkage_type")
                    pair_info = {
                        "reactant": {
                            "headgroup": hr,
                            "linkage": lr,
                            "chain_count": react_chain_count,
                            "carbons": react_carbons,
                            "double_bonds": react_dbs,
                            "total_carbons": sum(react_carbons) if react_carbons else None,
                            "total_double_bonds": sum(react_dbs) if react_dbs else None,
                        },
                        "product": {
                            "headgroup": hp,
                            "linkage": lp,
                            "chain_count": prod_chain_count,
                            "carbons": prod_carbons,
                            "double_bonds": prod_dbs,
                            "total_carbons": sum(prod_carbons) if prod_carbons else None,
                            "total_double_bonds": sum(prod_dbs) if prod_dbs else None,
                        },
                        "rule": {
                            "can_convert_to": rules.get("can_convert_to"),
                            "required_acyl_chains": rule.get("required_acyl_chains"),
                            "require_same_linkage": rule.get("require_same_linkage", True),
                            "is_molspecies": rule.get("is_molspecies"),
                        },
                    }
                    pairs_info.append(pair_info)
                except Exception:
                    # ignore detail assembly errors but don't fail evaluation
                    pass

        evaluation: Dict[str, Any] = {"possible": possible}
        if reasons:
            evaluation["explanation"] = "; ".join(reasons)
        else:
            # build a richer explanation summarizing the matched pairs
            if pairs_info:
                parts = []
                for p in pairs_info:
                    r = p["reactant"]
                    q = p["product"]
                    rule = p["rule"]
                    parts.append(
                        f"{r['headgroup']}->{q['headgroup']}: react(link={r['linkage']},C={r['total_carbons']},DB={r['total_double_bonds']})"
                        + f"|prod(link={q['linkage']},C={q['total_carbons']},DB={q['total_double_bonds']})"
                        + f"|require_same_linkage={rule.get('require_same_linkage')} can_convert_to={rule.get('can_convert_to')}"
                    )
                evaluation["explanation"] = "; ".join(parts)
            else:
                evaluation["explanation"] = "no_matched_pairs"

        evaluation["details"] = pairs_info
        return evaluation
