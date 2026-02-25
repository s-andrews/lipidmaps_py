import re
from typing import Any, Dict, List, Optional

lipid_reaction_rules: Dict[str, Any] = {
    "metadata": {
        "version": "0.1",
        "description": "Rules for lipid headgroup remodeling and fatty acid mass shifts for quantitation.",
        "units": "Daltons (Da)",
        "usage_note": "Conserve acyl tails (Total C:DB) during headgroup swaps.",
    },
    "headgroups": {
        "PA": {"name": "Phosphatidic Acid", "mass_shift": 79.9663, "can_convert_to": ["PC", "PE", "PS", "PI", "PG", "DAG"]},
        "PC": {"name": "Phosphatidylcholine", "mass_shift": 165.0555, "can_convert_to": ["PA", "LPC", "DAG", "SM"]},
        "PE": {"name": "Phosphatidylethanolamine", "mass_shift": 123.0085, "can_convert_to": ["PA", "PC", "PS", "LPE"]},
        "PS": {"name": "Phosphatidylserine", "mass_shift": 167.0222, "can_convert_to": ["PE", "PA"]},
        "PI": {"name": "Phosphatidylinositol", "mass_shift": 242.0192, "can_convert_to": ["PA", "PIP", "PIP2", "PIP3"]},
        "PG": {"name": "Phosphatidylglycerol", "mass_shift": 154.0031, "can_convert_to": ["CL", "PA"]},
        "DAG": {"name": "Diacylglycerol", "mass_shift": 17.0027, "can_convert_to": ["PA", "PC", "PE", "TAG"]},
        "TAG": {"name": "Triacylglycerol", "mass_shift": "VARIABLE_BY_TAIL", "can_convert_to": ["DAG"]},
    },
    "fatty_acids": {
        "12:0": {"name": "Lauric", "mass_delta": 182.1671},
        "14:0": {"name": "Myristic", "mass_delta": 210.1984},
        "14:1": {"name": "Myristoleic", "mass_delta": 208.1827},
        "16:0": {"name": "Palmitic", "mass_delta": 238.2297},
        "16:1": {"name": "Palmitoleic", "mass_delta": 236.2140},
        "18:0": {"name": "Stearic", "mass_delta": 266.2610},
        "18:1": {"name": "Oleic", "mass_delta": 264.2453},
        "18:2": {"name": "Linoleic", "mass_delta": 262.2297},
        "18:3": {"name": "Linolenic", "mass_delta": 260.2140},
        "20:0": {"name": "Arachidic", "mass_delta": 294.2923},
        "20:1": {"name": "Gadoleic", "mass_delta": 292.2766},
        "20:2": {"name": "Eicosadienoic", "mass_delta": 290.2610},
        "20:3": {"name": "DGLA", "mass_delta": 288.2453},
        "20:4": {"name": "Arachidonic", "mass_delta": 286.2297},
        "20:5": {"name": "EPA", "mass_delta": 284.2140},
        "22:0": {"name": "Behenic", "mass_delta": 322.3236},
        "22:1": {"name": "Erucic", "mass_delta": 320.3079},
        "22:4": {"name": "Adrenic", "mass_delta": 314.2610},
        "22:5": {"name": "DPA", "mass_delta": 312.2453},
        "22:6": {"name": "DHA", "mass_delta": 310.2297},
        "24:0": {"name": "Lignoceric", "mass_delta": 350.3549},
        "24:1": {"name": "Nervonic", "mass_delta": 348.3392},
    },
    "biochemical_deltas": {
        "methylation_CH2": 14.0156,
        "desaturation_H2_loss": -2.0156,
        "PE_to_PS_delta": 44.0137,
        "PE_to_PC_delta": 42.0470,
        "hydration_H2O": 18.0106,
    },
}


def _extract_headgroup_from_name(name: Optional[str]) -> Optional[str]:
    """Attempt to extract a headgroup prefix from a compound name string.

    Examples: 'PC(16:0/18:1)' -> 'PC', 'PE 16:0/18:1' -> 'PE'
    """
    if not name:
        return None
    m = re.match(r"^([A-Za-z0-9]+)", name.strip())
    if not m:
        return None
    hg = m.group(1)
    return hg if hg in lipid_reaction_rules.get("headgroups", {}) else None


def reactions_possible_in_dataset(dataset, reactions: List[Any]) -> List[Any]:
    """Return subset of `reactions` that are plausible given `dataset` and the headgroup rules.

    Heuristic used:
    - Determine which headgroups are present in the dataset by scanning lipid `input_name` prefixes.
    - For each reaction, extract headgroups from reactant and product component names.
    - A reaction is considered possible if at least one reactant headgroup is present in the dataset
      and at least one product headgroup is listed as a `can_convert_to` target for that reactant
      according to `lipid_reaction_rules`.

    This is intentionally conservative and rule-driven; it favors reactions that map known
    headgroup conversions and where the dataset actually contains the required starting headgroups.
    """
    present_hgs = set()
    for l in getattr(dataset, 'lipids', []):
        try:
            name = getattr(l, 'input_name', None)
            hg = _extract_headgroup_from_name(name)
            if hg:
                present_hgs.add(hg)
        except Exception:
            continue

    possible = []
    for rxn in reactions:
        # gather reactant and product headgroups
        reactant_hgs = set()
        product_hgs = set()
        for role in ('reactants', 'products'):
            items = getattr(rxn, role, None) or (rxn.get(role) if isinstance(rxn, dict) else None) or []
            for item in items:
                # item may be dict or object
                if isinstance(item, dict):
                    name = item.get('compound_name') or item.get('name') or item.get('compound')
                else:
                    name = getattr(item, 'compound_name', None) or getattr(item, 'name', None) or getattr(item, 'compound', None)
                hg = _extract_headgroup_from_name(name)
                if hg:
                    if role == 'reactants':
                        reactant_hgs.add(hg)
                    else:
                        product_hgs.add(hg)

        # If no reactant headgroups extracted, skip (not enough info)
        if not reactant_hgs:
            continue

        # Evaluate rule compatibility
        matched = False
        for rh in reactant_hgs:
            if rh not in lipid_reaction_rules.get('headgroups', {}):
                continue
            allowed = set(lipid_reaction_rules['headgroups'][rh].get('can_convert_to', []))
            # if any product headgroup matches allowed targets, and dataset contains the reactant headgroup,
            # mark reaction as possible
            if product_hgs and (allowed & product_hgs) and (rh in present_hgs):
                matched = True
                break
            # If products unknown, still consider plausible if dataset contains reactant
            if not product_hgs and (rh in present_hgs):
                matched = True
                break

        if matched:
            possible.append(rxn)

    return possible
