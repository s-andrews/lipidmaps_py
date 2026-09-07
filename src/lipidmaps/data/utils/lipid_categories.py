"""LIPID MAPS category for each BioPAN lipid class.

BioPAN's class and species filter trees group lipid classes under a top-level
heading. The legacy implementation read that heading from the PostgreSQL
``biopan_reaction.class`` column (``R/lib.r::get_lp_classification``), which held
BioPAN's own pathway groupings -- "Glycerolipids and Glycerophospholipids",
"Ether lipids", "Cholesterol biosynthesis", "TCA cycle" and so on.

This module replaces that with the LIPID MAPS classification proper -- the eight
categories of ``lm_main.core`` -- with two headings kept from BioPAN's own scheme:
the glycerolipids and glycerophospholipids are shown together, and ether lipids are
shown separately. It is shipped statically, like
:mod:`lipidmaps.data.utils.fa_reactions`, so the exporter keeps its independence
from the ``biopan_reaction`` table.

Every entry below was resolved against the LIPID MAPS database rather than
assigned by hand. A few are worth knowing about because they are not what
chemical intuition suggests:

* ``MA`` / ``M5P`` / ``M5PP`` (mevalonate and its phosphates) are **Fatty Acyls**
  in LIPID MAPS -- ``Mevalonic acid``, ``Mevalonate-P`` and ``Mevalonate-PP`` all
  carry ``core = 1`` -- not Prenol Lipids.
* The isoprenoid intermediates downstream of them (``IPP``, ``DMAPP``, ``GPP``,
  ``FPP``, ``SQ``) *are* Prenol Lipids.
* Ether and plasmalogen classes (``O-PC``, ``P-PE``, ``O-DG``, ...) are broken out
  under "Ether lipids". LIPID MAPS itself files them under Glycerophospholipids and
  Glycerolipids -- there is no ``core`` for them -- but BioPAN has always shown them
  as their own heading, and that is more useful here than burying them among the
  diacyl classes.
* LIPID MAPS cores 2 and 3 are shown together as "Glycerolipids and
  Glycerophospholipids", as the legacy ``biopan_reaction.class`` did.
  :data:`CORE_CATEGORIES` below still records the two separately, since it is
  reference data for ``lm_main.core``; only the display heading is merged.

To regenerate from the source database::

    -- the class names BioPAN uses
    SELECT DISTINCT trim(unnest(string_to_array(reaction, ','))) AS lp
      FROM biopan_reaction WHERE type = 'lipid';

    -- and the LIPID MAPS category for a representative member of each
    SELECT name, core FROM lm_main WHERE name ILIKE '<representative name>';
"""

from __future__ import annotations

from typing import Dict, List

# lm_main.core -> LIPID MAPS category name, in the canonical LIPID MAPS order.
CORE_CATEGORIES: Dict[int, str] = {
    1: "Fatty Acyls",
    2: "Glycerolipids",
    3: "Glycerophospholipids",
    4: "Sphingolipids",
    5: "Sterol Lipids",
    6: "Prenol Lipids",
    7: "Saccharolipids",
    8: "Polyketides",
}

# Headings kept from BioPAN's scheme rather than taken straight from lm_main.core:
# cores 2 and 3 shown together, and ether/plasmalogen classes broken out.
GLYCERO_LIPIDS = "Glycerolipids and Glycerophospholipids"
ETHER_LIPIDS = "Ether lipids"

# Display order for the filter tree; anything unmapped sorts last. Follows the
# canonical LIPID MAPS order, with Ether lipids after the glycerophospholipids
# they are derived from.
CATEGORY_ORDER: List[str] = [
    "Fatty Acyls",
    GLYCERO_LIPIDS,
    ETHER_LIPIDS,
    "Sphingolipids",
    "Sterol Lipids",
    "Prenol Lipids",
    "Saccharolipids",
    "Polyketides",
]

# Heading used for classes with no LIPID MAPS category.
UNCLASSIFIED = "Other"

LIPID_CATEGORY: Dict[str, str] = {
    # --- Fatty Acyls -------------------------------------------------------
    "FA": "Fatty Acyls",
    "FACOA": "Fatty Acyls",
    "Acetyl-CoA": "Fatty Acyls",
    "Acetoacetyl-CoA": "Fatty Acyls",
    "acetoacetyl-Co-A": "Fatty Acyls",
    "Butyryl-CoA": "Fatty Acyls",
    "Crotonoyl-CoA": "Fatty Acyls",
    "Propionyl-CoA": "Fatty Acyls",
    "Methylmalonyl-CoA": "Fatty Acyls",
    "Ethylmalonyl-CoA": "Fatty Acyls",
    "3-HMG-CoA": "Fatty Acyls",
    "MA": "Fatty Acyls",
    "M5P": "Fatty Acyls",
    "M5PP": "Fatty Acyls",

    # --- Glycerolipids and Glycerophospholipids (lm_main.core 2 and 3) -----
    "MG": GLYCERO_LIPIDS,
    "DG": GLYCERO_LIPIDS,
    "TG": GLYCERO_LIPIDS,
    "PA": GLYCERO_LIPIDS,
    "PC": GLYCERO_LIPIDS,
    "PE": GLYCERO_LIPIDS,
    "PG": GLYCERO_LIPIDS,
    "PI": GLYCERO_LIPIDS,
    "PS": GLYCERO_LIPIDS,
    "PIP": GLYCERO_LIPIDS,
    "PIP2": GLYCERO_LIPIDS,
    "PIP3": GLYCERO_LIPIDS,
    "LPA": GLYCERO_LIPIDS,
    "LPC": GLYCERO_LIPIDS,
    "LPE": GLYCERO_LIPIDS,
    "LPG": GLYCERO_LIPIDS,
    "LPI": GLYCERO_LIPIDS,
    "LPS": GLYCERO_LIPIDS,
    "CDP-DG": GLYCERO_LIPIDS,
    "CL": GLYCERO_LIPIDS,

    # --- Ether lipids ------------------------------------------------------
    # alkyl (O-) and plasmalogen (P-) forms; LIPID MAPS files these under
    # Glycerophospholipids and Glycerolipids, BioPAN shows them separately.
    "O-PC": ETHER_LIPIDS,
    "O-PE": ETHER_LIPIDS,
    "O-LPA": ETHER_LIPIDS,
    "O-LPC": ETHER_LIPIDS,
    "O-LPE": ETHER_LIPIDS,
    "O-DG": ETHER_LIPIDS,
    "P-PC": ETHER_LIPIDS,
    "P-PE": ETHER_LIPIDS,
    "P-LPA": ETHER_LIPIDS,
    "P-LPC": ETHER_LIPIDS,
    "P-LPE": ETHER_LIPIDS,

    # --- Sphingolipids -----------------------------------------------------
    "SPB": "Sphingolipids",
    "SPBP": "Sphingolipids",
    "SPC": "Sphingolipids",
    "dhSPB": "Sphingolipids",
    "dhSPBP": "Sphingolipids",
    "Cer": "Sphingolipids",
    "dhCer": "Sphingolipids",
    "Cer1P": "Sphingolipids",
    "dhCer1P": "Sphingolipids",
    "Cer-PE": "Sphingolipids",
    "SM": "Sphingolipids",
    "dhSM": "Sphingolipids",
    "Gal-Cer": "Sphingolipids",
    "Glc-Cer": "Sphingolipids",
    "Hex-Cer": "Sphingolipids",
    "Hex-dhCer": "Sphingolipids",
    "Lac-Cer": "Sphingolipids",
    "Lac-dhCer": "Sphingolipids",
    "Gal-SPB": "Sphingolipids",

    # --- Sterol Lipids -----------------------------------------------------
    "ST": "Sterol Lipids",
    "CE": "Sterol Lipids",
    "cholesterol": "Sterol Lipids",
    "desmosterol": "Sterol Lipids",
    "lanosterol": "Sterol Lipids",
    "lathosterol": "Sterol Lipids",
    "zymosterol": "Sterol Lipids",
    "zymostenol": "Sterol Lipids",
    "7-DHC": "Sterol Lipids",
    "8-DHC": "Sterol Lipids",
    "7-dehydrodesmosterol": "Sterol Lipids",
    "T-MAS": "Sterol Lipids",
    "4beta-HC": "Sterol Lipids",
    "7alpha-HC": "Sterol Lipids",
    "24-HC": "Sterol Lipids",
    "25-HC": "Sterol Lipids",
    "27-HC": "Sterol Lipids",
    "24_25-epoxycholesterol": "Sterol Lipids",
    "24_25-epoxylanosterol": "Sterol Lipids",

    # --- Prenol Lipids -----------------------------------------------------
    "IPP": "Prenol Lipids",
    "DMAPP": "Prenol Lipids",
    "GPP": "Prenol Lipids",
    "FPP": "Prenol Lipids",
    "SQ": "Prenol Lipids",
    "2_3-epoxySQ": "Prenol Lipids",
    "2_3;24_25-diepoxySQ": "Prenol Lipids",
}

# Case-insensitive index, built once.
_LOOKUP: Dict[str, str] = {key.lower(): value for key, value in LIPID_CATEGORY.items()}


def category_for(class_name: str) -> str:
    """LIPID MAPS category for a BioPAN lipid class.

    Falls back to :data:`UNCLASSIFIED` so a class added to the reaction network
    without a mapping still appears in the tree instead of disappearing.
    """
    if not class_name:
        return UNCLASSIFIED
    return _LOOKUP.get(class_name.strip().lower(), UNCLASSIFIED)


def category_sort_key(category: str) -> tuple:
    """Sort categories into canonical LIPID MAPS order, unmapped ones last."""
    try:
        return (CATEGORY_ORDER.index(category), category)
    except ValueError:
        return (len(CATEGORY_ORDER), category)
