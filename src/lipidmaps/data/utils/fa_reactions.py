"""Fatty-acid reaction network for BioPAN pathway analysis.

This is the curated FA-to-FA conversion network (elongation / desaturation)
that the legacy BioPAN served from the PostgreSQL ``biopan_reaction`` table
(rows where ``type = 'fa'``). It is shipped statically here so the new
``lipidmaps_py`` exporter can build the "Fatty acid" pathway view without a
dependency on that database.

Each entry is a single-step FA conversion ``reactant -> product`` with the
enzyme gene symbols associated with the step. ``sub_class`` is the reaction kind
shown in BioPAN ("Elongation" or "Desaturation").

To regenerate from the source database:

    SELECT r.reaction, r.sub_class,
           array_to_string(array_agg(g.gene_symbol), ',') AS genes
    FROM biopan_reaction r
    LEFT JOIN biopan_reaction_gene rg ON r.reaction_id = rg.reaction_id
    LEFT JOIN biopan_gene g ON g.gene_id = rg.gene_id
    WHERE r.type = 'fa'
    GROUP BY r.reaction_id, r.reaction, r.sub_class;
"""

from __future__ import annotations

from typing import Dict, List, Tuple

# (reactant, product, sub_class, genes)
_FA_REACTIONS: Tuple[Tuple[str, str, str, Tuple[str, ...]], ...] = (
    ("FA(16:0)", "FA(16:1)", "Desaturation", ("SCD3",)),
    ("FA(16:0)", "FA(18:0)", "Elongation", ("ELOVL6",)),
    ("FA(16:1)", "FA(16:2)", "Desaturation", ("FADS2",)),
    ("FA(16:1)", "FA(18:1)", "Elongation", ("ELOVL5", "ELOVL6")),
    ("FA(18:0)", "FA(18:1)", "Desaturation", ("SCD1",)),
    ("FA(18:0)", "FA(20:0)", "Elongation", ("ELOVL1", "ELOVL7", "ELOVL3")),
    ("FA(18:1)", "FA(18:2)", "Desaturation", ("FADS2",)),
    ("FA(18:1)", "FA(20:1)", "Elongation", ("ELOVL3",)),
    ("FA(18:2)", "FA(18:3)", "Desaturation", ("FADS2",)),
    ("FA(18:2)", "FA(20:2)", "Elongation", ("ELOVL5",)),
    ("FA(18:3)", "FA(18:4)", "Desaturation", ("FADS2",)),
    ("FA(18:3)", "FA(20:3)", "Elongation", ("ELOVL5", "ELOVL7")),
    ("FA(18:4)", "FA(20:4)", "Elongation", ("ELOVL5",)),
    ("FA(20:0)", "FA(22:0)", "Elongation", ("ELOVL1", "ELOVL7", "ELOVL3")),
    ("FA(20:1)", "FA(22:1)", "Elongation", ("ELOVL3",)),
    ("FA(20:2)", "FA(20:3)", "Desaturation", ("FADS1",)),
    ("FA(20:3)", "FA(20:4)", "Desaturation", ("FADS1",)),
    ("FA(20:4)", "FA(20:5)", "Desaturation", ("FADS1",)),
    ("FA(20:4)", "FA(22:4)", "Elongation", ("ELOVL5", "ELOVL2")),
    ("FA(20:5)", "FA(22:5)", "Elongation", ("ELOVL5", "ELOVL2")),
    ("FA(22:0)", "FA(24:0)", "Elongation", ("ELOVL1", "ELOVL7", "ELOVL3")),
    ("FA(22:1)", "FA(24:1)", "Elongation", ("ELOVL3",)),
    ("FA(22:4)", "FA(24:4)", "Elongation", ("ELOVL2",)),
    ("FA(22:5)", "FA(24:5)", "Elongation", ("ELOVL2",)),
    ("FA(24:0)", "FA(26:0)", "Elongation", ("ELOVL1", "ELOVL3")),
    ("FA(24:4)", "FA(24:5)", "Desaturation", ("FADS2",)),
    ("FA(24:5)", "FA(24:6)", "Desaturation", ("FADS2",)),
    ("FA(24:6)", "FA(34:6)", "Elongation", ("ELOVL2", "ELOVL4")),
    ("FA(26:0)", "FA(28:0)", "Elongation", ("ELOVL4",)),
    ("FA(28:0)", "FA(30:0)", "Elongation", ("ELOVL4",)),
)


def get_fa_reactions() -> List[Dict[str, object]]:
    """Return the FA reaction network as a list of dicts."""

    return [
        {
            "reactant": reactant,
            "product": product,
            "sub_class": sub_class,
            "genes": list(genes),
        }
        for reactant, product, sub_class, genes in _FA_REACTIONS
    ]
