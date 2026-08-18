"""Shunt / alternate-route pathway detection.

A "shunt" (or alternate route) is a pathway branch that produces the same or
related products via a different sequence of reactions -- the canonical example
being the "Cholesterol biosynthesis (Shunt)" branch of sterol biosynthesis.

Shunt is treated as a *general* concept, not a cholesterol-only special case.

Detection precedence (most authoritative first):
  1. A structured boolean flag on the pathway (``is_shunt`` / ``shunt``) -- this
     is what the LIPID MAPS API emits once the ``reaction_pathways.is_shunt``
     column is populated.
  2. A structured classification field containing a shunt marker
     (``pathway_class`` / ``reaction_class``).
  3. A text fallback: a shunt marker appearing in the pathway name
     (``pathway_name`` / ``name``) or type (``pathway_type`` / ``type``). This
     keeps detection working on legacy data that only carries the "(Shunt)"
     naming convention.

Precedence rule: any positive signal means shunt. An explicit structured
``is_shunt=False`` is authoritative only for *false* -- it suppresses the text
fallback so a curator can override a misleading name.
"""

from typing import Any, Iterable, Mapping, Optional

# Case-insensitive substrings that mark a pathway as a shunt / alternate route.
SHUNT_MARKERS: tuple = ("shunt", "alternate route", "alt route")


def _has_marker(value: Any) -> bool:
    """True if ``value`` (a string or iterable of strings) contains a marker."""
    if value is None:
        return False
    if isinstance(value, str):
        text = value.lower()
        return any(marker in text for marker in SHUNT_MARKERS)
    if isinstance(value, Mapping):
        return False
    if isinstance(value, Iterable):
        return any(_has_marker(item) for item in value)
    return False


def is_shunt_pathway(pathway: Optional[Mapping[str, Any]]) -> bool:
    """Return True if a single pathway dict denotes a shunt / alternate route."""
    if not isinstance(pathway, Mapping):
        return False

    # 1. Structured boolean flag is authoritative (for both True and False).
    for key in ("is_shunt", "shunt"):
        if key in pathway and pathway[key] is not None:
            return bool(pathway[key])

    # 2. Structured classification field.
    for key in ("pathway_class", "reaction_class"):
        if _has_marker(pathway.get(key)):
            return True

    # 3. Text fallback on name / type.
    for key in ("pathway_name", "name", "pathway_type", "type"):
        if _has_marker(pathway.get(key)):
            return True

    return False


def reaction_is_shunt(pathways: Optional[Iterable[Mapping[str, Any]]]) -> bool:
    """Return True if ANY pathway on the reaction is a shunt pathway."""
    if not pathways:
        return False
    return any(is_shunt_pathway(pathway) for pathway in pathways)
