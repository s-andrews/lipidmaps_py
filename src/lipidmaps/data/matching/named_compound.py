"""Named-compound reaction matcher.

Handles reactions between specific named molecular species that change
composition and carry no acyl chains -- primarily sterol biosynthesis and its
shunt / alternate-route branch (e.g. lanosterol -> cholesterol). Such reactions
cannot be paired by acyl composition (as SAME_STRUCTURE does), so they are
matched purely by compound identity/presence: if the reactant class and the
product class are both present in the dataset, the pair is emitted.
"""

from ..models.species_reaction import (
    ClassReaction,
    ReactionMatchResult,
    ReactionType,
)
from .base import MatcherContext, ReactionMatcher


class NamedCompoundMatcher(ReactionMatcher):
    """Matches named-species pairs by identity/presence, ignoring composition.

    Used for sterol / shunt reactions where the reactant and product are
    specific compounds rather than chain-varying lipid classes.
    """

    @property
    def reaction_type(self) -> ReactionType:
        return ReactionType.NAMED_COMPOUND

    def match(
        self,
        class_reaction: ClassReaction,
        context: MatcherContext,
    ) -> ReactionMatchResult:
        reactants = context.get_species(class_reaction.reactant_class)
        products = context.get_species(class_reaction.product_class)

        result = self._create_result(class_reaction, reactants, products)

        if not reactants or not products:
            return result

        # Presence-based pairing: composition is intentionally ignored. Normally
        # 1:1 since each named class maps to a single species, but pair the full
        # cross product to be robust against duplicates.
        for reactant in reactants:
            for product in products:
                result.pairs.append(
                    self._create_pair(reactant, product, class_reaction)
                )

        result.pairs_matched = len(result.pairs)
        return result
