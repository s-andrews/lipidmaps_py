
import logging
import requests
logger = logging.getLogger(__name__)

"""TODO: This will be upgraded to pydantic and improved for validation and error handling"""
class ReactionChecker:
    def __init__(self, base_url, endpoint='/api/reactions'):
        self.api_url = base_url.rstrip('/') + endpoint
        self.reactions = {}
        logger.info(f"Initialized ReactionChecker with URL: {self.api_url}")

    def check_reactions(self, lm_ids, search_type="lm_id"):
        payload = {
            "search_source": "biopan",
            "search_type": search_type,
            "lm_ids": lm_ids
        }
        try:
            logger.info(f"Sending reaction check request for LM IDs: {lm_ids}")
            response = requests.post(self.api_url, json=payload, timeout=10)
            response.raise_for_status()

            data = response.json()

            def has_lm_main(comps):
                return any((c.get('compound_type') == 'lm_main') for c in (comps or []))

            def names_list(comps):
                if not comps:
                    return ''
                parts = []
                for c in comps:
                    name = c.get('compound_name') or c.get('compound_lm_id') or c.get('compound_generic_id') or ''
                    if name:
                        parts.append(name)
                return '; '.join(parts)

            filtered = []
            for reaction in data:
                reactants = reaction.get('reactants', [])
                products = reaction.get('products', []) 

                # keep only lm_main components
                filtered_reactants = [c for c in (reactants or []) if c.get('compound_type') == 'lm_main']
                filtered_products = [c for c in (products or []) if c.get('compound_type') == 'lm_main']

                # skip reaction if neither side has lm_main
                if not (filtered_reactants or filtered_products):
                    continue

                # update reaction to contain only lm_main components and name from those
                reactant_names = names_list(filtered_reactants)
                product_names = names_list(filtered_products)
                r_name = f"{reactant_names} -> {product_names}"
                reaction['reactants'] = filtered_reactants
                reaction['products'] = filtered_products
                reaction['reaction_name'] = r_name
                filtered.append(reaction)

            self.reactions['reactions'] = filtered
        except requests.RequestException as e:
            logger.error(f"Reaction API call failed: {e}")
            self.reactions['reactions'] = []
            self.reactions['error'] = str(e)
        return self.reactions

