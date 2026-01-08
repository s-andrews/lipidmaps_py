#!/usr/bin/env python3

import lipidmaps

##################################
# Data Import and Regularisation #
##################################

# Go to the web API and reguarlise the lipid names in the file
#
# We might have auto-detect for certain standard formats from
# upstream programs (eg msdial)
lipid_data = lipidmaps.import_data("mydata.csv", lipid_col=1, sample_cols=[4, 5, 6, 7])

lipid_data = lipidmaps.import_msdial("mydata_msdial.csv")


###############################
# Properties of imported data #
###############################

# These worked
print(f"Imported Lipids: {lipid_data.successful_import_count()}")

# These didn't
print(f"Unrecognised Lipids: {lipid_data.failed_import_count()}")

# Maybe something about lipids where we understood the name but we didn't have an LMID for it?

# We can also get a list of the names which failed.
failed_names = lipid_data.failed_import_names()

for failed_name in failed_names:
    print(f"Failed to import {failed_name}")


##################################
# Data Access in imported lipids #
##################################

# Names of samples in the file
sample_names = lipid_data.samples()

# Object for lipids
lipids = lipid_data.lipids()

for lipid in lipids:
    print(f"Original name: {lipid.users_name}")
    print(f"Official name: {lipid.offical_name}")
    print(f"Specific LMID: {lipid.lmid}")
    print(
        f"Generic LMID: {lipid.generic_lmid}"
    )  # Do we need to allow multiple of these?

    # Access the users quantitated data for each lipid
    for sample in sample_names:
        print(
            f"Value for {lipid.official_name} in {sample} is {lipid_data.get_value_for_lipid(lipid, sample)}"
        )

#############
# Reactions #
#############

# Extract the set of reactions for the imported lipids via another
# lipidmaps web api call

# They need to say which species the data comes from, and whether they
# want only complete reactions (all lipid reactants and products were observed)
# or also include reactions where just some components were present.  We could
# possibly also include an option to retrieve all reactions for the species
# regardless of whether the lipids in our data were present.

reactions = lipid_data.get_reactions(species="human", complete=True)

# Access the lipid data in reactions

for reaction in reactions:
    print(f"Looking at reaction {reaction.name}")

    reactants = reaction.reactants()
    for reactant in reactants:
        if reactant.is_lipid():
            print(
                f"Looking at reactant {reactant.official_name} with id {reactant.lmid}"
            )

            # The lipids which are valid here are a function of the dataset
            # not the sample
            valid_lipids = lipid_data.get_lipids_for_reaction_component(reactant)

            for valid_lipid in valid_lipids:
                print(
                    f"Lipid {valid_lipid.official_name} is valid for component {reactant.name} in reaction {reaction.name}"
                )

            # For quantitation we shouldn't use the valid lipids but let the API do it so
            # we can't get it wrong.
            for sample in sample_names:
                reactant_quantitation = lipid_data.get_value_for_reaction_component(
                    reactant, method="sum"
                )
                print(
                    f"Value for {reactant.name} in {reaction.name} in {sample} is {reactant_quantitation}"
                )

    # We should have the same logic for products
    products = reaction.products()

    # Should we also have a method to pair up reactants and products?
    # We'll need this internally I guess? Is it more complex than this
    # when we're dealing with fatty acids?

    paired_components = reaction.paired_components()
    for i, paired_component in enumerate(paired_components):
        print(
            f"For pair {i} reactant is {paired_component[0].official_name} product is {paired_component[1].official_name}"
        )


###################
# Reaction Chains #
###################

reaction_chains = reactions.make_reaction_chains()

# Not quite sure of the best API here, but we need some way to iterate through the reactions
# so something like this should suffice for the API

for reaction_node in reaction_chains:
    if reaction_node.has_children():
        sub_reactions = reaction_node.children()

# The logic for the iterration could then be in the user's code.
