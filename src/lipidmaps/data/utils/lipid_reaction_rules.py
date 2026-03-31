import re
from typing import Any, Dict, List, Optional

lipid_reaction_rules: Dict[str, Any] = {
    "metadata": {
        "version": "0.2-delta-schema",
        "description": "Lipid headgroup and acyl-chain remodeling rules using delta-based schema.",
        "units": "Daltons (Da)",
        "usage_note": "Acyl-chain conservation is enforced via delta rules; linkage type is preserved unless stated."
    },

    "headgroup_reactions": {

        # ----------------------------------------------------------
        # PC / LPC (ester)
        # ----------------------------------------------------------
        "PC": {
            "name": "Phosphatidylcholine",
            "linkage_type": "ester",
            "acyl_chains": 2,
            "mass_shift": 165.0555,
            "can_convert_to": ["PA", "PS", "LPC", "DG", "SM"],

            "conversion_rules": {

                # PC → LPC (PLA2)
                "LPC": {
                    "reaction_requirements": {"external_compounds": []},
                    "acyl_chain_change": -1,
                    "require_same_linkage": True,
                    "reaction_type": "deacylation"
                },

                # PC → PS (base-exchange with serine)
                "PS": {
                    "reaction_requirements": {"external_compounds": ["serine"]},
                    "acyl_chain_change": 0,
                    "require_same_linkage": True,
                    "reaction_type": "headgroup_swap"
                },

                # PC → PA (phospholipase D)
                "PA": {
                    "reaction_requirements": {"external_compounds": ["serine"]},
                    "acyl_chain_change": 0,
                    "require_same_linkage": True,
                    "reaction_type": "headgroup_swap"
                },

                # ✅ PC → DG (phospholipase C)
                "DG": {
                    "reaction_requirements": {"external_compounds": []},
                    "acyl_chain_change": 0,
                    "require_same_linkage": True,
                    "reaction_type": "headgroup_cleavage"
                },

                # ✅ PC → SM (sphingomyelin synthase)
                "SM": {
                    "reaction_requirements": {"external_compounds": ["Cer"]},
                    "acyl_chain_change": 0,
                    "require_same_linkage": True,
                    "reaction_type": "headgroup_transfer"
                }
            }
        },

        "LPC": {
            "name": "Lyso-PC",
            "linkage_type": "ester",
            "acyl_chains": 1,
            "mass_shift": None,
            "can_convert_to": ["PC"],
            "conversion_rules": {
                "PC": {
                    "reaction_requirements": {"external_compounds": ["acylcoa"]},
                    "acyl_chain_change": 1,
                    "require_same_linkage": True,
                    "reaction_type": "acylation"
                }
            }
        },

        # ----------------------------------------------------------
        # Ether PC O-
        # ----------------------------------------------------------
        "PC O-": {
            "name": "Ether phosphatidylcholine",
            "linkage_type": "ether_alkyl",
            "acyl_chains": 1,
            "mass_shift": None,
            "can_convert_to": ["LPC O-", "DG O-"],
            "conversion_rules": {
                # Deacylation (PLA2-like)
                "LPC O-": {
                    "reaction_requirements": {"external_compounds": []},
                    "acyl_chain_change": 0,
                    "require_same_linkage": True,
                    "reaction_type": "deacylation"
                },
                
                # Headgroup cleavage (PLC-like) → ether DG
                "DG O-": {
                    "reaction_requirements": {"external_compounds": []},
                    "acyl_chain_change": 0,
                    "require_same_linkage": True,
                    "reaction_type": "headgroup_cleavage"
                }

            }
        },

        "DG O-": {
            "name": "Ether diacylglycerol",
            "linkage_type": "ether_alkyl",
            "acyl_chains": 1,   # 1 ether + 1 acyl chain total = 1 modifiable chain
            "mass_shift": None,
            "can_convert_to": ["TG O-", "LPE O-", "PC O-"],

            "conversion_rules": {

                # DG O- → TG O- (acylation)
                "TG O-": {
                    "reaction_requirements": {"external_compounds": ["acylcoa"]},
                    "acyl_chain_change": 1,
                    "require_same_linkage": True,
                    "reaction_type": "acylation"
                },

                # DG O- → PC O- (headgroup attachment)
                "PC O-": {
                    "reaction_requirements": {"external_compounds": ["CDP-choline"]},
                    "acyl_chain_change": 0,
                    "require_same_linkage": True,
                    "reaction_type": "headgroup_attachment"
                },

                # DG O- → LPE O- (deacylation)
                "LPE O-": {
                    "reaction_requirements": {"external_compounds": []},
                    "acyl_chain_change": -1,
                    "require_same_linkage": True,
                    "reaction_type": "deacylation"
                }
            }
        },

        "TG O-": {
            "name": "Ether triacylglycerol",
            "linkage_type": "ether_alkyl",
            "acyl_chains": 2,   # 1 ether + 2 acyl
            "mass_shift": "VARIABLE_BY_TAIL",
            "can_convert_to": ["DG O-"],

            "conversion_rules": {
                "DG O-": {
                    "reaction_requirements": {"external_compounds": []},
                    "acyl_chain_change": -1,
                    "require_same_linkage": True,
                    "reaction_type": "deacylation"
                }
            }
        },

        "TG P-": {
            "name": "Plasmalogen TG (vinyl ether TG)",
            "linkage_type": "ether_vinyl",
            "acyl_chains": 2,
            "mass_shift": "VARIABLE_BY_TAIL",
            "can_convert_to": ["DG P-"],

            "conversion_rules": {
                "DG P-": {
                    "reaction_requirements": {"external_compounds": []},
                    "acyl_chain_change": -1,
                    "require_same_linkage": True,
                    "reaction_type": "deacylation"
                }
            }
        },
        "LPC O-": {
            "name": "Ether lyso-PC",
            "linkage_type": "ether_alkyl",
            "acyl_chains": 1,
            "mass_shift": None,
            "can_convert_to": ["PC O-"],
            "conversion_rules": {
                "PC O-": {
                    "reaction_requirements": {"external_compounds": ["acylcoa"]},
                    "acyl_chain_change": 1,
                    "require_same_linkage": True,
                    "reaction_type": "acylation"
                }
            }
        },

        # ----------------------------------------------------------
        # Plasmalogen PC P-
        # ----------------------------------------------------------
        "PC P-": {
            "name": "Plasmalogen PC (vinyl ether)",
            "linkage_type": "ether_vinyl",
            "acyl_chains": 1,
            "mass_shift": None,
            "can_convert_to": ["LPC P-"],
            "conversion_rules": {
                "LPC P-": {
                    "reaction_requirements": {"external_compounds": []},
                    "acyl_chain_change": 0,
                    "require_same_linkage": True,
                    "reaction_type": "plasmalogen"
                }
            }
        },

        "LPC P-": {
            "name": "Lyso-plasmalogen PC",
            "linkage_type": "ether_vinyl",
            "acyl_chains": 1,
            "mass_shift": None,
            "can_convert_to": ["PC P-"],
            "conversion_rules": {
                "PC P-": {
                    "reaction_requirements": {"external_compounds": ["acylcoa"]},
                    "acyl_chain_change": 1,
                    "require_same_linkage": True,
                    "reaction_type": "acylation"
                }
            }
        },

        # ----------------------------------------------------------
        # PE / LPE (ester)
        # ----------------------------------------------------------
        "PE": {
            "name": "Phosphatidylethanolamine",
            "linkage_type": "ester",
            "acyl_chains": 2,
            "mass_shift": 123.0085,
            "can_convert_to": ["PA", "PC", "PS", "LPE"],

            "conversion_rules": {

                # ----------------------------------------------------------
                # PE → LPE (PLA2 deacylation)
                # ----------------------------------------------------------
                "LPE": {
                    "reaction_requirements": {"external_compounds": []},
                    "acyl_chain_change": -1,
                    "require_same_linkage": True,
                    "reaction_type": "deacylation"
                },

                # ----------------------------------------------------------
                # PE → PC (PEMT methylation pathway)
                # Requires 3 SAM molecules, but we only model externally as "SAM"
                # ----------------------------------------------------------
                "PC": {
                    "reaction_requirements": {"external_compounds": ["SAM"]},
                    "acyl_chain_change": 0,
                    "require_same_linkage": True,
                    "reaction_type": "methylation"
                },

                # ----------------------------------------------------------
                # PE → PS (base exchange, PSS1/2)
                # ----------------------------------------------------------
                "PS": {
                    "reaction_requirements": {"external_compounds": ["serine"]},
                    "acyl_chain_change": 0,
                    "require_same_linkage": True,
                    "reaction_type": "headgroup_swap"
                },

                # ----------------------------------------------------------
                # PE → PA (PLD-like phosphodiester cleavage)
                # ----------------------------------------------------------
                "PA": {
                    "reaction_requirements": {"external_compounds": []},
                    "acyl_chain_change": 0,
                    "require_same_linkage": True,
                    "reaction_type": "headgroup_cleavage"
                }

                # ----------------------------------------------------------
                # OPTIONAL: PE → DG (PLC-like cleavage)
                # Rare biologically but allowed for model symmetry with PC → DG.
                # Uncomment if desired:
                #
                # "DG": {
                #     "reaction_requirements": {"external_compounds": []},
                #     "acyl_chain_change": 0,
                #     "require_same_linkage": True,
                #     "reaction_type": "headgroup_cleavage"
                # }
                # ----------------------------------------------------------
            }
        },

        "LPE": {
            "name": "Lyso-PE",
            "linkage_type": "ester",
            "acyl_chains": 1,
            "mass_shift": None,
            "can_convert_to": ["PE"],
            "conversion_rules": {
                "PE": {
                    "reaction_requirements": {"external_compounds": ["acylcoa"]},
                    "acyl_chain_change": 1,
                    "require_same_linkage": True,
                    "reaction_type": "acylation"
                }
            }
        },

        "DG P-": {
            "name": "Plasmalogen DG (vinyl ether DG)",
            "linkage_type": "ether_vinyl",
            "acyl_chains": 1,
            "mass_shift": None,
            "can_convert_to": ["TG P-", "PC P-", "LPE P-"],

            "conversion_rules": {

                # DG P- → TG P- (acylation)
                "TG P-": {
                    "reaction_requirements": {"external_compounds": ["acylcoa"]},
                    "acyl_chain_change": 1,
                    "require_same_linkage": True,
                    "reaction_type": "acylation"
                },

                # DG P- → PC P- (CDP-choline attachment)
                "PC P-": {
                    "reaction_requirements": {"external_compounds": ["CDP-choline"]},
                    "acyl_chain_change": 0,
                    "require_same_linkage": True,
                    "reaction_type": "headgroup_attachment"
                },

                # DG P- → LPE P- (deacylation)
                "LPE P-": {
                    "reaction_requirements": {"external_compounds": []},
                    "acyl_chain_change": -1,
                    "require_same_linkage": True,
                    "reaction_type": "deacylation"
                }
            }
        },
        # ----------------------------------------------------------
        # Ether PE O-
        # ----------------------------------------------------------
        "PE O-": {
            "name": "Ether PE",
            "linkage_type": "ether_alkyl",
            "acyl_chains": 1,
            "mass_shift": None,
            "can_convert_to": ["LPE O-"],
            "conversion_rules": {
                "LPE O-": {
                    "reaction_requirements": {"external_compounds": []},
                    "acyl_chain_change": 0,
                    "require_same_linkage": True,
                    "reaction_type": "deacylation"
                }
            }
        },

        "LPE O-": {
            "name": "Ether LPE",
            "linkage_type": "ether_alkyl",
            "acyl_chains": 1,
            "mass_shift": None,
            "can_convert_to": ["PE O-"],
            "conversion_rules": {
                "PE O-": {
                    "reaction_requirements": {"external_compounds": ["acylcoa"]},
                    "acyl_chain_change": 1,
                    "require_same_linkage": True,
                    "reaction_type": "acylation"
                }
            }
        },

        # ----------------------------------------------------------
        # Plasmalogen PE P-
        # ----------------------------------------------------------
        "PE P-": {
            "name": "Plasmalogen PE",
            "linkage_type": "ether_vinyl",
            "acyl_chains": 1,
            "can_convert_to": ["LPE P-"],
            "conversion_rules": {
                "LPE P-": {
                    "reaction_requirements": {"external_compounds": []},
                    "acyl_chain_change": 0,
                    "require_same_linkage": True,
                    "reaction_type": "plasmalogen"
                }
            }
        },

        "LPE P-": {
            "name": "Lyso-plasmalogen PE",
            "linkage_type": "ether_vinyl",
            "acyl_chains": 1,
            "can_convert_to": ["PE P-"],
            "conversion_rules": {
                "PE P-": {
                    "reaction_requirements": {"external_compounds": ["acylcoa"]},
                    "acyl_chain_change": 1,
                    "require_same_linkage": True,
                    "reaction_type": "acylation"
                }
            }
        },

        # ----------------------------------------------------------
        # PG
        # ----------------------------------------------------------
        "PG": {
            "name": "Phosphatidylglycerol",
            "linkage_type": "ester",
            "acyl_chains": 2,
            "mass_shift": 154.0031,

            "can_convert_to": ["CL", "PA", "LPG"],

            "conversion_rules": {

                # ----------------------------------------------------------
                # PG → CL (cardiolipin synthesis via CLS)
                # ----------------------------------------------------------
                "CL": {
                    "reaction_requirements": {"external_compounds": ["CDP-DAG"]},
                    "acyl_chain_change": 0,
                    "require_same_linkage": True,
                    "reaction_type": "headgroup_condensation"
                },

                # ----------------------------------------------------------
                # PG → PA (phospholipase-D-like cleavage)
                # Removes glycerol headgroup.
                # ----------------------------------------------------------
                "PA": {
                    "reaction_requirements": {"external_compounds": []},
                    "acyl_chain_change": 0,
                    "require_same_linkage": True,
                    "reaction_type": "headgroup_cleavage"
                },

                # ----------------------------------------------------------
                # PG → LPG (PLA2 deacylation)
                # Mirrors PC→LPC and PE→LPE.
                # ----------------------------------------------------------
                "LPG": {
                    "reaction_requirements": {"external_compounds": []},
                    "acyl_chain_change": -1,
                    "require_same_linkage": True,
                    "reaction_type": "deacylation"
                }
            }
        },

        # ----------------------------------------------------------
        # LPA / PA
        # ----------------------------------------------------------
        "LPA": {
            "name": "Lyso-PA",
            "acyl_chains": 1,
            "linkage_type": "ester",
            "can_convert_to": ["PA"],
            "conversion_rules": {
                "PA": {
                    "reaction_requirements": {"external_compounds": ["acylcoa"]},
                    "acyl_chain_change": 1,
                    "require_same_linkage": True,
                    "reaction_type": "acylation"
                }
            }
        },

        "PA": {
            "name": "Phosphatidic acid",
            "linkage_type": "ester",
            "acyl_chains": 2,
            "mass_shift": 79.9663,
            "can_convert_to": ["PC", "PE", "PS", "PI", "PG", "DG", "LPA"],

            "conversion_rules": {

                # ----------------------------------------------------------
                # PA → DG (phosphatidic acid phosphatase / Lipin)
                # ----------------------------------------------------------
                "DG": {
                    "reaction_requirements": {"external_compounds": []},
                    "acyl_chain_change": 0,
                    "require_same_linkage": True,
                    "reaction_type": "dephosphorylation"
                },

                # ----------------------------------------------------------
                # PA → LPA (PLA2 deacylation)
                # ----------------------------------------------------------
                "LPA": {
                    "reaction_requirements": {"external_compounds": []},
                    "acyl_chain_change": -1,
                    "require_same_linkage": True,
                    "reaction_type": "deacylation"
                },

                # ----------------------------------------------------------
                # PA → PC (Kennedy pathway: CDP-choline + DAG → PC)
                # Modeled directly as headgroup addition to PA.
                # ----------------------------------------------------------
                "PC": {
                    "reaction_requirements": {"external_compounds": ["CDP-choline"]},
                    "acyl_chain_change": 0,
                    "require_same_linkage": True,
                    "reaction_type": "headgroup_attachment"
                },

                # ----------------------------------------------------------
                # PA → PE (Kennedy pathway: CDP-ethanolamine)
                # ----------------------------------------------------------
                "PE": {
                    "reaction_requirements": {"external_compounds": ["CDP-ethanolamine"]},
                    "acyl_chain_change": 0,
                    "require_same_linkage": True,
                    "reaction_type": "headgroup_attachment"
                },

                # ----------------------------------------------------------
                # PA → PS (via CDP-DAG → PS synthase)
                # ----------------------------------------------------------
                "PS": {
                    "reaction_requirements": {"external_compounds": ["CDP-DAG", "serine"]},
                    "acyl_chain_change": 0,
                    "require_same_linkage": True,
                    "reaction_type": "headgroup_attachment"
                },

                # ----------------------------------------------------------
                # PA → PI (via CDP-DAG → PI synthase)
                # ----------------------------------------------------------
                "PI": {
                    "reaction_requirements": {"external_compounds": ["CDP-DAG", "inositol"]},
                    "acyl_chain_change": 0,
                    "require_same_linkage": True,
                    "reaction_type": "headgroup_attachment"
                },

                # ----------------------------------------------------------
                # PA → PG (via CDP-DAG → PGP → PG)
                # ----------------------------------------------------------
                "PG": {
                    "reaction_requirements": {"external_compounds": ["CDP-DAG", "glycerol-3-phosphate"]},
                    "acyl_chain_change": 0,
                    "require_same_linkage": True,
                    "reaction_type": "headgroup_attachment"
                }
            }
        },

        # ----------------------------------------------------------
        # Sphingolipids (always delta = 0)
        # ----------------------------------------------------------
        "Cer": {
            "name": "Ceramide",
            "linkage_type": "amide",
            "has_sphingoid": True,
            "acyl_chains": 1,
            "mass_shift": None,
            "can_convert_to": ["SM"],
            "conversion_rules": {
                "SM": {
                    "reaction_requirements": {"external_compounds": ["PC"]},
                    "acyl_chain_change": 0,
                    "require_same_linkage": True,
                    "reaction_type": "sphingolipid"
                }
            }
        },

        "dhCer": {
            "name": "Dihydroceramide",
            "linkage_type": "amide",
            "has_sphingoid": True,
            "acyl_chains": 1,
            "mass_shift": None,
            "can_convert_to": ["Cer","dhSM"],
            "conversion_rules": {
                "Cer": {
                    "reaction_requirements": {"external_compounds": []},
                    "acyl_chain_change": 0,
                    "require_same_linkage": True,
                    "reaction_type": "sphingolipid"
                },
                "dhSM": {
                    "reaction_requirements": {"external_compounds": ["PC"]},
                    "acyl_chain_change": 0,
                    "require_same_linkage": True,
                    "reaction_type": "sphingolipid"
                }
            }
        },

        "SM": {
            "name": "Sphingomyelin",
            "linkage_type": "amide",
            "has_sphingoid": True,
            "acyl_chains": 1,
            "mass_shift": None,
            "can_convert_to": ["Cer"],
            "conversion_rules": {
                "Cer": {
                    "reaction_requirements": {"external_compounds": []},
                    "acyl_chain_change": 0,
                    "require_same_linkage": True,
                    "reaction_type": "sphingolipid"
                }
            }
        },

        # ----------------------------------------------------------
        # PS / LPS
        # ----------------------------------------------------------
        "LPS": {
            "name": "Lyso-PS",
            "linkage_type": "ester",
            "acyl_chains": 1,
            "mass_shift": None,
            "can_convert_to": ["PS"],
            "conversion_rules": {
                "PS": {
                    "reaction_requirements": {"external_compounds": ["acylcoa"]},
                    "acyl_chain_change": 1,
                    "require_same_linkage": True,
                    "reaction_type": "acylation"
                }
            }
        },

        "PS": {
            "name": "Phosphatidylserine",
            "linkage_type": "ester",
            "acyl_chains": 2,
            "mass_shift": 167.0222,
            "can_convert_to": ["PE","PA","LPS"],
            "conversion_rules": {
                "LPS": {
                    "reaction_requirements": {"external_compounds": []},
                    "acyl_chain_change": -1,
                    "require_same_linkage": True,
                    "reaction_type": "deacylation"
                }
            }
        },

        # ----------------------------------------------------------
        # PI / LPI
        # ----------------------------------------------------------
        "LPI": {
            "name": "Lyso-PI",
            "linkage_type": "ester",
            "acyl_chains": 1,
            "mass_shift": None,
            "can_convert_to": ["PI"],
            "conversion_rules": {
                "PI": {
                    "reaction_requirements": {"external_compounds": ["acylcoa"]},
                    "acyl_chain_change": 1,
                    "require_same_linkage": True,
                    "reaction_type": "acylation"
                }
            }
        },

        "PI": {
            "name": "Phosphatidylinositol",
            "linkage_type": "ester",
            "acyl_chains": 2,
            "mass_shift": 242.0192,
            "can_convert_to": ["PA","PIP","LPI"]
        },

        "MG": {
            "name": "Monoacylglycerol",
            "linkage_type": "ester",
            "acyl_chains": 1,
            "mass_shift": None,  # optional, MG mass is tail + glycerol backbone
            "can_convert_to": ["DG"],

            "conversion_rules": {
                "DG": {
                    "reaction_requirements": {"external_compounds": ["acylcoa"]},
                    "acyl_chain_change": 1,
                    "require_same_linkage": True,
                    "reaction_type": "acylation"
                }
            }
        },
        # ----------------------------------------------------------
        # DG / TG
        # ----------------------------------------------------------
        "DG": {
            "name": "Diacylglycerol",
            "linkage_type": "ester",
            "acyl_chains": 2,
            "mass_shift": 17.0027,
            "can_convert_to": ["PA","PC","PE","TG", "MG"],
            "conversion_rules": {          

                # DG → PA (phosphorylation)
                "PA": {
                    "reaction_requirements": {"external_compounds": ["ATP"]},
                    "acyl_chain_change": 0,
                    "require_same_linkage": True,
                    "reaction_type": "phosphorylation"
                },

                # DG → PC (Kennedy pathway)
                "PC": {
                    "reaction_requirements": {"external_compounds": ["CDP-choline"]},
                    "acyl_chain_change": 0,
                    "require_same_linkage": True,
                    "reaction_type": "headgroup_attachment"
                },

                # DG → PE (Kennedy pathway)
                "PE": {
                    "reaction_requirements": {"external_compounds": ["CDP-ethanolamine"]},
                    "acyl_chain_change": 0,
                    "require_same_linkage": True,
                    "reaction_type": "headgroup_attachment"
                },

                # DG → TG (acylation)
                "TG": {
                    "reaction_requirements": {"external_compounds": ["acylcoa"]},
                    "acyl_chain_change": 1,
                    "require_same_linkage": True,
                    "reaction_type": "acylation"
                },
                
                # ✅ DG → MG (deacylation)
                "MG": {
                    "reaction_requirements": {"external_compounds": []},
                    "acyl_chain_change": -1,
                    "require_same_linkage": True,
                    "reaction_type": "deacylation"
                }


            }
        },

        "TG": {
            "name": "Triacylglycerol",
            "linkage_type": "ester",
            "acyl_chains": 3,
            "mass_shift": "VARIABLE_BY_TAIL",
            "can_convert_to": ["DG"],
            "conversion_rules": {
                "DG": {
                    "reaction_requirements": {"external_compounds": []},
                    "acyl_chain_change": -1,
                    "require_same_linkage": True,
                    "reaction_type": "deacylation"
                }
            }
        }
    },

    # ----------------------------------------------------------
    # Fatty acid masses (unchanged from original)
    # ----------------------------------------------------------
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
        "24:1": {"name": "Nervonic", "mass_delta": 348.3392}
    },

    # ----------------------------------------------------------
    # Biochemical deltas (unchanged)
    # ----------------------------------------------------------
    "biochemical_deltas": {
        "methylation_CH2": 14.0156,
        "desaturation_H2_loss": -2.0156,
        "PE_to_PS_delta": 44.0137,
        "PE_to_PC_delta": 42.0470,
        "hydration_H2O": 18.0106
    }
}