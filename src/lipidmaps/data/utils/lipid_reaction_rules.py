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
            "can_convert_to": ["PA","PS","LPC","DAG","SM"],
            "conversion_rules": {
                "LPC": {
                    "reaction_requirements": {"external_compounds": []},
                    "acyl_chain_change": -1,
                    "require_same_linkage": True,
                    "reaction_type": "deacylation"
                },
                "PS": {
                    "reaction_requirements": {"external_compounds": ["serine"]},
                    "acyl_chain_change": 0,
                    "require_same_linkage": True,
                    "reaction_type": "headgroup_swap"
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
            "can_convert_to": ["LPC O-"],
            "conversion_rules": {
                "LPC O-": {
                    "reaction_requirements": {"external_compounds": []},
                    "acyl_chain_change": 0,
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
            "can_convert_to": ["PA","PC","PS","LPE"],
            "conversion_rules": {
                "LPE": {
                    "reaction_requirements": {"external_compounds": []},
                    "acyl_chain_change": -1,
                    "require_same_linkage": True,
                    "reaction_type": "deacylation"
                }
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
            "can_convert_to": ["CL","PA"]
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
            "can_convert_to": ["PC","PE","PS","PI","PG","DAG"]
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

        # ----------------------------------------------------------
        # DAG / TAG
        # ----------------------------------------------------------
        "DAG": {
            "name": "Diacylglycerol",
            "linkage_type": "ester",
            "acyl_chains": 2,
            "mass_shift": 17.0027,
            "can_convert_to": ["PA","PC","PE","TAG"],
            "conversion_rules": {
                "TAG": {
                    "reaction_requirements": {"external_compounds": ["acylcoa"]},
                    "acyl_chain_change": 1,
                    "require_same_linkage": True,
                    "reaction_type": "acylation"
                }
            }
        },

        "TAG": {
            "name": "Triacylglycerol",
            "linkage_type": "ester",
            "acyl_chains": 3,
            "mass_shift": "VARIABLE_BY_TAIL",
            "can_convert_to": ["DAG"],
            "conversion_rules": {
                "DAG": {
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