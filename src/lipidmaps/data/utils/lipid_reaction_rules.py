import re
from typing import Any, Dict, List, Optional

lipid_reaction_rules: Dict[str, Any] = {
    "metadata": {
        "version": "0.1",
        "description": "Rules for lipid headgroup remodeling and fatty acid mass shifts for quantitation.",
        "units": "Daltons (Da)",
        "usage_note": "Conserve acyl tails (Total C:DB) during headgroup swaps.",
    },
    "headgroup_reactions": {
      "PC": {
        "name": "Phosphatidylcholine",
        "linkage_type": "ester",
        "acyl_chains": 2,
        "mass_shift": 165.0555,
        "can_convert_to": ["PA", "LPC", "DAG", "SM"],
        "conversion_rules": {
          "LPC": {"required_compound": "fa", "required_acyl_chains": 1, "is_molspecies": False, "require_same_linkage": True}
        }
      },
      "LPC": {
        "name": "Monoacylglycerophosphocholine",
        "linkage_type": "ester",
        "acyl_chains": 1,
        "mass_shift": None,
        "can_convert_to": ["PC"],
        "conversion_rules": {
          "PC": {"required_compound": "facoa", "required_acyl_chains": 2, "is_molspecies": False, "require_same_linkage": True}
        }
      },
      "PC O-": {
        "name": "Monoalkylglycerophosphocholine (ether)",
        "linkage_type": "ether_alkyl",
        "acyl_chains": 1,
        "mass_shift": None,
        "can_convert_to": ["LPC O-"],
        "conversion_rules": {
          "LPC O-": {"required_compound": "fa", "required_acyl_chains": 1, "is_molspecies": False, "require_same_linkage": True}
        }
      },
      "LPC O-": {
        "name": "Monoalkyl-lyso-phosphocholine (ether)",
        "linkage_type": "ether_alkyl",
        "acyl_chains": 1,
        "mass_shift": None,
        "can_convert_to": ["PC O-"],
        "conversion_rules": {
          "PC O-": {"required_compound": "facoa", "required_acyl_chains": 2, "is_molspecies": False, "require_same_linkage": True}
        }
      },
      "PC P-": {
        "name": "Plasmenyl-phosphocholine (plasmalogen)",
        "linkage_type": "ether_vinyl",
        "acyl_chains": 1,
        "mass_shift": None,
        "can_convert_to": ["LPC P-"],
        "conversion_rules": {
          "LPC P-": {"required_compound": "fa", "required_acyl_chains": 1, "is_molspecies": False, "require_same_linkage": True}
        }
      },
      "LPC P-": {
        "name": "Plasmenyl-lyso-phosphocholine (plasmalogen)",
        "linkage_type": "ether_vinyl",
        "acyl_chains": 1,
        "mass_shift": None,
        "can_convert_to": ["PC P-"],
        "conversion_rules": {
          "PC P-": {"required_compound": "facoa", "required_acyl_chains": 2, "is_molspecies": False, "require_same_linkage": True}
        }
      },
      "PE": {
        "name": "Phosphatidylethanolamine",
        "linkage_type": "ester",
        "acyl_chains": 2,
        "mass_shift": 123.0085,
        "can_convert_to": ["PA", "PC", "PS", "LPE"],
        "conversion_rules": {
          "LPE": {"required_compound": "fa", "required_acyl_chains": 1, "is_molspecies": False, "require_same_linkage": True}
        }
      },
      "LPE": {
        "name": "Monoacylglycerophosphoethanolamine",
        "linkage_type": "ester",
        "acyl_chains": 1,
        "mass_shift": None,
        "can_convert_to": ["PE"],
        "conversion_rules": {
          "PE": {"required_compound": "facoa", "required_acyl_chains": 2, "is_molspecies": False, "require_same_linkage": True}
        }
      },
      "PI": {
        "name": "Phosphatidylinositol",
        "acyl_chains": 2,
        "mass_shift": 242.0192,
        "can_convert_to": ["PA", "PIP", "PIP2", "PIP3"]
      },
      "PE O-": {
        "name": "Monoalkylglycerophosphoethanolamine (ether)",
        "linkage_type": "ether_alkyl",
        "acyl_chains": 1,
        "mass_shift": None,
        "can_convert_to": ["LPE O-"],
        "conversion_rules": {
          "LPE O-": {"required_compound": "fa", "required_acyl_chains": 1, "is_molspecies": False, "require_same_linkage": True}
        }
      },
      "LPE O-": {
        "name": "Monoalkyl-lyso-phosphoethanolamine (ether)",
        "linkage_type": "ether_alkyl",
        "acyl_chains": 1,
        "mass_shift": None,
        "can_convert_to": ["PE O-"],
        "conversion_rules": {
          "PE O-": {"required_compound": "facoa", "required_acyl_chains": 2, "is_molspecies": False, "require_same_linkage": True}
        }
      },
      "PE P-": {
        "name": "Plasmenyl-phosphoethanolamine (plasmalogen)",
        "linkage_type": "ether_vinyl",
        "acyl_chains": 1,
        "mass_shift": None,
        "can_convert_to": ["LPE P-"],
        "conversion_rules": {
          "LPE P-": {"required_compound": "fa", "required_acyl_chains": 1, "is_molspecies": False, "require_same_linkage": True}
        }
      },
      "LPE P-": {
        "name": "Plasmenyl-lyso-phosphoethanolamine (plasmalogen)",
        "linkage_type": "ether_vinyl",
        "acyl_chains": 1,
        "mass_shift": None,
        "can_convert_to": ["PE P-"],
        "conversion_rules": {
          "PE P-": {"required_compound": "facoa", "required_acyl_chains": 2, "is_molspecies": False, "require_same_linkage": True}
        }
      },
      "PG": {
        "name": "Phosphatidylglycerol",
        "linkage_type": "ester",
        "acyl_chains": 2,
        "mass_shift": 154.0031,
        "can_convert_to": ["CL", "PA"]
      },
      "LPA": {
        "name": "Lysophosphatidic acid",
        "linkage_type": "ester",
        "acyl_chains": 1,
        "mass_shift": None,
        "can_convert_to": ["PA"],
        "conversion_rules": {
          "PA": {"required_compound": "facoa", "required_acyl_chains": 2, "is_molspecies": False, "require_same_linkage": True}
        }
      },
      "PA": {
        "name": "Phosphatidic Acid (augment)",
        "linkage_type": "ester",
        "acyl_chains": 2,
        "mass_shift": 79.9663,
        "can_convert_to": ["PC", "PE", "PS", "PI", "PG", "DAG"]
      },
      "Cer": {
        "name": "Ceramide",
        "linkage_type": "amide",
        "has_sphingoid": True,
        "sphingoid_chains": 1,
        "acyl_chains": 1,
        "mass_shift": None,
        "can_convert_to": ["SM"],
        "conversion_rules": {
          "SM": {"required_compound": None, "required_acyl_chains": 1, "is_molspecies": False, "require_same_linkage": True}
        }
      },
      "dhCer": {
        "name": "Dihydroceramide",
        "linkage_type": "amide",
        "has_sphingoid": True,
        "sphingoid_chains": 1,
        "acyl_chains": 1,
        "mass_shift": None,
        "can_convert_to": ["Cer", "dhSM"],
        "conversion_rules": {
          "Cer": {"required_compound": None, "required_acyl_chains": 1, "is_molspecies": False, "require_same_linkage": True},
          "dhSM": {"required_compound": None, "required_acyl_chains": 1, "is_molspecies": False, "require_same_linkage": True}
        }
      },
      "SM": {
        "name": "Sphingomyelin",
        "linkage_type": "amide",
        "has_sphingoid": True,
        "sphingoid_chains": 1,
        "acyl_chains": 1,
        "mass_shift": None,
        "can_convert_to": ["Cer"],
        "conversion_rules": {
          "Cer": {"required_compound": None, "required_acyl_chains": 1, "is_molspecies": False, "require_same_linkage": True}
        }
      },
      "LPS": {
        "name": "Lyso-Phosphatidylserine",
        "linkage_type": "ester",
        "acyl_chains": 1,
        "mass_shift": None,
        "can_convert_to": ["PS"],
        "conversion_rules": {
          "PS": {"required_compound": "facoa", "required_acyl_chains": 2, "is_molspecies": False, "require_same_linkage": True}
        }
      },
      "PS": {
        "name": "Phosphatidylserine (augment)",
        "linkage_type": "ester",
        "acyl_chains": 2,
        "mass_shift": 167.0222,
        "can_convert_to": ["PE", "PA", "LPS"],
        "conversion_rules": {
          "LPS": {"required_compound": "fa", "required_acyl_chains": 1, "is_molspecies": False, "require_same_linkage": True}
        }
      },
      "LPI": {
        "name": "Lyso-Phosphatidylinositol",
        "linkage_type": "ester",
        "acyl_chains": 1,
        "mass_shift": None,
        "can_convert_to": ["PI"],
        "conversion_rules": {
          "PI": {"required_compound": "facoa", "required_acyl_chains": 2, "is_molspecies": False, "require_same_linkage": True}
        }
      },
      "PI": {
        "name": "Phosphatidylinositol (augment)",
        "linkage_type": "ester",
        "acyl_chains": 2,
        "mass_shift": 242.0192,
        "can_convert_to": ["PA", "PIP", "LPI"],
      },
      "DAG": {
        "name": "Diacylglycerol (augment)",
        "linkage_type": "ester",
        "mass_shift": 17.0027,
        "can_convert_to": ["PA", "PC", "PE", "TAG"],
        "conversion_rules": {
          "TAG": {"required_compound": "facoa", "required_acyl_chains": 3, "is_molspecies": False, "require_same_linkage": True}
        }
      },
      "DAG": {
        "name": "Diacylglycerol",
        "linkage_type": "ester",
        "mass_shift": 17.0027,
        "can_convert_to": ["PA", "PC", "PE", "TAG"]
      },
      "TAG": {
        "name": "Triacylglycerol",
        "linkage_type": "ester",
        "mass_shift": "VARIABLE_BY_TAIL",
        "can_convert_to": ["DAG"]
      }
  },
  
  "fatty_acids": {
    "12:0": { "name": "Lauric", "mass_delta": 182.1671 },
    "14:0": { "name": "Myristic", "mass_delta": 210.1984 },
    "14:1": { "name": "Myristoleic", "mass_delta": 208.1827 },
    "16:0": { "name": "Palmitic", "mass_delta": 238.2297 },
    "16:1": { "name": "Palmitoleic", "mass_delta": 236.2140 },
    "18:0": { "name": "Stearic", "mass_delta": 266.2610 },
    "18:1": { "name": "Oleic", "mass_delta": 264.2453 },
    "18:2": { "name": "Linoleic", "mass_delta": 262.2297 },
    "18:3": { "name": "Linolenic", "mass_delta": 260.2140 },
    "20:0": { "name": "Arachidic", "mass_delta": 294.2923 },
    "20:1": { "name": "Gadoleic", "mass_delta": 292.2766 },
    "20:2": { "name": "Eicosadienoic", "mass_delta": 290.2610 },
    "20:3": { "name": "DGLA", "mass_delta": 288.2453 },
    "20:4": { "name": "Arachidonic", "mass_delta": 286.2297 },
    "20:5": { "name": "EPA", "mass_delta": 284.2140 },
    "22:0": { "name": "Behenic", "mass_delta": 322.3236 },
    "22:1": { "name": "Erucic", "mass_delta": 320.3079 },
    "22:4": { "name": "Adrenic", "mass_delta": 314.2610 },
    "22:5": { "name": "DPA", "mass_delta": 312.2453 },
    "22:6": { "name": "DHA", "mass_delta": 310.2297 },
    "24:0": { "name": "Lignoceric", "mass_delta": 350.3549 },
    "24:1": { "name": "Nervonic", "mass_delta": 348.3392 }
  },
  "biochemical_deltas": {
    "methylation_CH2": 14.0156,
    "desaturation_H2_loss": -2.0156,
    "PE_to_PS_delta": 44.0137,
    "PE_to_PC_delta": 42.0470,
    "hydration_H2O": 18.0106
  }
}
