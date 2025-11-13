"""
Configuration settings for BioPAN tool
UI should allow users to modify these settings as needed.
"""
LEVEL = "MS" # Class, Species, Molecular Species
LOG_FILE = "biopan.log"
MAX_RETRIES = 3
TIMEOUT = 30  # seconds
CONDITIONS = {
    "refmet": True,
    "lipidlynxx": False,
    "conditions": ["control", "condition_1", "condition_2", "condition_3"],
    "labels": ["replicate_1", "replicate_2", "replicate_3"]
}