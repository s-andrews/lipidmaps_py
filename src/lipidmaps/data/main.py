import os
import logging

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[logging.FileHandler("lipidmaps_py.log"), logging.StreamHandler()],
)


def main():
    # Basic functionality for CLI
    logger = logging.getLogger(__name__)
    logger.info("LIPID MAPS Python API suite initialized.")


if __name__ == "__main__":
    main()
