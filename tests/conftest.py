from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "src"

if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

# Needed so tests can import the repo-level `scripts` package under
# --import-mode=importlib (which does not add the rootdir to sys.path).
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(1, str(REPO_ROOT))