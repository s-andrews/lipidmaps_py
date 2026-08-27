import os

# Base URL for the ReactionChecker / LMSD reactions service.
# Can be overridden with the environment variable `LMSD_REACTIONS_BASE_URL`.
LMSD_REACTIONS_BASE_URL = os.getenv(
	"LMSD_REACTIONS_BASE_URL", "https://dev.lipidmaps.org"
)

# User-Agent sent on all outbound HTTP requests. Identifies this package to
# remote services; some hosts block the default `python-requests` UA via bot
# protection, so we always send an explicit identifier.
USER_AGENT = "lipidmaps_py"
DEFAULT_HEADERS = {"User-Agent": USER_AGENT}
