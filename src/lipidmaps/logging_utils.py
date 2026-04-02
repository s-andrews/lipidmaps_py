from __future__ import annotations

import logging
import os
from logging.handlers import RotatingFileHandler
from pathlib import Path


LOG_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
DEFAULT_APP_NAME = "lipidmaps_py"
_HANDLER_NAMES = {
    "info": f"{DEFAULT_APP_NAME}.info_file",
    "error": f"{DEFAULT_APP_NAME}.error_file",
    "stream": f"{DEFAULT_APP_NAME}.stream",
}


def _resolve_log_dir(log_dir: str | Path | None = None) -> Path:
    if log_dir is not None:
        return Path(log_dir).expanduser().resolve()

    env_log_dir = os.getenv("LIPIDMAPS_LOG_DIR")
    if env_log_dir:
        return Path(env_log_dir).expanduser().resolve()

    package_root = Path(__file__).resolve().parents[2]
    if (package_root / "setup.py").exists() and (package_root / "src").exists():
        return package_root / "logs"

    return (Path.cwd() / "logs").resolve()


def get_log_paths(log_dir: str | Path | None = None, app_name: str = DEFAULT_APP_NAME) -> dict[str, Path]:
    resolved_dir = _resolve_log_dir(log_dir)
    return {
        "directory": resolved_dir,
        "info": resolved_dir / f"{app_name}.log",
        "error": resolved_dir / f"{app_name}.error.log",
    }


def _build_handler(path: Path, level: int, name: str) -> RotatingFileHandler:
    handler = RotatingFileHandler(path, maxBytes=5 * 1024 * 1024, backupCount=5)
    handler.setLevel(level)
    handler.set_name(name)
    handler.setFormatter(logging.Formatter(LOG_FORMAT))
    return handler


def _remove_managed_handlers(root_logger: logging.Logger) -> None:
    managed_names = set(_HANDLER_NAMES.values())
    for handler in list(root_logger.handlers):
        if handler.get_name() in managed_names:
            root_logger.removeHandler(handler)
            handler.close()


def configure_logging(
    log_dir: str | Path | None = None,
    *,
    app_name: str = DEFAULT_APP_NAME,
    level: int = logging.INFO,
    console: bool = True,
    force: bool = False,
) -> Path:
    paths = get_log_paths(log_dir=log_dir, app_name=app_name)
    paths["directory"].mkdir(parents=True, exist_ok=True)

    root_logger = logging.getLogger()
    root_logger.setLevel(level)

    if force:
        _remove_managed_handlers(root_logger)

    existing_names = {handler.get_name() for handler in root_logger.handlers}

    if _HANDLER_NAMES["info"] not in existing_names:
        root_logger.addHandler(_build_handler(paths["info"], level, _HANDLER_NAMES["info"]))

    if _HANDLER_NAMES["error"] not in existing_names:
        root_logger.addHandler(
            _build_handler(paths["error"], logging.ERROR, _HANDLER_NAMES["error"])
        )

    if console and _HANDLER_NAMES["stream"] not in existing_names:
        stream_handler = logging.StreamHandler()
        stream_handler.setLevel(level)
        stream_handler.set_name(_HANDLER_NAMES["stream"])
        stream_handler.setFormatter(logging.Formatter(LOG_FORMAT))
        root_logger.addHandler(stream_handler)

    logging.captureWarnings(True)
    return paths["directory"]