import logging
import sys
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lipidmaps.logging_utils import configure_logging, get_log_paths


def test_configure_logging_creates_log_files(tmp_path):
    log_dir = tmp_path / "logs"
    resolved_dir = configure_logging(log_dir=log_dir, force=True)
    logger = logging.getLogger("lipidmaps.tests.logging")

    logger.info("info message from test")
    logger.error("error message from test")

    for handler in logging.getLogger().handlers:
        flush = getattr(handler, "flush", None)
        if callable(flush):
            flush()

    paths = get_log_paths(log_dir=log_dir)
    assert resolved_dir == log_dir.resolve()
    assert paths["info"].exists()
    assert paths["error"].exists()
    assert "info message from test" in paths["info"].read_text()
    assert "error message from test" in paths["error"].read_text()


def test_configure_logging_does_not_duplicate_handlers(tmp_path):
    log_dir = tmp_path / "logs"
    configure_logging(log_dir=log_dir, force=True)
    root_logger = logging.getLogger()
    managed_before = sorted(
        handler.get_name()
        for handler in root_logger.handlers
        if handler.get_name() and handler.get_name().startswith("lipidmaps_py.")
    )

    configure_logging(log_dir=log_dir)
    managed_after = sorted(
        handler.get_name()
        for handler in root_logger.handlers
        if handler.get_name() and handler.get_name().startswith("lipidmaps_py.")
    )

    assert managed_after == managed_before