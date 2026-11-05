"""Central logging bootstrap for trackcell.

Library code only *emits* records; it never configures handlers.  Users opt in::

    import logging
    logging.basicConfig(level=logging.INFO)
    logging.getLogger("trackcell").setLevel(logging.INFO)

Importing :mod:`trackcell` attaches a :class:`logging.NullHandler` to the
``trackcell`` logger, so library consumers never see
"No handlers could be found" warnings, and progress output stays silent until
they explicitly configure logging.
"""

from __future__ import annotations

import logging

LOGGER_NAME = "trackcell"

# Attach a NullHandler once, so importing trackcell never triggers
# "No handlers could be found" warnings for library consumers.
_root = logging.getLogger(LOGGER_NAME)
if not any(isinstance(h, logging.NullHandler) for h in _root.handlers):
    _root.addHandler(logging.NullHandler())


def get_logger(name: str) -> logging.Logger:
    """Return a child logger under the ``trackcell`` namespace.

    ``name`` may be a bare module name (``"read_sto"``) or a fully-qualified
    module name (``"trackcell.io.read_sto"``); both resolve under ``trackcell``.
    """
    if name == LOGGER_NAME or name.startswith(LOGGER_NAME + "."):
        return logging.getLogger(name)
    return logging.getLogger(f"{LOGGER_NAME}.{name}")
