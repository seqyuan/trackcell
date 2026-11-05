Logging
=======

TrackCell emits progress and diagnostic messages through the standard
:mod:`logging` module under the ``trackcell`` logger namespace.  Library code
only *emits* records; it never configures handlers, so importing ``trackcell``
is silent and never prints to stdout or stderr on its own (a
:class:`logging.NullHandler` is attached so you never see
"No handlers could be found" warnings).

Reading spatial data can take minutes for large datasets.  TrackCell used to
report progress with ``print()``; as of the logging refactor those messages are
``INFO``-level records that you can silence, capture, or route freely.

Enabling progress output
------------------------

Configure logging explicitly, for example in a notebook or script::

    import logging
    logging.basicConfig(level=logging.INFO)
    logging.getLogger("trackcell").setLevel(logging.INFO)

    import trackcell as tcl
    adata = tcl.io.read_sto_cellbin("B00001.cellbin.gef")
    # INFO progress such as:
    #   [STO cellbin] 12,345 cells [0.4s]
    #   [STO cellbin] Borders: 12,300 valid polygons, 45 fallback [0.9s]

Capturing progress into a list (e.g. for a GUI log pane)::

    import logging

    records: list[str] = []

    class _Collect(logging.Handler):
        def emit(self, record: logging.LogRecord) -> None:
            records.append(record.getMessage())

    logger = logging.getLogger("trackcell")
    logger.addHandler(_Collect())
    logger.setLevel(logging.INFO)

Silencing a noisy reader
------------------------

The default is already silent.  If you enabled ``logging.basicConfig``
globally but want TrackCell quiet, raise its level::

    logging.getLogger("trackcell").setLevel(logging.WARNING)

Non-fatal data problems are warnings
------------------------------------

Issues that do not stop a read but degrade the result — for example a missing
tissue image or ``scalefactors_json.json`` — are reported with
``warnings.warn(..., RuntimeWarning)`` so they integrate with
``pytest.warns`` and ``warnings.filterwarnings``::

    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        adata = tcl.io.read_hd_cellseg(datapath="outs/segmented_outputs")

User-facing report output that is intentionally *displayed* (for example
``tcl.io.atera_info(..., verbose=True)``) is still printed, because it is
opted into by the caller.
