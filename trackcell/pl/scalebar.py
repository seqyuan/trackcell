"""
Physical scale bars for spatial transcriptomics plots.

A scale bar is a short line of known physical length (e.g. 50 µm) drawn on top of
a spatial plot. Its physical length must be *derived from the data*, not hard-coded,
because the same image can be displayed at very different zoom levels: a whole
section (~4 mm) needs a 500 µm bar, whereas a 200 µm zoom-in needs a 10 µm bar.

Coordinates in ``adata.obsm['spatial']`` are not always in micrometers:

* Xenium / Atera / MERSCOPE centroid coordinates are already in µm,
* Visium / Visium HD ``obsm['spatial']`` are in pixels of the full-resolution image,
  so ``scalefactors['microns_per_pixel']`` must be applied.

:func:`resolve_microns_per_unit` inspects the AnnData metadata to decide, and
:func:`add_scale_bar` draws the bar in *axes* coordinates so it stays valid after
cropping, axis inversion (``invert_y``), and tight bounding-box saving.

Examples
--------
>>> import trackcell as tcl
>>> tcl.pl.spatial_cell(adata, color="celltype", scalebar=True)      # auto length
>>> tcl.pl.spatial_cell(adata, color="celltype", scalebar="50um")    # fixed length
>>> tcl.pl.spatial_cell(                                        # doctest: +SKIP
...     adata,
...     color="celltype",
...     scalebar=dict(length="auto", loc="lower left", color="white",
...                   frame=True, fontsize=12),
... )
>>> tcl.pl.spatial_cell(adata, color="celltype", crop_coord=(1000, 1200, 800, 1000),
...                     scalebar=True)                          # zooms in -> shorter bar
"""

from __future__ import annotations

import re
import warnings
from typing import Any, Dict, Iterator, Optional, Tuple, Union, cast

import numpy as np
from matplotlib.artist import Artist
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from matplotlib.transforms import Bbox
from matplotlib.transforms import offset_copy

__all__ = [
    "add_scale_bar",
    "resolve_microns_per_unit",
    "nice_scale_length",
]

#: Unit strings that mean "micrometer".
_UM_ALIASES = {
    "um",
    "µm",
    "μm",
    "micron",
    "microns",
    "micrometer",
    "micrometers",
    "micrometre",
    "micrometres",
}
#: Unit strings that mean "pixel / arbitrary data unit".
#:
#: The empty string is deliberately absent: a *missing* unit is not the same as a
#: declared pixel unit. While ``""`` was a member, ``_parse_length("50")`` parsed
#: the absent suffix as pixels while ``_parse_length(50)`` returned ``None`` and
#: fell back to µm, so ``length=50`` and ``length="50"`` drew bars differing by
#: the µm-per-pixel factor (3x on Visium HD) with nothing but a pair of quotes to
#: explain it.
_PX_ALIASES = {"px", "pixel", "pixels", "data", "unit", "units", "none"}

#: Human readable label suffix for each unit family.
_UNIT_LABEL = {"um": "µm", "px": "px"}

#: Anchor position of the bar inside the axes, as (x, y) fractions.
_LOC_FACTORS = {
    "lower right": (1.0, 0.0),
    "lower left": (0.0, 0.0),
    "lower center": (0.5, 0.0),
    "upper right": (1.0, 1.0),
    "upper left": (0.0, 1.0),
    "upper center": (0.5, 1.0),
    "center right": (1.0, 0.5),
    "center left": (0.0, 0.5),
    "center": (0.5, 0.5),
    "right": (1.0, 0.5),
    "left": (0.0, 0.5),
    "bottom": (0.5, 0.0),
    "top": (0.5, 1.0),
}


# ─────────────────────────────────────────────────────────────────────────────
# Unit resolution
# ─────────────────────────────────────────────────────────────────────────────


def _iter_unit_hints(
    adata: object, library_id: Optional[str] = None
) -> Iterator[Optional[str]]:
    """Yield candidate coordinate-unit strings found in AnnData metadata."""
    uns = getattr(adata, "uns", None)
    if uns is None:
        return

    # Atera readers store a canonical record in uns['atera']
    atera = uns.get("atera") or {}
    if isinstance(atera, dict):
        yield atera.get("coord_unit")
        yield atera.get("spatial_units")

    spatial = uns.get("spatial") or {}
    if library_id is not None:
        spatial = {library_id: spatial.get(library_id, {})} if spatial else spatial
    for info in spatial.values():
        if not isinstance(info, dict):
            continue
        yield info.get("coord_unit")
        metadata = info.get("metadata") or {}
        if isinstance(metadata, dict):
            yield metadata.get("coord_unit")
            yield metadata.get("spatial_units")


def _normalize_unit(unit: object) -> Optional[str]:
    """Map a free-form unit string onto ``'um'`` or ``'px'``; ``None`` if unknown."""
    if unit is None:
        return None
    text = str(unit).strip().lower()
    if text in _UM_ALIASES:
        return "um"
    if text in _PX_ALIASES:
        return "px"
    return None


def _microns_per_pixel(
    adata: object, library_id: Optional[str] = None
) -> Optional[float]:
    """Read ``microns_per_pixel`` from scalefactors, if present and sensible."""
    spatial = getattr(adata, "uns", {}).get("spatial") or {}
    infos = (
        [spatial.get(library_id, {})]
        if library_id is not None
        else list(spatial.values())
    )
    for info in infos:
        if not isinstance(info, dict):
            continue
        for key in ("microns_per_pixel", "micron_per_pixel", "pixel_size"):
            value = (info.get("scalefactors") or {}).get(key)
            if value is None:
                continue
            try:
                value = float(value)
            except (TypeError, ValueError):
                continue
            if np.isfinite(value) and value > 0:
                return value
    return None


def resolve_microns_per_unit(
    adata: object = None,
    library_id: Optional[str] = None,
    *,
    microns_per_unit: Optional[float] = None,
    coord_unit: Optional[str] = None,
) -> Tuple[float, str]:
    """Resolve how many micrometers one unit of ``obsm['spatial']`` represents.

    Parameters
    ----------
    adata
        AnnData with spatial metadata. May be ``None`` when
        ``microns_per_unit`` is given explicitly.
    library_id
        Spatial library key used to look up per-sample metadata.
    microns_per_unit
        Explicit override (µm per data unit). Takes precedence over everything.
    coord_unit
        Explicit hint for the coordinate unit (``"um"`` or ``"px"``).

    Returns
    -------
    (microns_per_unit, unit_family)
        ``microns_per_unit`` is the conversion factor from data units to µm;
        ``unit_family`` is ``"um"`` (coordinates are physical) or ``"px"``
        (coordinates are pixels, factor derived from scalefactors).

    Notes
    -----
    Resolution order:

    1. explicit ``microns_per_unit``;
    2. a ``coord_unit`` / ``spatial_units`` metadata hint (Atera, trackcell readers);
    3. ``scalefactors['microns_per_pixel']`` (Visium / Visium HD pixels);
    4. fallback ``1.0`` — coordinates are assumed to be in µm, which holds for
       Xenium, Atera and MERSCOPE outputs. Pass ``microns_per_unit`` explicitly
       to override if your object uses pixel coordinates without scalefactors.
    """
    if microns_per_unit is not None:
        factor = float(microns_per_unit)
        if not np.isfinite(factor) or factor <= 0:
            raise ValueError(
                f"`microns_per_unit` must be a positive finite number, got {microns_per_unit!r}."
            )
        return factor, "um"

    def _pixel_factor() -> float:
        found = _microns_per_pixel(adata, library_id)
        if found is None:
            raise ValueError(
                "Coordinates were declared as pixels but no "
                "`scalefactors['microns_per_pixel']` was found. Pass "
                "`microns_per_unit` explicitly (e.g. "
                "`scalebar=dict(microns_per_unit=0.2125)`)."
            )
        return found

    hint = _normalize_unit(coord_unit)
    if hint == "px":
        return _pixel_factor(), "px"
    if hint == "um":
        return 1.0, "um"

    if adata is not None:
        declared_px = False
        for raw in _iter_unit_hints(adata, library_id):
            resolved = _normalize_unit(raw)
            if resolved == "um":
                return 1.0, "um"
            if resolved == "px":
                declared_px = True
                break
        known_factor = _microns_per_pixel(adata, library_id)
        if known_factor is not None:
            return known_factor, "px"
        if declared_px:
            return _pixel_factor(), "px"

    return 1.0, "um"


# ─────────────────────────────────────────────────────────────────────────────
# Length selection
# ─────────────────────────────────────────────────────────────────────────────


def nice_scale_length(target: float) -> float:
    """Round ``target`` down to the nearest "nice" length (1, 2 or 5 × 10ⁿ).

    Parameters
    ----------
    target : float
        Desired length in the same unit as the return value.

    Returns
    -------
    float
        The largest value of the form ``{1, 2, 5} × 10ⁿ`` that is ``<= target``.
    """
    target = float(target)
    if not np.isfinite(target) or target <= 0:
        raise ValueError(f"`target` must be a positive finite number, got {target!r}.")
    exponent = int(np.floor(np.log10(target)))
    for multiplier in (5.0, 2.0, 1.0):
        candidate = multiplier * 10.0**exponent
        if candidate <= target * (1 + 1e-9):
            return float(candidate)
    return float(10.0 ** (exponent - 1))


def _format_length(length: float) -> str:
    if length >= 1:
        return f"{length:.0f}" if abs(length - round(length)) < 1e-9 else f"{length:g}"
    for digits in (1, 2, 3):
        text = f"{length:.{digits}f}".rstrip("0").rstrip(".")
        if text not in ("", "0"):
            return text
    return f"{length:g}"


_LENGTH_RE = re.compile(r"^\s*([0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)\s*([^\s0-9]*)\s*$")


def _parse_length(length: object) -> Tuple[Optional[float], Optional[str]]:
    """Parse ``5``, ``"50um"``, ``"100 µm"``, ``"20px"``, ``"auto"``, ``None``."""
    if length is None:
        return None, None
    if isinstance(length, (int, float, np.integer, np.floating)) and not isinstance(
        length, bool
    ):
        return float(length), None
    text = str(length).strip()
    if text.lower() in ("auto", "true", "none", ""):
        return None, None
    match = _LENGTH_RE.match(text)
    if not match:
        raise ValueError(
            f"Could not parse scale bar length {length!r}. Use a number (µm), "
            "a string such as '50um' / '100um' / '20px', or 'auto'."
        )
    value = float(match.group(1))
    raw_unit = match.group(2)
    if raw_unit == "":
        # No suffix given: keep the caller's default unit (µm) rather than
        # guessing pixels. See the note on ``_PX_ALIASES``.
        return value, None
    resolved = _normalize_unit(raw_unit)
    if resolved is None:
        # Falling back to µm here turned ``length="20mm"`` into a bar labelled
        # "20 µm" -- a silent 1000x error that renders perfectly, which is the
        # one failure mode this module exists to avoid.
        raise ValueError(
            f"Unrecognised unit {raw_unit!r} in scale bar length {length!r}. "
            "Use µm/um (also 'micron'/'micrometre'), px/pixel/data, or pass a "
            "bare number plus `unit=` / `microns_per_unit=`."
        )
    return value, resolved


# ─────────────────────────────────────────────────────────────────────────────
# Drawing
# ─────────────────────────────────────────────────────────────────────────────


def _axes_points(ax: Axes) -> Tuple[float, float]:
    """Axes box size in points (width, height).

    ``ax.apply_aspect()`` is called first because ``ax.bbox`` is only final once
    matplotlib has applied the aspect ratio (``set_aspect("equal")``, which every
    spatial plot uses, shrinks the box inside its subplot slot). Measuring before
    that reports a box wider than the one actually drawn, and every clamped
    margin comes out too small — which is how a bar ends up on top of the spine.
    """
    fig = ax.figure
    dpi = fig.dpi if fig is not None else 100.0
    try:
        ax.apply_aspect()
    except Exception:  # noqa: BLE001 - headless axes without a figure manager
        pass
    bbox = ax.bbox
    return bbox.width * 72.0 / dpi, bbox.height * 72.0 / dpi


class _ScaleBarPlacer(Artist):
    """Invisible artist that re-places the scale bar right before every draw."""

    def __init__(self, place, *, zorder: float):
        super().__init__()
        self._place = place
        self.set_zorder(zorder)

    def draw(self, renderer):
        self._place(renderer)

    def get_window_extent(self, renderer=None):  # pragma: no cover - empty artist
        return Bbox.null()


def _text_extent_points(
    ax, text: str, fontsize: float, renderer=None
) -> Tuple[float, float]:
    """Measure the rendered size of ``text`` in points.

    Falls back to a width estimate when no renderer is available (interactive
    backends before the first draw), which keeps the clamp conservative rather
    than letting the label escape the axes.
    """
    fig = ax.figure
    if fig is not None:
        if renderer is None:
            try:
                renderer = fig.canvas.get_renderer()
            except Exception:  # noqa: BLE001 - backend without a renderer yet
                renderer = None
        if renderer is not None:
            dpi = fig.dpi if fig.dpi else 100.0
            probe = fig.text(0, 0, text, fontsize=fontsize)
            try:
                bbox = probe.get_window_extent(renderer=renderer)
                width = bbox.width * 72.0 / dpi
                height = bbox.height * 72.0 / dpi
            finally:
                probe.remove()
            if np.isfinite(width) and np.isfinite(height) and width > 0:
                return float(width), float(height)
    # Rough fallback: average glyph width ~0.62 em, line height ~1.2 em.
    return float(0.62 * fontsize * max(len(text), 1)), float(1.2 * fontsize)


def add_scale_bar(
    ax: Axes,
    adata: object = None,
    *,
    length: Union[float, str, None] = "auto",
    unit: Optional[str] = None,
    microns_per_unit: Optional[float] = None,
    library_id: Optional[str] = None,
    basis: str = "spatial",
    loc: str = "lower right",
    pad: float = 0.03,
    offset: Tuple[float, float] = (0.0, 0.0),
    frac: float = 0.18,
    color: str = "black",
    linewidth: float = 2.5,
    alpha: float = 1.0,
    fontsize: float = 10,
    caps: bool = True,
    cap_size: Optional[float] = None,
    label: Optional[str] = None,
    label_position: str = "auto",
    label_pad: float = 3.0,
    frame: bool = False,
    frame_color: str = "white",
    frame_alpha: float = 0.85,
    frame_pad: float = 4.0,
    zorder: float = 20,
) -> Dict[str, Any]:
    """Draw a physical scale bar on a spatial :class:`matplotlib.axes.Axes`

    The bar is positioned and sized in *axes* coordinates, so it survives axis
    inversion (``invert_y``), cropping, and ``bbox_inches='tight'`` saving. The
    drawn length is converted from µm to data units through
    :func:`resolve_microns_per_unit`.

    Parameters
    ----------
    ax
        Target axes (must already have its final data limits set).
    adata
        AnnData used to resolve the coordinate unit. Optional when
        ``microns_per_unit`` is provided.
    length
        Bar length. ``"auto"``/``None`` picks a "nice" length (1/2/5 × 10ⁿ) equal
        to roughly ``frac`` of the visible width. A number is interpreted in
        ``unit`` (µm by default); strings such as ``"50um"``, ``"100 µm"`` or
        ``"20px"`` are accepted.
    unit
        Unit for a numeric ``length``. ``"um"`` (default) or ``"px"``/``"data"``
        to use raw data units.
    microns_per_unit
        Explicit µm per data unit. Overrides metadata inference.
    library_id
        Spatial library key used for unit inference.
    basis
        Only used for the informative warning when limits are unusable.
    loc
        Anchor of the bar inside the axes: ``"lower right"`` (default),
        ``"lower left"``, ``"upper right"``, ``"upper left"``, ``"lower center"``,
        ``"upper center"``, ``"center"``, ``"center right"``, ``"center left"``,
        ``"right"``, ``"left"``, ``"top"``, ``"bottom"``.
    pad
        Inset from the axes border, as a fraction of the axes size.
    offset
        Extra ``(dx, dy)`` shift in points, applied to the bar and its label.
    frac
        Fraction of the visible width the auto length should occupy
        (default ``0.18``).
    color
        Bar and label color.
    linewidth
        Bar thickness in points (also the cap thickness).
    alpha
        Opacity of bar, caps, and label.
    fontsize
        Label font size in points.
    caps
        Whether to draw end caps (making the bar look like ``|——|``).
    cap_size
        Cap height in points. Defaults to ``4 × linewidth``.
    label
        Custom label text. By default ``"<length> µm"``.
    label_position
        ``"auto"`` (above the bar, below for ``upper*`` anchors), ``"above"``,
        ``"below"`` or ``"none"``.
    label_pad
        Gap between bar and label, in points (default 3).
    frame
        Draw a translucent rounded box behind the bar and label (useful on busy
        H&E backgrounds).
    frame_color, frame_alpha, frame_pad
        Frame face color, opacity, and padding in points.
    zorder
        Matplotlib z-order of the bar (label ``+0.1``, frame ``-0.1``).

    Returns
    -------
    dict
        Resolved properties: ``length`` (in µm), ``label``, ``microns_per_unit``,
        ``unit``, ``axes_fraction`` (bar length as a fraction of the x range).

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> import trackcell as tcl
    >>> fig, ax = plt.subplots()
    >>> ax.scatter(*adata.obsm["spatial"].T, s=1)                # doctest: +SKIP
    >>> ax.set_aspect("equal")                                    # doctest: +SKIP
    >>> tcl.pl.add_scale_bar(ax, adata, length="100um", loc="lower left")
    """
    loc_key = str(loc).strip().lower().replace("_", " ")
    if loc_key not in _LOC_FACTORS:
        raise ValueError(f"Unknown `loc` {loc!r}. Choose from: {sorted(_LOC_FACTORS)}.")
    label_position = str(label_position).strip().lower()
    if label_position not in ("auto", "above", "below", "none"):
        raise ValueError(
            "`label_position` must be 'auto', 'above', 'below' or 'none', "
            f"got {label_position!r}."
        )

    raw_length, unit_from_length = _parse_length(length)
    resolved_unit = _normalize_unit(unit) if unit is not None else None
    if unit is not None and resolved_unit is None:
        # ``loc`` and ``label_position`` are validated above; ``unit`` used to be
        # the one argument where a typo (``"micronz"``) was silently read as µm.
        raise ValueError(
            f"Unrecognised `unit` {unit!r}. Choose from µm/um, micron/micrometre, "
            "or px/pixel/data."
        )
    if unit_from_length is not None and resolved_unit not in (None, unit_from_length):
        raise ValueError(
            f"Conflicting units: `length={length!r}` implies {unit_from_length!r} "
            f"but `unit={unit!r}` was given."
        )
    resolved_unit = unit_from_length or resolved_unit or "um"

    if resolved_unit == "px":
        mpp, unit_family = 1.0, "px"
    else:
        mpp, unit_family = resolve_microns_per_unit(
            adata,
            library_id,
            microns_per_unit=microns_per_unit,
        )
    # The bar length is expressed in µm unless the user asked for raw data units,
    # even when the coordinates themselves are pixels (unit_family == 'px').
    length_unit = "px" if resolved_unit == "px" else "um"

    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    x_span = abs(float(xlim[1]) - float(xlim[0]))
    y_span = abs(float(ylim[1]) - float(ylim[0]))
    if not np.isfinite(x_span) or x_span <= 0:
        raise ValueError(
            "Scale bar needs finite, non-zero axis limits. Set the data limits on "
            "`ax` (e.g. via `spatial_cell(..., crop_coord=...)`) before calling "
            "`add_scale_bar`."
        )
    if not np.isfinite(y_span) or y_span <= 0:
        y_span = x_span

    physical_span = x_span * mpp  # visible width in µm
    if raw_length is None:
        bar_length = nice_scale_length(physical_span * float(frac))
    else:
        bar_length = float(raw_length)
    if not np.isfinite(bar_length) or bar_length <= 0:
        raise ValueError(
            f"Scale bar length must be positive and finite, got {bar_length!r}."
        )

    length_data = bar_length / mpp
    axes_fraction = length_data / x_span
    if axes_fraction > 0.95:
        raise ValueError(
            f"Scale bar of {_format_length(bar_length)} "
            f"{_UNIT_LABEL[length_unit]} spans {axes_fraction:.0%} of the visible "
            "width. Use a shorter `length` or disable the scale bar."
        )
    if axes_fraction < 0.02:
        warnings.warn(
            f"Scale bar of {_format_length(bar_length)} {_UNIT_LABEL[length_unit]} "
            f"covers only {axes_fraction:.1%} of the visible width and may be hard "
            "to read. Consider a longer `length` or `frac`.",
            stacklevel=3,
        )

    if label is None:
        label = f"{_format_length(bar_length)} {_UNIT_LABEL[length_unit]}"

    # ── artists are created once, then *placed* on every draw ───────────────
    #
    # Placement depends on the axes box in points and on the rendered label size.
    # Neither is final when this function returns: ``set_aspect("equal")``,
    # ``tight_layout()``, ``savefig(bbox_inches="tight")`` and window resizes all
    # move the box afterwards, and a bar clamped against the old box ends up on
    # top of the spine. So the geometry is recomputed inside a placer artist that
    # draws just before the bar/label/frame on every render pass.
    default_cap = max(4.0, 4.0 * linewidth)
    bar_pt = max(axes_fraction * _axes_points(ax)[0], 1.0)
    cap_pt = (
        float(cap_size)
        if cap_size is not None
        else float(min(default_cap, max(2.0, 0.6 * bar_pt)))
    )

    trans = offset_copy(
        ax.transAxes,
        fig=cast("Figure", ax.figure),
        x=0.0,
        y=0.0,
        units="points",
    )

    (bar_line,) = ax.plot(
        [0.0, 0.0],
        [0.0, 0.0],
        transform=trans,
        color=color,
        linewidth=linewidth,
        alpha=alpha,
        solid_capstyle="butt",
        clip_on=False,
        zorder=zorder,
        marker="|" if caps else None,
        markersize=cap_pt,
        markeredgewidth=linewidth,
    )

    label_artist = None
    if label and label_position != "none":
        label_artist = ax.annotate(
            label,
            xy=(0.5, 0.0),
            xycoords=trans,
            xytext=(0.0, 0.0),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=fontsize,
            color=color,
            alpha=alpha,
            clip_on=False,
            zorder=zorder + 0.1,
        )

    frame_patch = None
    if frame:
        frame_patch = Rectangle(
            (0.0, 0.0),
            0.0,
            0.0,
            transform=trans,
            facecolor=frame_color,
            alpha=frame_alpha,
            edgecolor="none",
            clip_on=False,
            zorder=zorder - 0.1,
        )
        ax.add_patch(frame_patch)

    warned = {"fit": False}

    def _place(renderer=None):
        """Compute the clamped geometry in axes fractions and apply it."""
        ax_w_pt, ax_h_pt = _axes_points(ax)
        ax_w_pt = ax_w_pt if ax_w_pt > 0 else 1.0
        ax_h_pt = ax_h_pt if ax_h_pt > 0 else 1.0

        # Caps are ``|`` markers: length ``cap_pt`` vertically, ``linewidth``
        # horizontally, so they bound the group on both axes.
        stroke_half = linewidth / 2.0
        breathing = 1.0  # points of margin so the stroke never sits on the spine
        frame_extra = float(frame_pad) if frame else 0.0
        margin_x_pt = stroke_half + breathing + frame_extra
        margin_y_pt = max(stroke_half, cap_pt / 2.0) + breathing + frame_extra

        # An offset larger than the axes can never keep the group inside, so it
        # is clamped first: a requested shift is a hint, staying in frame is the rule.
        max_off_x = max(ax_w_pt / 2.0 - margin_x_pt, 0.0)
        max_off_y = max(ax_h_pt / 2.0 - margin_y_pt, 0.0)
        off_x = float(np.clip(float(offset[0]), -max_off_x, max_off_x))
        off_y = float(np.clip(float(offset[1]), -max_off_y, max_off_y))

        x_lo = (margin_x_pt - off_x) / ax_w_pt
        x_hi = 1.0 - (margin_x_pt + off_x) / ax_w_pt
        y_lo = (margin_y_pt - off_y) / ax_h_pt
        y_hi = 1.0 - (margin_y_pt + off_y) / ax_h_pt
        if x_hi < x_lo:  # request larger than the box: fall back to the padding
            x_lo, x_hi = pad, 1.0 - pad
        if y_hi < y_lo:
            y_lo, y_hi = pad, 1.0 - pad

        placement = offset_copy(
            ax.transAxes,
            fig=cast("Figure", ax.figure),
            x=off_x,
            y=off_y,
            units="points",
        )
        bar_line.set_transform(placement)
        if label_artist is not None:
            label_artist.xycoords = placement
        if frame_patch is not None:
            frame_patch.set_transform(placement)

        x_anchor, y_anchor = _LOC_FACTORS[loc_key]
        if y_anchor <= 0.0:
            y_bar = y_lo
        elif y_anchor >= 1.0:
            y_bar = y_hi
        else:
            y_bar = float(np.clip(y_anchor, y_lo, y_hi))

        if axes_fraction <= max(x_hi - x_lo, 0.0):
            if x_anchor >= 1.0:
                x_end = x_hi
                x_start = x_end - axes_fraction
            elif x_anchor <= 0.0:
                x_start = x_lo
                x_end = x_start + axes_fraction
            else:
                center = float(
                    np.clip(
                        x_anchor, x_lo + axes_fraction / 2.0, x_hi - axes_fraction / 2.0
                    )
                )
                x_start = center - axes_fraction / 2.0
                x_end = center + axes_fraction / 2.0
        else:
            # Shortening the bar would misreport the scale, so keep the true
            # length and say so instead of silently drawing a wrong one. Warn
            # once, at call time: repeating this on every draw would be noise.
            if not warned["fit"] and renderer is None:
                warned["fit"] = True
                warnings.warn(
                    f"A {_format_length(bar_length)} {_UNIT_LABEL[length_unit]} "
                    "scale bar does not fit inside the axes together with its caps "
                    "and label; it will overflow. Use a shorter `length`, a smaller "
                    "`frac`, or `caps=False`/`frame=False`.",
                    stacklevel=3,
                )
            center = (x_lo + x_hi) / 2.0
            x_start = center - axes_fraction / 2.0
            x_end = center + axes_fraction / 2.0

        bar_line.set_xdata([x_start, x_end])
        bar_line.set_ydata([y_bar, y_bar])

        position = label_position
        x_center = (x_start + x_end) / 2.0
        ha = "center"
        text_w_pt = text_h_pt = 0.0
        if label_artist is not None and position != "none":
            text_w_pt, text_h_pt = _text_extent_points(ax, label, fontsize, renderer)
            if position == "auto":
                position = "below" if y_anchor >= 1.0 else "above"
            band = (float(label_pad) + text_h_pt + margin_y_pt - breathing) / ax_h_pt
            if position == "above" and y_bar + band > y_hi:
                position = "below"
            elif position == "below" and y_bar - band < y_lo:
                position = "above"

            text_w = text_w_pt / ax_w_pt
            if text_w <= max(axes_fraction, 1e-9):
                x_center = float(
                    np.clip(x_center, x_lo + text_w / 2.0, x_hi - text_w / 2.0)
                )
            else:
                # Wider label than the bar: keep it attached to the bar's inner
                # end instead of floating in the middle of the axes.
                if x_anchor >= 1.0:
                    x_center, ha = x_end, "right"
                elif x_anchor <= 0.0:
                    x_center, ha = x_start, "left"
                if ha == "right":
                    x_center = max(x_center, x_lo + text_w)
                elif ha == "left":
                    x_center = min(x_center, x_hi - text_w)
                else:
                    x_center = float(
                        np.clip(x_center, x_lo + text_w / 2.0, x_hi - text_w / 2.0)
                    )

            sign = 1.0 if position == "above" else -1.0
            label_artist.xy = (x_center, y_bar)
            label_artist.xyann = (0.0, sign * float(label_pad))
            label_artist.set_ha(ha)
            label_artist.set_va("bottom" if sign > 0 else "top")
        elif label_artist is not None:
            label_artist.xy = (x_center, y_bar)

        if frame_patch is not None:
            pad_x = float(frame_pad) / ax_w_pt
            pad_y = float(frame_pad) / ax_h_pt
            if label_artist is not None and position != "none":
                text_extra = (float(label_pad) + text_h_pt) / ax_h_pt
            else:
                text_extra = 0.0
            if position == "below":
                y_low = y_bar - text_extra - pad_y
                y_high = y_bar + pad_y
            else:
                y_low = y_bar - pad_y
                y_high = y_bar + text_extra + pad_y
            frame_patch.set_xy((min(x_start, x_end) - pad_x, y_low))
            frame_patch.set_width(abs(x_end - x_start) + 2 * pad_x)
            frame_patch.set_height(max(y_high - y_low, 1e-6))

    placer = _ScaleBarPlacer(_place, zorder=zorder - 0.2)
    ax.add_artist(placer)
    _place(None)  # so the artists are sane before the first draw

    return {
        "length": bar_length,
        "label": label,
        "microns_per_unit": mpp,
        "unit": unit_family,
        "length_unit": length_unit,
        "axes_fraction": axes_fraction,
    }


def _maybe_add_scale_bar(
    ax: Axes, scalebar: object, adata: object, library_id: Optional[str] = None
) -> Optional[Dict[str, Any]]:
    """Interpret the ``scalebar=`` argument of the plotting functions.

    ``scalebar`` may be ``None``/``False`` (skip), ``True`` (auto length), a
    number, a string such as ``"50um"``, or a dict of :func:`add_scale_bar`
    keyword arguments.
    """
    if scalebar is None or scalebar is False:
        return None
    if isinstance(scalebar, dict):
        options: Dict[str, Any] = dict(scalebar)
    else:
        options = {"length": scalebar}
    options.setdefault("library_id", library_id)
    return add_scale_bar(ax, adata, **options)
