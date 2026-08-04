"""Shared plot dimensions and memory-safe headless rendering helpers."""
from __future__ import annotations

import gc
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable


@dataclass(frozen=True)
class PlotProfile:
    """Pixel-exact output profile used by static scientific figures."""

    name: str
    width_px: int
    height_px: int
    dpi: int

    @property
    def figsize(self) -> tuple[float, float]:
        return self.width_px / self.dpi, self.height_px / self.dpi


SCREEN_PROFILE = PlotProfile("screen", 2480, 3508, 180)
REPORT_PROFILE = PlotProfile("report", 2480, 3508, 600)
PLOT_PROFILES = {
    SCREEN_PROFILE.name: SCREEN_PROFILE,
    REPORT_PROFILE.name: REPORT_PROFILE,
}


@dataclass
class RenderJob:
    """One isolated rendering request.

    ``build_plot`` must return ``(figure, axes)``. ``axes`` may be a single
    Matplotlib axis or an iterable of axes whose labels need normalization.
    """

    output: Path
    build_plot: Callable[[tuple[float, float], int], tuple[object, object]]
    x_label_count: int = 0
    y_label_count: int = 0


def adaptive_label_size(item_count: int, *, minimum: float = 4.0, maximum: float = 9.0) -> float:
    """Return a readable tick size without hiding scientific categories."""
    count = max(1, int(item_count or 0))
    return max(minimum, min(maximum, 11.0 - (count ** 0.5) * 0.55))


def configure_heatmap_axes(ax, x_count: int = 0, y_count: int = 0) -> None:
    """Apply the common heatmap label policy required by both profiles."""
    from matplotlib import pyplot as plt

    ax.tick_params(
        axis="x",
        labelsize=adaptive_label_size(x_count),
        pad=2,
    )
    ax.tick_params(
        axis="y",
        labelsize=adaptive_label_size(y_count),
        pad=2,
    )
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")


def enable_constrained_layout(fig) -> None:
    """Enable constrained layout across the supported Matplotlib versions."""
    import matplotlib

    parts = str(matplotlib.__version__).split(".")
    try:
        version = int(parts[0]), int(parts[1])
    except (IndexError, ValueError):
        version = (0, 0)
    if version >= (3, 6) and hasattr(fig, "set_layout_engine"):
        fig.set_layout_engine("constrained")
    else:
        fig.set_constrained_layout(True)


def _axes_list(value: object) -> list[object]:
    if value is None:
        return []
    if hasattr(value, "tick_params"):
        return [value]
    if isinstance(value, Iterable):
        return [item for item in value if hasattr(item, "tick_params")]
    return []


def render_plot(job: RenderJob, profile: PlotProfile) -> Path:
    """Render one plot with Agg and release every Matplotlib allocation."""
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    output = Path(job.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig = None
    try:
        fig, axes = job.build_plot(profile.figsize, profile.dpi)
        for ax in _axes_list(axes):
            configure_heatmap_axes(ax, job.x_label_count, job.y_label_count)
        try:
            enable_constrained_layout(fig)
        except (AttributeError, ValueError):
            pass
        # Matplotlib truncates an exact floating-point product on some older
        # releases (3508 / 180 * 180 becomes 3507.999...). The subpixel guard
        # preserves the requested physical size while keeping the PNG exact.
        guarded_size = (
            (profile.width_px + 0.01) / profile.dpi,
            (profile.height_px + 0.01) / profile.dpi,
        )
        fig.set_size_inches(*guarded_size, forward=False)
        fig.savefig(output, dpi=profile.dpi, facecolor=fig.get_facecolor())
        return output
    finally:
        if fig is not None:
            fig.clear()
            plt.close(fig)
        plt.close("all")
        gc.collect()
