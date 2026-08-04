from __future__ import annotations

import tempfile
from pathlib import Path

from PIL import Image

from luna_gui.core.plot_profiles import (
    REPORT_PROFILE,
    SCREEN_PROFILE,
    PlotProfile,
    RenderJob,
    reference_tick_indices,
    render_plot,
)


def _builder(figsize, dpi):
    from matplotlib import pyplot as plt

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi, constrained_layout=True)
    ax.imshow([[0, 1], [1, 0]], aspect="auto")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["long label A", "long label B"])
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["row A", "row B"])
    return fig, ax


def test_screen_and_report_profiles_are_pixel_exact_and_release_figures() -> None:
    from matplotlib import pyplot as plt

    with tempfile.TemporaryDirectory() as tmp:
        for profile in (SCREEN_PROFILE, REPORT_PROFILE):
            output = Path(tmp) / f"{profile.name}.png"
            result = render_plot(
                RenderJob(output, _builder, x_label_count=2, y_label_count=2),
                profile,
            )
            with Image.open(result) as image:
                assert image.size == (2480, 3508)
            assert plt.get_fignums() == []
    assert SCREEN_PROFILE.dpi == 180
    assert REPORT_PROFILE.dpi == 300


def test_heatmap_reference_ticks_use_nine_frames_or_two_values() -> None:
    trajectory = reference_tick_indices(101, trajectory=True)
    values = reference_tick_indices(101, trajectory=False)

    assert len(trajectory) == 9
    assert trajectory[0] == 0
    assert trajectory[4] == 50
    assert trajectory[-1] == 100
    assert all(left + right == 100 for left, right in zip(trajectory, reversed(trajectory)))
    assert values == [0, 100]


def test_large_synthetic_heatmap_releases_all_matplotlib_figures() -> None:
    import numpy as np
    from matplotlib import pyplot as plt

    def large_builder(figsize, dpi):
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi, constrained_layout=True)
        ax.imshow(np.zeros((400, 120), dtype=np.uint8), aspect="auto")
        return fig, ax

    with tempfile.TemporaryDirectory() as tmp:
        render_plot(
            RenderJob(Path(tmp) / "large.png", large_builder, 120, 400),
            PlotProfile("test", 600, 800, 100),
        )

    assert plt.get_fignums() == []
