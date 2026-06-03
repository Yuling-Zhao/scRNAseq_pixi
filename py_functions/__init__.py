"""Callable Python helpers for the scRNA-seq workflow."""

from __future__ import annotations

from typing import Any


def semi_automated_marker_annotation(*args: Any, **kwargs: Any) -> Any:
    """Run semi-automated marker annotation.

    Heavy scientific dependencies are imported only when the function is called.
    """
    from .semi_manual_annotation import semi_automated_marker_annotation as _impl

    return _impl(*args, **kwargs)

__all__ = ["semi_automated_marker_annotation"]
