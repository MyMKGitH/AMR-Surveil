"""AMR-Surveil: reproducible genomic and epidemiological AMR surveillance."""

from amr_surveil._version import __version__
from amr_surveil.pipeline import run_pipeline

__all__ = ["__version__", "run_pipeline"]
