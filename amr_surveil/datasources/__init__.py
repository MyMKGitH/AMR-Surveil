"""Public-data hub: bounded, auditable access to openly available AMR data.

See :mod:`amr_surveil.datasources.registry` for the supported sources and their terms.
"""

from amr_surveil.datasources.bvbrc import fetch_bvbrc, parse_bvbrc, parse_bvbrc_export
from amr_surveil.datasources.common import PublicDataset, write_inputs
from amr_surveil.datasources.ecdc import parse_ecdc_atlas
from amr_surveil.datasources.fetch import FetchError, FetchResult, default_cache_dir, fetch
from amr_surveil.datasources.geo import country_to_iso3, iso3_to_name
from amr_surveil.datasources.ncbi_pathogen import parse_pathogen_table
from amr_surveil.datasources.registry import SOURCES, SourceInfo, get_source

__all__ = [
    "SOURCES",
    "FetchError",
    "FetchResult",
    "PublicDataset",
    "SourceInfo",
    "country_to_iso3",
    "default_cache_dir",
    "fetch",
    "fetch_bvbrc",
    "get_source",
    "iso3_to_name",
    "parse_bvbrc",
    "parse_bvbrc_export",
    "parse_ecdc_atlas",
    "parse_pathogen_table",
    "write_inputs",
]
