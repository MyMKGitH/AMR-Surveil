"""Catalogue of supported public sources with licence/terms facts and how each one is accessed.

``terms_checked`` records when (and from which page) the terms below were read by the maintainers.
A source's terms can change: the app shows this text next to every dataset and writes it into the
sample table's ``data_license`` column, so the run provenance carries it. AMR-Surveil never bundles
or redistributes data from these sources; it fetches on the user's request to the user's machine, or
reads a file the user downloaded.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class SourceInfo:
    source_id: str
    name: str
    provider: str
    access: str  # "live_api" | "file_import" | "link_only"
    homepage: str
    hosts: tuple[str, ...]
    scope: str
    terms: str
    terms_url: str
    terms_checked: str
    attribution: str
    caveats: tuple[str, ...]
    live_verified: bool


SOURCES: dict[str, SourceInfo] = {
    "bvbrc": SourceInfo(
        source_id="bvbrc",
        name="BV-BRC AMR phenotypes",
        provider="Bacterial and Viral Bioinformatics Resource Center (NIAID-funded)",
        access="live_api",
        homepage="https://www.bv-brc.org/",
        hosts=("www.bv-brc.org",),
        scope="Isolate-level lab AST (S/I/R, MIC) joined to genome metadata (country, year, host).",
        terms=(
            "BV-BRC states that all public data are freely available without restriction on their use "
            "and contain no PHI/PII. AMR-Surveil does not redistribute the data."
        ),
        terms_url="https://p3c.theseed.org/docs/system_documentation/data_management_sharing.html",
        terms_checked="2026-11-12",
        attribution="Cite BV-BRC and the original data submitters/publications of the genomes you use.",
        caveats=(
            "Genomes in public repositories are not a random sample: outbreak and resistant isolates "
            "are over-represented, so proportions are NOT prevalence.",
            "No patient identifier is available, so patient-level deduplication is impossible.",
            "Query and field names follow the published BV-BRC API reference but could not be tested "
            "against the live server by the maintainers; report any mismatch as an issue.",
        ),
        live_verified=False,
    ),
    "ncbi_pathogen": SourceInfo(
        source_id="ncbi_pathogen",
        name="NCBI Pathogen Detection isolates (TSV import)",
        provider="NCBI / NLM",
        access="file_import",
        homepage="https://www.ncbi.nlm.nih.gov/pathogens/",
        hosts=(),
        scope="Isolates Browser / metadata TSV with AST phenotypes (BioSample antibiograms) and AMR genotypes.",
        terms=(
            "NCBI places no restrictions on use or distribution of its data, but individual submitters "
            "may claim rights NCBI cannot verify; NCBI cannot grant permission on their behalf. "
            "Therefore AMR-Surveil only reads files you download yourself and never bundles them."
        ),
        terms_url="https://www.ncbi.nlm.nih.gov/home/about/policies/",
        terms_checked="2026-11-12",
        attribution="Cite NCBI Pathogen Detection (accessed date) and the BioSample submitters.",
        caveats=(
            "Antibiogram method and breakpoint standard are usually not in the table; records are "
            "labelled 'unspecified' and comparability flags will say so.",
            "Pathogen Detection accessions are not archival; save the file and record the download date.",
            "Genotypes are not imported by default: an empty AMR_genotypes cell cannot distinguish "
            "'screened, nothing found' from 'not screened'.",
            "Column layout follows NCBI's help pages; inspect the detected-columns preview before use.",
        ),
        live_verified=False,
    ),
    "ecdc_atlas": SourceInfo(
        source_id="ecdc_atlas",
        name="ECDC Surveillance Atlas (EARS-Net) CSV import",
        provider="European Centre for Disease Prevention and Control",
        access="file_import",
        homepage="https://atlas.ecdc.europa.eu/public/index.aspx",
        hosts=(),
        scope="Aggregate national resistance counts by organism, antimicrobial group, country and year.",
        terms=(
            "ECDC open-data policy: free of charge under a standardised open-access licence; users must "
            "acknowledge the data sources and ECDC in any resulting product. Aggregate data only; "
            "provided as is."
        ),
        terms_url="https://www.ecdc.europa.eu/en/about-ecdc/who-we-are/key-documents/eueea-routine-surveillance-open-data-policy",
        terms_checked="2026-11-12",
        attribution="Source: ECDC Surveillance Atlas of Infectious Diseases / EARS-Net, data retrieved by the user.",
        caveats=(
            "Aggregate counts only; they are kept as a separate stream and never joined to isolates.",
            "Atlas export column layout varies by indicator; the importer detects columns and shows a "
            "preview. Rows lacking counts (percent only) are skipped and reported.",
            "Antimicrobial entries are drug groups (e.g. third-generation cephalosporins), not single drugs.",
        ),
        live_verified=False,
    ),
    "who_glass": SourceInfo(
        source_id="who_glass",
        name="WHO GLASS dashboard (link-out)",
        provider="World Health Organization",
        access="link_only",
        homepage="https://data.who.int/dashboards/amr/overview",
        hosts=(),
        scope="Global AMR indicators; WHO publishes machine-readable datasets and APIs.",
        terms="Datasets are published under CC BY 4.0 with WHO's added terms (no misrepresentation, no endorsement).",
        terms_url="https://data.who.int/dashboards/amr/overview",
        terms_checked="2026-11-12",
        attribution="WHO. Global Antimicrobial Resistance and Use Surveillance System (GLASS).",
        caveats=(
            "Not yet imported; open the dashboard and upload a downloaded aggregate table as a CSV.",
        ),
        live_verified=False,
    ),
}


def get_source(source_id: str) -> SourceInfo:
    try:
        return SOURCES[source_id]
    except KeyError as exc:
        raise KeyError(f"unknown source {source_id!r}; known: {sorted(SOURCES)}") from exc


def license_text(source_id: str, retrieved_at: str) -> str:
    info = get_source(source_id)
    return f"{info.name}: {info.terms} Retrieved {retrieved_at}. Terms: {info.terms_url}"
