"""Canonical table contracts and controlled vocabularies."""

ISSUE_COLUMNS = ["table", "severity", "code", "message", "row_count", "examples"]
SEVERITIES = ("ERROR", "WARNING", "INFO")
VALID_AST_CATEGORIES = {"S", "I", "R"}
VALID_SECTORS = {
    "human",
    "animal",
    "veterinary",
    "food",
    "environment",
    "aquatic",
    "aquaculture",
    "wastewater",
    "unknown",
}
VALID_CARE_SECTORS = {"community", "hospital", "total", "other"}
VALID_CONSUMPTION_SOURCES = {"sales", "distribution", "dispensing", "reimbursement", "other"}

#: How trustworthy the *detection method* behind an AMR call is.
#: ``curated_database``  -- alignment/HMM/point-mutation search against a curated, versioned
#:                          resistance database (e.g. AMRFinderPlus, ResFinder, CARD/RGI).
#: ``homology_search``   -- generic sequence-similarity hit (e.g. remote BLAST against nt);
#:                          a similar sequence is not a validated resistance-gene call.
#: ``heuristic_motif``   -- short peptide/regex motif scan with no curated reference.
#: ``unspecified``       -- tier not stated by the exporter (legacy files).
EVIDENCE_TIERS = ("curated_database", "homology_search", "heuristic_motif", "unspecified")
DEFAULT_EXCLUDED_EVIDENCE_TIERS = ("heuristic_motif", "homology_search")

DEDUPLICATION_RULES = (
    "none",
    "first_per_patient_species_year",
    "first_per_patient_species_specimen_year",
)

SAMPLE_COLUMNS = [
    "sample_id",
    "organism",
    "collection_date",
    "country",
    "region",
    "source_sector",
    "specimen",
    "patient_id",
    "host",
    "sequencing_accession",
    "sequencing_type",
    "source_study",
    "data_license",
]
AMR_CALL_COLUMNS = [
    "sample_id",
    "organism",
    "gene",
    "identity_pct",
    "coverage_pct",
    "contig",
    "start",
    "end",
    "method",
    "database",
    "database_version",
    "source_system",
    "element_type",
    "determinant_class",
    "hit_method",
    "evidence_tier",
]
SCREENING_COLUMNS = [
    "sample_id",
    "method",
    "tool_version",
    "database",
    "database_version",
    "identity_threshold_pct",
    "coverage_threshold_pct",
    "source_system",
]
AST_COLUMNS = [
    "sample_id",
    "organism",
    "antimicrobial",
    "result",
    "mic",
    "mic_unit",
    "zone_mm",
    "method",
    "breakpoint_standard",
    "breakpoint_version",
    "interpretation_date",
    "testing_lab",
]
BREAKPOINT_COLUMNS = [
    "organism",
    "antimicrobial",
    "standard",
    "version",
    "mic_s_le",
    "mic_r_gt",
    "zone_s_ge",
    "zone_r_lt",
    "source_url",
]
AGGREGATE_COLUMNS = [
    "year",
    "country",
    "source_sector",
    "organism",
    "antimicrobial",
    "resistant_n",
    "tested_n",
    "source_name",
    "metric_definition",
    "source_url",
]
PREDICTION_COLUMNS = [
    "sample_id",
    "organism",
    "antimicrobial",
    "predicted_mic",
    "predicted_mic_unit",
    "predicted_class",
    "prediction_score",
    "model_name",
    "model_version",
    "checkpoint_id",
    "training_data_version",
    "split_strategy",
    "source_system",
    "prediction_unit",
    "protein_cluster",
    "decision_threshold",
    "probability_temperature",
    "evaluation_split",
]
PLASMID_COLUMNS = [
    "sample_id",
    "replicon",
    "mobility_prediction",
    "method",
    "database",
    "database_version",
    "source_system",
    "evidence_tier",
]
INCIDENCE_COLUMNS = [
    "year",
    "country",
    "source_sector",
    "organism",
    "antimicrobial",
    "resistant_infection_n",
    "population_n",
    "population_multiplier",
    "source_name",
    "case_definition",
    "source_url",
]
CONSUMPTION_COLUMNS = [
    "year",
    "country",
    "care_sector",
    "antimicrobial",
    "atc_code",
    "ddd_total",
    "population_n",
    "period_days",
    "aware_group",
    "source_name",
    "data_source_type",
    "source_url",
]
GENE_MAP_COLUMNS = ["organism", "gene", "antimicrobial", "mapping_source", "match_type"]
