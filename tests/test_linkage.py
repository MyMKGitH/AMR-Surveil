import pandas as pd

from amr_surveil.linkage import build_matched_isolates, build_prediction_context


def test_linkage_uses_exact_sample_id_not_country_year():
    samples = pd.DataFrame(
        {
            "sample_id": ["ISO-1", "ISO-2"],
            "organism": ["E. coli", "E. coli"],
            "country": ["X", "X"],
            "source_sector": ["human", "human"],
            "collection_date": ["2024-01-01", "2024-01-02"],
        }
    )
    calls = pd.DataFrame({"sample_id": ["ISO-1"], "gene": ["blaTEM-1"]})
    ast = pd.DataFrame(
        {
            "sample_id": ["ISO-2"],
            "organism": ["E. coli"],
            "antimicrobial": ["DrugA"],
            "result": ["R"],
            "method": ["BMD"],
            "breakpoint_standard": ["EUCAST"],
            "breakpoint_version": ["16.1"],
        }
    )
    paired = build_matched_isolates(samples, calls, ast)
    assert paired.empty


def test_matched_isolate_contains_only_exact_pair_and_caveat():
    samples = pd.DataFrame(
        {
            "sample_id": ["ISO-1"],
            "organism": ["E. coli"],
            "country": ["X"],
            "source_sector": ["human"],
            "collection_date": ["2024-01-01"],
        }
    )
    calls = pd.DataFrame({"sample_id": ["ISO-1"], "gene": ["blaTEM-1"]})
    ast = pd.DataFrame(
        {
            "sample_id": ["ISO-1"],
            "organism": ["E. coli"],
            "antimicrobial": ["DrugA"],
            "result": ["R"],
            "method": ["BMD"],
            "breakpoint_standard": ["EUCAST"],
            "breakpoint_version": ["16.1"],
        }
    )
    paired = build_matched_isolates(samples, calls, ast)
    assert len(paired) == 1
    assert paired.loc[0, "n_detected_determinants"] == 1
    assert paired.loc[0, "breakpoint_context_complete"]
    assert "not inferred" in paired.loc[0, "interpretation_boundary"]


def test_prediction_context_is_not_ast():
    samples = pd.DataFrame({"sample_id": ["ISO-1"], "organism": ["E. coli"]})
    predictions = pd.DataFrame(
        {"sample_id": ["ISO-1"], "antimicrobial": ["DrugA"], "predicted_mic": ["2"]}
    )
    result = build_prediction_context(samples, predictions)
    assert "predicted_mic" in result.columns
    assert "result" not in result.columns


def test_prediction_organism_conflict_is_excluded():
    samples = pd.DataFrame({"sample_id": ["ISO-1"], "organism": ["E. coli"]})
    predictions = pd.DataFrame(
        {
            "sample_id": ["ISO-1"],
            "organism": ["K. pneumoniae"],
            "antimicrobial": ["DrugA"],
            "predicted_mic": ["2"],
        }
    )
    result = build_prediction_context(samples, predictions)
    assert result.empty


def test_ambiguous_manifest_ids_are_not_linked():
    samples = pd.DataFrame({"sample_id": ["ISO-1", "ISO-1"], "organism": ["E. coli", "E. coli"]})
    calls = pd.DataFrame({"sample_id": ["ISO-1"], "gene": ["blaTEM-1"]})
    ast = pd.DataFrame({"sample_id": ["ISO-1"], "antimicrobial": ["DrugA"], "result": ["R"]})
    assert build_matched_isolates(samples, calls, ast).empty


def test_curated_gene_drug_stratification_enriches_missing_organism_from_manifest():
    from amr_surveil.linkage import gene_drug_stratification

    samples = pd.DataFrame({"sample_id": ["ISO-1", "ISO-2"], "organism": ["E. coli", "E. coli"]})
    calls = pd.DataFrame(
        {
            "sample_id": ["ISO-1", "ISO-2"],
            "gene": ["blaTEM-1", "blaTEM-1"],
            "organism": [pd.NA, pd.NA],
        }
    )
    ast = pd.DataFrame(
        {
            "sample_id": ["ISO-1", "ISO-2"],
            "organism": [pd.NA, pd.NA],
            "antimicrobial": ["Amoxicillin", "Amoxicillin"],
            "result": ["R", "S"],
            "method": ["BMD", "BMD"],
            "breakpoint_standard": ["EUCAST", "EUCAST"],
            "breakpoint_version": ["16.1", "16.1"],
        }
    )
    mapping = pd.DataFrame(
        {
            "organism": ["E. coli"],
            "gene": ["blaTEM-1"],
            "antimicrobial": ["Amoxicillin"],
            "mapping_source": ["curated-test"],
        }
    )
    result = gene_drug_stratification(calls, ast, mapping, samples=samples)
    assert len(result) == 1
    assert result.loc[0, "n_gene_positive_with_ast"] == 2
    assert result.loc[0, "n_R"] == 1
    assert result.loc[0, "R_pct_among_gene_positive_with_ast"] == 50.0


def test_gene_drug_stratification_excludes_organism_conflicts():
    from amr_surveil.linkage import gene_drug_stratification

    samples = pd.DataFrame({"sample_id": ["ISO-1"], "organism": ["E. coli"]})
    calls = pd.DataFrame({"sample_id": ["ISO-1"], "gene": ["blaTEM-1"], "organism": ["E. coli"]})
    ast = pd.DataFrame(
        {
            "sample_id": ["ISO-1"],
            "organism": ["K. pneumoniae"],
            "antimicrobial": ["Amoxicillin"],
            "result": ["R"],
        }
    )
    mapping = pd.DataFrame(
        {
            "organism": ["E. coli"],
            "gene": ["blaTEM-1"],
            "antimicrobial": ["Amoxicillin"],
            "mapping_source": ["curated-test"],
        }
    )
    assert gene_drug_stratification(calls, ast, mapping, samples=samples).empty


def test_matched_linkage_handles_non_string_organism_labels_without_crashing():
    samples = pd.DataFrame({"sample_id": ["ISO-1"], "organism": [123]})
    calls = pd.DataFrame({"sample_id": ["ISO-1"], "gene": ["geneA"], "organism": ["123"]})
    ast = pd.DataFrame(
        {
            "sample_id": ["ISO-1"],
            "organism": [123],
            "antimicrobial": ["DrugA"],
            "result": ["R"],
            "method": ["BMD"],
            "breakpoint_standard": ["EUCAST"],
            "breakpoint_version": ["16.1"],
        }
    )
    paired = build_matched_isolates(samples, calls, ast)
    assert len(paired) == 1
