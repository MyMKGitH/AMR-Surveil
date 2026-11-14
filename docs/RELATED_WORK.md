# Related software and data sources

Reviewed 12 November 2026 from the projects' own pages and papers. Statements marked *(title only)* come from a search-result title and were **not** read in depth; treat them as pointers, not claims.

## Where AMR-Surveil sits

| Tool | What it is | Overlap | What AMR-Surveil adds / does not do |
|---|---|---|---|
| **AMRnet** | Free open-source (GPLv3) dashboard of genome-derived AMR for *Salmonella* Typhi, *Klebsiella pneumoniae*, *N. gonorrhoeae*, *E. coli*, *Shigella* and non-typhoidal *Salmonella*; per-country maps (countries with fewer than 20 samples greyed), annual trends, genotype views, REST API; data from Pathogenwatch and Enterobase. States its data are not curated for sampling bias. | Country map with minimum-n suppression, trend views, sampling-bias caveat | AMR-Surveil analyses **your** data with an audit trail (QA gate, comparability flags, provenance manifest, genotype–phenotype go/no-go). It does **not** host curated global genome collections. |
| **WHONET / BacLink** | WHO Collaborating Centre desktop software for laboratory-based surveillance: data entry, import from lab systems, antibiograms, alerts *(title/summary only)* | Antibiogram | AMR-Surveil is not a laboratory data-management system. Importing WHONET files is on the [roadmap](ROADMAP.md). |
| **Pathogenwatch / Microreact** | Genomic surveillance web platforms from the Centre for Genomic Pathogen Surveillance: AMR prediction, typing, phylogeny and mapping *(title only)* | Maps, AMR determinants | No sequence analysis, typing or phylogeny here; AMR-Surveil consumes already-called determinants. |
| **AMRmap** | Interactive web platform for AMR surveillance data in Russia *(title only: PMC7994358)* | Interactive resistance views | Country-specific; not reviewed in depth. |
| **NCBI Pathogen Detection** | Isolate metadata, AST phenotypes (from BioSample antibiograms) and AMRFinderPlus genotypes for surveillance isolates | Data source | Imported from a table you download. |
| **BV-BRC** | Public resource with lab AST phenotypes joined to genomes and an open API | Data source | Queried live, bounded. |
| **ECDC Surveillance Atlas / EARS-Net** | Official EU/EEA aggregate resistance counts | Data source | Imported as a separate aggregate stream. |
| **WHO GLASS dashboard** | Global indicators with machine-readable datasets and APIs (CC BY 4.0 plus WHO conditions) | Data source | Link-out only so far. |

## Features taken from this review

* **Map with suppression** (AMRnet greys countries under 20 samples): default minimum 20 isolates per country, suppressed rows keep counts but get no colour.
* **Cumulative antibiogram** in the style of CLSI M39 (S/(S+I+R), minimum isolate count, contexts never merged).
* **Per-source terms panel**: each public dataset shows its licence/terms, required attribution and caveats, and the same text is written to the sample table's `data_license` column so it travels into the run manifest.
* **Sampling-bias warnings** on every map and public dataset, as AMRnet does.

## Considered and not (yet) done

| Idea | Why not now |
|---|---|
| Bundling EUCAST/CLSI breakpoint tables | Version-sensitive and licence-restricted; the user supplies the table and its hash is recorded. |
| MDR / co-resistance profiles, trend tests | Valuable but need a stated definition (e.g. an ECDC/CDC MDR standard) and validation first. See [roadmap](ROADMAP.md). |
| Phylogeny, SNP clusters, transmission | Out of scope; relatedness is not transmission. |
| Bundled public datasets | Never: terms differ per submitter (NCBI says it cannot grant permission on submitters' behalf). Data are fetched or imported on the user's machine instead. |
| WHO GLASS import | The dashboard says datasets and APIs exist, but the exact file layout was not verified. |

## Sources

* AMRnet documentation: <https://amrnet.readthedocs.io/>
* ECDC EU/EEA routine surveillance open data policy: <https://www.ecdc.europa.eu/en/about-ecdc/who-we-are/key-documents/eueea-routine-surveillance-open-data-policy>
* WHO GLASS dashboard and terms: <https://data.who.int/dashboards/amr/overview>
* BV-BRC data management and sharing: <https://p3c.theseed.org/docs/system_documentation/data_management_sharing.html>; API: <https://www.patricbrc.org/api/doc/>
* NCBI policies: <https://www.ncbi.nlm.nih.gov/home/about/policies/>; Pathogen Detection help: <https://www.ncbi.nlm.nih.gov/pathogens/pathogens_help/>
* AMRmap: <https://pmc.ncbi.nlm.nih.gov/articles/PMC7994358>
