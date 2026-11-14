# Reference library

These are the public sources that informed the design (reference date 31 October 2026). They provide scientific and operational context; they are not embedded standards or automatically executed data sources. Check the latest version and effective date before applying a protocol or breakpoint.

## Surveillance and reporting

1. WHO. [Global Antimicrobial Resistance and Use Surveillance System (GLASS)](https://www.who.int/initiatives/glass).
2. WHO. [Global antibiotic resistance surveillance report 2025](https://www.who.int/publications/i/item/9789240116337) — analyses 2023 data; do not label the report's data as 2025 surveillance observations.
3. WHO. [GLASS manual for AMR surveillance in common bacterial infections](https://www.who.int/publications/i/item/9789240076600).
4. ECDC. [EARS-Net data and surveillance scope](https://www.ecdc.europa.eu/en/about-us/networks/disease-networks-and-laboratory-networks/ears-net-data).
5. ECDC. [EARS-Net annual epidemiological reports](https://www.ecdc.europa.eu/en/antimicrobial-resistance/surveillance-and-disease-data/annual-epidemiological-reports).
6. ECDC. [Antimicrobial consumption surveillance (ESAC-Net)](https://www.ecdc.europa.eu/en/antimicrobial-consumption).
7. ECDC. [Antimicrobial resistance surveillance and disease data](https://www.ecdc.europa.eu/en/antimicrobial-resistance/surveillance-and-disease-data).
8. Norwegian NORM/NORM-VET. [Annual reports](https://www.vetinst.no/en/surveillance-programmes/norm-norm-vet-report).
9. WHO Regional Office for Europe. [CAESAR programme](https://www.who.int/europe/initiatives/central-asian-and-european-surveillance-of-antimicrobial-resistance-caesar).
10. ECDC/EFSA. [EU summary reports on AMR in zoonotic and indicator bacteria](https://www.ecdc.europa.eu/en/publications-and-data/monitoring-antimicrobial-resistance-zoonotic-and-indicator-bacteria).
11. EMA. [European surveillance of veterinary antimicrobial consumption (ESVAC) and ESUAvet transition](https://www.ema.europa.eu/en/veterinary-regulatory-overview/antimicrobial-resistance-veterinary-medicine/european-surveillance-veterinary-antimicrobial-consumption-esvac-2009-2023).
12. [DANMAP annual reports](https://www.danmap.org/reports).
13. Swedish Veterinary Agency. [Swedres-Svarm reports](https://www.sva.se/en/what-we-do/antibiotics/svarm-resistance-monitoring/swedres-svarm-reports/).

## Laboratory interpretation and data quality

14. EUCAST. [Clinical breakpoint tables](https://www.eucast.org/clinical_breakpoints/). Always capture the table version and applicable date range.
15. EUCAST. [Definitions of S, I and R](https://www.eucast.org/bacteria/clinical-breakpoints-and-interpretation/definition-of-s-i-and-r/). EUCAST `I` means susceptible, increased exposure.
16. ECDC. [External quality assessment of EARS-Net laboratories](https://www.ecdc.europa.eu/en/publications-data/antimicrobial-resistance-external-quality-assessment-laboratories-ears-net-2024).

## Microbial genomics and AMR bioinformatics

17. Galaxy Training Network. [AMR gene-detection learning pathway](https://training.galaxyproject.org/training-material/learning-pathways/amr-gene-detection.html).
18. Galaxy Training Network. [Identification of AMR genes in an assembled bacterial genome](https://training.galaxyproject.org/training-material/topics/genome-annotation/tutorials/amr-gene-detection/tutorial.html).
19. EMBL-EBI Training. [Investigating antimicrobial resistance with EMBL-EBI resources](https://www.ebi.ac.uk/training/online/courses/investigating-antimicrobial-resistance/).
20. NCBI. [Pathogen Detection help and processing documentation](https://www.ncbi.nlm.nih.gov/pathogens/pathogens_help/).
21. NCBI. [AMRFinderPlus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/).
22. [CARD — Comprehensive Antibiotic Resistance Database](https://card.mcmaster.ca/).
23. [ResFinder](https://genomicepidemiology.org/services/).
24. [MOB-suite](https://github.com/phac-nml/mob-suite).
25. [Bakta](https://github.com/oschwengers/bakta).
26. [Nextflow learning resources](https://www.nextflow.io/learn/).
27. [nf-core/bacass](https://github.com/nf-core/bacass).

## Statistics, software and governance

28. [R for Data Science (2e)](https://r4ds.hadley.nz/).
29. Software Carpentry. [The Unix Shell](https://swcarpentry.github.io/shell-novice/).
30. GitHub. [README files and repository documentation](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-readmes).
31. Hugging Face. [Docker Spaces](https://huggingface.co/docs/hub/spaces-sdks-docker).

## Methods and data-source citations

32. Wilson EB. Probable inference, the law of succession, and statistical inference. *Journal of the American Statistical Association*. 1927;22(158):209–212. doi:[10.1080/01621459.1927.10502953](https://doi.org/10.1080/01621459.1927.10502953). — Source of the score interval used for every proportion in this software.
33. CLSI. *M39: Analysis and Presentation of Cumulative Antimicrobial Susceptibility Test Data*, 5th edition (2022). [clsi.org](https://clsi.org/standards/products/microbiology/documents/m39). — Basis for the antibiogram conventions (% susceptible, first isolate per patient per analysis period, minimum isolate count). The standard is copyrighted and not embedded; check the current edition.
34. Feldgarden M, Brover V, Gonzalez-Escalona N, et al. AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. *Scientific Reports*. 2021;11. doi:[10.1038/s41598-021-91456-0](https://doi.org/10.1038/s41598-021-91456-0). — Cite when AMRFinderPlus output is ingested.
35. Olson RD, Assaf R, Brettin T, et al. Introducing the Bacterial and Viral Bioinformatics Resource Center (BV-BRC): a resource combining PATRIC, IRD and ViPR. *Nucleic Acids Research*. 2023. doi:[10.1093/nar/gkac1003](https://doi.org/10.1093/nar/gkac1003). — Cite when BV-BRC AMR data are retrieved.

## Use of references

- Standards and databases can change; do not rely on the versions in an old report without checking official archives.
- This software does not download or apply EUCAST/CLSI breakpoints and does not replace the applicable interpretive standard.
- Upstream software/database citations should be added to each study report, including release, parameters, database snapshot and access date.
- Dataset licence and access conditions are separate from the software's MIT licence.
