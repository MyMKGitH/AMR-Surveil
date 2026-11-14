# External WGS workflow boundary

AMR-Surveil does not run raw read processing. It accepts canonical, versioned outputs from specialist workflows. The recommended adoption sequence is:

1. Learn read QC, species/contamination checks, assembly, assembly QC and annotation separately in Galaxy Training Network.
2. Repeat a small organism-appropriate example from Linux command line.
3. Record software/database versions and parameters in the exported result tables and project notes.
4. Only then wrap a selected and tested organism-specific workflow in Nextflow or Snakemake (for example, evaluate a current nf-core bacterial assembly workflow).
5. Export AMR calls, plasmid results, typing and metadata in the AMR-Surveil contracts.

Before adding an executable Nextflow/Snakemake pipeline, choose one bacterial species, one sequencing strategy, a pinned toolchain, database versions, reference inputs and a versioned expected-output test set. Avoid a tool zoo with no coherent validated scope.
