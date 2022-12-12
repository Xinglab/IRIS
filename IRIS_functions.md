# IRIS modules

[back to IRIS quick guide](README.md)

For questions about input file format, see example folder.

For test runs with files in example folder, users will need to modify directories of `fin_matrices`, etc.

## `format`

When starting from standard output of [rMATS](https://github.com/Xinglab/rmats-turbo), users should use this step to 1) reformat splice junction counts into a PSI (percent-spliced-in) value matrix, and 2) index and 3) move the PSI matrix for IRIS screening (when -d is enabled).
```
usage: IRIS format [-h] -t {SE,RI,A3SS,A5SS} -n DATA_NAME -s {1,2}
                   [-c COV_CUTOFF] [-i] [-e] [-d IRIS_DB_PATH] [--novelSS]
                   [--gtf GTF]
                   rmats_mat_path_manifest rmats_sample_order

required arguments:
  rmats_mat_path_manifest
                        txt manifest of path(s) to rMATS output folder(s)
  rmats_sample_order    TXT file manifest of corresponding rMATS input sample
                        order file(s). Required input for rMATS
  -t {SE,RI,A3SS,A5SS}, --splicing-event-type {SE,RI,A3SS,A5SS}
                        String of splicing event types based on rMATS
                        definition (SE,RI,A3SS,A5SS).Used to name output file
  -n DATA_NAME, --data-name DATA_NAME
                        Defines dataset name (disease state, study name, group
                        name etc.). Used during IRIS screening
  -s {1,2}, --sample-name-field {1,2}
                        Specifies sample name field (1- SJ count file name, 2-
                        SJ count folder name), for each sample the name should
                        match their name in "rmats_sample_order"

optional arguments:
  -h, --help            show this help message and exit
  -c COV_CUTOFF, --cov-cutoff COV_CUTOFF
                        Average coverage filter for merged matrix (Default is
                        10)
  -i, --sample-based-filter
                        Coverage filter by individual sample not by entire
                        input group. (Default is disabled)
  -e, --merge-events-only
                        Do not perform matrix merge, only merge events list
  -d IRIS_DB_PATH, --iris-db-path IRIS_DB_PATH
                        Path to store the formatted/indexed AS matrix.
                        Strongly recommend to store the AS matrix to the IRIS
                        db by setting the path to the directory containing
                        folders of pre-index AS reference
                        ("full_path/IRIS_data.vX/db"). Default is current
                        location.
  --novelSS             Enable formatting events with splice junctions
                        containing novelSS. (Different and a subset of rMATS
                        novelSS definition. Default is False)
  --gtf GTF             Path to the Genome annotation GTF file. Required input
                        when novelSS is enabled.
```

## `screen`

This step takes a user-defined screening parameter file ([example/NEPC_test.para](example/NEPC_test.para)), and performs comparisons against reference databases, and returns tumor-associated, tumor-recurrent, and tumor-specific AS events based on user-defined criteria.

When the -t option is enabled, the screening step translates identified tumor AS events into peptide sequences that can be used in the prediction step.
```
usage: IRIS screen [-h] -p PARAMETER_FIN
                   [--splicing-event-type {SE,RI,A3SS,A5SS}] -o OUTDIR [-t]
                   [-g GTF] [--all-orf] [--ignore-annotation]
                   [--remove-early-stop] [--min-sample-count MIN_SAMPLE_COUNT]
                   [--use-existing-test-result]

required arguments:
  -p PARAMETER_FIN, --parameter-fin PARAMETER_FIN
                        File of 'IRIS screen' parameters
  --splicing-event-type {SE,RI,A3SS,A5SS}
                        String of splicing event types based on rMATS
                        definition (SE,RI,A3SS,A5SS).Used to name output file.
                        (Default is SE event)
  -o OUTDIR, --outdir OUTDIR
                        Directory of IRIS screening results

optional arguments:
  -h, --help            show this help message and exit
  -t, --translating     Translates IRIS-screened tumor splice junctions into
                        peptides
  -g GTF, --gtf GTF     The Genome annotation GTF file. Required by IRIS
                        translate option.
  --all-orf             Perform the 3 ORF translation. ORF known in the
                        UniProtKB will be labeled as uniprotFrame in the bed
                        file (Default is to use the known ORF ONLY)
  --ignore-annotation   Perform 3 ORF translation without annotating known ORF
                        from the UniProtKB (Default is disabled)
  --remove-early-stop   Discard the peptide if containing early stop codon
                        (Default is keep the truncated peptide)
  --min-sample-count MIN_SAMPLE_COUNT
                        The minimum number of non-missing sample in the input
                        group for an event to be considered for testiing. Once
                        specified, removed events will be written to "notest"
                        file. (Default is no minimum)
  --use-existing-test-result
                        Skip testing and use existing testing result (Default
                        is run full testing steps)
```
Additionally, `screen_sjc` can be performed as part of the 'tumor-specificity' screen. 
```
usage: IRIS screen_sjc [-h] -p PARAMETER_FIN
                       --splicing-event-type {SE,RI,A3SS,A5SS}
                       -e EVENT_LIST_FILE -o OUTDIR
                       [--use-existing-test-result]
                       [--tumor-read-cov-cutoff TUMOR_READ_COV_CUTOFF]
                       [--normal-read-cov-cutoff NORMAL_READ_COV_CUTOFF]
```
Optionally, `screen_cpm` can be performed to as a higher stringent 'tumor-association' screen or less stringent 'tumor-specificity' by using normalized splice junction counts (in CPM).
```
usage: IRIS screen_cpm [-h] -p PARAMETER_FIN
                       --splicing-event-type {SE,RI,A3SS,A5SS}
                       -e EVENT_LIST_FILE -o OUTDIR
                       [--use-existing-test-result]
```

## `predict`

This step takes the screening result and performs annotation of extracellular and HLA-binding epitope predictions to discover immunotherapy targets.

IRIS prediction of HLA-binding epitopes is a massive prediction job that can utilize a compute cluster. The `prediction` step will create scripts to perform subtasks. If properly configured, those subtask scripts can be executed concurrently by snakemake.
```
usage: IRIS predict [-h] --task-dir TASK_DIR -p PARAMETER_FIN
                    [-t {SE,RI,A3SS,A5SS}] [--iedb-local IEDB_LOCAL]
                    [-m MHC_LIST] [--extracellular-only] [--tier3-only]
                    [--gene-exp-matrix GENE_EXP_MATRIX] [-c DELTAPSI_COLUMN]
                    [-d DELTAPSI_CUT_OFF] [-e EPITOPE_LEN_LIST] [--all-orf]
                    [--extracellular-anno-by-junction]
                    IRIS_screening_result_path

required arguments:
  IRIS_screening_result_path
                        Directory of IRIS screening results
  --task-dir TASK_DIR   Directory to write individual task scripts
  -p PARAMETER_FIN, --parameter-fin PARAMETER_FIN
                        File of parameters used in 'IRIS screen'
  -t {SE,RI,A3SS,A5SS}, --splicing-event-type {SE,RI,A3SS,A5SS}
                        String of splicing event types based on rMATS
                        definition (SE,RI,A3SS,A5SS).Used to name output file.
                        (Default is SE event)

optional arguments:
  -h, --help            show this help message and exit
  --iedb-local IEDB_LOCAL
                        Specify local IEDB location (if installed)
  -m MHC_LIST, --mhc-list MHC_LIST
                        List of HLA/MHC types among samples. HLA type follows
                        seq2HLA format
  --extracellular-only  Only predict CAR-T Targets. Will not predict HLA
                        binding.
  --tier3-only          To only run predict on events passing all screen
                        tiers, which is the tier3 output. Will be much faster
                        when both the tier1 and tier3 were used.
  --gene-exp-matrix GENE_EXP_MATRIX
                        Tab-delimited matrix of gene expression vs. samples
  -c DELTAPSI_COLUMN, --deltaPSI-column DELTAPSI_COLUMN
                        Column of deltaPSI value in matrix, 1-based (Default
                        is 5th column)
  -d DELTAPSI_CUT_OFF, --deltaPSI-cut-off DELTAPSI_CUT_OFF
                        Defines cutoff of deltaPSI (or other metric) to select
                        tumor-enriched splice form (Default is 0)
  -e EPITOPE_LEN_LIST, --epitope-len-list EPITOPE_LEN_LIST
                        Epitope length for prediction (Default is 9,10,11)
  --all-orf             Perform prediction based on 3 ORF translation
                        peptides. Enable this if translation/screening used
                        this option (Default is False)
  --extracellular-anno-by-junction
                        By default, CAR-T targets are annotated by association
                        of event with extracellular domain. This option
                        annotates target based on a junction (not recommended)
```

## `epitope_post`

```
usage: IRIS epitope_post [-h] -p PARAMETER_FIN -o OUTDIR
                         [-t {SE,RI,A3SS,A5SS}] -m MHC_BY_SAMPLE
                         [-e GENE_EXP_MATRIX] [--tier3-only] [--keep-exist]
                         [--epitope-len-list EPITOPE_LEN_LIST]
                         [--no-match-to-canonical-proteome]
                         [--no-uniqueness-annotation]
                         [--ic50-cut-off IC50_CUT_OFF]

required arguments:
  -p PARAMETER_FIN, --parameter-fin PARAMETER_FIN
                        File of parameters used in IRIS screen
  -o OUTDIR, --outdir OUTDIR
                        Directory of IRIS screening results
  -t {SE,RI,A3SS,A5SS}, --splicing-event-type {SE,RI,A3SS,A5SS}
                        String of splicing event types based on rMATS
                        definition (SE,RI,A3SS,A5SS).Used to name output file
                        (Default is SE event)
  -m MHC_BY_SAMPLE, --mhc-by-sample MHC_BY_SAMPLE
                        Tab-delimited matrix of HLA/MHC type vs. samples. HLA
                        type follows seq2HLA format
  -e GENE_EXP_MATRIX, --gene-exp-matrix GENE_EXP_MATRIX
                        Tab-delimited matrix of gene expression vs. samples

optional arguments:
  -h, --help            show this help message and exit
  --tier3-only          Only predict tier3 events. Will be much faster.
  --keep-exist          Do not rewrite a new postive prediction file when the
                        file existed. Default is False
  --epitope-len-list EPITOPE_LEN_LIST
                        Epitope length for prediction (Default is 9,10,11)
  --no-match-to-canonical-proteome
                        Matches epitopes to UniProt canonical protein
                        sequences as an annotation.
  --no-uniqueness-annotation
                        Matches epitopes to all IRIS translated junction
                        peptides in the same analysis as an annotation.
  --ic50-cut-off IC50_CUT_OFF
                        Specifies IC50 cut-off to define HLA-binding epitopes
                        (Default is 500)
```

## `process_rnaseq`

When starting from a fastq file, users should use this step to perform RNA-Seq alignment and quantification. This module uses STAR and cufflinks. This module only takes one sample (can be multiple fastq files) for each run. Users are recommended to run this module in parallel (use `makesubsh_mapping` for snakemake).
```
usage: IRIS process_rnaseq [-h] --starGenomeDir STARGENOMEDIR --gtf GTF -p
                           SAMPLEID_OUTDIR [--db-length DB_LENGTH] [--mapping]
                           [--quant] [--sort]
                           readsFilesRNA

required arguments:
  --starGenomeDir STARGENOMEDIR
                        The path to the STAR indexed reference genome. Pass to
                        the "genomeDir" parameter in STAR
  --gtf GTF             Path to the Genome annotation GTF file
  -p SAMPLEID_OUTDIR, --sampleID-outdir SAMPLEID_OUTDIR
                        Output directory where sample ID will be used as the
                        output folder name
  --db-length DB_LENGTH
                        Pass to the "sjdbOverhang" parameter in STAR. Default
                        is 100
  readsFilesRNA         Specify the path to the paired-end FASTQ files for the
                        sample. Files are seperated eperated by ",".

optional arguments:
  -h, --help            show this help message and exit
  --mapping             Only perform reads mapping
  --quant               Only perform gene expression and AS quantification
  --sort                Only perform BAM file sorting
```

## `makesubsh_mapping`

Run `process_rnaseq` jobs in parallel on HPC or cloud based on snakemake.
```
usage: IRIS makesubsh_mapping [-h] [--fastq-folder-dir FASTQ_FOLDER_DIR]
                              --starGenomeDir STARGENOMEDIR --gtf GTF
                              --data-name DATA_NAME --outdir OUTDIR
                              --label-string LABEL_STRING --task-dir TASK_DIR

required arguments:
  --fastq-folder-dir FASTQ_FOLDER_DIR
                        Specify the path to the higher level of all folders
                        containing FASTQ files
  --starGenomeDir STARGENOMEDIR
                        The path to the STAR indexed reference genome. Pass to
                        the "genomeDir" parameter in STAR
  --gtf GTF             Path to the Genome annotation GTF file
  --data-name DATA_NAME
                        Data set name used to name submission shell scripts
                        files.
  --outdir OUTDIR       Output directory for folders of aligned BAM files
  --label-string LABEL_STRING
                        String in the fastq file name between the reads pair
                        number and "fastq/fq". This is used to recognize
                        paired-end reads. e.g. For FASTQ_file_L1_R2.fastq.gz,
                        the label string is the "." between "2" and "fastq".
  --task-dir TASK_DIR   Directory to write individual task scripts

optional arguments:
  -h, --help            show this help message and exit
```

## `makesubsh_rmats`

After running `process_rnaseq`, this step can be used to prepare files to run rMATS-turbo in parallel.
```
usage: IRIS makesubsh_rmats [-h] --rMATS-path RMATS_PATH --bam-dir BAM_DIR
                            [--bam-prefix BAM_PREFIX] --gtf GTF --data-name
                            DATA_NAME --task-dir TASK_DIR [--novelSS]
                            [--read-length READ_LENGTH]

required arguments:
  --rMATS-path RMATS_PATH
                        Path to the rMATS-turbo script.
  --bam-dir BAM_DIR     The path one level higher to folders containing BAM
                        file generated by "process_rnaseq".
  --bam-prefix BAM_PREFIX
                        BAM file prefix (Default is
                        "Aligned.sortedByCoord.out")
  --gtf GTF             Path to the Genome annotation GTF file
  --data-name DATA_NAME
                        Data set name used to name submission shell scripts
  --task-dir TASK_DIR   Directory to write individual task scripts

optional arguments:
  -h, --help            show this help message and exit
  --novelSS             Enable rMATS novelSS option to include novel splice
                        site detected from the RNA-seq data (Default is False)
  --read-length READ_LENGTH
                        User defined read length instead of using STAR maaping
                        log file to define automatically.
```

## `makesubsh_rmatspost`

After running `makesubsh_rmats`, this step can be used to merge files to generate final rMATS-turbo results.
```
usage: IRIS makesubsh_rmatspost [-h] --rMATS-path RMATS_PATH --bam-dir BAM_DIR
                                --gtf GTF --data-name DATA_NAME [--novelSS]
                                --task-dir TASK_DIR

required arguments:
  --rMATS-path RMATS_PATH
                        Path to the rMATS-turbo scripte
  --bam-dir BAM_DIR     The path one level higher to folders containing BAM
                        file generated by "process_rnaseq".
  --gtf GTF             Path to the Genome annotation GTF file
  --data-name DATA_NAME
                        Data set name used to name submission shell scripts
  --task-dir TASK_DIR   Directory to write individual task scripts

optional arguments:
  -h, --help            show this help message and exit
  --novelSS             Enable rMATS novelSS option to include novel splice
                        site detected from the RNA-seq data (Default is False)
```

## `exp_matrix`

After running `process_rnaseq`, if samples of interest are all processed, users can use this script to generate a gene expression matrix, which will be used as annotations in downstream IRIS prediction and/or proteomics reports.
```
usage: IRIS exp_matrix [-h] [--exp-cutoff EXP_CUTOFF] -o OUTDIR -n DATA_NAME
                       gene_exp_file_list

required arguments:
  gene_exp_file_list    A txt manifest of path(s) of cufflinks gene expression
                        output(s).
  -n DATA_NAME, --data-name DATA_NAME
                        Name of the dataset (disease state, study name, group
                        name etc.).

optional arguments:
  -h, --help            show this help message and exit
  --exp-cutoff EXP_CUTOFF
                        Gene expression cut-off based on FPKM (Default is 1)
  -o OUTDIR, --outdir OUTDIR
                        Output directory for IRIS exp_matrix
```

## `index`

This step is incorporated by formatting. For users who already have a matrix of AS PSI values (generated by rMATS or another tool), this command could finish the indexing and other steps to prepare for IRIS screening.
```
usage: IRIS index [-h] -t {SE,RI,A3SS,A5SS} -n DATA_NAME [-c COV_CUTOFF]
                  [-o OUTDIR]
                  splicing_matrix

required arguments:
  splicing_matrix       Tab-delimited matrix of splicing events (row) vs.
                        sample IDs (col)
  -t {SE,RI,A3SS,A5SS}, --splicing-event-type {SE,RI,A3SS,A5SS}
                        String of splicing event types based on rMATS
                        definition (SE,RI,A3SS,A5SS).Used to name output file
  -n DATA_NAME, --data-name DATA_NAME
                        Name of data matrix (disease state, study name, group
                        name, etc.) being indexed. Used by IRIS during
                        screening

optional arguments:
  -h, --help            show this help message and exit
  -c COV_CUTOFF, --cov-cutoff COV_CUTOFF
                        For the naming purpose, Input average coverage cutoff
                        used when generating the PSI matrix (Default is 10)
  -o OUTDIR, --outdir OUTDIR
                        Output directory for IRIS database
```

## `translate`

```
usage: IRIS translate [-h] -g REF_GENOME -t {SE,RI,A3SS,A5SS} --gtf GTF -o
                      OUTDIR [--all-orf] [--ignore-annotation]
                      [--remove-early-stop] [-c DELTAPSI_COLUMN]
                      [-d DELTAPSI_CUT_OFF] [--no-tumor-form-selection]
                      [--check-novel]
                      as_input

required arguments:
  as_input              Inputs AS event coordinates and delta PSI values
  -g REF_GENOME, --ref-genome REF_GENOME
                        Specifies reference genome (FASTA format) location
  -t {SE,RI,A3SS,A5SS}, --splicing-event-type {SE,RI,A3SS,A5SS}
                        String of splicing event types based on rMATS
                        definition (SE,RI,A3SS,A5SS).Used to name output file
  --gtf GTF             Path to the Genome annotation GTF file. Used to define
                        exon ends for microexons
  -o OUTDIR, --outdir OUTDIR
                        Defines IRIS translation output directory

optional arguments:
  -h, --help            show this help message and exit
  --all-orf             Perform the 3 ORF translation. ORF known in the
                        UniProtKB will be labeled as uniprotFrame in the bed
                        file (Default is to use the known ORF ONLY)
  --ignore-annotation   Perform 3 ORF translation without annotating known ORF
                        from the UniProtKB (Default is disabled)
  --remove-early-stop   Discard the peptide if containing early stop codon
                        (Default is keep the truncated peptide)
  -c DELTAPSI_COLUMN, --deltaPSI-column DELTAPSI_COLUMN
                        Column of deltaPSI value in matrix, 1-based (Default
                        is 5th column)
  -d DELTAPSI_CUT_OFF, --deltaPSI-cut-off DELTAPSI_CUT_OFF
                        Defines cutoff of deltaPSI (or other metric) used to
                        select tumor-enriched splice form (Default is 0)
  --no-tumor-form-selection
                        Translates splicing junctions derived from both
                        skipping and inclusion forms (Default is False)
  --check-novel         Translates splicing junctions derived from novel
                        splice sites only using information passed from
                        screen_novelss (Default is False)
```

## `makesubsh_hla`

This step uses the RNA-Seq FASTQ file to infer the HLA type of a sample.
```
usage: IRIS makesubsh_hla [-h] [--fastq-folder-dir FASTQ_FOLDER_DIR]
                          --data-name DATA_NAME -o OUTDIR --label-string
                          LABEL_STRING --task-dir TASK_DIR

required arguments:
  --fastq-folder-dir FASTQ_FOLDER_DIR
                        Specify the path to the higher level of all folders
                        containing FASTQ files
  --data-name DATA_NAME
                        Data set name used to name submission shell scripts.
  -o OUTDIR, --outdir OUTDIR
                        Output directory for folders of seq2hla result
  --label-string LABEL_STRING
                        String in the fastq file name between the reads pair
                        number and "fastq/fq". This is used to recognize
                        paired-end reads. e.g. For FASTQ_file_L1_R2.fastq.gz,
                        the label string is the "." between "2" and "fastq".
  --task-dir TASK_DIR   Directory to write individual task scripts

optional arguments:
  -h, --help            show this help message and exit
```

## `pep2epitope`

This module is a wrapper of prediction tools (IEDB) for predicting peptide-HLA binding. The `prediction` and `epitope_post` modules can generate scripts to run this module in parallel and summarize the result into one TCR target report.
```
usage: IRIS pep2epitope [-h] [-e EPITOPE_LEN_LIST] [-a HLA_ALLELE_LIST] -o
                        OUTDIR [--iedb-local IEDB_LOCAL]
                        [--ic50-cut-off IC50_CUT_OFF]
                        junction_pep_input

required arguments:
  junction_pep_input    Inputs junction peptides
  -e EPITOPE_LEN_LIST, --epitope-len-list EPITOPE_LEN_LIST
                        Epitope length for prediction (Default is 9,10,11)
  -a HLA_ALLELE_LIST, --hla-allele-list HLA_ALLELE_LIST
                        List of HLA types (Default is HLA-A*01:01,
                        HLA-B*08:01, HLA-C*07:01)
  -o OUTDIR, --outdir OUTDIR
                        Define output directory of pep2epitope
  --iedb-local IEDB_LOCAL
                        Specify local IEDB location (if installed)
  --ic50-cut-off IC50_CUT_OFF
                        Cut-off based on median value of consensus-predicted
                        IC50 values (Default is 500)

optional arguments:
  -h, --help            show this help message and exit
```
