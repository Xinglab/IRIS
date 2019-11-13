[back to IRIS quick guide](README.md)

For questions about input file format, see example folder. 

For test runs with files in example folder, users will need to modify directories of 'fin_matrices', etc.

##### formatting
When starting from standard output of [rMATS](http://rnaseq-mats.sourceforge.net), users should use this step to 1) reformat splice junction counts into a PSI (percent-spliced-in) value matrix, and 2) index and 3) move the PSI matrix for IRIS screening (when -d is enabled). 
```bash
IRIS formatting -h
usage: IRIS formatting [-h] -t {SE,RI,A3,A5} -n DATA_NAME -s {1,2}
                       [-c COV_CUTOFF] [-e] [-d IRIS_DB_PATH]
                       rmats_mat_path_manifest rmats_sample_order

required arguments:
  rmats_mat_path_manifest
                        txt manifest of path(s) to rMATS output folder(s)
  rmats_sample_order    txt manifest of corresponding rMATS input sample order file(s) 
                        Required input for rMATS
  -t {SE,RI,A3,A5}, --splicing_event_type {SE,RI,A3,A5}
                        String of splicing event types based on rMATS definition (SE,RI,A3,A5) 
                        Used to name output file
  -n DATA_NAME, --data-name DATA_NAME
                        Defines dataset name (disease state, study name, group name etc.) 
                        Used during IRIS screening
  -s {1,2}, --sample-name-field {1,2}
                        Specifies sample name field for each sample in sample order file(s) 
                        listed by "rmats_sample_order" (1- BAM file name, 2- BAM folder name)

optional arguments:
  -h, --help            Shows help message and exits
  -c COV_CUTOFF, --cov-cutoff COV_CUTOFF
                        Average coverage filter for merged matrix (default is 10)
  -e, --merge-events-only
                        Do not perform matrix merge, only merge events list
  -d IRIS_DB_PATH, --iris-db-path IRIS_DB_PATH
                        Path to IRIS database 
                        Formatted/indexed AS matrices are stored here and used during IRIS screening
```

##### screening
This step takes a user-defined screening parameter file ([example](example/Test.para)), which performs comparisons against reference databases, and returns tumor-associated, tumor-recurrent, and tumor-specific AS events based on user-defined criteria. 

When the -t option is enabled, the screening step translates identified tumor AS events into peptide sequences that can be used in the prediction step. 
```bash
IRIS screening -h
usage: IRIS screening [-h] [-o OUTDIR] [-t] parameter_fin

required arguments:
  parameter_fin         File of IRIS screening parameters
  -o OUTDIR, --outdir OUTDIR
                        Directory of IRIS screening results

optional arguments:
  -h, --help            Shows help message and exits
  -t, --translating     Translates IRIS-screened tumor splice junctions into peptides
```

##### prediction 
This step takes the screening result and performs annotation of extracellular and HLA-binding epitope predictions to discover immunotherapy targets. 

IRIS prediction of HLA-binding epitopes is a massive prediction job that requires access to computing clusters with the SGE system for completion. The 'prediction' step will create qsub scripts for job array submission.

###### Perform extracellular (CAR-T) target annotation & prepare epitope (TCR) target prediction
 ```
IRIS prediction -h
usage: IRIS prediction [-h] [-p PARAMETER_FIN] [--iedb-local IEDB_LOCAL]
                       [-c DELTAPSI_COLUMN] [-d DELTAPSI_CUT_OFF] -m MHC_LIST
                       [--extracellular-anno-by-junction]
                       IRIS_screening_result_path

required arguments:
  IRIS_screening_result_path
                        Input AS event coordinates and PSI values
  -p PARAMETER_FIN, --parameter-fin PARAMETER_FIN
                        File of parameters used in IRIS screening
  --iedb-local IEDB_LOCAL
                        Specify local IEDB location (if installed)
  -m MHC_LIST, --mhc-list MHC_LIST
                        List of HLA/MHC types among samples 
                        HLA type follows seq2HLA format

optional arguments:
  -h, --help            Shows help message and exits
  -c DELTAPSI_COLUMN, --deltaPSI-column DELTAPSI_COLUMN
                        Column of deltaPSI value in matrix, 1-based (default is 5th column)
  -d DELTAPSI_CUT_OFF, --deltaPSI-cut-off DELTAPSI_CUT_OFF
                        Defines cutoff of deltaPSI (or other metric) to select tumor-enriched 
                        splice form (default is 0)
  --extracellular-anno-by-junction
                        By default, CAR-T targets are annotated by association of event 
                        with extracellular domain 
                        This option annotates target based on a junction (not recommended)
 ```

###### Epitope (TCR) target prediction (requires SGE system)
 ```
IRIS epitope_post -h
usage: IRIS epitope_post [-h] -p PARAMETER_FIN -o OUTDIR -m MHC_BY_SAMPLE
                         [-e GENE_EXP_MATRIX] [--ic50-cut-off IC50_CUT_OFF]

required arguments:
  -p PARAMETER_FIN, --parameter_fin PARAMETER_FIN
                        File of parameters used in IRIS screening
  -o OUTDIR, --outdir OUTDIR
                        Directory of IRIS screening results
  -m MHC_BY_SAMPLE, --mhc-by-sample MHC_BY_SAMPLE
                        Tab-delimited matrix of HLA/MHC type vs. samples
                        HLA type follows seq2HLA format
  -e GENE_EXP_MATRIX, --gene-exp-matrix GENE_EXP_MATRIX
                        Tab-delimited matrix of gene expression vs. samples

optional arguments:
  -h, --help            Shows help message and exits
  --ic50-cut-off IC50_CUT_OFF
                        Specifies IC50 cut-off to define HLA-binding epitopes (default is 500)
 ```
 
##### process_rnaseq
When starting from a FASTQ file, users should use this step to perform RNA-Seq alignment and quantification. This module uses STAR and cufflinks. This module only takes one sample (can be multiple FASTQ files) for each run. Users are recommended to run this module in parallel in the SGE system. 
```
IRIS process_rnaseq -h
usage: IRIS process_rnaseq [-h] --starGenomeDir STARGENOMEDIR --gtf GTF -p
                           SAMPLEID_OUTDIR [--db-length DB_LENGTH] [--mapping]
                           [--quant] [--sort]
                           readsFilesRNA

required arguments:
  --starGenomeDir STARGENOMEDIR
                        Path to STAR-indexed reference genome 
                        Passes to "genomeDir" parameter in STAR
  --gtf GTF             Genome annotation file.
  -p SAMPLEID_OUTDIR, --sampleID-outdir SAMPLEID_OUTDIR
                        Output directory, where sample ID will be used as output folder name
  --db-length DB_LENGTH
                        Passes to "sjdbOverhang" parameter in STAR (default is 100)
  readsFilesRNA         Specifies path to paired-end FASTQ files for sample 
                        Files separated by ","

optional arguments:
  -h, --help            Shows help message and exits
  --mapping             Only perform reads mapping
  --quant               Only perform gene expression and AS quantification
  --sort                Only perform BAM file sorting
 ```
 
##### makeqsub_rmats (requires SGE system)
After running 'process_rnaseq', this step can be used to prepare files to run rMATS-turbo in parallel in the SGE system.
```
IRIS makeqsub_rmats -h
usage: IRIS makeqsub_rmats [-h] --rMATS-path RMATS_PATH --bam-dir BAM_DIR
                           --gtf GTF --read-length READ_LENGTH

required arguments:
  --rMATS-path RMATS_PATH
                        Path to rMATS-turbo script
  --bam-dir BAM_DIR     Path one level higher to folders containing BAM file generated by 'process_rnaseq'
  --gtf GTF             Genome annotation file
  --read-length READ_LENGTH
                        Passes to "readLength" parameter in rMATS-turbo

optional arguments:
  -h, --help            Shows help message and exits 
 ```

##### exp_matrix
After running 'process_rnaseq', if samples of interest are all processed, users can use this script to generate a gene expression matrix, which will be used as annotations in downstream IRIS prediction and/or proteomics reports.
```
IRIS exp_matrix -h
usage: IRIS exp_matrix [-h] [--exp-cutoff EXP_CUTOFF] -o OUTDIR -n DATA_NAME
                       gene_exp_file_list

required arguments:
  gene_exp_file_list    txt manifest of path(s) of cufflinks gene expression output(s)
  -n DATA_NAME, --data-name DATA_NAME
                        Name of dataset (disease state, study name, group name, etc.)

optional arguments:
  -h, --help            Shows help message and exits
  --exp-cutoff EXP_CUTOFF
                        Gene expression cut-off based on FPKM (default is 1)
  -o OUTDIR, --outdir OUTDIR
                        Output directory for IRIS exp_matrix
```

##### indexing
This step is incorporated by formatting. For users who already have a matrix of AS PSI values (generated by rMATS or another tool), this command could finish the indexing and other steps to prepare for IRIS screening.
```bash 
IRIS indexing -h
usage: IRIS indexing [-h] -n DATA_NAME [-d DB_DIR] splicing_matrix

required arguments:
  splicing_matrix       Tab-delimited matrix of splicing events (row) vs. sample IDs (col)
  -n DATA_NAME, --data-name DATA_NAME
                        Name of data matrix (disease state, study name, group name, etc.) being 
                        formatted & indexed 
                        Used by IRIS during screening

optional arguments:
  -h, --help            Shows help message and exits
  -d DB_DIR, --db-dir DB_DIR
                        Directory of IRIS database 
                        Program creates a folder in this directory for IRIS to recognize
```

##### translation
```bash 
IRIS translation -h
usage: IRIS translation [-h] -g REF_GENOME -o OUTDIR [-c DELTAPSI_COLUMN]
                        [-d DELTAPSI_CUT_OFF] [--no-tumor-form-selection]
                        as_input

required arguments:
  as_input              Inputs AS event coordinates and PSI values
  -g REF_GENOME, --ref-genome REF_GENOME
                        Specifies reference genome (FASTA format) location
  -o OUTDIR, --outdir OUTDIR
                        Defines IRIS translation output directory

optional arguments:
  -h, --help            Show help message and exits
  -c DELTAPSI_COLUMN, --deltaPSI-column DELTAPSI_COLUMN
                        Column of deltaPSI value in matrix, 1-based (default is 5th column)
  -d DELTAPSI_CUT_OFF, --deltaPSI-cut-off DELTAPSI_CUT_OFF
                        Defines cutoff of deltaPSI (or other metric) used to select tumor-enriched 
                        splice form (default is 0)
  --no-tumor-form-selection
                        Translates splicing junctions derived from both skipping and inclusion forms
 ```
 
##### seq2hla 
This step uses the RNA-Seq FASTQ file to infer the HLA type of a sample.
```
IRIS seq2hla -h
usage: IRIS seq2hla [-h] -b SEQ2HLA_PATH -p SAMPLEID_OUTDIR readsFilesCaseRNA

required arguments:
  -b SEQ2HLA_PATH, --seq2hla-path SEQ2HLA_PATH
                        Path to seq2hla folder
  -p SAMPLEID_OUTDIR, --sampleID-outdir SAMPLEID_OUTDIR
                        Output directory, where sample ID will be used as output folder name
  readsFilesCaseRNA     Tumor sample paired-end fastq files, separated by ","

optional arguments:
  -h, --help            Shows help message and exits
 ```
 
##### pep2epitope
This module is a wrapper of prediction tools (IEDB) for predicting peptide-HLA binding. The 'prediction' and 'epitope_post' modules can make qsub submissions to run this module in parallel and summarize the result into one TCR target report.
```IRIS pep2epitope -h
usage: IRIS pep2epitope [-h] [-e EPITOPE_LEN_LIST] [-a HLA_ALLELE_LIST] -o
                        OUTDIR [--iedb-local IEDB_LOCAL]
                        [--ic50-cut-off IC50_CUT_OFF]
                        junction_pep_input

required arguments:
  junction_pep_input    Inputs AS event coordinates and PSI values
  -e EPITOPE_LEN_LIST, --epitope-len-list EPITOPE_LEN_LIST
                        Epitope length for prediction (default is 9,10,11)
  -a HLA_ALLELE_LIST, --hla-allele-list HLA_ALLELE_LIST
                        List of HLA types (default is HLA-A*01:01, HLA-B*08:01, HLA-C*07:01)
  -o OUTDIR, --outdir OUTDIR
                        Define output directory of pep2epitope
  --iedb-local IEDB_LOCAL
                        Specify local IEDB location (if installed)
  --ic50-cut-off IC50_CUT_OFF
                        Cut-off based on median value of consensus-predicted IC50 values (default is 500)
                   ```

optional arguments:
  -h, --help            Shows help message and exits
