# IRIS: Isoform peptides from RNA splicing for Immunotherapy target Screening



### Quick guide
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Usage](#usage)
  - [Usage - individual modules (for customized pipelines)](#individual-modules)
  - [Usage - streamlined major modules (for common use)](#streamlined-major-modules)
- [Example](#example)
- [Output](#example-output)
- [Contact](#contact)
- [Publication](#citation)



### Dependencies

#### Core dependencies (required for major IRIS modules/steps - formatting, screening, and prediction):
- python 2.7.x (numpy, scipy, seaborn, pyBigWig, etc.)
- [IEDB stand-alone 20130222 2.15.5 (2.22.1 is not fully tested)](http://tools.iedb.org/main/download/)
- [bedtools 2.29.0](https://bedtools.readthedocs.io/en/latest/)

#### Other dependencies (required for processing raw RNA-Seq and MS data)
- [STAR 2.5.3](https://github.com/alexdobin/STAR/releases/tag/2.5.3a): required for IRIS RNA-seq processing
- [samtools 1.3](https://sourceforge.net/projects/samtools/files/samtools/): required for IRIS RNA-seq processing
- [rMATS-turbo](http://rnaseq-mats.sourceforge.net): required for IRIS RNA-seq processing
- [Cufflinks 2.2.1](http://cole-trapnell-lab.github.io/cufflinks/install/): required for IRIS RNA-seq processing
- [seq2HLA](https://bitbucket.org/sebastian_boegel/seq2hla/src/default/): required for HLA typing; requires [bowtie](http://bowtie-bio.sourceforge.net/index.shtml)
- [MS GF+ (v2018.07.17)](https://github.com/MSGFPlus/msgfplus): required for MS search; requiring [Java](https://www.java.com/en/download/)



### Installation
Two steps to set up IRIS:
#### 1. Download
##### 1.1 Download IRIS program
The IRIS program can be downloaded directly from the repository, as shown below:
```
git clone https://github.com/Xinglab/IRIS.git
cd IRIS
```
__For full functionality, IRIS requires use of the SGE system. For users who want to use functions involving SGE (see [Usage](#usage) for details), please check IRIS/config.py to ensure qsub parameters are correct before moving to the next step.__
##### 1.2 Download IRIS db 
IRIS loads a big-data reference database of splicing events and other genomic annotations. \
These data are included in [IRIS_data.tgz](https://drive.google.com/file/d/1TaswpWPnEd4TXst46jsa9XSMzLsbzjOQ/view?usp=sharing) (a Google Drive link; size ~10 GB). Users need to move this file to the IRIS folder for streamlined installation. 
##### 1.3 Download IEDB MHC I prediction tools
Download IEDB_MHC_I-X.XX.X.tar.gz from IEDB website (see [Dependencies](#dependencies)). Create a folder named 'IEDB' in the IRIS folder, then move the downloaded gz file to the 'IEDB' folder.
#### 2. Install and configure
Under the IRIS folder, do:
```
./install.sh.bak 
```
Follow instructions to finish the installation of conda, python and its dependencies, bedtools, the downloaded IEDB package, and the IRIS data and packages. 

### Usage
- For streamlined AS-derived target discovery, please follow [major modules](#streamlined-major-modules) and run the corresponding toy example.
- For customized pipeline development, please check [all modules](#individual-modules) of IRIS.

#### Individual modules
IRIS provides individual modules/steps, allowing users to build pipelines for their customized needs.\
For a description of each [module/step](IRIS_modules.md), including RNA-seq preprocessing, HLA typing, proteo-transcriptomic MS searching, visualization, etc., please click [here](IRIS_modules.md) or the subheader above.
```
usage: IRIS [-h] [--version]
            {formatting,screening,prediction,epitope_post,process_rnaseq,makeqsub_rmats,exp_matrix,indexing,translation,pep2epitope,screening_plot,seq2hla,parse_hla,ms_makedb,ms_search,ms_parse}
            ...

IRIS -- IRIS

positional arguments:
  {formatting,screening,prediction,epitope_post,process_rnaseq,makeqsub_rmats,exp_matrix,indexing,translation,pep2epitope,screening_plot,seq2hla,parse_hla,ms_makedb,ms_search,ms_parse}
    formatting          Formats AS matrices from rMATS, followed by indexing for IRIS
    screening           Screens AS-derived tumor antigens using big-data reference
    prediction          Predicts and annotates AS-derived TCR (pre-prediction) and CAR-T targets
    epitope_post        Post-prediction step to summarize predicted TCR targets
    process_rnaseq      Processes RNA-Seq FASTQ files to quantify gene expression and AS
    makeqsub_rmats      Makes qsub files for running rMATS-turbo 'prep' step
    exp_matrix          Makes a merged gene expression matrix from multiple cufflinks results
    indexing            Indexes AS matrices for IRIS
    translation         Translates AS junctions into junction peptides
    pep2epitope         Wrapper to run IEDB for peptide-HLA binding prediction
    screening_plot      Makes stacked/individual violin plots for list of AS events
    seq2hla             Wrapper to run seq2HLA for HLA typing using RNA-Seq
    parse_hla           Summarizes seq2HLA results of all input samples into matrices for IRIS use
    ms_makedb           Generates proteo-transcriptomic database for MS search
    ms_search           Wrapper to run MSGF+ for MS search
    ms_parse            Parses MS search results to generate tables of identified peptides

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit

For command line options of each sub-command, type: IRIS COMMAND -h
```


#### Streamlined major modules
The common use of IRIS immunotherapy target discovery comprises three major steps. For a quick test, see [Example](#example), in which a shell script is provided for a streamlined example run:
- Step 1. IRIS formatting (& indexing) 
```
usage: IRIS formatting [-h] -t {SE,RI,A3,A5} -n DATA_NAME -s {1,2}
                       [-c COV_CUTOFF] [-e] [-d IRIS_DB_PATH]
                       rmats_mat_path_manifest rmats_sample_order
```

- Step 2. IRIS screening (& translation) 
Here is a [description of the parameter file](example/paramenter_file_description.txt) and an [example file](example/Test.para).
```
usage: IRIS screening [-h] [-o OUTDIR] [-t] parameter_fin
```

- Step 3. IRIS prediction (predicts both extracellular targets and epitopes; __requires SGE system__)
```
usage: IRIS prediction [-h] [-p PARAMETER_FIN] [--iedb-local IEDB_LOCAL]
                       [-c DELTAPSI_COLUMN] [-d DELTAPSI_CUT_OFF] -m MHC_LIST
                       [--extracellular-anno-by-junction]
                       IRIS_screening_result_path
                          
usage: IRIS epitope_post [-h] -p PARAMETER_FIN -o OUTDIR -m MHC_BY_SAMPLE
                         [-e GENE_EXP_MATRIX] [--ic50-cut-off IC50_CUT_OFF]
```


### Example
We provide a wrapper ([run_example.sh](run_example.sh)) to run the above [IRIS streamlined major modules](#streamlined-major-modules) using [example files](example), included in the IRIS package. For customized pipeline development, we recommend that users use this script as a reference. Under the IRIS folder, do:
```
./run_example.sh
```
__As mentioned in [Usage](#usage), this example run will involve submitting the job array to the SGE system.__ It will take < 5 min for the formatting and screening steps and usually < 15 min for the prediction step (SGE job arrays).\
A successful test run will generate the following result files (row numbers are displayed before each file name):
```
      0 Glioma_test.notest.txt
     13 Glioma_test.primary.txt
      3 Glioma_test.primary.txt.ExtraCellularAS.txt
     11 Glioma_test.prioritized.txt
      3 Glioma_test.prioritized.txt.ExtraCellularAS.txt
     13 Glioma_test.test.all.txt
     13 primary/epitope_summary.junction-based.txt
     74 primary/epitope_summary.peptide-based.txt
    148 primary/pred_filtered.score500.txt 
     11 prioritized/epitope_summary.junction-based.txt
     45 prioritized/epitope_summary.peptide-based.txt
     84 prioritized/pred_filtered.score500.txt
```
__Users can refer to relative paths in the parameter file Test.para, the file manifest matrice.txt, and the file samples.txt. These relative paths were made for the example run. Users will need to change the path for their own analyses.__

### Example output
Final reports are shown in __bold__ font.

#### Screening results
[TASK/DATA_NAME].test.all.txt: All AS events tested by IRIS screening

[TASK/DATA_NAME].notest.txt: During screening, AS events skipped due to no variance or no available comparisons

[TASK/DATA_NAME].primary.txt: Tumor AS events after comparison to tissue-matched normal panel ('primary' events)

[TASK/DATA_NAME].prioritized.txt: Tumor AS events after comparison to tissue-matched normal panel, tumor panel, and normal tissue panel ('prioritized' AS events)

#### CAR-T annotation reports
__[TASK/DATA_NAME].primary.txt.ExtraCellularAS.txt__: Tumor AS events in 'primary' set that are associated with protein extracellular annotation and may be used for CAR-T targets

__[TASK/DATA_NAME].prioritized.txt.ExtraCellularAS.txt__: Tumor AS events in 'prioritized' set that are associated with protein extracellular annotation and may be used for CAR-T targets

#### TCR prediction reports
primary/pred_filtered.score500.txt: IEDB prediction outputs for AS junction peptides from 'primary' set with HLA-peptide binding IC50 values passing user-defined cut-off

__primary/epitope_summary.peptide-based.txt__: AS-derived epitopes from 'primary' set that are predicted to bind user-defined HLA type

__primary/epitope_summary.junction-based.txt__: Epitope-producing AS junctions from 'primary' set that are predicted to bind user-defined HLA type

prioritized/pred_filtered.score500.txt: IEDB prediction outputs for AS junction peptides from 'prioritized' set with HLA-peptide binding IC50 value passing user-defined cut-off

__prioritized/epitope_summary.peptide-based.txt__: AS-derived epitopes from 'prioritized' set that are predicted to bind user-defined HLA type

__prioritized/epitope_summary.junction-based.txt__: Epitope-producing AS junctions from 'prioritized' set that are predicted to bind user-defined HLA type 





 ### Contact
Yang Pan <panyang@ucla.edu>

Yi Xing <yxing@ucla.edu>
                     


### Publication
Manuscript in submission
 
