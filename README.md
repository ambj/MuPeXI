# MuPeXI #
------------------------------------------------


MuPeXI - Mutant peptide extractor and Informer

* Authors: Anne-Mette Bjerregaard and Aron C. Eklund 
* Date: September 2016
* version: 1.1 



## Introduction

Personalization of immunotherapies such as cancer vaccines and adoptive T cell therapy will depend on identification of patient-specific neo-epitopes that can be specifically targeted. MuPeXI, the Mutant Peptide Extractor and Informer, is a program that uses somatic mutation data, HLA binding prediction, self-similarity comparison and gene expression to identify and assess potential neo-epitopes as targets for immunotherapy. 

MuPeXI.py extracts peptides of user-defined lengths around missense variant mutations, indels and frameshifts. Information from each mutation is annotated together with the mutant and normal peptides in the file output.


### Algorithmic Flow Chart  

![](/doc/Mupexi_flow_chart.png)



## Dependencies  


To run MuPeXI the following software and packages must be installed:

#### Required software:
* [Python 2.7](https://www.python.org/download/releases/2.7/)
* [netMHCpan 2.8](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCpan)
* [Varaint Effect Predictor (VEP)](http://www.ensembl.org/info/docs/tools/vep/index.html) 

#### Python modules:
* [Biopython](http://biopython.org/wiki/Download)
* [numpy](http://www.numpy.org/)
* [pandas](http://pandas.pydata.org/)

These packages are included if downloading python through [Anaconda](https://www.continuum.io/downloads); the full anaconda package list can be found [here](https://docs.continuum.io/anaconda/pkg-docs). The Python 2.7 standard library list can be found [here](https://docs.python.org/2.7/library/).




## Installation  


1. Install all software listed above. 
2. Install Biopython, pandas and nympy with the following command:

        pip install biopython
        pip install numpy
        pip install pandas


3. Download or clone the MuPeXI repository to your local system

        git clone https://github.com/ambj/MuPeXI.git

4. If not already on your system, reference files from GRCh38 should be obtained. These include cDNA, peptide and cosmic references, look under references in the user manual for a detailed description.

5. Fill the config.ini file  
    * State the full path to netMHCpan 3.0 and VEP.
    * State the full path to the reference files:
        - cDNA
        - peptide
        - cosmic

   Additional peptide references and liftover paths can be stated in the config.ini file, see the user manual for detailed information. Instructions on how to fill the config.ini file is found within the file. `config.ini` is automatically detected in the same directory as `MuPeXI.py` script but can also be placed elsewhere and referred to by the `-c` option. 





## Usage  

After installation MuPeXI is called as follows. Here is an example where the reference files are stated in the config.ini file and netMHCpan-3.0 is run for the HLA types HLA-A01:01 and HLA-B08:01. The config file is specified using the `-c` option and the expression file with the `-e` option.

    path/to/MuPeXI.py -v mutation_call_file.vcf -a HLA-A01:01,HLA-B08:01 -c path/to/config.ini -e expression_file.tsv


MuPeXI can be used for both peptide extraction, giving immunogenicity information for peptide selection, and for printing of fasta file mutant-peptide-library for optimal mass spectrometry peptide detection. 

All options can be explored using the usage information with the `-h` option:   

    > path/to/MuPeXI.py -h

        MuPeXI - Mutant Peptide Extractor and Informer
        version 2016-08-11

        The current version of this program is available from
        https://github.com/ambj/MuPeXI

        MuPeXI.py accepts a VCF file describing somatic mutations as input, and from this 
        derives a set of mutated peptides of specified length(s). These mutated peptides 
        are returned in a table along with various annotations that may be useful for 
        predicting the immunogenicity of each peptide.

        Usage: ./MuPeXI.py -v <VCF-file> [options]

                                                                                    DEFAULT
        Required arguments:
        -v, --vcf-file <file>   VCF file of variant calls, preferably from
                                MuTect (only SNVs) or MuTect2 (SNVs and indels)

        Recommended arguments:
        -a, --alleles           HLA alleles, comma separated.                       HLA-A02:01
        -l, --length            Peptide length, given as single number,             9
                                range (9-11) or comma separated (9,10,11).
        -e, --expression-file   Expression file, tab separated
                                ((ENST*/ENSG*)   mean)

        Optional arguments affecting output files:
        -o, --output-file       Output file name.                                   <VEP-file>.mupexi
        -d, --out-dir           Output directory - full path.                       current directory
        -p, --prefix            Prefix for output files - will be applied           <VCF-file>
                                to all (.mupexi, .log, .fasta) unless specified 
                                by -o or -L.
        -L, --log-file          Logfile name.                                       <VCF-file>.log
        -m, --mismatch-number   Maximum number of mismatches to search for in       4
                                normal peptide match. 

        Other options (these do not take values)
        -f, --make-fasta        Create FASTA file with long peptides 
                                - mutation in the middle
        -c, --config-file       Path to the config.ini file                         current directory
        -t, --keep-temp         Retain all temporary files
        -M, --mismatch-only     Print only mismatches in normal peptide sequence 
                                and otherwise use dots (...AA.....)
        -w, --webface           Run in webserver mode
        -g, --hg19              Perform liftover HG19 to GRCh38.
                                Requires local picard installation with paths
                                stated in the config file
        -E, --expression-type   Setting if the expression values in the expression  transcript
                                files are determined on transcript or gene level.
                                (transcript/gene)
        -h, --help              Print this help information

        REMEMBER to state references in the config.ini file



## User Manual 
For detailed information about usage, input and output files, test examples and data preparation read the [MuPeXI User Manual](/doc/MuPeXI_User_Manual.md)


## Contact   

Anne-Mette Bjerregaard  
ambj@cbs.dtu.dk

or 

Aron Charles Eklund  
eklund@cbs.dtu.dk


Department of Bio and Health Informatics  
Technical University of Denmark  
Kemitorvet 208, 2800 Lyngby, Denmark  