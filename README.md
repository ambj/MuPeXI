# MuPeXI: Mutant Peptide eXtractor and Informer #

Given a list of somatic mutations (VCF file) as input, MuPeXI returns a table containing
all mutated peptides (neo-peptides) of user-defined lengths, along with several pieces
of information relevant for identifying which of these neo-peptides are likely to serve as
neo-epitopes. 

NEW: MuPeXI is now tested and compatible for suquencing data of murine origin.

#### Authors: 
Anne-Mette Bjerregaard and Aron C. Eklund 

#### License: 
To be determined, but certainly free for academic use.

#### Citation:
Bjerregaard AM, Nielsen M, Hadrup SR, Szallasi Z, Eklund AC.  
MuPeXI: Prediction of neo-epitopes from tumor sequencing data.  
Cancer Immunol Immunother. 2017 Apr 20. doi: 10.1007/s00262-017-2001-3. [Epub ahead of print]  
PubMed ID: [28429069](https://www.ncbi.nlm.nih.gov/pubmed/28429069)  
You can read the paper here: http://rdcu.be/rwVP

#### Web servers:
For limited data, MuPeXI can be run on our
[web server](http://www.cbs.dtu.dk/services/MuPeXI/)

A mouse specific [web server](http://www.cbs.dtu.dk/services/MuPeXI-mouse/)

## Dependencies

#### Hardware:
MuPeXI currently runs only on x86_64 machines running Linux or Darwin.

#### Required software:
* [Python 2.7](https://www.python.org/download/releases/2.7/)
* [NetMHCpan 4.0](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCpan)
* [Variant Effect Predictor (VEP)](http://www.ensembl.org/info/docs/tools/vep/index.html) 

#### Required Python modules:
* [Biopython](http://biopython.org/wiki/Download)
* [numpy](http://www.numpy.org/)
* [pandas](http://pandas.pydata.org/)

Note: These modules are already included if using Python through
[Anaconda](https://www.continuum.io/downloads).

#### Optional software, required only for liftover from HG19
* [Picard tools](https://broadinstitute.github.io/picard/)
* [Java 8](https://java.com/en/download/help/linux_x64rpm_install.xml)


## Installation  

1. Install all software listed above.

2. Download or clone the MuPeXI repository to your local system

        git clone https://github.com/ambj/MuPeXI.git

3. Obtain the reference files from GRCh38. These include cDNA, peptide and COSMIC
files; see the References section in the [user manual](/doc/MuPeXI_User_Manual.md)
for a detailed description.

4. Fill in the config.ini file  
    * Provide the full path to NetMHCpan and VEP.
    * Provide the full path to the reference files:
        - cDNA
        - peptide
        - COSMIC

   Additional peptide references and liftover paths can be provided in the config.ini
   file; see the user manual for detailed information. Instructions on how to fill in 
   the config.ini file are found within the file. `config.ini` is automatically found if 
   it is in the same directory as `MuPeXI.py` script, but it can also be placed elsewhere
   and specified by the `-c` option. 


## Usage  

Here is a simple example in which somatic mutation calls and gene expression data are
provided, and MHC binding is predicted for HLA types HLA-A01:01 and HLA-B08:01. 

    path/to/MuPeXI.py -v mutations.vcf -a HLA-A01:01,HLA-B08:01 -e expression.tsv

MuPeXI can be used for both peptide extraction, giving immunogenicity information for
peptide selection (the default), and for generation of a FASTA-formatted mutant-peptide
file suitable for input to mass spectrometry peptide search software (with the `-f` 
option). 

All options can be displayed using the usage information with the `-h` option:   

    path/to/MuPeXI.py -h


## User Manual 
For detailed information about usage, input and output files, test examples and data
preparation read the [MuPeXI User Manual](/doc/MuPeXI_User_Manual.md)


## FAQs

* We pronounce it moo-PECKS-ee


## Contact   

Anne-Mette Bjerregaard  
ambj@bioinformatics.dtu.dk

or 

Aron Charles Eklund  
eklund@bioinformatics.dtu.dk


Department of Bio and Health Informatics  
Technical University of Denmark  
http://www.bioinformatics.dtu.dk/english


## Algorithmic Flow Chart  

![](/doc/Mupexi_flow_chart.png)
