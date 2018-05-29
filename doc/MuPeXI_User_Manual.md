# MuPeXI User Manual

Anne-Mette Bjerregaard  
Department of Bio and Health Informatics  
Technical University of Denmark  
Lyngby, Denmark.


## Table of Contents
1. [General Description](#general-description)  
2. [Dependencies](#dependencies)  
    - [Required software](#required-software)  
    - [Python packages](#python-packages)  
3. [Installation](#installation)  
4. [Usage](#usage)  
5. [Input Files](#input-files)  
    - [VCF file](#vcf-file)  
    - [Expression file](#expression-file)  
    - [References](#references)  
6. [Output Files](#output-files)  
    - [Column explanation](#column-explanation)  
7. [Test Example](#test-example)  
    - [Data preparation](#data-preparation)  
        * [Recommended preprocessing of next generation sequencing (NGS)
data](#recommended-preprocessing-of-next-generation-sequencing-(ngs)-data)  
        * [Data cleanup](#data-cleanup)
        * [WXS data](#wxs-data)
        * [RNAseq](#rnaseq)
        * [HLA typing](#hla-typing)  




## General Description

Personalization of immunotherapies such as cancer vaccines and adoptive T cell therapy
will depend on identification of patient-specific neo-epitopes that can be specifically
targeted. MuPeXI, the Mutant Peptide Extractor and Informer, is a program that uses
somatic mutation data, HLA binding prediction, self-similarity comparison and gene
expression to identify and assess potential neo-epitopes as targets for immunotherapy. 
`MuPeXI.py` extracts peptides of user-defined lengths around missense variant mutations,
indels and frameshifts. Information from each mutation is annotated together with the
mutant and normal peptides in the file output.


## Dependencies  


To run MuPeXI the following software and packages must be installed:

#### Required software:
* [Python 2.7](https://www.python.org/download/releases/2.7/)
* [netMHCpan 4.0](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCpan)
* [Varaint Effect Predictor (VEP)](http://www.ensembl.org/info/docs/tools/vep/index.html) 

#### Python modules:
* [Biopython](http://biopython.org/wiki/Download)
* [numpy](http://www.numpy.org/)
* [pandas](http://pandas.pydata.org/)

These packages are included if downloading python through
[Anaconda](https://www.continuum.io/downloads); the full anaconda package list can be
found [here](https://docs.continuum.io/anaconda/pkg-docs). The Python 2.7 standard
library list can be found [here](https://docs.python.org/2.7/library/). If installing
Anaconda it should not be necessary to do the second step of the installation guide. 



## Installation  


1. Install all software listed above. 
2. Install Biopython, pandas and nympy with the following command:

        pip install biopython
        pip install numpy
        pip install pandas


3. Download or clone the MuPeXI repository to your local system

        git clone https://github.com/ambj/MuPeXI.git

4. If not already on your system, reference files from GRCh38 should be obtained. These
include cDNA, peptide and cosmic references, look under References (below)
for a detailed description.

5. Fill in the config.ini file  
    Instructions on how to fill in the file are found within the file. The `config.ini` file
    is automatically detected in the same directory as `MuPeXI.py` script but can also be
    placed elsewhere and referred to by the `-c` option.

    State the path to netMHCpan and VEP in the config.ini file.

        [netMHC]
        # Please specify the netMHCpan binary path
        MHC = your/path/to/netMHCpan-4.0/netMHCpan

        [EnsemblVEP]
        # Please specify binary path to Ensembls Variant effect predictor (VEP), and the directory of containing the cache database. 
        VEP = your/path/to/ensembl-tools-release-85/scripts/variant_effect_predictor/variant_effect_predictor.pl
        VEPdir = your/path/to/ensembl-tools-release-85/scripts/variant_effect_predictor

    You should be aware that the version of VEP library you use should match the
    references used (peptide and cDNA).  E.g. in the example above is used version/release
    85 of GRCh38. 
    The path to reference files should also be stated in the config.ini file, together
    with FULL path to the small C script named "pepmatch_db" downloaded with MuPeXI. 
    All reference peptides of the defined length, extracted from the proteome reference,
    can be defined as a file in the config.ini file. This will save time, avoiding
    chop-up of the proteome reference for each run. As many files as peptides length
    should be stated, if not stated MuPeXI runs the chop-up for the peptide lengths where
    no peptide reference file is provided. 

        [References]
        # Please specify the binary path to the references used (optional)
        cDNA = your/path/to/human_GRCh38/cDNA/Homo_sapiens.GRCh38.78.cdna.all.fa
        pep = your/path/to/human_GRCh38/pep/Homo_sapiens.GRCh38.78.pep.all.fa
        cosmic = your/path/to/cosmic/Census_allWed_Feb_17_09-33-40_2016.tsv
        pep9 = your/path/to/reference_peptide_9.txt

        [PeptideMatch]
        PM = your/path/to/MuPeXI/apps/pepmatch_db

    If the user wishes for MuPeXI to do a lift-over from GRCh37 / HG19 to GRCh38 this is
    possible if a local version of picard tools 2.5 is installed. The paths to java8 and
    picard-tools should also be stated in the config.ini file.

        [LiftOver]
        fasta = your/path/to/references/GRCh38.fa
        chain = your/path/to/references/HG19toGRCh38.chain
        java8 = your/path/to/java8
        picard = your/path/to/picard-tools-2.5.0




## Usage  

After installation MuPeXI is called as follows. Here is an example where the reference
files are stated in the config.ini file and netMHCpan is run for the HLA types
HLA-A01:01 and HLA-B08:01. The config file is specified using the `-c` option and the
expression file with the `-e` option.

    path/to/MuPeXI.py -v mutation_call_file.vcf -a HLA-A01:01,HLA-B08:01 -c path/to/config.ini -e expression_file.tsv


MuPeXI can be used for both peptide extraction, giving immunogenicity information for
peptide selection, and for printing of fasta file mutant-peptide-library for optimal mass
spectrometry peptide detection. 

All options can be explored using the usage information with the `-h` option:   

    > path/to/MuPeXI.py -h

        MuPeXI - Mutant Peptide Extractor and Informer

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
        -E, --expression-type   Are the expression values in the expression         transcript
                                files determined on transcript or gene level.
                                (transcript/gene)

        Optional arguments affecting output files:
        -o, --output-file       Output file name.                                   <VEP-file>.mupexi
        -d, --out-dir           Output directory - full path.                       current directory
        -p, --prefix            Prefix for output files - will be applied           <VCF-file>
                                to all (.mupexi, .log, .fasta) unless specified 
                                by -o or -L.
        -L, --log-file          Logfile name.                                       <VCF-file>.log
        -m, --mismatch-number   Maximum number of mismatches to search for in       4
                                normal peptide match.
        -A, --assembly          The assembly version to run VEP.                    GRCh38
        -s, --species           Species to analyze (human / mouse / mouse_black6 /  human
                                mouse_balbc)
                                If mouse is set default assembly is GRCm38.
                                Remember to download the correct VEP cache
                                and state species corresponding MHC alleles.

        Optional arguments affecting computational process:
        -F, --fork              Number of processors running VEP.                   2
        
	Other options (these do not take values)
        -f, --make-fasta        Create FASTA file with long peptides 
                                - mutation in the middle
        -c, --config-file       Path to the config.ini file                         current directory
        -t, --keep-temp         Retain all temporary files
        -M, --mismatch-only     Print only mismatches in normal peptide sequence 
                                and otherwise use dots (...AA.....)
        -w, --webface           Run in webserver mode
        -g, --liftover          Perform liftover HG19 to GRCh38.
                                Requires local picard installation with paths
                                stated in the config file
        -n, --netmhc-full-anal  Run NetMHCpan 4.0 with the full analysis including 
                                both eluted ligand (EL) and binding affinity (BA) 
                                prediction output (priority score calculated from 
                                EL rank score)
        -h, --help              Print this help information

        REMEMBER to state references in the config.ini file



## Input Files 

MuPeXI accepts a VCF file of somatic mutation calls optimally obtained from either MuTect
or MuTect2. The VCF files do not need to be modified; the "raw" output VCF file can
be put directly into MuPeXI. 

### VCF file 

Compact example of a VCF file:

        ##fileformat=VCFv4.2
        ##GATKCommandLine.MuTect2=<ID=MuTect2,Version=3.5-0-g36282e4
        ##SAMPLE=<ID=NORMAL,SampleName=TCGA-XV-A9W5_N
        ##SAMPLE=<ID=TUMOR,SampleName=TCGA-XV-A9W5_T
        ##reference=file:///home/projects/pr_46630/data/references/human_GRCh38/GCA_000001405.15_GRCh38_full_analysis_set.fa
        #CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT  TUMOR   NORMAL
        chr1    948711  .   C   G   .   germline_risk . . . .
        chr1    1657358 .   T   TA  .   alt_allele_in_normal . . . .
        chr1    1986752 rs4233028   A   G   .   germline_risk . . . .   
        chr1    3431704 rs2493274   G   C   .   germline_risk . . . .  
        chr1    3631978 rs2244942   T   C   .   germline_risk . . . .  
        chr1    3839305 rs1891940   T   C   .   clustered_events;germline_risk  . . . .

A full example of a VCF file can be found on the MuPeXI webserver
[here](http://www.cbs.dtu.dk/services/MuPeXI/files/example.vcf). 

### Expression file

It is optional, but preferable, to provide a file with expression values as input to add
the expression of each transcript where a mutated peptide was extracted.
The expression files used for testing MuPeXI were generated from raw RNA-seq data using
Kallisto. The files should be tab separated and include Ensembl transcript ID (ENST) and
mean expression. 

        ENST00000456328.2   0.868567715
        ENST00000450305.2   0
        ENST00000488147.1   2.72373575
        ENST00000619216.1   0
        ENST00000473358.1   0
        ENST00000469289.1   0
        ENST00000607096.1   0

A full example of an expression file can be found on the MuPeXI webserver
[here](http://www.cbs.dtu.dk/services/MuPeXI/files/example_expression.tsv).

It should be noted that MuPeXI takes both expression values determined on transcript and
gene level, though transcript is preferable. If gene level is used (ENSG...) the `-E gene`
option should be used. 

### References 
The following references are required for MuPeXI to run:
* Peptide  
    The peptide reference is a FASTA file containing all peptides of the human proteome.
* cDNA  
    The cDNA reference is a FASTA file containing all cDNA sequences of the human
    proteome.  

These references can be acquired from the
[Ensembl website](http://www.ensembl.org/Homo_sapiens/Info/Index).  
The most recent release is found under Gene annotation > Download genes, cDNAs, ncRNA,
protein (FASTA) > pep (for peptide reference) and > cdna (for cDNA reference).  
It should be emphasized that it is of very high importance that the references and VEP
match in release version (e.g. release-85).


The following reference are optional but preferable:
* Cosmic  
        TSV file containing known cancer driver genes. The cancer gene census can be
        downloaded from the [COSMIC](http://cancer.sanger.ac.uk/census) website.  

#### Murine Specific References 
We used the following murine specific references for our data analysis of BALBc and C57BL/6 mice strains. 

| Reference for 	| Species     	| Description | Link |
| -----------           | ----------- 	| ----------- | ----------- |
| NGS Analysis 		| C57BL/6     	| SNP from Sanger mouse project | ftp://ftp-mouse.sanger.ac.uk/current_snps/strain_specific_vcfs/C57BL_6NJ.mgp.v5.snps.dbSNP142.vcf.gz |
| NGS Analysis  	| BALBc		| SNP from Sanger mouse project | ftp://ftp-mouse.sanger.ac.uk/current_snps/strain_specific_vcfs/BALB_cJ.mgp.v5.snps.dbSNP142.vcf.gz |
| NGS Analysis          | C57BL/6       | Indels from Sanger mouse project | ftp://ftp-mouse.sanger.ac.uk/current_snps/strain_specific_vcfs/C57BL_10J.mgp.v5.indels.dbSNP142.normed.vcf.gz |
| NGS Analysis          | BALBc       	| Indels from Sanger mouse project | ftp://ftp-mouse.sanger.ac.uk/current_snps/strain_specific_vcfs/BALB_cJ.mgp.v5.indels.dbSNP142.normed.vcf.gz |
| NGS Analysis          | C57BL/6       | Ensembl genome DNA reference | ftp://ftp.ensembl.org/pub/release-89/fasta/mus_musculus_c57bl6nj/dna/ |
| NGS Analysis          | BALBc       | Ensembl genome DNA reference | ftp://ftp.ensembl.org/pub/release-89/fasta/mus_musculus_balbcj/dna/ |
| MuPeXI 		| C57BL/6 	| cDNA from Ensembl | ftp://ftp.ensembl.org/pub/release-89/fasta/mus_musculus_c57bl6nj/cdna/ |
| MuPeXI		| BALBc 	| cDNA from Ensembl | ftp://ftp.ensembl.org/pub/release-89/fasta/mus_musculus_balbcj/cdna/ |
| MuPeXI                | C57BL/6       | Peptide from Ensembl | ftp://ftp.ensembl.org/pub/release-89/fasta/mus_musculus_c57bl6nj/pep/ |
| MuPeXI                | BALBc         | Peptide from Ensembl | ftp://ftp.ensembl.org/pub/release-89/fasta/mus_musculus_balbcj/pep/ |


## Output Files 
MuPeXI can output up to three files.

The output files are the following: 
1.  .mupexi 
    The main output file containing a TSV file with the extracted mutated peptides and
    all the information needed to choose the wanted peptides. 
2.  . log
    The log file containing information about the run. 
3.  .fasta
    The FASTA file with the long peptides, before they are chopped up into user defined
    length, with the mutation centered. (use the option `-f`)

### Column explanation

The prediction output (.mupexi) for each peptide pair consists of the following columns:

| Column Name           | Description |
| -----------           | ----------- |
| HLA allele            | Allele name |
| Normal peptide        | Peptide from reference corresponding to the mutant peptide. |
| Normal MHC affinity   | Predicted binding affinity of normal peptide in nanoMolar units. |
| Normal MHC % rank     | %Rank of prediction score for nomal peptides. |
| Mutant peptide        | The extracted mutant peptide. |
| Mutant MHC affinity   | Predicted binding affinity of mutant peptide in nanoMolar units. |
| Mutant MHC % rank     | %Rank of prediction score for mutant peptides. |
| Gene ID               | Ensembl gene ID |
| Transcript ID         | Ensembl transcript ID |
| Amino acid change     | Amino acid change annotated in VEP file. |
| Allele Frequency      | Genomic allele frequency detected by MuTect2. |
| Mismatches            | Mismatches between normal and mutant peptide. |
| Peptide position      | Position of amino acid change in the peptide. Can be a range in the case of insertions and frameshifts. |
| Chr                   | Chromosome position annotated in the VEP file. |
| Genomic position      | Genome nucleotide position annotated in the VEP file. |
| Protein position      | Amino acid position annotated in the VEP file. |
| Mutation cons.        | The consequence annotated in the VEP file translated into single letter abbreviations: M; Missense variant, I; In-frame insertion, D; In-frame deletion, F; Frameshift variant |
| Gene symbol           | HUGO symbol corresponding to the Ensembl transcript id. |
| Cancer driver gene    | Yes if the HUGO symbol is in the cosmic reference list, No if it is not. |
| Expression Level      | Expression of the transcript which the mutant peptide was extracted from. |
| Mutant affinity score | Calculated binding affinity score of the mutant peptide, based on a negative logistic function of the mutant MHC %Rank score. This is used to calculate the final prioritization score. |
| Normal affinity score | Calculated binding affinity score of the normal peptide, based on a positive logistic function of the normal MHC %Rank score. This is used to calculate the final prioritization score. |
| Expression score      | Calculated Expression score of the transcript expression level. This is used to calculate the final prioritization score. |
| Priority score        | Calculated prioritization dependent on HLA binding, gene expression, normal and mutant peptide binding ratio and allele frequency. |

NetMHCpan output:
%Rank of prediction score to a set of 200.000 random natural 9mer peptides. For more
information go to NetMHCpan
[output format](http://www.cbs.dtu.dk/services/NetMHCpan/output.php).



## Test example 

To run the provided test files with MuPeXI the following command can be run: 

        path/to/MuPeXI.py -v test.vcf -c path/to/config.ini -e expression_test.tsv

For additional fasta file output:

        path/to/MuPeXI.py -v test.vcf -c path/to/config.ini -e expression_test.tsv -f

Print only the mismatch amino acid for the normal peptide:

        path/to/MuPeXI.py -v test.vcf -c path/to/config.ini -e expression_test.tsv -m

## Data preparation

Raw sequencing data should be obtained and the following steps is a recommendation on how
to processes this data prior to using MuPeXI. 

### Recommended preprocessing of next generation sequencing (NGS) data  
The following data should be obtained:

* WXS (or WGS) of tumor and corresponding normal sample 
* RNAseq of tumor sample  

Raw sequencing data is processed to extract variant calls, and optionally gene expression
values.  Therefore, either whole exome (WXS) or whole genome (WGS) sequencing data is
required from both tumor and matching normal sample. For optimal utilization of MuPeXI,
RNA sequencing (RNAseq) data should be obtained as well from the tumor sample. 

#### Data cleanup  
First, it is important to check the quality of both WXS and RNAseq data. This can be done
by running the wrapper tool
[Trim Galore](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/), which
combines the functions of [Cutadapt](http://cutadapt.readthedocs.org/en/stable/index.html)
and [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/): trimming the
reads below an average phred score of 20 (default value), cutting out standard adaptors,
such as those from Illumina, and in the end running FastQC to display the FastQC report.
This enables the user to check if the quality is satisfactory and if non-standard
adapters are seen in the overrepresented sequences. 

#### WXS data  
For the WXS pipeline the preprocessing steps for variant calling should be done following
the 
[GATK best practice guidelines](https://software.broadinstitute.org/gatk/best-practices/).
The reads are aligned using default [BWA](http://bio-bwa.sourceforge.net/) mem options;
here it is important to provide the Read Group for the following steps to be possible.
Then indel realignment and base recalibration are done to reduce false positive variant
calls; a detailed description of these steps can be found on the GATK website.  
Single nucleotide variant (SNV) and indel calling should be done using
[MuTect2](https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_cancer_m2_MuTect2.php).
MuTect2 is designed to call variants, both SNVs and indels, from tumor samples evaluating
reads from matched tumor and normal samples.  
The standard output format from MuTect2 is a VCF file.  
If MuTect or Mutect2 is not used for variant calling, the genomic allele frequency will
not be taken into account in the MuPeXI priority score. 

#### RNAseq  
After quality check the expression scores are obtained by pseudo-aligning the raw reads
using [Kallisto](https://pachterlab.github.io/kallisto/). Since this is very fast, the
runs are bootstrapped (e.g. 500 times) giving a variance of each expression score
obtained. The expression score is normalized by Kallisto in transcripts per million (TPM).

#### HLA typing  
If the HLA alleles for the sample(s) have not been typed by other means, they can be
inferred from sequencing data by several different algorithms. Among these are
[ATHLATES](https://www.broadinstitute.org/scientific-community/science/projects/viral-genomics/athlates),
[Polysolver](https://www.broadinstitute.org/cancer/cga/polysolver),
[OptiType](https://github.com/FRED-2/OptiType) and
[PHLAT](https://sites.google.com/site/phlatfortype/). 
We used OptiType with default settings, first filtering the reads with
[RazerS3](http://www.seqan.de/projects/razers/).  
The outputs files from the NGS analysis; VCF file, HLA alleles and expression file, are
now ready to be processed by MuPeXI.



![](/doc/NGS_analysis.png)  
**Pre-processing of raw sequencing data**. Flowchart of an example NGS pipeline to
process raw sequencing data in preparation to predict neo-epitopes with MuPeXI.
BWA: Burrows-Wheeler Aligner. GATK: Genome Analysis Toolkit. WXS: Whole exome sequencing.  

