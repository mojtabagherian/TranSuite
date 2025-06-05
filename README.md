# TranSuite

## Original Version
A software suite for accurate identification, annotation, translation, and feature 
characterization of annotate transcripts.
Reference: 
* Juan C Entizne, Wenbin Guo, Cristiane P.G. Calixto, Mark Spensley, Nikoleta 
Tzioutziou, Runxuan Zhang, and John W. S. Brown; [*TranSuite: a software suite for 
accurate translation and characterization of transcripts*](https://www.biorxiv.org/
content/10.1101/2020.12.15.422989v1) 

Original TranSuite repository: https://github.com/anonconda/TranSuite

----------------------------
## Modified TranSuite
----------------------------
A modified version of TranSuite with additional features for comprehensive transcript analysis.


Modified by:
- Mojtaba Bagherian
- James Lloyd

----------------------------
Overview of Features
----------------------------
This modified TranSuite builds upon the original TranSuite software with the following changes:

### Bug Fixes
- **Chimeric Gene Handling**: Fixed an issue with the `--chimeric` parameter encountered when using the Araport 11 gtf file handling across Findlorf and Transfeat modules. The parameter is now consistently handled with a default value of `None`, preventing AttributeError exceptions when the parameter is not provided. This ensures more robust processing of chimeric genes across all TranSuite modules.

## Modified Splice Junction Analysis
Precise measurements of splice junction distances are essential for understanding nonsense-mediated decay (NMD) mechanisms and transcript regulation. These calculations measure exonic distances between stop codons and splice junctions to predict the fate of transcripts.

### Core Functions

#### 1. Exon Distance Calculation
Calculates the precise distance between genomic positions along the mature transcript, excluding intronic sequences. This biological accuracy is critical for NMD analysis, where distances are measured on processed mRNA.

#### 2. Upstream Splice Junction Distance (UpstreamEJ)
Measures the exonic distance from the nearest upstream exon-exon junction to the stop codon. Identifies the closest splice junction before the stop codon and calculates the distance along the mature transcript. This measurement is important for understanding exon junction complex (EJC) positioning and its regulatory effects.

#### 3. Downstream Splice Junction Distance (DownstreamEJ)
Calculates the exonic distance from the stop codon to the farthest downstream junction in the 3' UTR. This is the primary NMD determinant - transcripts with stop codons >50 nucleotides upstream of a downstream exon junction are typically NMD substrates. The measurement focuses on the most distant junction for conservative NMD prediction.

### Technical Implementation
The modular design includes:

- `calculate_exon_distance()`: Core distance calculation excluding introns
- `calculate_utr_lengths()`: Accurate 3' UTR length determination
- `calculate_splice_junction_distances()`: Comprehensive USJ/DSJ analysis
- `write_splice_junction_data()`: Structured output generation

All functions are strand-aware and handle edge cases robustly.

### Output Enhancement
New file: `<outname>_splice_junctions.csv`
Contains: T_ID, Strand, UpstreamEJ, DownstreamEJ, 3UTRlength, PTC_dEJ

The PTC_dEJ column provides binary NMD classification (Yes/No for ≥50nt threshold), enabling both automated analysis and custom threshold studies.

### Biological Impact
These measurements enable:

- Accurate NMD prediction based on the established 50-nucleotide rule
- Alternative splicing impact assessment on transcript stability
- Transcript quality evaluation for identifying aberrant isoforms
- Post-transcriptional regulation studies with quantitative metrics

The modified splice junction analysis makes TranSuite a comprehensive tool for understanding how transcript structure influences gene expression regulation through NMD and related quality control mechanisms.


----------------------------
# Table of Contents
----------------------------

   * [Overview](#overview)
   * [Installation](#installation)
   * [Modules](#modules)
   * [Input files](#input-files)
   * [FindLORF](#findlorf)
      * [Command and options](#command-and-options)
      * [Output files](#output-files)
   * [TransFix](#transfix)
      * [Command and options](#command-and-options-1)
      * [Output files](#output-files-1)
   * [TransFeat](#transfeat)
      * [Command and options](#command-and-options-2)
      * [Output files](#output-files-2)
   * [Auto](#auto)
      * [Command and options](#command-and-options-3)
      * [Output files](#output-files-3)
   * [Contact](#contact)
   * [Citation](#citation)
   * [License](#license)


----------------------------
# Overview
----------------------------

TranSuite is a software for identifying coding sequences of transcripts, selecting translation start sites at gene-level, generating accurate translations of transcript isoforms, and identifying and characterizing multiple coding related features, such as: coding potential, similar-translation features, and multiple NMD-related signals. TranSuite consists of three independent modules, FindLORF, TransFix and TransFeat. Each module can be run independently or as a [pipeline with a single command](#auto).

- FindLORFS - finds and annotates the longest ORF of each transcript
- TransFix - "fixes" the same translation start codon AUG in all the transcripts in a gene and re-annotates the resulting ORF of the transcripts
- TransFeat - identifies structural features of transcripts, coding potential and NMD signals
- Auto - executes the whole pipeline (FindLORFS, TransFix, and TransFeat) in tandem


----------------------------
# Installation
----------------------------

TranSuite has been developed in Python 3.6 

TranSuite requires the following packages:
- [BioPython v1.78](https://anaconda.org/anaconda/biopython)


Packages installations commands:
```
conda install -c anaconda biopython
```

TranSuite is ready to use. The compressed file can be directly downloaded from the [GitHub repository](https://github.com/anonconda). Once uncompressed, TranSuite can be used directly from the command line by specifying the path to the main executable `transuite.py`

Additionally, in the future it will also be possible to install TranSuite through popular Python installation managers PyPI and Anaconda:


----------------------------
# Modules
----------------------------

TranSuite works with a module / module-options structure:

```
transuite.py module options
```
where the modules are:

- **FindLORF**    : Finds and annotates the longest ORF of each transcript.
- **TransFix**       : Fix the same translation start codon AUG in all the transcripts in a gene and re-annotates the resulting ORF of the transcripts.
- **TransFeat**       : Identify structural features of transcripts, coding potential and NMD signals.
- **Auto**        : Execute the whole pipeline (FindLORFS, TransFix, and TransFeat) in tandem.

Each module can be executed as follow:

```
python /path/to/transuite.py module options
```

For example, to observe the help documentation of Auto module:

```
python /path/to/transuite.py Auto --help
```


----------------------------
# Input files
----------------------------

All of TranSuite modules (FindLORF, TransFix and TransFeat) use the same input files format:

- The transcriptome annotation to analyze, in [GTF format](https://www.ensembl.org/info/website/upload/gff.html)
- The transcripts nucleotide sequences, in [FASTA format](http://bioinformatics.org/annhyb/examples/seq_fasta.html)

Please note that the programs assume that nucleotide sequences in the FASTA file represent the exonic region of the transcript. Please beware that errors will happen if the user provide a FASTA file of the coding-region (CDS) sequences instead. 


Additional notes:

- FindLORF will parse only the "exon" [feature](https://www.ensembl.org/info/website/upload/gff.html) information from the GTF file
- TransFix and TransFeat will extract "CDS" [feature](https://www.ensembl.org/info/website/upload/gff.html) information from the GTF file
- Transcripts without annotated "CDS" features will remain unprocessed by TransFix, and they will be identified as "No ORF" by TransFeat
- When executing the whole pipeline (Auto module), TranSuite will automatically forward the appropiate GTF file to each module. That is: FindLORF output GTF will be use as input by TransFix, and TransFix output GTF will be TransFeat input


----------------------------
**FindLORF**
==============

----------------------------

FindLORF identifies and annotate ORF information in newly transcriptome annotations. Firstly, FindLORF translates each transcript sequence in its three frames of translation according to its annotated strand and stores the relative start and stop codon positions of all the resulting ORFs. Secondly, FindLORF selects the longest ORF for each transcript as its putative CDS region. Finally, FindLORF annotates the CDS using the genomic information contained in the transcriptome annotation to convert the relative ORF start-stop codon positions in the transcript sequence into genomic co-ordinates. The FindLORF module takes as input the transcriptome annotation to be curated (GTF format) and the transcripts exon sequences (FASTA format). See [Input files](#input-files) above.

## Command and options
Command to run FindLORF:
```
python transuite.py FindLORF [options]
```
```
python transuite.py FindLORF --gtf <input-gtf.gtf> --fasta <input-fasta.fa> --outpath </path/for/output-folder> --outname <outname> --cds <30>
```

List of options available:
- **--gtf**: Transcriptome annotation file in GTF format
- **--fasta**: Transcripts nucleotide fasta file
- **--outpath**: Path of the output folder
- **--outname**: Prefix for the output files
- **--cds**: Minimum number of amino-acids an ORF must have to be considered as a potential CDS. Default: 30 AA

Example:
```
python transuite.py FindLORF --gtf ./test_dataset/subset_AtRTD2.gtf --fasta ./test_dataset/subset_AtRTD2_transcripts.fa --outpath ./test_dataset/test_output --outname  --cds 30
```

## Output files
FindLORF automatically generates a subfolder to store the output files:
/**&lt;outpath&gt;**/**&lt;outname&gt;**\_longorf/

FindLORF generates the following output files:
1. *GTF* file with the longest ORF in the transcripts annotated as its CDS
2. *FASTA* files of the transcripts CDS regions (nucleotide, and peptide)
3. Log *CSV* files reporting transcripts that could not be annotated, for example: due to lack of an AUG
4. A *JSON* file containing the transcripts ORF relative coordinates


-------------------
**TransFix**
==============
-------------------

TransFix provides more biologically accurate translations by selecting the authentic translation start site for a gene, "fixing" this location and using it to translate the gene transcripts and annotating the resulting CDS of the translations. We define the authentic translation start site as the site used to produce the full-length protein of the gene. In detail, TransFix firstly extracts the CDS co-ordinates of the transcripts from the transcriptome annotation and groups the transcripts according to their gene of origin. Then, TransFix selects the start codon of the longest annotated CDS in the gene as the representative translation start site and translates all of the transcripts in the gene from the "fixed" translation start site. Finally, TransFix annotates the genomic co-ordinates of the resulting stop codons. In some cases, transcript isoforms do not contain the "fixed" translation start site due to an AS event or an alternative transcription start site. To account for this, TransFix tracks those transcripts that are not translated during the first fix AUG/translation cycle and they are then processed through a second fix AUG/translation cycle to determine and annotate their valid translation start-sites.

### **Command and options** ###
Command to run TransFix:
```
python transuite.py TransFix [options]
```
```
python transuite.py TransFix --gtf <input-gtf.gtf> --fasta <input-fasta.fa> --outpath </path/for/output-folder> --outname <outname> --iter <5>
```

List of options available:
- **--gtf**: Transcriptome annotation file in GTF format
- **--fasta**: Transcripts nucleotide fasta file
- **--outpath**: Path of the output folder
- **--outname**: Prefix for the output files
- **--iter**: Maximum number of 'start-fixing & translation' cycles to identify alternative start sites. Default: 5
- **--chimeric**: Table indicating chimeric genes in the annotation (Optional) 

Example:
```
python transuite.py TransFix --gtf ./test_dataset/test_output/test_run_longorf/test_run_longorf.gtf --fasta ./test_dataset/subset_AtRTD2_transcripts.fa --outpath ./test_dataset/test_output --outname test_run --iter 5
```

## Output files
TransFix automatically generates a subfolder to store the output files:
/**&lt;outpath&gt;**/**&lt;outname&gt;**\_transfix/

TransFix generates the following output files:
1. *GTF* file with the fixed CDS coordinates at the gene-level
2. *FASTA* files of the transcripts CDS regions (nucleotide, and peptide)
3. Multiple log *CSV* files: a) log files reporting transcripts that could not be annotated, for example for lack of an annotated CDS; and b) logfile tracking the *fixing* cycle at which the AUG was annotated


-------------------
**TransFeat**
==============
-------------------

TransFeat extracts and processes the transcripts CDS information contained in transcriptome annotations to infer multiple characteristics of the genes, transcripts and their coding potential, and it to reports theis information in an easily accessible format.

### **Command and options** ###
Command to run TransFeat:
```
python transuite.py TransFeat [options]
```
```
python transuite.py TransFeat --gtf <input-gtf.gtf> --fasta <input-fasta.fa> --outpath </path/for/output-folder> --outname <outname> --pep <30> --ptc <70>
```
List of options available:
- **--gtf**: Transcriptome annotation file in GTF format
- **--fasta**: Transcripts nucleotide fasta file
- **--outpath**: Path of the output folder
- **--outname**: Prefix for the output files
- **--pep**: Minimum number of amino-acids a translation must have to be consider a peptide. Default: 100 AA
- **--ptc**: Minimum CDS length percentage below which a transcript is considered prematurely terminated (PTC). Default: 70%

Example:
```
python transuite.py TransFeat --gtf ./test_dataset/test_output/test_run_transfix/test_run_transfix.gtf --fasta ./test_dataset/subset_AtRTD2_transcripts.fa --outpath ./test_dataset/test_output --outname test_run --pep 30 --ptc 70
```

**Note:** When running the analysis on the above *test_dataset* you will get a *WARNING* message regarding a number of features not present in the TransFeat table. This is expected given the small number of transcripts in the test dataset.

### **Output files** ###
TransFeat automatically generates a subfolder to store the output files:
/**&lt;outpath&gt;**/**&lt;outname&gt;**\_longorf/

TransFeat generates the following output files:
1. A main *CSV* table reporting the transcripts coding features
2. *FASTA* files of: (**1**) transcripts with CDS (both protein-coding and unproductive), (**2**) transcripts classified by coding-potentiality (protein-coding transcripts, non-coding genes), (**3**) transcripts alternative ORFs (uORF, ldORF) 
3. Multiple *CSV* tables reporting number of transcripts and transcripts-features subdivided by gene categories
4. Multiple *CSV* tables reporting co-ordinates and sequences of transcripts subdivided by feature categories (ldORF, NMD)
5. A *JSON* file containing the transcripts ORF relative coordinates


-------------------
**Auto**
==============
-------------------

This module performs FindLORF, TransFix, and TransFeat analysis in tandem with a single command. 

### **Command and options** ###
Command to run Auto:
```
python transuite.py Auto [options]
```
```
python transuite.py Auto --gtf <input-gtf.gtf> --fasta <input-fasta.fa> --outpath </path/for/output-folder> --outname <outname> --cds <30> --iter <5> --pep <100> --ptc <70>
```
List of options available:
- **--gtf**: Transcriptome annotation file in GTF format
- **--fasta**: Transcripts nucleotide fasta file
- **--outpath**: Path of the output folder
- **--outname**: Prefix for the output files
- **--cds**: Minimum number of amino-acids an ORF must have to be considered as a potential CDS. Default: 30 AA
- **--iter**: Maximum number of 'start-fixing & translation' cycles to identify alternative start sites. Default: 5
- **--pep**: Minimum number of amino-acids a translation must have to be consider a peptide. Default: 100 AA
- **--ptc**: Minimum CDS length percentage below which a transcript is considered prematurely terminated (PTC). Default: 70%
- **--chimeric**: Table indicating chimeric genes in the annotation (Optional)

Example:
```
python transuite.py Auto --gtf ./test_dataset/subset_AtRTD2.gtf --fasta ./test_dataset/subset_AtRTD2_transcripts.fa --outpath ./test_dataset/test_output --outname test_run --cds 30 --iter 5 --pep 100 --ptc 70
```

### **Output files** ###
The Auto module automatically generates all of the modules subfolders and their respective output files:
/**&lt;outpath&gt;**/**&lt;outname&gt;**\_longorf/
/**&lt;outpath&gt;**/**&lt;outname&gt;**\_transfix/
/**&lt;outpath&gt;**/**&lt;outname&gt;**\_transfeat/

The main output files of TranSuite pipeline are:
1. The *GTF* file generate by TransFix
2. The main *CSV* feature table generate by TransFeat
3. The *FASTA* files generate by TransFix and/or as classified by TransFeat (Coding transcripts, Non-coding genes)
4. Any of the multiple log files generated during the analysis


----------------------------
# Contact
----------------------------

For questions about the modified features, please contact: mojtaba.bagherian@uwa.edu.au

For questions about the original TranSuite, please contact: e.entizne@dundee.ac.uk

----------------------------
Citation
----------------------------

When using this TranSuite version, please cite both:

The original TranSuite paper
This modified version

----------------------------
# License
----------------------------
This TranSuite version is released under the same [MIT license](https://opensource.org/licenses/MIT) as the original TranSuite.
