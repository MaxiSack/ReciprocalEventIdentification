# Idenify alternative splicing events from differentially expressed genes

This script was used to identify alternative splicing event types for DEAD-dependent differentially expressed genes.

Created by Maximilian Sack, Bioinformatics Group, Department of Computer Science and Interdisciplinary Centre for Bioinformatics, Leipzig University

Part of the Study:

A structured RNA balances DEAD-box RNA helicase function in plant alternative splicing control (2026)
Rica Burgardt, Julia Bauer, Maren Reinhardt, Natalie Rupp, Christoph Engel, Lukas Hellmann, Maximilian Sack, Zasha Weinberg, and Andreas Wachter

Cite that paper if you use this script. It also contains a detailed description of the logic applied by this script.

## Usage

Example command line:
```
python3 findSpliceEvents.py 
	AtRTDv2_QUASI_19April2016.gtf 
	DTU_transcripts.csv 
	--output_all SpliceEvents_all.tab 
	--output_regulated SpliceEvents_regulated.tab 
	--output_reciprocal SpliceEvents_reciprocal.tab
```
Ensure that the input files have the correct format, like this:

Transcript table from AtRTD:
```
Chr1    Araport11       exon    108946  111699  .       +       .       transcript_id "AT1G01260_P1"; gene_id "AT1G01260"; Note "basic";
Chr1    Araport11       exon    199527  199763  .       +       .       transcript_id "AT1G01550_P1"; gene_id "AT1G01550"; Note "BPS1-like";
Chr1    Araport11       exon    199890  199959  .       +       .       transcript_id "AT1G01550_P1"; gene_id "AT1G01550"; Note "BPS1-like";
Chr1    Araport11       exon    200511  201775  .       +       .       transcript_id "AT1G01550_P1"; gene_id "AT1G01550"; Note "BPS1-like";
```
See also `AtRTDv2_QUASI_19April2016.gtf` on [AtRTDv2](https://ics.hutton.ac.uk/atRTD/).

Differential gene expression table from RNAseq:
```
AT1G01260_P1,X03h.EST-X03h.mock,0.00395097538878148,0.171207076302256,0.369400232338061,up-regulated
AT1G01260_P3,X03h.EST-X03h.mock,0.00395097538878148,-0.171207076302256,-1.19466698377057,down-regulated
```

## Demo
This repository contains a small demo dataset. You can use it to test the scripts functionality by running the following command:
```
python3 findSpliceEvents.py demo/demo_AtRTDv2.gtf demo/demo_DTU_transcripts.csv --output_all demo/demo_SpliceEvents_all.tab --output_regulated demo/demo_SpliceEvents_regulated.tab --output_reciprocal demo/demo_SpliceEvents_reciprocal.tab
```
It takes 0 seconds to run.

## System Requirements
The script requries Python 3 to run and as such should be operating system independent.

We tested it with the following software configurations:

| Operation System     | Python version |
| -------------------- | -------------- |
| Ubuntu 22.04         | Python 3.10.12 |
| Ubuntu 24.04         | Python 3.12.3  |
| Windows 11 Education | Python 3.13.5  |

The runtime for the full dataset from the study was between 1.4 and 6 seconds, depending on hardware.

## Installation
[Python 3](https://www.python.org/downloads/) needs to be installed.
The script itself is then executed through python.

