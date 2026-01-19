
Script to identify event types for DEAD-dependent differentially expressed genes

Created by Maximilian Sack, Bioinformatics Group, Department of Computer Science and Interdisciplinary Centre for Bioinformatics, Leipzig University

Part of the Study:

A structured RNA balances DEAD-box RNA helicase function in plant alternative splicing control (2026)

Rica Burgardt, Julia Bauer, Maren Reinhardt, Natalie Rupp, Christoph Engel, Lukas Hellmann, Maximilian Sack, Zasha Weinberg, and Andreas Wachter

Cite that paper if you use this script.


Example command line:
```
python3 findSpliceEvents.py 
	AtRTDv2_QUASI_19April2016.gtf 
	DTU_transcripts_cpm4_sample1_3h_6h_all_transcripts_3h.tab 
	--output_all SpliceEvents_all.tab 
	--output_regulated SpliceEvents_regulated.tab 
	--output_reciprocal SpliceEvents_reciprocal.tab
```
Ensure that the input files have the correct format, like this:

Differential gene expression table:
```
AT1G01260_P1,X03h.EST-X03h.mock,0.00395097538878148,0.171207076302256,0.369400232338061,up-regulated
AT1G01260_P3,X03h.EST-X03h.mock,0.00395097538878148,-0.171207076302256,-1.19466698377057,down-regulated
```
Transcript table from AtRTD:
```
Chr1    Araport11       exon    108946  111699  .       +       .       transcript_id "AT1G01260_P1"; gene_id "AT1G01260"; Note "basic";
Chr1    Araport11       exon    199527  199763  .       +       .       transcript_id "AT1G01550_P1"; gene_id "AT1G01550"; Note "BPS1-like";
Chr1    Araport11       exon    199890  199959  .       +       .       transcript_id "AT1G01550_P1"; gene_id "AT1G01550"; Note "BPS1-like";
Chr1    Araport11       exon    200511  201775  .       +       .       transcript_id "AT1G01550_P1"; gene_id "AT1G01550"; Note "BPS1-like";
```
