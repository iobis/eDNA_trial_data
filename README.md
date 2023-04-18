# eDNA_trial_data
Analysis of eDNA trial samples


# eDNA expeditions trial samples

This document shows the overall results of the eDNA expeditions trial sample analysis. 

The analysis consisted of a mixture of samples of 

- eDNA filters
- control filters
- pcr controls

Which were first amplified in 3 different types of reactions:

- Mifish multiplex (all mifish and mimammal primers combined in one reaction)
- All multiplex (all primers combined in one reaction)
- All singleplex (all primers amplified separately, then combined for sequencing)

The samples are from the eDNA expeditions project, and are from three locations: 

- Wadden Sea 
- Brazil
- Australia

Random positive samples (PCR bands detected) from each of the countries were combined together for the analysis. 

The primers that were used for the analysis were:

- Mifish-UE, Mimammal-UEB
- COI (Leray fragment)
- 16S rRNA for vertebrates
- Teleo 12S rRNA 

The purpose of this analysis was to test the sample processing pipeline and to understand which reaction type works best with the samples. 

All sequences were analysed with the [PacMAN pipeline](https://github.com/iobis/PacMAN-pipeline), with each primer set analysed separately. Cutadapt was used to search for the primer sequences and discard any reads were the primer sequences were not found. For the MifishUE and the Mimammal, the different primer sets in each category were combined as a consensus sequence. Still there is a large overlap within the seqeunces found by Mimammal and Mifish primers, as the primers are very similar. The overlap in ASVs is removed in the analysis stage. 

The sequences were obtained with a paired run of a 150 bp Novaseq kit. Only paired sequences from the mifish and mimammal analyses were merged in the analysis. The COI and 16S rRNA markers were too long to be merged (313 bp and 250 bp respectively) as the length of good quality sequences for the 16S rRNA was also about 125 bp. For Teleo, it is a short read, but there were large differences in the amount of primers recognized in forward reads or reverse reads by cutadapt (possibly something to still check why). Cutadapt removes reads if the other pair is not found in paired-end mode, therefore the reads for teleo were run separately. In this case, each read (both forward and reverse), cover the full range of the marker. 

A reference database built with the [reference database creator](https://github.com/gjeunen/reference_database_creator) either from ncbi or from the mitofish database (for the mifish primers), and cut to the specific primer region that was used. Remaining unknown sequences were searched against the full ncbi-nt database. 

The sequence annotation in this first trial was not strict, with an 85% identity cutoff for both the targeted database file and blast. Therefore the annotation results need to be considered carefully. 
