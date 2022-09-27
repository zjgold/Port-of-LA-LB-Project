# A Manager’s Guide to Using eDNA Metabarcoding in Marine Ecosystems

Zachary Gold<sup>1</sup>, Adam R. Wall<sup>2</sup>, Teia M. Schweizer<sup>3</sup>, N. Dean Pentcheff<sup>2</sup>, Emily E. Curd<sup>4</sup>, Paul H. Barber<sup>1</sup>, Rachel S. Meyer<sup>1,5</sup>, Robert Wayne<sup>1</sup>, Kevin Stolzenbach<sup>6</sup>, Kat Prickett<sup>7</sup>, Justin Luedy<sup>8</sup>, Regina Wetzer<sup>2</sup>


<sup>1</sup>Department of Ecology and Evolutionary Biology, University of California Los Angeles, Los Angeles, CA, USA <br/>
<sup>2</sup>Diversity Initiative for the Southern California Ocean (DISCO), Natural History Museum of Los Angeles County, Los Angeles, CA, USA<br/>
<sup>3</sup>Department of Fish and Wildlife Conservation Biology, Colorado State University, Fort Collins, CO, USA<br/>
<sup>4</sup>Department of Natural Sciences, Landmark College, Putney, VT, USA<br/>
<sup>5</sup>Department of Ecology and Evolutionary Biology, University of California Los Angeles, Los Angeles, CA, USA<br/>
<sup>6</sup>Wood Environment and Infrastructure, Inc., San Diego, CA, USA<br/>
<sup>7</sup>Port of Los Angeles, San Pedro, CA, USA<br/>
<sup>8</sup>Port of Long Beach, Long Beach, CA, USA<br/>


## Description
This page is dedicated to hosting data and code generated for the manuscript. <br/>
Pre-print is available here:
2018 Report to Port of Los Angeles and Port of Long Beach is available here: *eDNA Report Final.pdf*

Included on this page is
1. Scripts used to Conduct Analyses

    * /analysis

      * *ports_data_analysis_20220619.Rmd* This script does the main analyses in the paper.

  * /data/fish

      * *20210420_ports_decontamination.R* This script processes the "raw" *Anacapa Toolkit* output and runs decontamination to remove poorly sequenced samples and contaminant ASVs. Used on the MiFish *12S* Universal Teleost and Elasmobranch data.

  * /data/invert

      * *220220130_ports_CO1_16S_decontamination.R* This script processes the "raw" *Anacapa Toolkit* output and runs decontamination to remove poorly sequenced samples and contaminant ASVs.Used on the *CO1* and *16S* metazoan data.

  * /data/Anacapa_output_old

      * *decontam.R* This script processes the "raw" *Anacapa Toolkit* output and runs decontamination to remove poorly sequenced samples and contaminant ASVs. Used on the MiFish *12S* Universal Teleost and Elasmobranch data that was processed with an older reference database that lacked hundreds of CA species barcodes.

2. Data

  * *abundance_trawl* Trawl data in counts of individual fish.
  * *biomass_trawl* Trawl data in biomass of fish.
  * *la_ports_meta_data.txt* Metadata of eDNA and trawl data associated with stations/bottles sampled.

  * /fish

    1. *miu_c19_fishcard_taxonomy_tables* *Anacapa Toolkit* Output of MiFish *12S* Universal Teleost data using the *CRUX* Global Reference database from [Gold et al. 2021](https://onlinelibrary.wiley.com/doi/epdf/10.1111/1755-0998.13450)       
    2. *miu_fishcard_12S_all_taxonomy_tables* *Anacapa Toolkit* Output of MiFish *12S* data using the *CRUX* Local Reference database from Gold et al. 2021. [See GitHub](https://github.com/zjgold/FishCARD) & [See Gold et al. 2021 Dryad](https://doi.org/10.5068/D1H963)
    3. *elas_c19_fishcard_taxonomy_tables* *Anacapa Toolkit* Output of MiFish *12S* Elasmobranch data using the *CRUX* Global Reference database.   
    4. *elas_fishcard_12S_all_taxonomy_tables* *Anacapa Toolkit* Output of MiFish 12S data using the *CRUX* Local Reference database.
    5. *decontamination/* Directory with decontamination outputs.
    6. *CA_fish_list_20210216.csv* List of California fishes obtained from Gold et al. 2021

  * /Anacapa_output_old

    1. *12S_U* *Anacapa Toolkit* Output of MiFish *12S* Universal Teleost data using an old *CRUX* Reference database with no additional CA species barcodes added.[See Anacapa Dryad](https://doi.org/10.1111/2041-210x.13214)       
    2. *12S_E* *Anacapa Toolkit* Output of MiFish *12S* Elasmobranch data using an old *CRUX* Reference database with no additional CA species barcodes added.
    3. *decontamination/* Directory with decontamination outputs.

  * /invert

    1. *CO1_taxonomy_tables* *Anacapa Toolkit* Output of  *CO1* metazoan data using a *CRUX* Reference database with additional marine barcodes added from the Mo'orea BIOCODE project.       
    2. *Metazoa_16s_taxonomy_tables* *Anacapa Toolkit* Output of MiFish *16S* Elasmobranch data using an old *CRUX* Reference database with no additional species barcodes added.
    3. *decontamination/* Directory with decontamination outputs.


**Note** *A few files needed for analyses are not immediately included on GitHub due to size, please see [See Dryad]INSERT LINK*

**Raw sequence inputs, reference databases, and Anacapa scripts  are made available on Dryad**

# Dryad readme with more details:


## GENERAL INFORMATION

1. Date of eDNA sample collection: 20–21 August 2018

2. Geographic location of data collection: Ports of Los Angeles and Long Beach, California, USA

3. Funding sources that supported the collection of data: The Port of Los Angeles and Port of Long Beach, UC Office of the President Catalyst Program (CA-16-376437), the Howard Hughes Medical Institute Grant (GT10483), and NSF-GRFP 2015204395 (to Zachary Gold).

4. Recommended citation for this dataset: Gold, Zachary et al. (2021), Data from:A Manager’s Guide to Using eDNA Metabarcoding in Marine Ecosystems, Dryad, dataset

## DATA & FILE OVERVIEW

1. Description of the data set

These data were generated from midwater trawls and eDNA metabarcoding of seawater samples collected within the Ports of Los Angeles and Long Beach, California, USA. For midwater trawls, we trawled the station for 5 minutes with a 7.6 m semi-balloon otter trawl with 2.5-cm side mesh and 1.3-cm mesh cod-end deployed from the R/V Early Bird II (Figure 2D), following standard protocols of the Ports Biosurvey program, as detailed in the 2018 Biosurvey report (Wood, 2021). All fish were identified to species level and standard length was recorded. See Table S2 for details of trawl locations and depths. Trawls were conducted within 15–90 minutes of the eDNA seawater sampling. Resulting data included counts and biomass of fish observed in trawls along with associated metadata from each trawl.
For eDNA metabarcoding samples, we sampled three replicate 1 L seawater samples at four locations 2 m above the sea floor at each of the 7 sampling stations along an approximately 150 m transect aligned with the trawl track, using the R/V Yellowfin as a sampling platform. Seawater was gravity filtered onto 0.22µm Sterivex filter cartridges and stored on dry ice before being transferred to a -20˚C freezer. We exacted eDNA from the frozen Sterivex filters within one week of collection using the DNeasy Blood and Tissue Qiagen Kit (Spens et al., 2017). We amplified eDNA with four primer sets (described below) modified with Illumina Nextera XT adapters (Illumina, San Diego, CA) in triplicate PCR replicates. We prepared sequencing libraries following the methods of (Curd et al., 2019) which involved a second PCR indexing step to tag libraries so they can later be distinguished from each other (See S1 Appendix 1). Resulting indexed libraries were bead cleaned to remove fragments <200bp, were quantified with the Qubit BR assay (Thermofisher, Waltham, MA, USA), then sequenced on an Illumina MiSeq (Illumina Inc., La Jolla, CA, USA) with V3 kit to obtain paired-end 2x 300 bp reads at the Technology Center for Genomics & Bioinformatics (UCLA, CA, USA) (See S1 appendix 1 for further details). We used the Anacapa Toolkit (Curd et al., 2019) for amplicon sequence variant (ASV) parsing, quality control, and taxonomic assignment (Curd et al., 2019; Gold et al., 2021a) (See S1 Appendix 1 for full description). Importantly, we note that we employed stricter sequence alignment parameters within the Anacapa classifier for the CO1 and 16S locus (95% identity and query coverage) than the 12S loci (80% identity and query coverage) given the lack of complete reference databases for the CO1 and 16S loci (Curd et al. 2019). Separate reference databases were generated for the 12S, COI, and 16S from all publicly available sequencing data in GenBank from October 2019 using CRUX with default parameters (Curd et al. 2019). A more comprehensive 12S reference database was created by supplementing the above database with barcodes from 252 native fish species as detailed in Gold et al. 2021a. Resulting data included raw DNA sequences (fastq files), intermediate processed data resulting from the processing of raw sequences through the Anacapa Toolkit, metadata associated with each bottle collection, and reference databases used to assign taxonomy to derived amplicon sequence variants.

2. File List:


    Directory 1: CO1/

        Directory 1.a: CO1_fasta_and_taxonomy/

          File 1.a.1 Name: CO1_.fasta
          File 1.a.1 Description: Leray CO1 reference database fasta file with accession numbers and associated DNA sequence

          File 1.a.2 Name: CO1_taxonomy.txt
          File 1.a.2 Description: Leray CO1 reference database taxonomy table with accession numbers and full taxonomic paths (Domain; Phylum;Class;Order;Family;Genus;Species)


    Directory 2: Metazoa_16s/

        Directory 2.a: Metazoa_16s_fasta_and_taxonomy/

          File 2.a.1 Name: Metazoa_16s_.fasta
          File 2.a.1 Description: Kelly 16S reference database fasta file with accession numbers and associated DNA sequence

          File 2.a.2 Name: Metazoa_16s_taxonomy.txt
          File 2.a.2 Description: Kelly 16S reference database taxonomy table with accession numbers and full taxonomic paths (Domain; Phylum;Class;Order;Family;Genus;Species)


    Directory 3: data/

      File 3.a.1 Name: abundance_trawl.txt
      File 3.a.1 Description: Data table of fish counts found in each trawl (Species x Site) including scientific names, common names, and counts of fish found in each trawl.

      File 3.a.2 Name: biomass_trawl.txt
      File 3.a.2 Description: Data table of fish biomass (kg) found in each trawl (Species x Site) including scientific names, common names, and counts of fish found in each trawl.

      File 3.a.3 Name: eDNA Report Final.pdf
      File 3.a.3 Description: Grey literature 2018 Report to Port of Los Angeles and Port of Long Beach. Preliminary analyses and findings.

      File 3.a.4 Name: la_ports_meta_data.txt
      File 3.a.4 Description: Metadata table for eDNA samples. Columns include Seq_number (name of sequence fastq file),	New_name (Updated internal name in Site_Biorep_Tech_rep),	Site (Port of Los Angeles and Long Beach trawl site name),	Bio_rep (location within station),	Site_station (station_location combined),	Station (station sampled),	Tech_rep (bottle replicate within a location),	Replicate_per_station (unique bottle number within a station),	k (number of bottles within a station collected),	k_station (number of bottles within a locaiton collected),	Time_passed_h (length of time samples were gravity filtered),	Start_volume_mL (starting volume of sea water collected from enteral feeding pouch),	End_volume_mL (final volume of sea water collected from enteral feeding pouch),	Volume_filtered_mL (total volume filtered),	Latitude,	Longitude,	Salinity (PSU),	DO_mg_l (dissolved oxygen mg/L),	pH,	CTD_depth (maximum depth of site),	Sample_control (demarcation of samples and controls), Control_Type (labeling of control type, here all Blanks or negative controls)

      File 3.a.5 Name: Supplemental_tables_20220619.xlsx
      File 3.a.5 Description: Supplemental Tables from the associated manuscript. Title of each tab provides description of table. See manuscript and supplemental materials for full details of each table.

      Directory 3.a: Anacapa_output_old/

        Directory 3.a.1: 12S_E/

          Directory 3.a.1.a: 12S_taxonomy_tables/
          Description Directory 3.a.1.a:  Anacapa Toolkit Output files for eDNA metabarcoding data generated from MiFish elasmobranch primers and processed with the CRUX generated 12S reference database generated in Curd et al. 2019 (See https://doi.org/10.5061/dryad.mf0126f). The detailed and brief taxonomy files are described in detail in Curd et al. 2019 and associated Anacapa Toolkit software pages. The most relevant file here is ~/Summary_by_percent_confidence/60/12S_ASV_sum_by_taxonomy_60.txt which is the final ASV table with taxonomy assigned with a confidence cutoff score of 60 used for our analyses within the manuscript.

            File 3.a.1.a.1 Name: 12S_ASV_taxonomy_detailed.txt
            File 3.a.1.a.1 Description: Detailed output file from Anacapa Toolkit. Columns are as follows: 12S_seq_number (ASV name),	sequence (DNA sequence of ASV),	sequencesF (forward sequence of unmerged ASVs, otherwise empty),	sequencesR (reverse sequence of unmerged ASVs, otherwise empty),	forward_12S_seq_number (column of zeros, not important),	columns of Sample data [file naming convention is barcode_samplename e.g. X12S_LA11.1.R1.S25.L001] (counts of sequences for each ASV),	merged_12S_seq_number column of zeros, not important),	reverse_12S_seq_numbercolumn of zeros, not important),	unmerged_12S_seq_numbercolumn of zeros, not important),	single_or_multiple_hit (single or multiple sequence alignment hits from bowtie2 sequence alignment),	end_to_end_or_local(end to end or local sequence alignment hits from bowtie2 sequence alignment),	max_percent_id (maximum percent sequence alignment ID from bowtie2 sequence alignment),	input_sequence_length (ASV sequence length), 	taxonomy (assigned taxonomy from Anacapa Toolkit's Bayesian Lowest Common Ancestor classifier, full taxonomic path : Domain; Phylum;Class;Order;Family;Genus;Species), 	taxonomy_confidence (Taxonomic path and Bayesian cutoff score (bcc) for each taxonomic level Domain;Domain bcc; Phylum;Phylum bcc;Class;Class bcc;Order;Order bcc;Family;Family bcc;Genus;Genus bcc;Species;Species bcc;),	accessions (accessions of top aligned sequences).

            File 3.a.1.a.2 Name: 12S_ASV_taxonomy_brief.txt
            File 3.a.1.a.2 Description: Brief output file from Anacapa Toolkit. Columns are as follows: 12S_seq_number (ASV name),	columns of Sample data [file naming convention is barcode_samplename e.g. X12S_LA11.1.R1.S25.L001] (counts of sequences for each ASV),	taxonomy (assigned taxonomy from Anacapa Toolkit's Bayesian Lowest Common Ancestor classifier, full taxonomic path : Domain; Phylum;Class;Order;Family;Genus;Species), 	taxonomy_confidence (Taxonomic path and Bayesian cutoff score (bcc) for each taxonomic level Domain;Domain bcc; Phylum;Phylum bcc;Class;Class bcc;Order;Order bcc;Family;Family bcc;Genus;Genus bcc;Species;Species bcc;),	accessions (accessions of top aligned sequences).

            Directory 3.a.1.a.1: Summary_by_percent_confidence/
            Description Directory 3.a.1.a.1: Taxonomy tables filtered by Bayesian cutoff score. For example, ~/60/ has only retained taxonomic assignments with bayesian cutoff scores > 60. Directories include 40,50,60,70,80,90,95,100 each which is filtered at directory labeled cutoff score. Only the 60 directory is described in detail below but this applies to all subfolders within this directory.
                Directory 3.a.1.a.1.a: 60/
                    File 3.a.1.a.1.a.1 Name: 12S_ASV_raw_taxonomy_60.txt
                    File 3.a.1.a.1.a.1 Description: 12S_seq_number (ASV name),	columns of Sample data [file naming convention is barcode_samplename e.g. X12S_LA11.1.R1.S25.L001] (counts of sequences for each ASV),	taxonomy (assigned taxonomy from Anacapa Toolkit's Bayesian Lowest Common Ancestor classifier, full taxonomic path : Domain; Phylum;Class;Order;Family;Genus;Species).

                    File 3.a.1.a.1.a.2 Name: 12S_ASV_sum_by_taxonomy_60.txt
                    File 3.a.1.a.1.a.2 Description: This table  has grouped and summarized reads for all ASVs with the same taxonomic path. For example, if there were 4 ASVs assigned to the same species, the reads of each of those species at each site were summed together. Column names: sum.taxonomy (assigned taxonomy from Anacapa Toolkit's Bayesian Lowest Common Ancestor classifier, full taxonomic path : Domain; Phylum;Class;Order;Family;Genus;Species).,	columns of Sample data [file naming convention is barcode_samplename e.g. X12S_LA11.1.R1.S25.L001] (counts of sequences for each ASV).


          Directory 3.a.2: 12S_U/

            Directory 3.a.2.a: 12S_taxonomy_tables/
            Description Directory 3.a.2.a:  Anacapa Toolkit Output files for eDNA metabarcoding data generated from MiFish Universal Teleost primers and processed with the CRUX generated 12S reference database generated in Curd et al. 2019 (See https://doi.org/10.5061/dryad.mf0126f). The detailed and brief taxonomy files are described in detail in Curd et al. 2019 and associated Anacapa Toolkit software pages. The most relevant file here is ~/Summary_by_percent_confidence/60/12S_ASV_sum_by_taxonomy_60.txt which is the final ASV table with taxonomy assigned with a confidence cutoff score of 60 used for our analyses within the manuscript.

                File 3.a.2.a.1 Name: 12S_ASV_taxonomy_detailed.txt
                File 3.a.2.a.1 Description: Detailed output file from Anacapa Toolkit. Columns are as follows: 12S_seq_number (ASV name),	sequence (DNA sequence of ASV),	sequencesF (forward sequence of unmerged ASVs, otherwise empty),	sequencesR (reverse sequence of unmerged ASVs, otherwise empty),	forward_12S_seq_number (column of zeros, not important),	columns of Sample data [file naming convention is barcode_samplename e.g. X12S_LA11.1.R1.S25.L001] (counts of sequences for each ASV),	merged_12S_seq_number column of zeros, not important),	reverse_12S_seq_numbercolumn of zeros, not important),	unmerged_12S_seq_numbercolumn of zeros, not important),	single_or_multiple_hit (single or multiple sequence alignment hits from bowtie2 sequence alignment),	end_to_end_or_local(end to end or local sequence alignment hits from bowtie2 sequence alignment),	max_percent_id (maximum percent sequence alignment ID from bowtie2 sequence alignment),	input_sequence_length (ASV sequence length), 	taxonomy (assigned taxonomy from Anacapa Toolkit's Bayesian Lowest Common Ancestor classifier, full taxonomic path : Domain; Phylum;Class;Order;Family;Genus;Species), 	taxonomy_confidence (Taxonomic path and Bayesian cutoff score (bcc) for each taxonomic level Domain;Domain bcc; Phylum;Phylum bcc;Class;Class bcc;Order;Order bcc;Family;Family bcc;Genus;Genus bcc;Species;Species bcc;),	accessions (accessions of top aligned sequences).

                File 3.a.2.a.2 Name: 12S_ASV_taxonomy_brief.txt
                File 3.a.2.a.2 Description: Brief output file from Anacapa Toolkit. Columns are as follows: 12S_seq_number (ASV name),	columns of Sample data [file naming convention is barcode_samplename e.g. X12S_LA11.1.R1.S25.L001] (counts of sequences for each ASV),	taxonomy (assigned taxonomy from Anacapa Toolkit's Bayesian Lowest Common Ancestor classifier, full taxonomic path : Domain; Phylum;Class;Order;Family;Genus;Species), 	taxonomy_confidence (Taxonomic path and Bayesian cutoff score (bcc) for each taxonomic level Domain;Domain bcc; Phylum;Phylum bcc;Class;Class bcc;Order;Order bcc;Family;Family bcc;Genus;Genus bcc;Species;Species bcc;),	accessions (accessions of top aligned sequences).

                Directory 3.a.2.a.1: Summary_by_percent_confidence/
                Description Directory 3.a.2.a.1: Taxonomy tables filtered by Bayesian cutoff score. For example, ~/60/ has only retained taxonomic assignments with bayesian cutoff scores > 60. Directories include 40,50,60,70,80,90,95,100 each which is filtered at directory labeled cutoff score. Only the 60 directory is described in detail below but this applies to all subfolders within this directory.
                    Directory 3.a.2.a.1.a: 60/

                        File 3.a.2.a.1.a.1 Name: 12S_ASV_raw_taxonomy_60.txt
                        File 3.a.2.a.1.a.1 Description: 12S_seq_number (ASV name),	columns of Sample data [file naming convention is barcode_samplename e.g. X12S_LA11.1.R1.S25.L001] (counts of sequences for each ASV),	taxonomy (assigned taxonomy from Anacapa Toolkit's Bayesian Lowest Common Ancestor classifier, full taxonomic path : Domain; Phylum;Class;Order;Family;Genus;Species).

                        File 3.a.2.a.1.a.2 Name: 12S_ASV_sum_by_taxonomy_60.txt
                        File 3.a.2.a.1.a.2 Description: This table  has grouped and summarized reads for all ASVs with the same taxonomic path. For example, if there were 4 ASVs assigned to the same species, the reads of each of those species at each site were summed together. Column names: sum.taxonomy (assigned taxonomy from Anacapa Toolkit's Bayesian Lowest Common Ancestor classifier, full taxonomic path : Domain; Phylum;Class;Order;Family;Genus;Species).,	columns of Sample data [file naming convention is barcode_samplename e.g. X12S_LA11.1.R1.S25.L001] (counts of sequences for each ASV).


      Directory 3.b: fish/

        File 3.b.1 Name: CA_fish_list_20210216.csv
        File 3.b.1 Description: List of California Current Large Marine Ecosystem fishes from Gold et al. 2021 (See https://doi.org/10.5068/D1H963).

          Directory 3.b.1: elas_c19_fishcard_taxonomy_tables/
          Description Directory 3.b.1:  Anacapa Toolkit Output files for eDNA metabarcoding data generated from MiFish elasmobranch primers and processed with the global CRUX generated 12S reference database generated in from Gold et al. 2021(See https://doi.org/10.5068/D1H963). The most relevant file here is ~/Summary_by_percent_confidence/60/c19_fishcard_ASV_raw_taxonomy_60.txt which is the final ASV table with taxonomy assigned with a confidence cutoff score of 60 used for our analyses within the manuscript.

              File 3.b.1.a Name: c19_fishcard_ASV_taxonomy_detailed.txt
              File 3.b.1.a Description: Detailed output file from Anacapa Toolkit. Columns are as follows: c19_fishcard_seq_number (ASV name),	sequence (DNA sequence of ASV),	sequencesF (forward sequence of unmerged ASVs, otherwise empty),	sequencesR (reverse sequence of unmerged ASVs, otherwise empty),	forward_c19_fishcard_seq_number (column of zeros, not important),	columns of Sample data [file naming convention is barcode_samplename e.g. c19_fishcard_LA11.1.R1.S25.L001] (counts of sequences for each ASV),	merged_c19_fishcard_seq_number column of zeros, not important),	reverse_c19_fishcard_seq_numbercolumn of zeros, not important),	unmerged_c19_fishcard_seq_numbercolumn of zeros, not important),	single_or_multiple_hit (single or multiple sequence alignment hits from bowtie2 sequence alignment),	end_to_end_or_local(end to end or local sequence alignment hits from bowtie2 sequence alignment),	max_percent_id (maximum percent sequence alignment ID from bowtie2 sequence alignment),	input_sequence_length (ASV sequence length), 	taxonomy (assigned taxonomy from Anacapa Toolkit's Bayesian Lowest Common Ancestor classifier, full taxonomic path : Domain; Phylum;Class;Order;Family;Genus;Species), 	taxonomy_confidence (Taxonomic path and Bayesian cutoff score (bcc) for each taxonomic level Domain;Domain bcc; Phylum;Phylum bcc;Class;Class bcc;Order;Order bcc;Family;Family bcc;Genus;Genus bcc;Species;Species bcc;),	accessions (accessions of top aligned sequences).

              File 3.b.1.b Name: c19_fishcard_ASV_taxonomy_brief.txt
              File 3.b.1.b Description: Brief output file from Anacapa Toolkit. Columns are as follows: c19_fishcard_seq_number (ASV name),	columns of Sample data [file naming convention is barcode_samplename e.g. c19_fishcard_LA11.1.R1.S25.L001] (counts of sequences for each ASV),	taxonomy (assigned taxonomy from Anacapa Toolkit's Bayesian Lowest Common Ancestor classifier, full taxonomic path : Domain; Phylum;Class;Order;Family;Genus;Species), 	taxonomy_confidence (Taxonomic path and Bayesian cutoff score (bcc) for each taxonomic level Domain;Domain bcc; Phylum;Phylum bcc;Class;Class bcc;Order;Order bcc;Family;Family bcc;Genus;Genus bcc;Species;Species bcc;),	accessions (accessions of top aligned sequences).

              Directory 3.b.1.a: Summary_by_percent_confidence/
              Description Directory 3.b.1.a: Taxonomy tables filtered by Bayesian cutoff score. For example, ~/60/ has only retained taxonomic assignments with bayesian cutoff scores > 60. Directories include 40,50,60,70,80,90,95,100 each which is filtered at directory labeled cutoff score. Only the 60 directory is described in detail below but this applies to all subfolders within this directory.

                  Directory 3.b.1.a.1: 60/

                        File 3.b.1.a.1.a Name: c19_fishcard_ASV_raw_taxonomy_60.txt
                        File 3.b.1.a.1.a.1 Description: c19_fishcard_seq_number (ASV name),	columns of Sample data [file naming convention is barcode_samplename e.g. c19_fishcard_LA11.1.R1.S25.L001] (counts of sequences for each ASV),	taxonomy (assigned taxonomy from Anacapa Toolkit's Bayesian Lowest Common Ancestor classifier, full taxonomic path : Domain; Phylum;Class;Order;Family;Genus;Species).

                        File 3.b.1.a.1.b Name: c19_fishcard_ASV_sum_by_taxonomy_60.txt
                        File 3.b.1.a.1.b Description: This table  has grouped and summarized reads for all ASVs with the same taxonomic path. For example, if there were 4 ASVs assigned to the same species, the reads of each of those species at each site were summed together. Column names: sum.taxonomy (assigned taxonomy from Anacapa Toolkit's Bayesian Lowest Common Ancestor classifier, full taxonomic path : Domain; Phylum;Class;Order;Family;Genus;Species).,	columns of Sample data [file naming convention is barcode_samplename e.g. c19_fishcard_LA11.1.R1.S25.L001] (counts of sequences for each ASV).

          Directory 3.b.2: elas_fishcard_12S_all_taxonomy_tables/
          Description Directory 3.b.1:  Anacapa Toolkit Output files for eDNA metabarcoding data generated from MiFish elasmobranch primers and processed with the regional CRUX generated 12S reference database generated in from Gold et al. 2021(See https://doi.org/10.5068/D1H963). The most relevant file here is ~/Summary_by_percent_confidence/60/fishcard_12S_all_ASV_raw_taxonomy_60.txt which is the final ASV table with taxonomy assigned with a confidence cutoff score of 60 used for our analyses within the manuscript.

              File 3.b.2.a Name: fishcard_12S_ASV_taxonomy_detailed.txt
              File 3.b.2.a Description: Detailed output file from Anacapa Toolkit. Columns are as follows: fishcard_12S_seq_number (ASV name),	sequence (DNA sequence of ASV),	sequencesF (forward sequence of unmerged ASVs, otherwise empty),	sequencesR (reverse sequence of unmerged ASVs, otherwise empty),	forward_fishcard_12S_seq_number (column of zeros, not important),	columns of Sample data [file naming convention is barcode_samplename e.g. fishcard_12S_LA11.1.R1.S25.L001] (counts of sequences for each ASV),	merged_fishcard_12S_seq_number column of zeros, not important),	reverse_fishcard_12S_seq_numbercolumn of zeros, not important),	unmerged_fishcard_12S_seq_numbercolumn of zeros, not important),	single_or_multiple_hit (single or multiple sequence alignment hits from bowtie2 sequence alignment),	end_to_end_or_local(end to end or local sequence alignment hits from bowtie2 sequence alignment),	max_percent_id (maximum percent sequence alignment ID from bowtie2 sequence alignment),	input_sequence_length (ASV sequence length), 	taxonomy (assigned taxonomy from Anacapa Toolkit's Bayesian Lowest Common Ancestor classifier, full taxonomic path : Domain; Phylum;Class;Order;Family;Genus;Species), 	taxonomy_confidence (Taxonomic path and Bayesian cutoff score (bcc) for each taxonomic level Domain;Domain bcc; Phylum;Phylum bcc;Class;Class bcc;Order;Order bcc;Family;Family bcc;Genus;Genus bcc;Species;Species bcc;),	accessions (accessions of top aligned sequences).

              File 3.b.2.b Name: fishcard_12S_ASV_taxonomy_brief.txt
              File 3.b.2.b Description: Brief output file from Anacapa Toolkit. Columns are as follows: fishcard_12S_seq_number (ASV name),	columns of Sample data [file naming convention is barcode_samplename e.g. fishcard_12S_LA11.1.R1.S25.L001] (counts of sequences for each ASV),	taxonomy (assigned taxonomy from Anacapa Toolkit's Bayesian Lowest Common Ancestor classifier, full taxonomic path : Domain; Phylum;Class;Order;Family;Genus;Species), 	taxonomy_confidence (Taxonomic path and Bayesian cutoff score (bcc) for each taxonomic level Domain;Domain bcc; Phylum;Phylum bcc;Class;Class bcc;Order;Order bcc;Family;Family bcc;Genus;Genus bcc;Species;Species bcc;),	accessions (accessions of top aligned sequences).

              Directory 3.b.2.a: Summary_by_percent_confidence/
                  Description Directory 3.b.2.a: Taxonomy tables filtered by Bayesian cutoff score. For example, ~/60/ has only retained taxonomic assignments with bayesian cutoff scores > 60. Directories include 40,50,60,70,80,90,95,100 each which is filtered at directory labeled cutoff score. Only the 60 directory is described in detail below but this applies to all subfolders within this directory.

                    Directory 3.b.2.a.1: 60/

                        File 3.b.2.a.1.a Name: fishcard_12S_ASV_raw_taxonomy_60.txt
                        File 3.b.2.a.1.a.1 Description: fishcard_12S_seq_number (ASV name),	columns of Sample data [file naming convention is barcode_samplename e.g. fishcard_12S_LA11.1.R1.S25.L001] (counts of sequences for each ASV),	taxonomy (assigned taxonomy from Anacapa Toolkit's Bayesian Lowest Common Ancestor classifier, full taxonomic path : Domain; Phylum;Class;Order;Family;Genus;Species).

                        File 3.b.2.a.1.b Name: fishcard_12S_ASV_sum_by_taxonomy_60.txt
                        File 3.b.2.a.1.b Description: This table  has grouped and summarized reads for all ASVs with the same taxonomic path. For example, if there were 4 ASVs assigned to the same species, the reads of each of those species at each site were summed together. Column names: sum.taxonomy (assigned taxonomy from Anacapa Toolkit's Bayesian Lowest Common Ancestor classifier, full taxonomic path : Domain; Phylum;Class;Order;Family;Genus;Species).,	columns of Sample data [file naming convention is barcode_samplename e.g. fishcard_12S_LA11.1.R1.S25.L001] (counts of sequences for each ASV).

          Directory 3.b.3: miu_c19_fishcard_taxonomy_tables/
          Description Directory 3.b.1:  Anacapa Toolkit Output files for eDNA metabarcoding data generated from MiFish Universal Teleost primers and processed with the global CRUX generated 12S reference database generated in from Gold et al. 2021(See https://doi.org/10.5068/D1H963). The most relevant file here is ~/Summary_by_percent_confidence/60/c19_fishcard_ASV_raw_taxonomy_60.txt which is the final ASV table with taxonomy assigned with a confidence cutoff score of 60 used for our analyses within the manuscript.

              File 3.b.3.a Name: c19_fishcard_ASV_taxonomy_detailed.txt
              File 3.b.3.a Description: Detailed output file from Anacapa Toolkit. Columns are as follows: c19_fishcard_seq_number (ASV name),	sequence (DNA sequence of ASV),	sequencesF (forward sequence of unmerged ASVs, otherwise empty),	sequencesR (reverse sequence of unmerged ASVs, otherwise empty),	forward_c19_fishcard_seq_number (column of zeros, not important),	columns of Sample data [file naming convention is barcode_samplename e.g. c19_fishcard_LA11.1.R1.S25.L001] (counts of sequences for each ASV),	merged_c19_fishcard_seq_number column of zeros, not important),	reverse_c19_fishcard_seq_numbercolumn of zeros, not important),	unmerged_c19_fishcard_seq_numbercolumn of zeros, not important),	single_or_multiple_hit (single or multiple sequence alignment hits from bowtie2 sequence alignment),	end_to_end_or_local(end to end or local sequence alignment hits from bowtie2 sequence alignment),	max_percent_id (maximum percent sequence alignment ID from bowtie2 sequence alignment),	input_sequence_length (ASV sequence length), 	taxonomy (assigned taxonomy from Anacapa Toolkit's Bayesian Lowest Common Ancestor classifier, full taxonomic path : Domain; Phylum;Class;Order;Family;Genus;Species), 	taxonomy_confidence (Taxonomic path and Bayesian cutoff score (bcc) for each taxonomic level Domain;Domain bcc; Phylum;Phylum bcc;Class;Class bcc;Order;Order bcc;Family;Family bcc;Genus;Genus bcc;Species;Species bcc;),	accessions (accessions of top aligned sequences).

              File 3.b.3.b Name: c19_fishcard_ASV_taxonomy_brief.txt
              File 3.b.3.b Description: Brief output file from Anacapa Toolkit. Columns are as follows: c19_fishcard_seq_number (ASV name),	columns of Sample data [file naming convention is barcode_samplename e.g. c19_fishcard_LA11.1.R1.S25.L001] (counts of sequences for each ASV),	taxonomy (assigned taxonomy from Anacapa Toolkit's Bayesian Lowest Common Ancestor classifier, full taxonomic path : Domain; Phylum;Class;Order;Family;Genus;Species), 	taxonomy_confidence (Taxonomic path and Bayesian cutoff score (bcc) for each taxonomic level Domain;Domain bcc; Phylum;Phylum bcc;Class;Class bcc;Order;Order bcc;Family;Family bcc;Genus;Genus bcc;Species;Species bcc;),	accessions (accessions of top aligned sequences).

              Directory 3.b.3.a: Summary_by_percent_confidence/
              Description Directory 3.b.3.a: Taxonomy tables filtered by Bayesian cutoff score. For example, ~/60/ has only retained taxonomic assignments with bayesian cutoff scores > 60. Directories include 40,50,60,70,80,90,95,100 each which is filtered at directory labeled cutoff score. Only the 60 directory is described in detail below but this applies to all subfolders within this directory.

                  Directory 3.b.3.a.1: 60/

                    File 3.b.3.a.1.a Name: c19_fishcard_ASV_raw_taxonomy_60.txt
                    File 3.b.3.a.1.a.1 Description: c19_fishcard_seq_number (ASV name),	columns of Sample data [file naming convention is barcode_samplename e.g. c19_fishcard_LA11.1.R1.S25.L001] (counts of sequences for each ASV),	taxonomy (assigned taxonomy from Anacapa Toolkit's Bayesian Lowest Common Ancestor classifier, full taxonomic path : Domain; Phylum;Class;Order;Family;Genus;Species).

                    File 3.b.3.a.1.b Name: c19_fishcard_ASV_sum_by_taxonomy_60.txt
                    File 3.b.3.a.1.b Description: This table  has grouped and summarized reads for all ASVs with the same taxonomic path. For example, if there were 4 ASVs assigned to the same species, the reads of each of those species at each site were summed together. Column names: sum.taxonomy (assigned taxonomy from Anacapa Toolkit's Bayesian Lowest Common Ancestor classifier, full taxonomic path : Domain; Phylum;Class;Order;Family;Genus;Species).,	columns of Sample data [file naming convention is barcode_samplename e.g. c19_fishcard_LA11.1.R1.S25.L001] (counts of sequences for each ASV).

          Directory 3.b.4: miu_fishcard_12S_all_taxonomy_tables/
          Description Directory 3.b.1:  Anacapa Toolkit Output files for eDNA metabarcoding data generated from MiFish Universal Teleost primers and processed with the regional CRUX generated 12S reference database generated in from Gold et al. 2021(See https://doi.org/10.5068/D1H963). The most relevant file here is ~/Summary_by_percent_confidence/60/fishcard_12S_all_ASV_raw_taxonomy_60.txt which is the final ASV table with taxonomy assigned with a confidence cutoff score of 60 used for our analyses within the manuscript.

              File 3.b.4.a Name: fishcard_12S_ASV_taxonomy_detailed.txt
              File 3.b.4.a Description: Detailed output file from Anacapa Toolkit. Columns are as follows: fishcard_12S_seq_number (ASV name),	sequence (DNA sequence of ASV),	sequencesF (forward sequence of unmerged ASVs, otherwise empty),	sequencesR (reverse sequence of unmerged ASVs, otherwise empty),	forward_fishcard_12S_seq_number (column of zeros, not important),	columns of Sample data [file naming convention is barcode_samplename e.g. fishcard_12S_LA11.1.R1.S25.L001] (counts of sequences for each ASV),	merged_fishcard_12S_seq_number column of zeros, not important),	reverse_fishcard_12S_seq_numbercolumn of zeros, not important),	unmerged_fishcard_12S_seq_numbercolumn of zeros, not important),	single_or_multiple_hit (single or multiple sequence alignment hits from bowtie2 sequence alignment),	end_to_end_or_local(end to end or local sequence alignment hits from bowtie2 sequence alignment),	max_percent_id (maximum percent sequence alignment ID from bowtie2 sequence alignment),	input_sequence_length (ASV sequence length), 	taxonomy (assigned taxonomy from Anacapa Toolkit's Bayesian Lowest Common Ancestor classifier, full taxonomic path : Domain; Phylum;Class;Order;Family;Genus;Species), 	taxonomy_confidence (Taxonomic path and Bayesian cutoff score (bcc) for each taxonomic level Domain;Domain bcc; Phylum;Phylum bcc;Class;Class bcc;Order;Order bcc;Family;Family bcc;Genus;Genus bcc;Species;Species bcc;),	accessions (accessions of top aligned sequences).

              File 3.b.4.b Name: fishcard_12S_ASV_taxonomy_brief.txt
              File 3.b.4.b Description: Brief output file from Anacapa Toolkit. Columns are as follows: fishcard_12S_seq_number (ASV name),	columns of Sample data [file naming convention is barcode_samplename e.g. fishcard_12S_LA11.1.R1.S25.L001] (counts of sequences for each ASV),	taxonomy (assigned taxonomy from Anacapa Toolkit's Bayesian Lowest Common Ancestor classifier, full taxonomic path : Domain; Phylum;Class;Order;Family;Genus;Species), 	taxonomy_confidence (Taxonomic path and Bayesian cutoff score (bcc) for each taxonomic level Domain;Domain bcc; Phylum;Phylum bcc;Class;Class bcc;Order;Order bcc;Family;Family bcc;Genus;Genus bcc;Species;Species bcc;),	accessions (accessions of top aligned sequences).

              Directory 3.b.4.a: Summary_by_percent_confidence/
              Description Directory 3.b.4.a: Taxonomy tables filtered by Bayesian cutoff score. For example, ~/60/ has only retained taxonomic assignments with bayesian cutoff scores > 60. Directories include 40,50,60,70,80,90,95,100 each which is filtered at directory labeled cutoff score. Only the 60 directory is described in detail below but this applies to all subfolders within this directory.

                Directory 3.b.4.a.1: 60/

                    File 3.b.4.a.1.a Name: fishcard_12S_ASV_raw_taxonomy_60.txt
                    File 3.b.4.a.1.a.1 Description: fishcard_12S_seq_number (ASV name),	columns of Sample data [file naming convention is barcode_samplename e.g. fishcard_12S_LA11.1.R1.S25.L001] (counts of sequences for each ASV),	taxonomy (assigned taxonomy from Anacapa Toolkit's Bayesian Lowest Common Ancestor classifier, full taxonomic path : Domain; Phylum;Class;Order;Family;Genus;Species).

                    File 3.b.4.a.1.b Name: fishcard_12S_ASV_sum_by_taxonomy_60.txt
                    File 3.b.4.a.1.b Description: This table  has grouped and summarized reads for all ASVs with the same taxonomic path. For example, if there were 4 ASVs assigned to the same species, the reads of each of those species at each site were summed together. Column names: sum.taxonomy (assigned taxonomy from Anacapa Toolkit's Bayesian Lowest Common Ancestor classifier, full taxonomic path : Domain; Phylum;Class;Order;Family;Genus;Species).,	columns of Sample data [file naming convention is barcode_samplename e.g. fishcard_12S_LA11.1.R1.S25.L001] (counts of sequences for each ASV).

      Directory 3.c: fish/

        File 3.c.1 Name: cal_nemo_10102019.txt
        File 3.c.1 Description: List of invasive species in California as listed in the CalNEMO database from October 2019.

          Directory 3.c.1: CO1_taxonomy_tables/
          Description Directory 3.c.1:  Anacapa Toolkit Output files for eDNA metabarcoding data generated from Leray CO1 primers and processed with the CRUX generated Co1 reference database described above. The most relevant file here is ~/Summary_by_percent_confidence/60/CO1_ASV_raw_taxonomy_60.txt which is the final ASV table with taxonomy assigned with a confidence cutoff score of 60 used for our analyses within the manuscript.

            File 3.c.1.a Name: CO1_ASV_taxonomy_detailed.txt
            File 3.c.1.a Description: Detailed output file from Anacapa Toolkit. Columns are as follows: CO1_seq_number (ASV name),	sequence (DNA sequence of ASV),	sequencesF (forward sequence of unmerged ASVs, otherwise empty),	sequencesR (reverse sequence of unmerged ASVs, otherwise empty),	forward_CO1_seq_number (column of zeros, not important),	columns of Sample data [file naming convention is barcode_samplename e.g. CO1_LA11.1.R1.S25.L001] (counts of sequences for each ASV),	merged_CO1_seq_number column of zeros, not important),	reverse_CO1_seq_numbercolumn of zeros, not important),	unmerged_CO1_seq_numbercolumn of zeros, not important),	single_or_multiple_hit (single or multiple sequence alignment hits from bowtie2 sequence alignment),	end_to_end_or_local(end to end or local sequence alignment hits from bowtie2 sequence alignment),	max_percent_id (maximum percent sequence alignment ID from bowtie2 sequence alignment),	input_sequence_length (ASV sequence length), 	taxonomy (assigned taxonomy from Anacapa Toolkit's Bayesian Lowest Common Ancestor classifier, full taxonomic path : Domain; Phylum;Class;Order;Family;Genus;Species), 	taxonomy_confidence (Taxonomic path and Bayesian cutoff score (bcc) for each taxonomic level Domain;Domain bcc; Phylum;Phylum bcc;Class;Class bcc;Order;Order bcc;Family;Family bcc;Genus;Genus bcc;Species;Species bcc;),	accessions (accessions of top aligned sequences).

            File 3.c.2.b Name: CO1_ASV_taxonomy_brief.txt
            File 3.c.2.b Description: Brief output file from Anacapa Toolkit. Columns are as follows: CO1_seq_number (ASV name),	columns of Sample data [file naming convention is barcode_samplename e.g. CO1_LA11.1.R1.S25.L001] (counts of sequences for each ASV),	taxonomy (assigned taxonomy from Anacapa Toolkit's Bayesian Lowest Common Ancestor classifier, full taxonomic path : Domain; Phylum;Class;Order;Family;Genus;Species), 	taxonomy_confidence (Taxonomic path and Bayesian cutoff score (bcc) for each taxonomic level Domain;Domain bcc; Phylum;Phylum bcc;Class;Class bcc;Order;Order bcc;Family;Family bcc;Genus;Genus bcc;Species;Species bcc;),	accessions (accessions of top aligned sequences).

            Directory 3.c.1.a: Summary_by_percent_confidence/
            Description Directory 3.c.1.a: Taxonomy tables filtered by Bayesian cutoff score. For example, ~/60/ has only retained taxonomic assignments with bayesian cutoff scores > 60. Directories include 40,50,60,70,80,90,95,100 each which is filtered at directory labeled cutoff score. Only the 60 directory is described in detail below but this applies to all subfolders within this directory.

              Directory 3.c.1.a.1: 60/

                File 3.c.1.a.1 Name: CO1_ASV_raw_taxonomy_60.txt
                File 3.c.1.a.1 Description: CO1_seq_number (ASV name),	columns of Sample data [file naming convention is barcode_samplename e.g. CO1_LA11.1.R1.S25.L001] (counts of sequences for each ASV),	taxonomy (assigned taxonomy from Anacapa Toolkit's Bayesian Lowest Common Ancestor classifier, full taxonomic path : Domain; Phylum;Class;Order;Family;Genus;Species).

                File 3.c.1.a.2 Name: CO1_ASV_sum_by_taxonomy_60.txt
                File 3.c.1.a.2 Description: This table  has grouped and summarized reads for all ASVs with the same taxonomic path. For example, if there were 4 ASVs assigned to the same species, the reads of each of those species at each site were summed together. Column names: sum.taxonomy (assigned taxonomy from Anacapa Toolkit's Bayesian Lowest Common Ancestor classifier, full taxonomic path : Domain; Phylum;Class;Order;Family;Genus;Species).,	columns of Sample data [file naming convention is barcode_samplename e.g. CO1_LA11.1.R1.S25.L001] (counts of sequences for each ASV).

          Directory 3.c.2: Metazoa_16s_taxonomy_tables/
          Description Directory 3.c.1:  Anacapa Toolkit Output files for eDNA metabarcoding data generated from Kelly 16S primers and processed with the  CRUX generated Metazoa_16S reference database described above. The most relevant file here is ~/Summary_by_percent_confidence/60/Metazoa_16s_ASV_raw_taxonomy_60.txt which is the final ASV table with taxonomy assigned with a confidence cutoff score of 60 used for our analyses within the manuscript.

            File 3.c.2.a Name: Metazoa_16s_ASV_taxonomy_detailed.txt
            File 3.c.2.a Description: Detailed output file from Anacapa Toolkit. Columns are as follows: Metazoa_16s_seq_number (ASV name),	sequence (DNA sequence of ASV),	sequencesF (forward sequence of unmerged ASVs, otherwise empty),	sequencesR (reverse sequence of unmerged ASVs, otherwise empty),	forward_Metazoa_16s_seq_number (column of zeros, not important),	columns of Sample data [file naming convention is barcode_samplename e.g. Metazoa_16s_LA11.1.R1.S25.L001] (counts of sequences for each ASV),	merged_Metazoa_16s_seq_number column of zeros, not important),	reverse_Metazoa_16s_seq_numbercolumn of zeros, not important),	unmerged_Metazoa_16s_seq_numbercolumn of zeros, not important),	single_or_multiple_hit (single or multiple sequence alignment hits from bowtie2 sequence alignment),	end_to_end_or_local(end to end or local sequence alignment hits from bowtie2 sequence alignment),	max_percent_id (maximum percent sequence alignment ID from bowtie2 sequence alignment),	input_sequence_length (ASV sequence length), 	taxonomy (assigned taxonomy from Anacapa Toolkit's Bayesian Lowest Common Ancestor classifier, full taxonomic path : Domain; Phylum;Class;Order;Family;Genus;Species), 	taxonomy_confidence (Taxonomic path and Bayesian cutoff score (bcc) for each taxonomic level Domain;Domain bcc; Phylum;Phylum bcc;Class;Class bcc;Order;Order bcc;Family;Family bcc;Genus;Genus bcc;Species;Species bcc;),	accessions (accessions of top aligned sequences).

            File 3.c.2.b Name: Metazoa_16s_ASV_taxonomy_brief.txt
            File 3.c.2.b Description: Brief output file from Anacapa Toolkit. Columns are as follows: Metazoa_16s_seq_number (ASV name),	columns of Sample data [file naming convention is barcode_samplename e.g. Metazoa_16s_LA11.1.R1.S25.L001] (counts of sequences for each ASV),	taxonomy (assigned taxonomy from Anacapa Toolkit's Bayesian Lowest Common Ancestor classifier, full taxonomic path : Domain; Phylum;Class;Order;Family;Genus;Species), 	taxonomy_confidence (Taxonomic path and Bayesian cutoff score (bcc) for each taxonomic level Domain;Domain bcc; Phylum;Phylum bcc;Class;Class bcc;Order;Order bcc;Family;Family bcc;Genus;Genus bcc;Species;Species bcc;),	accessions (accessions of top aligned sequences).

            Directory 3.c.2.a: Summary_by_percent_confidence/
            Description Directory 3.c.2.a: Taxonomy tables filtered by Bayesian cutoff score. For example, ~/60/ has only retained taxonomic assignments with bayesian cutoff scores > 60. Directories include 40,50,60,70,80,90,95,100 each which is filtered at directory labeled cutoff score. Only the 60 directory is described in detail below but this applies to all subfolders within this directory.

            Directory 3.c.2.a.1: 60/

              File 3.c.2.a.1 Name: Metazoa_16s_ASV_raw_taxonomy_60.txt
              File 3.c.2.a.1 Description: Metazoa_16s_seq_number (ASV name),	columns of Sample data [file naming convention is barcode_samplename e.g. Metazoa_16s_LA11.1.R1.S25.L001] (counts of sequences for each ASV),	taxonomy (assigned taxonomy from Anacapa Toolkit's Bayesian Lowest Common Ancestor classifier, full taxonomic path : Domain; Phylum;Class;Order;Family;Genus;Species).

              File 3.c.2.a.2 Name: Metazoa_16s_ASV_sum_by_taxonomy_60.txt
              File 3.c.2.a.2 Description: This table  has grouped and summarized reads for all ASVs with the same taxonomic path. For example, if there were 4 ASVs assigned to the same species, the reads of each of those species at each site were summed together. Column names: sum.taxonomy (assigned taxonomy from Anacapa Toolkit's Bayesian Lowest Common Ancestor classifier, full taxonomic path : Domain; Phylum;Class;Order;Family;Genus;Species).,	columns of Sample data [file naming convention is barcode_samplename e.g. Metazoa_16s_LA11.1.R1.S25.L001] (counts of sequences for each ASV).

### How to convert data bases into Anacapa Formatted databases
    Directory 1.b: CO1_bowtie2_database/
    Directory 1.b Description: This is a derived intermediate data product, bowtie2 formatted reference database, that is too large to store on GitHub. This bowtie2(https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)formatted reference database (.bt2) is used by the Anacapa Toolkit for fast sequence alignments during the Anacapa Classifier step. To visualize these files first install bowtie2 (using homebrew, brew install bowtie2). Files are then directly derived from Files 1.a.1 & 1.a.2 using the following script: bowtie2-build -f ~/CO1_.fasta ~/CO1/CO1_bowtie2_database/CO1_bowtie2_index . The following code allows for inspection of the raw fasta file: bowtie2-inspect ~/CO1/CO1_bowtie2_database/CO1_bowtie2_index . Warning this is a large file and such inspection may use a lot of computer memory.

    Directory 2.b: Metazoa_16s_bowtie2_database/
    Directory 2.b Description:his is a derived intermediate data product, bowtie2 formatted reference database, that is too large to store on GitHub. This bowtie2(https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)formatted reference database (.bt2) is used by the Anacapa Toolkit for fast sequence alignments during the Anacapa Classifier step. To visualize these files first install bowtie2 (using homebrew, brew install bowtie2). Files are then directly derived from Files 2.a.1 & 2.a.2 using the following script: bowtie2-build -f ~/Metazoa_16s_.fasta ~/Metazoa_16s/Metazoa_16s_bowtie2_database/Metazoa_16s_bowtie2_index . The following code allows for inspection of the raw fasta file: bowtie2-inspect ~/Metazoa_16s/Metazoa_16s_bowtie2_database/Metazoa_16s_bowtie2_index . Warning this is a large file and such inspection may use a lot of computer memory.
