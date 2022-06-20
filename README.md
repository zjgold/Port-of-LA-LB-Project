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

    /analysis

      * *ports_data_analysis_20220619.Rmd* This script does the main analyses in the paper.

  /data/fish

      * *20210420_ports_decontamination.R* This script processes the "raw" *Anacapa Toolkit* output and runs decontamination to remove poorly sequenced samples and contaminant ASVs. Used on the MiFish *12S* Universal Teleost and Elasmobranch data.

  /data/invert
      * *220220130_ports_CO1_16S_decontamination.R* This script processes the "raw" *Anacapa Toolkit* output and runs decontamination to remove poorly sequenced samples and contaminant ASVs.Used on the *CO1* and *16S* metazoan data.

    /data/Anacapa_output_old
      * *decontam.R* This script processes the "raw" *Anacapa Toolkit* output and runs decontamination to remove poorly sequenced samples and contaminant ASVs. Used on the MiFish *12S* Universal Teleost and Elasmobranch data that was processed with an older reference database that lacked hundreds of CA species barcodes.

2. Data

  1. *abundance_trawl* Trawl data in counts of individual fish.
  2. *biomass_trawl* Trawl data in biomass of fish.
  3. *la_ports_meta_data.txt* Metadata of eDNA and trawl data associated with stations/bottles sampled.

  /fish

    1. *miu_c19_fishcard_taxonomy_tables* *Anacapa Toolkit* Output of MiFish *12S* Universal Teleost data using the *CRUX* Global Reference database from [Gold et al. 2021](https://onlinelibrary.wiley.com/doi/epdf/10.1111/1755-0998.13450)       
    2. *miu_fishcard_12S_all_taxonomy_tables* *Anacapa Toolkit* Output of MiFish *12S* data using the *CRUX* Local Reference database from Gold et al. 2021. [See GitHub](https://github.com/zjgold/FishCARD) & [See Gold et al. 2021 Dryad](https://doi.org/10.5068/D1H963)
    3. *elas_c19_fishcard_taxonomy_tables* *Anacapa Toolkit* Output of MiFish *12S* Elasmobranch data using the *CRUX* Global Reference database.   
    4. *elas_fishcard_12S_all_taxonomy_tables* *Anacapa Toolkit* Output of MiFish 12S data using the *CRUX* Local Reference database.
    5. *decontamination/* Directory with decontamination outputs.
    6. *CA_fish_list_20210216.csv* List of California fishes obtained from Gold et al. 2021

  /Anacapa_output_old

  1. *12S_U* *Anacapa Toolkit* Output of MiFish *12S* Universal Teleost data using an old *CRUX* Reference database with no additional CA species barcodes added.[See Anacapa Dryad](https://doi.org/10.1111/2041-210x.13214)       
  2. *12S_E* *Anacapa Toolkit* Output of MiFish *12S* Elasmobranch data using an old *CRUX* Reference database with no additional CA species barcodes added.
  3. *decontamination/* Directory with decontamination outputs.

  /invert

  1. *CO1_taxonomy_tables* *Anacapa Toolkit* Output of  *CO1* metazoan data using a *CRUX* Reference database with additional marine barcodes added from the Mo'orea BIOCODE project.       
  2. *Metazoa_16s_taxonomy_tables* *Anacapa Toolkit* Output of MiFish *16S* Elasmobranch data using an old *CRUX* Reference database with no additional species barcodes added.
  3. *decontamination/* Directory with decontamination outputs.


**Note** *A few files needed for analyses are not immediately included on GitHub due to size, please see [See Dryad]INSERT LINK*

**Raw sequence inputs, reference databases, and Anacapa scripts  are made available on Dryad**
