# Single cell quantification of protein regulation in primary tissue


* [Project Website](https://decode.slavovlab.net) &nbsp; 
* [Preprint article](https://doi.org/10.1101/2024.08.26.609665)

The code in this repository was used to quantify the regulation of protein abundance via transcription, translation, and protein clearance across cell types and single cells via integrated analysis of single cell in vivo metabolic pulse proteomics and singel cell mRNA sequencing. 


## Processed data availible via Zenodo

Outputs from each stage of analysis and processed gene by single cell and cell type specific data matricies for each modality are located in the "dat" folder on [zenodo](google.com). A summary of the folder's contents can be found also on [zenodo](google.com).

## Reproducing analysis

### Overview of data processing scripts (numbered as in code directory)

#### Preprocessing
1. Preprocessing of single cell metabolic pulse data
2. Integrating proteomics data and annotating cell types
3. Calculating protein clearance rates accounting for amio acid recycling
4. Accounting for missing data in computing cell type specific protein abundance
5. Integrating mRNA and proteomics data

#### Downstream analysis

6. Analysis of factors explaining absolute protein concentrations within cell type and the influence of cell growth
7. Analyis of influence of mRNA abundance, translation, and protein clearance on relative protein abundance across cell types.
8. Single cell covariation analysis within cell type
    a. mRNA-Protein comparison
    b. Across cell type comparison


### Instructions for reproducing results

#### If you want to reproduce only the downstream analysis 

1. Download the code this github repository.

2. Downoad the data folder from Zenodo, unzip and add the contents to the "dat" folder.

3. Run scripts 07, 08, and 09, order does not matter. Relevant functions for each scripts analysis are found at the beginning of each script.



#### If you want to reproduce the preprocessing

1. Download the code this github repository.

2. Downoad the data folder from Zenodo, unzip and add the contents to the "dat" folder.

3. Download the DIANN processed raw data from the "searched" section of MassIVE [MSV000093494](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=ac44b779d8a04ca285a263616796c3b8). Unzip and add the contents to the "raw_data" folder.


4. Preprocessing of the single cell data requires installing the QuantQC R package. For instructions installing go [here](google.com). 


5. Run Scripts in order starting with 4 individual preprocessing scripts in folder 01 and save data outputs throughout each script

6. After getting to 




-------------

