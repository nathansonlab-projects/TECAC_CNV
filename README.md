# TECAC_CNV
This project aims at identifying germline Copy Number Variations (CNVs) asscociated with TCGT.

## Project Flow Chart ##
![CNV_Project_Flow](https://user-images.githubusercontent.com/58447038/113348280-65ae3500-9304-11eb-966c-75169ded59f3.png)

## Processing Steps ##
* **1. Run PennCNV** to generate .rawcnv files <br />
* **2. Run format_rawcnv.sh** to format .rawcnv file for further processing <br />
* **3. Run subsetRawcnv.sh** to split .rawcnv file into duplication and deletion files <br />
* **4. Run runCNVRassoc.R** to run association testing and gene/cnvr mapping <br />
* **5. Run plotGeneInCNVR.R** to visualize CNVR/gene overlap <br />
<br />

## Tools ##
* Shell <br />
* **PennCNV** <br /> 
* **GenoCN** <br />
* **ParseCNV** <br /> 

## Supporting Files ##
* **TECAC_CNV_SUBLIST.txt**: The file containing BID and SSID for all samples that will be retained for CNV Calling <br />
<br />

* **TECAC_CNV_PHENO.txt**: Contains BID, SSID, and phenotype for all TECAC samples <br /> 
<br />

* **tecac.evec**: principle components for TECAC data, used as covariates <br />
<br />

* **header_all.txt**: The file containing all column names of the large, unsplit raw signal intensity file. <br />
<br />
