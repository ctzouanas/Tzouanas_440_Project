Overview: This repo contains the data and code used by Constantine Tzouanas and Jeff Hsiao as part of the project in MIT's 20.440 course, Spring 2020. More specifically, we are focusing on immune cells captured in Smillie et al.'s single-cell RNA-seq (scRNA-seq) datasets from their paper "Intra- and Inter-cellular Rewiring of the Human Colon during Ulcerative Colitis". To understand shifts in cellular abundances and cell states associated with the development of ulcerative colitis, this analysis specifically focuses on B cells. At present, the code is sufficient to conduct quality control of the scRNA-seq dataset, dimensionality reduction for summarization using PCA, clustering through k-nearest neighbors, and dimensionality reduction for visualization using UMAP, then to split the dataset based on the patient's health status (i.e., health control or ulcerative colitis). As part of our code (implemented in R), we use Seurat v3 (courtesy of Stuart et al., 2019). 
	Citations: 
		Smillie, C.S., et al. "Intra- and Inter-cellular Rewiring of the Human Colon during Ulcerative Colitis," Cell, 178(3), 714-730 (2019). DOI: 10.1016/j.cell.2019.06.029.
		Stuart, T., et al. "Comprehensive Integration of Single-Cell Data," Cell, 177(7), 1888-1902 (2019). DOI: 10.1016/j.cell.2019.05.031

Data: Data was generated by taking 68 colon biopsies from 30 patients, split across 18 patients with ulcerative colitis and 12 healthy controls. Single-cell RNA-sequencing was conducted using the commercially-available 10X platform, and the immune subset of their data was accessed through the Broad Institute's Single-Cell Portal, via the "Download" tab of the following link (https://singlecell.broadinstitute.org/single_cell/study/SCP259/intra-and-inter-cellular-rewiring-of-the-human-colon-during-ulcerative-colitis#study-summary). Due to limitations on Github's allowable file sizes, "gene_sorted-Imm.matrix.mtx" must be downloaded into the "RawData_Folder". "RawData_Folder" should be in the same folder/save location as the code file "Simillie_UC_v1.R" (mirroring the file structure found on this repo). 
	Citations:
		Smillie, C.S., et al. "Intra- and Inter-cellular Rewiring of the Human Colon during Ulcerative Colitis," Cell, 178(3), 714-730 (2019). DOI: 10.1016/j.cell.2019.06.029.

Folder Structure: The main folder of this repo contains three hidden files, this README, a .R code file, and three subfolders. 
	Hidden files: .gitignore specifies types of files that should be ignored by git in the process of pushing edits/changes to/from the central repository. .gitattributes establishes specific protocols for handling particular file types. .DS_Store is a vestigial file, hidden on the local folder used to generate this repo. 
	README: Description file used to convey the various components of this repo. 
	.R file: Code needed to reproduce the analysis conducted for the 20.440 project completed by Tzouanas and Hsiao. Please note the description in the "Data" section of this README, which specifies that "gene_sorted-Imm.matrix.mtx" must be downloaded into the "RawData_Folder", which should exist in the same folder as this .R code file. Additionally, for proper saving of graphs and R data objects, the subfolders "Plot_Folder" and "Object_Folder" must also exist in the same folder as this .R code file. 
	Subfolders: "RawData_Folder" contains the raw data downloaded from the Broad Institute's Single-Cell Portal, as described in the "Data" section of this README. "Plot_Folder" is used as the save location for output graphs from the .R code file, while "ProcessedData_Folder" is used as the save location for output RData and Seurat objects from the .R code file. 


Installation: To run this code, please download the .R code file from this repo. At your desired save location for this .R code file, please also download "Plot_Folder", "Object_Folder", and "RawData_Folder". Please download "gene_sorted-Imm.matrix.mtx" from https://singlecell.broadinstitute.org/single_cell/study/SCP259/intra-and-inter-cellular-rewiring-of-the-human-colon-during-ulcerative-colitis#study-summary, and save it into "RawData_Folder". In your preferred method for running R code, please ensure that the following libraries are installed: Seurat, data.table, ggplot2, gridExtra, gtable, plotly, scales, spam, dplyr, and patchwork. If any of them are not already installed, they can be installed by running "install.packages("<the package's name>")" in an R session command line. This code is based on Seurat v3.1.4. 