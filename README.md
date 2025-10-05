This repository contains the data and code to reproduce the results in the paper: 

Title: "Temperature-dependent supply and demand drive shrinking body sizes across levels of ecological organisation"
Author(s): Mikael Pontarp, Jason Griffiths, Pedro Rosero, Florian Altermatt, Yves Choffat, John DeLong, Aurelie Garnier, Suzanne Greene, Thomas Massie, Gian-Marco Palmara, Mathew Seymour, Owen L Petchey, and Frank Pennekamp

Instructions:
- The Data folder contains the raw data as a single (compressed) csv file (Individual_cell_size_data.csv.gz), as well as a csv file with the metadata (metadata.csv). 
- We used the R package drake to create a reproducible workflow for the statistical analysis. Use the plan.R file to run the pipeline. All functions are in the R folder (funs.R).
- Main and supplementary figures and tables can be found in the output folder.

<a href="https://doi.org/10.5281/zenodo.17273578"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.17273578.svg" alt="DOI"></a>
