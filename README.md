# Flow cytometry analysis using GatingSet in R

This code is an R script that can be used to analyze flow cytometry data for both the Broad and Koch Institute flow cytometers. It performs a series of gating steps to identify cells of interest and generates plots/spreadsheets to quantify the reporter protein expression.

## Prerequisites

To run the code, you need to have R and the following libraries installed: flowCore, flowWorkspace, ggcyto, knitr, dplyr, and broom.

## Running the code

1. Open the R console and set the working directory to the location of the script.
2. Modify the directory path in the code to specify the location of the fcs files.
3. Run the script.
4. The script will generate several plots, which will be saved as pdf files in the same directory.
5. The final plot will show the expression of the constitutive reporter gene, BFP, in each sample. This plot can be used to compare the reporter gene expression across different samples.

## Output

The script generates several plots and a CSV file containing fluorescence and gating statistics, with the number and percentage of cells in relevant gates. The plots show gates and the expression of the constitutive reporter gene, BFP, in each sample.

## Notes

The coordinates of the gating polygons and the reporter gene expression range should be adjusted based on the data being analyzed. These parameters are specified in the comments in the code and should be modified according to the user's needs.