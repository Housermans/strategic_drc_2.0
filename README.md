## ENGLISH
# Convert drug screen overview so that it works in this script:
Changes to drug screen overviews

You take a drug screen report produced by TECAN 300D This is an .xml file that is produced with the execution of the drugprint experiment.

In this .xml file (open in excel or spreadsheet program of your choice) you use the tabular tab.

Delete all rows below 386 (so that you read one plate at a time)

From this you copy the following columns to a new excel file: C: “Dispensed well” D: “Dispensed row”, E: “Dispensed col” From I: “Fluid name” all columns starting with "Conc. " The last 4 columns: “Volume (nL) DMSO normalization”, “Volume (nL) a+Tw normalization”, “Total well volume (nL)”, “DMSO \%”

add the following columns for Dispensedwell: “Organoid”, “Timepoint”,
Place an X for Organoid (so that it is clear that this will contain a string)
Fill each cell of Timepoint with D5
Add the following columns for condition: “Value”, “value_corr”, “GR”,
place the value 0.01 in each cell of the first row of the three columns so that it is clear that this will contain a numerical (float) value
Add the following column to after Concentration: “conc_condition”.
Change the name of column “fluid name” to “condition”,
Change the name of column “Conc. (Um) <drug>” to “conc_<drug>”. 
(use excel formula: =SUBSTITUTE(M1;CHAR(13);"") to get rid of the new line character)
Change the name of column “DMSO \%” to “DMSO_pct”

All conditions changed using a the R script 'Script0'

Now add the data of your experiment to the screening_overview Excel file following the example of my prefilled file. 

When you did, add the spectramax readout files to the first folder, add a folder with your experiment ID, and make sure the control file is labeled: ctrl, and the 
experimental files are labeled: d5_<organoid_name>