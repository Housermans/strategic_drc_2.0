## ENGLISH
# Convert drug screen overview so that it works in this script:
Changes to drug screen overviews

You take a drug screen report produced by TECAN 300D This is an .xml file that is produced with the execution of the drugprint experiment.

In this .xml file (open in excel or spreadsheet program of your choice) you use the tabular tab.

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

All conditions changed using an R script as below:

#create variable with fluid names for combination treatment
D$condition <- ifelse(D$condition == "2 Fluids" & !is.na(D$`conc_SN-38`) & !is.na(D$conc_CHEK1), "SN38_CHEK1", D$condition)
D$condition <- ifelse(D$condition == "2 Fluids" & !is.na(D$conc_Lapatinib) & !is.na(D$conc_Binimetinib), "binimetinib_lapatinib", D$condition) 
D$condition <- ifelse(D$condition == "2 Fluids" & !is.na(D$conc_Alpelisib) & !is.na(D$conc_Binimetinib), "binimetinib_alpelisib", D$condition) 
D$condition <- ifelse(D$condition == "2 Fluids" & !is.na(D$conc_Lapatinib) & !is.na(D$conc_Alpelisib), "alpelisib_lapatinib", D$condition) 
D$condition <- ifelse(D$condition == "2 Fluids" & !is.na(D$conc_Vinorelbine) & !is.na(D$conc_Navitoclax), "navitoclax_vinorelbine", D$condition)
D$condition <- ifelse(D$condition == "3 Fluids", "vi_bi_la", D$condition) 

The conc_condition you make as follows in excel: =IF(NOT(ISBLANK(K4));K4;IF(I4=“vi_bi_la”;T4;MAX(K4:T4))) Where K4: Concentration I4: condition name T4: conc_vinorelbine So, in the case of vi_bi_la (two concentrations the same, 1 varies) - the varying Otherwise, the maximum of all.

Don’t forget that if you add something to the structure of the table, you change this in all other tables and in the ‘standard_file’ otherwise you can’t connect them.

## DUTCH
# Drug screen overzicht omzetten zodat het werkt in dit script:
Wijzigingen aan drug screen overzichten

Je neemt een drug screen report geproduceerd door TECAN 300D
Dit is een .xml file dat geproduceerd wordt met het uitvoeren van het drugprint experiment.

In dit .xml file (openen in excel of spreadsheet programma naar keuze) gebruik je de tabular tab. 

Hiervan kopieer je de volgende kolommen naar een nieuw excel bestand: 
C: “Dispensed well”
D: “Dispensed row”,
E: “Dispensed col”
Vanaf I: "Fluid name" alle kolommen die starten met "Conc. "
De laatste 4 kolommen: "Volume (nL) DMSO normalization", 	"Volume (nL) a+Tw normalization",  "Total well volume (nL)", "DMSO \%" 

1. voeg de volgende kolommen toe voor Dispensedwell: “Organoid”, “Timepoint”, 
	- Plaats een X voor Organoid (zodat het duidelijk is dat deze een string gaat bevatten) 
	- Vul elk vakje van Timepoint met D5
2. Voeg de volgende kolommen toe voor condition: “Value”, “value_corr”, “GR”, 
	- plaats de value 0,01 in elke cel van de eerste rij van de drie kolommen zodat het duidelijk is dat deze een numerieke (float) waarde gaat bevatten
3. Voeg de volgende kolom to toe na Concentration: “conc_condition”. 
4. Verander de naam van kolom “fluid name” naar “condition”, 
5. Verander de naam van kolom “Conc. (Um) <drug>”  naar “conc_<drug>”.
6. Verander de naam van kolom "DMSO \%" naar "DMSO_pct"

Alle conditions veranderd adhv een R script zoals hieronder: 

#create variable with fluid names for combination treatment
D$condition <- ifelse(D$condition == "2 Fluids" & !is.na(D$`conc_SN-38`) & !is.na(D$conc_CHEK1), "SN38_CHEK1", D$condition)
D$condition <- ifelse(D$condition == "2 Fluids" & !is.na(D$conc_Lapatinib) & !is.na(D$conc_Binimetinib), "binimetinib_lapatinib", D$condition) 
D$condition <- ifelse(D$condition == "2 Fluids" & !is.na(D$conc_Alpelisib) & !is.na(D$conc_Binimetinib), "binimetinib_alpelisib", D$condition) 
D$condition <- ifelse(D$condition == "2 Fluids" & !is.na(D$conc_Lapatinib) & !is.na(D$conc_Alpelisib), "alpelisib_lapatinib", D$condition) 
D$condition <- ifelse(D$condition == "2 Fluids" & !is.na(D$conc_Vinorelbine) & !is.na(D$conc_Navitoclax), "navitoclax_vinorelbine", D$condition)
D$condition <- ifelse(D$condition == "3 Fluids", "vi_bi_la", D$condition) 

De conc_condition die maak je als volgt in excel: 
=IF(NOT(ISBLANK(K4));K4;IF(I4="vi_bi_la";T4;MAX(K4:T4)))
Waarin K4: Concentration
I4: condition naam
T4: conc_vinorelbine
Dus, in het geval van vi_bi_la (twee concentraties hetzelfde, 1 varieert) - de variërende
Anders, de maximum van allemaal. 

Niet vergeten dat je, als je wat toevoegt aan de structuur van de tabel, dat je dit wijzigt in alle andere tabellen en in de ‘standard_file’ anders kan je ze niet aan elkaar knopen. 
