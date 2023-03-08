# github_guide
Wijzigingen aan drug screen overzichten

Je neemt een drug screen report geproduceerd door TECAN 300D
Dit is een .xml file dat geproduceerd wordt met het uitvoeren van het drugprint experiment.

In dit .xml file (openen in excel of spreadsheet programma naar keuze) gebruik je de tabular tab. 

Hiervan neem je de volgende kolommen: 
(C: “Dispensed well”, D: “Dispensed row”, E: “Dispensed col”, I:”Fluid name”, J-?: “Concentration, en alle losse concentraties”, dan van volume DMSO normalization t/m DMSO %. Tussenliggende volumes laat je voor wat het is). 
Daaraan voeg je toe: “Organoid”, “Timepoint”, “Value”, “value_corr”, “GR”, “conc_condition”. 
Daaraan verander je de “fluid name” naar “condition”, “Conc. (Um) … drug” verander je in “conc_drug”. 

Alle conditions veranderd adhv script zoals hieronder: 

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
