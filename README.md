## Internship_M2
#Stores the relevant files and codes used during my M2 internship in the Musée de l'Homme, MNHN.
#Files and codes to get files done to run analyses:
- CapeVerde_PeoplingHistory_ABCpriors_RefTable_FM_jan2019.xlsx (file storing the demographic data we gathered and used to infer priors for our simulations with PLGS. All the data and references I personally found are in pale yellow cells)

- make_files_for_Arlequin.R (Make files in a proper format to be used in Arlequin for AMOVA and Fst claculations)

- make_gen_dist_files.R (To create ASD')


#Codes to get results in the order of "results" part in my report:
- R_mantel_age.R (correlations and Mantel tests WSD/age, ASD/age)

- R_Mantel_birthVSling.R (correlations and Mantel tests WSD/birthplace, WSD/residence place. Partial Mantel tests between WSD, birthplace and residence place)

- R_ling_partial_mantel_birthVSparent_birth.R (correlations and Mantel tests between WSD and each parental birthplace. Partial mantel tests betwen WSD and birthplace, adjusting for each parental birthplace)

- R_gen_ling_mantel_genVS_words.R (correlation and Mantel test ASD/WSD. Partial Mantel test adjusting for the birthplace)

- R_mantel_genVSbirthVSlivingplace.R (correlations and Mantel tests ASD/birthplace, WSD/residence place. Partial Mantel tests between ASD, birthplace and residence place)

- R_mantel_genVSparentbirth.R (correlations and Mantel tests between ASD and each parental birthplace. Partial mantel tests betwen ASD and birthplace, adjusting for each parental birthplace)

- R_MCA.R (MCA analyses and plot on linguistic data)

- R_lingNJ.R (to get the NJ tree of linguistic Fst between populations by birth island)

- R_FstByIslands.R (Graphical representation of linguistic Fst between populations by birth island)

- PLGS_param_maio.cpp (externalisation of PLGS parameters)

- PLGS_Scenario1_Split_Ling.cpp (PLGS: Simulation of Scenario 1 (split between Maio and Santiago) for linguistic data)

- PLGS_Scenario2_Merge_Ling.cpp (PLGS: Simulation of Scenario 2 (Maio and Santiago are indistinguishable populations) for linguistic data)

- R_stats_simul_ling.R (Results of ABC analyses for the simulations of linguistic history)

- PLGS_Scenario1_Split.cpp (PLGS: Simulation of Scenario 1 (split between Maio and Santiago) for genetic and linguistic data)

- PLGS_Scenario2_Merge.cpp (PLGS: Simulation of Scenario 2 (Maio and Santiago are indistinguishable populations) for genetic and linguistic data)

- R_stats_simul_gen_and_ling.R ((Results of ABC analyses for the simulations of genetic and linguistic histories)
