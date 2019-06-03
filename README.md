# Internship_M2
## Stores the relevant files and codes used during my M2 internship in the Musée de l'Homme, MNHN.
### Files and codes to get files done to run analyses:
- <b>CapeVerde_PeoplingHistory_ABCpriors_RefTable_FM_jan2019.xlsx</b> (file storing the demographic data we gathered and used to infer priors for our simulations with PLGS. All the data and references I personally found are in pale yellow cells)

- <b>make_files_for_Arlequin.R</b> (Make files in a proper format to be used in Arlequin for AMOVA and Fst claculations)

- <b>make_gen_dist_files.R</b> (To create ASD')


### Codes to get results in the order of the "results" part of my report:
- <b>R_mantel_age.R</b> (correlations and Mantel tests WSD/age, ASD/age)

- <b>R_Mantel_birthVSling.R</b> (correlations and Mantel tests WSD/birthplace, WSD/residence place. Partial Mantel tests between WSD, birthplace and residence place)

- <b>R_ling_partial_mantel_birthVSparent_birth.R</b> (correlations and Mantel tests between WSD and each parental birthplace. Partial mantel tests betwen WSD and birthplace, adjusting for each parental birthplace)

- <b>R_gen_ling_mantel_genVS_words.R</b> (correlation and Mantel test ASD/WSD. Partial Mantel test adjusting for the birthplace)

- <b>R_mantel_genVSbirthVSlivingplace.R</b> (correlations and Mantel tests ASD/birthplace, WSD/residence place. Partial Mantel tests between ASD, birthplace and residence place)

- <b>R_mantel_genVSparentbirth.R</b> (correlations and Mantel tests between ASD and each parental birthplace. Partial mantel tests betwen ASD and birthplace, adjusting for each parental birthplace)

- <b>R_MCA.R</b> (MCA analyses and plot on linguistic data)

- <b>R_lingNJ.R</b> (to get the NJ tree of linguistic Fst between populations by birth island)

- <b>R_FstByIslands.R</b> (Graphical representation of linguistic Fst between populations by birth island)

- <b>PLGS_param_maio.cpp</b> (externalisation of PLGS parameters)

- <b>PLGS_Scenario1_Split_Ling.cpp</b> (PLGS: Simulation of Scenario 1 (split between Maio and Santiago) for linguistic data)

- <b>PLGS_Scenario2_Merge_Ling.cpp</b> (PLGS: Simulation of Scenario 2 (Maio and Santiago are indistinguishable populations) for linguistic data)

- <b>R_stats_simul_ling.R</b> (Results of ABC analyses for the simulations of linguistic history)

- <b>PLGS_Scenario1_Split.cpp</b> (PLGS: Simulation of Scenario 1 (split between Maio and Santiago) for genetic and linguistic data)

- <b>PLGS_Scenario2_Merge.cpp</b> (PLGS: Simulation of Scenario 2 (Maio and Santiago are indistinguishable populations) for genetic and linguistic data)

- <b>R_stats_simul_gen_and_ling.R</b> ((Results of ABC analyses for the simulations of genetic and linguistic histories)
