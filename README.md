# IMLS_woodland_restoration

A repository that contains data and R scripts for running analyses performed in the paper "High-diversity seed additions promote herb-layer recovery 
during restoration of degraded oak woodland" written by Kaul et al., published in Ecological Solutions and Evidence (currently in review).

The files in this repository consists of 2 data sets and 1 R script. The files are described as follows:

IMLS_analysis_code.r - An R script which performs all analyses and produces all figures presented in Kaul et al. It uses the following data:

IMLS_all_species_raw.csv - A dataset containing relative abundance cover data for plant species across 5 samples of 10 transects, with 5 quadrats per transect. Species are rows. The spreadsheet contains the following columns:

`
sp_id	- species identificaiton number. Integers 1-356. Each taxon has a unique sp_id including those not identified to species level.\
family - taxonomic designation of plant family based on Ladd and Thomas 2015 for Missouri.
seeded_1_2	- factor with levels "yes", or "no" indicating if a species was added to the seed mix for transects 1 through 2.\
seeded_3_6	- factor with levels "yes", or "no" indicating if a species was added to the seed mix for transects 3 through 6.\
seeded_7_10	- factor with levels "yes", or "no" indicating if a species was added to the seed mix for transects 7 through 10.\
native	- factor with levels "yes", or "no" indicating if a species is native to the study region in MO, USA.\
duration	- factor with levels "annual","biennial", "perennial", "woody","shrub", or "tree" indicating duration of life form.\
fg - Functional group.	Factor with levels "forb","grass", "legume", "vine","tree", "fern", "sedge","vine"\
c_score	- conservatism score for the state of Missouri based on Ladd and Thomas 2015\
conservatism	- factor with levels "ruderal", "matrix", or "conservative". Ruderal (C=0-3), matrix (C=4-6), conservative (C=7-10).\
latin	- latin couplet\
taxon - unique 6-letter code generated for each taxon in the data set.\
syn_ladd_thomas - notes some common latin synonyms for some taxa.\
01-A.2017.Spring - all remaining columns denote exact percent cover estimates for each species. Column naming convention is transect-quadrat.year.season.\
    The first column thus denotes cover in transect 01, quadrat A (of posible A-E), sampled in 2017, in the spring.\
`

IMLS_seed_mixes.csv - A dataset containing compositional data for plant species sown in each of three mixes developed for restoraiton of herb diversity in a midwestern US oak woodland.

`
family - taxonomic designation of plant family based on Ladd and Thomas 2015 for Missouri.\
latin	- latin couplet\
taxon - unique 6-letter code generated for each taxon in the data set.\
t_1_2_oz - bulk mass in ounces (seed plus some other inert material such as stems) of seed added to the mix for transects 1-2.\
t_3_6_oz - bulk mass in ounces of seed added to the mix for transects 3-6.\
t_7_10_oz - bulk mass in ounces of seed added to the mix for transects 7-10.\
transects_seeded - integer value. 1-5. how many of the 5 seeded transects received each species.\
quadrats_seeded - integer value. 5,10,15,20,or 25. how many of the 25 seeded quadrats received each species.\
seeds_oz - estimate of the number of seeds in an ounce.\
source - factor with levels "stored","fresh", or "combo". Denotes the source and storage of seeds for each species. See manuscript for details.\
`
