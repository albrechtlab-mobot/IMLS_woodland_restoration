#### data processing and analysis for IMLS woodland restoration project ####
# last updated 12/02/2022


#### packages ####

getwd()

library(tidyverse) # for data wrangling
library(ggplot2) # for figure generation
library(vegan) # for calculating diversity indices
library(nlme) # for glmms
library(emmeans) # for extracting model predictions
library(DHARMa) # for analysis of model fit
library(scales) # for working with colors  
library(ggplotify) # for making a legend into its own ggplot
library(mvnormtest) # for testing multivariate normality 
library(reshape) # for melt() function
library(rstatix) # for function get_summary_stats()
library(ggpubr) # for function ggarrange()
library(goeveg) # for dimcheckMDS() function
library(coxed) # for bca() function used for alt. method of bootstrapping perMANOVA
library(car) # for leveneTest() for homogeneity of variance across groups
# library(glmmTMB) # for betabinomial glmms

#### color palette ####

#create color blind friendly color palette
cb_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#999999","#000000")
#scales::show_col(cb_palette)

cb_palette_2 <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7")
#scales::show_col(cb_palette_2)

cb_palette_3 <- c("#D55E00","#56B4E9", "#999999")
#scales::show_col(cb_palette_3)

cb_palette_4 <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7","#330099","#999999","#000000" )
#scales::show_col(cb_palette_4)


#### import data ####

imls_full_raw<- read.csv("G:/Postdoc Projects/Project 1 - IMLS/IMLS Analyses/IMLS_Analysis/IMLS_all_species_raw.csv",
                         header = T, 
                         stringsAsFactors=TRUE)

imls_seed_mixes<- read.csv("G:/Postdoc Projects/Project 1 - IMLS/IMLS Analyses/IMLS_Analysis/IMLS_seed_mixes.csv",
                         header = T, 
                         ) %>%
  select(-latin, -family,) %>%
  mutate(seeds_oz = as.integer(seeds_oz))


#### data processing 1 - separating seeded and sampled species ####

# 355 rows represent taxa sampled at any point or any species that were sown in any transect
# colnames(imls_full_raw)
temp1 <- select(imls_full_raw, -latin, -syn_ladd_thomas) %>% 
    # create column summing across all quadrats
    mutate(total = select(., X01.A.2017.Spring:X10.E.2019.Fall) %>% 
    rowSums(na.rm = TRUE))

not_sampled_spp <- temp1$total==0
sum(not_sampled_spp)
# 98 species were sown and never sampled

sampled_spp <- filter(temp1, total>0)
# 257 taxa were recorded across 5 samples over 3 years
rm(temp1)
rm(not_sampled_spp)

# extract only sown species
seeded_spp <- filter(imls_full_raw, 
                     seeded_1_2 == "yes" |
                     seeded_3_6 == "yes" |
                     seeded_7_10 == "yes") %>%
              select(-latin, -syn_ladd_thomas) %>%
              filter(taxon != "elysp") 
# note, we remove the Elymus sp. group which is treated as seeded for sampled community analyses
  
nrow(seeded_spp) # 169 species were included in at least one of the three mixes

#### data processing 2 - sampled community data ####

# create long-form data
sampled_long <- melt(sampled_spp,) %>% 
                filter(value>0) %>% # remove rows for zero counts
  
                # help the melt function by removing extra rows it melted
                filter(variable != "sp_id" & variable != "total" & variable != "c_score") %>% 
               
                # melt made two new columns called variable and value
                mutate(plot = variable, cover = value) %>% 
                select(-variable, -value) %>%
                separate(col=plot, into= c('transect','quadrat','year','season'))

# remove the "X" in transect that showed up upon data import
sampled_long$transect <- substring(sampled_long$transect, 2)      


# calculate max cover for each species in each quadrat within year
year_max <- group_by(sampled_long, year,transect, quadrat, taxon, native) %>% 
  summarise(max_cover = max(cover)
            ) %>%
  ungroup()


#calculate the proportional abundance of unknowns species to ensure it is low in all samples

unknowns <- group_by(year_max, year,transect)%>%
  summarise(sum_max = sum(max_cover),
            natives = sum(max_cover[native == "yes"],na.rm=T),
            exotics = sum(max_cover[native == "no"],na.rm=T),
            unknowns = sum(max_cover[is.na(native)]),
            knowns = sum_max - unknowns,
            prop_unknown = unknowns / sum_max
            )%>%
  ungroup()
mean(unknowns$prop_unknown) # average is < 1% in each transect in each sample


# create column for number of quadrats per transect. This accounts for 9D missing in 2018 and 2019
# remove rows for all unknown taxa (not identified to species level)
# dropped from 3031 to 2935 rows
year_max <- year_max %>%
  mutate(num_quadrat = as.numeric(ifelse(
  transect == "09" & year == "2018"| transect == "09" & year=="2019", "4","5"))) %>%
  drop_na(native)


year_max_wide <- year_max %>% 
  select(-native) %>%
  pivot_wider(names_from = taxon, values_from = max_cover)

# calculate mean total, and native quadrat (plot) richness within each transect in each year
quadrat_s <- group_by(year_max, year, transect) %>%
  summarise(
  num_quadrat = length(unique(quadrat)),
  mean_quad_s = length(taxon) / length(unique(quadrat)),
  mean_nat_quad_s = (length(na.omit(taxon)[native =="yes"])) / mean(num_quadrat),
  #mean_quad_s = length(taxon) / mean(num_quadrat) #alternate way of calculating this
  ) %>%
  ungroup()

# calculate average cover for each species in each transect in each year
transect_mean <- group_by(year_max, year, transect, taxon) %>%
  summarise(
  mean_t_cover = sum(max_cover) / mean(num_quadrat), # accounts for only 4 quadrats in transect 9 in years 2018, 2019
  ) %>%
  ungroup() 
transect_mean$transect <-as.factor(transect_mean$transect)

transect_mean_native <- group_by(year_max, year, transect, native, taxon) %>%
  summarise(
    mean_t_cover = sum(max_cover) / mean(num_quadrat), # accounts for only 4 quadrats in transect 9 in years 2018, 2019
  ) %>%
  ungroup() %>%
  filter(native == "yes") %>%
  select(-native)
transect_mean_native$transect <-as.factor(transect_mean_native$transect)


total_cover <- group_by(transect_mean,year,transect) %>%
  summarise(sum_cover = sum(mean_t_cover)
  ) %>%
  ungroup()


# use vegan to calculate Simpson's effective species number for each transect for all species, or natives only
transect_mean_wide <- pivot_wider(transect_mean, names_from = taxon , values_from = mean_t_cover )
transect_mean_wide[is.na(transect_mean_wide)] <- 0 # convert NA cells to zeros 
transect_vegan<- unite(transect_mean_wide,"trans_sample", c("year", "transect"), remove=T) %>%
  data.frame( row.names=1) #combine two columns to one, and then assign it as the row names

s <- as.numeric(apply(transect_vegan>0,1,sum))
simp_eff_spp_num <- as.numeric(diversity(transect_vegan, index = "invsimpson"))
pileou_j <- as.numeric(diversity(transect_vegan, index = "shannon")/log(s))
shan_h <-as.numeric(diversity(transect_vegan, index = "shannon"))
shan_eff_spp_num <- exp(shan_h)

transect_mean_wide_native <- pivot_wider(transect_mean_native, names_from = taxon , values_from = mean_t_cover )
transect_mean_wide_native[is.na(transect_mean_wide_native)] <- 0 # convert NA cells in year_max_wide to zeros 
transect_vegan_native<- unite(transect_mean_wide_native,"trans_sample", c("year", "transect"), remove=T) %>%
  data.frame( row.names=1) #combine two columns to one, and then assign it as the row names 
 
s_native <- as.numeric(apply(transect_vegan_native>0,1,sum))
simp_eff_spp_num_native <- as.numeric(diversity(transect_vegan_native, index = "invsimpson"))
pileou_j_native <- as.numeric(diversity(transect_vegan_native, index = "shannon")/log(s))
shan_h_native <- as.numeric(diversity(transect_vegan_native, index = "shannon"))
shan_eff_spp_num_native <- exp(shan_h_native)

# separate data into 3 groups based on management unit / seed mixes used

transects_1_2 <- filter(transect_mean, transect == "01" | transect == "02" )
nrow(transects_1_2) # 325

transects_3_6 <- filter(transect_mean, transect == "03" | transect == "04" | transect == "05" | transect == "06")
nrow(transects_3_6) # 567

transects_7_10 <- filter(transect_mean, transect == "07" | transect == "08" | transect == "09" | transect == "10")
nrow(transects_7_10) # 462


# merge species traits into our 3 data frames corresponding to management areas
transects_1_2_merged <- merge(x=transects_1_2, y=imls_full_raw[,c( "taxon",
                        "seeded_1_2", "native", "duration", "fg", "c_score",
                        "conservatism")], by="taxon", all.x =T) %>%
  merge(y=total_cover, by.x=c("year","transect"), by.y=c("year","transect")) %>%
  mutate(prop = mean_t_cover / sum_cover,
          weighted_c = prop * c_score )

transects_3_6_merged <- merge(x=transects_3_6, y=imls_full_raw[,c( "taxon",
                        "seeded_3_6", "native", "duration", "fg", "c_score",
                        "conservatism")], by="taxon", all.x =T) %>%
  merge(y=total_cover, by.x=c("year","transect"), by.y=c("year","transect")) %>%
  mutate(prop = mean_t_cover / sum_cover,
         weighted_c = prop * c_score )

transects_7_10_merged <- merge(x=transects_7_10, y=imls_full_raw[,c( "taxon",
                        "seeded_7_10", "native", "duration", "fg", "c_score",
                        "conservatism")], by="taxon", all.x =T) %>%
  merge(y=total_cover, by.x=c("year","transect"), by.y=c("year","transect")) %>%
  mutate(prop = mean_t_cover / sum_cover,
         weighted_c = prop * c_score )

# check that relative abundances sum to 1
group_by(transects_1_2_merged, year, transect)%>%
  summarise(
   sum_prop =sum(prop) 
  )%>%
  ungroup()

group_by(transects_3_6_merged, year, transect)%>%
  summarise(
    sum_prop =sum(prop) 
  )%>%
  ungroup()

group_by(transects_7_10_merged, year, transect)%>%
  summarise(
    sum_prop =sum(prop) 
  )%>%
  ungroup()

# within each transect in each year, calculate community-level responses
# NOTE: when using length() function you must use na.omit() within.
# you cannot use "na.rm=T" and R WILL count NA rows even if you specify to only count rows with Col val == ""

transects_1_2_summary <- group_by(transects_1_2_merged, year, transect) %>%
  summarise(
    total_s = length(taxon),
    native_s = length(na.omit(taxon[native == "yes"])),
    exotic_s = length(na.omit(taxon[native == "no"])),
    unknown_s = length(taxon[is.na(native)]),    
    known_s = total_s - unknown_s,
        
    total = sum(mean_t_cover),
    natives = sum(mean_t_cover[native == "yes"],na.rm=T),
    exotics = sum(mean_t_cover[native == "no"],na.rm=T),
    unknowns = sum(mean_t_cover[is.na(native)]),
    knowns = total - unknowns,
    # note: total cover = native + exotic + unknown
    
    prop_unknown = unknowns / total,
 
    # unknown species are excluded from all analyses   
    prop_native = natives / knowns,
    prop_exotic = exotics / knowns,
    
    grass = sum(mean_t_cover[fg == "grass"], na.rm=T),
    sedge = sum(mean_t_cover[fg == "sedge"], na.rm=T),
    forb = sum(mean_t_cover[fg == "forb"], na.rm=T),
    fern = sum(mean_t_cover[fg == "fern"], na.rm=T),
    legume = sum(mean_t_cover[fg == "legume"], na.rm=T),
    vine = sum(mean_t_cover[fg == "vine"], na.rm=T),
    shrub = sum(mean_t_cover[fg == "shrub"], na.rm=T),
    tree = sum(mean_t_cover[fg == "tree"], na.rm=T),
    
    grass_s = length(na.omit(taxon [fg == "grass"])),
    sedge_s = length(na.omit(taxon [fg == "sedge"])),
    forb_s = length(na.omit(taxon [fg == "forb"])),
    fern_s = length(na.omit(taxon [fg == "fern"])),
    legume_s = length(na.omit(taxon [fg == "legume"])),
    vine_s = length(na.omit(taxon [fg == "vine"])),
    shrub_s = length(na.omit(taxon [fg == "shrub"])),
    tree_s = length(na.omit(taxon [fg == "tree"])),
    
    prop_grass = grass / knowns,
    prop_sedge = sedge / knowns,
    prop_forb = forb / knowns,
    prop_fern = fern / knowns,
    prop_legume = legume / knowns,
    prop_vine = vine / knowns,
    prop_shrub = shrub / knowns,
    prop_tree = tree / knowns,
    sum_fg = prop_grass+prop_sedge+prop_forb+prop_fern+prop_legume+prop_vine+prop_shrub+prop_tree,
    
    mean_c = mean(c_score,na.rm=T), # this is native mean_c; exotics have NA for c_score
    cwmc = (sum(weighted_c[!is.na(weighted_c)])/native_s),
    fqi = mean_c * sqrt(native_s),
    weighted_fqi = cwmc * sqrt(native_s),
    
    ruderal = sum(mean_t_cover[conservatism == "ruderal"], na.rm=T),
    matrix = sum(mean_t_cover[conservatism == "matrix"], na.rm=T),
    conservative = sum(mean_t_cover[conservatism == "conservative"], na.rm=T),
    
    ruderal_s = length(na.omit(taxon[conservatism == "ruderal"])),
    matrix_s = length(na.omit(taxon[conservatism == "matrix"])),
    conservative_s = length(na.omit(taxon[conservatism == "conservative"])),
    
    prop_ruderal = ruderal / natives,
    prop_matrix = matrix / natives,
    prop_conservative = conservative / natives,
    sum_conservatism = prop_ruderal+prop_matrix+prop_conservative,
    
    seeded = sum(mean_t_cover[seeded_1_2 == "yes"]),
    prop_seeded = seeded / knowns,
    seeded_s = length(na.omit(taxon[seeded_1_2 == "yes"])),
  ) %>%
  ungroup()


transects_3_6_summary <- group_by(transects_3_6_merged, year, transect) %>%
  summarise(
    total_s = length(taxon),
    native_s = length(na.omit(taxon[native == "yes"])),
    exotic_s = length(na.omit(taxon[native == "no"])),
    unknown_s = length(taxon[is.na(native)]),    
    known_s = total_s - unknown_s,
    
    total = sum(mean_t_cover),
    natives = sum(mean_t_cover[native == "yes"],na.rm=T),
    exotics = sum(mean_t_cover[native == "no"],na.rm=T),
    unknowns = sum(mean_t_cover[is.na(native)]),
    knowns = total - unknowns,
    # note: total cover = native + exotic + unknown
    
    prop_unknown = unknowns / total,
    
    # unknown species are excluded from all analyses   
    prop_native = natives / knowns,
    prop_exotic = exotics / knowns,
    
    grass = sum(mean_t_cover[fg == "grass"], na.rm=T),
    sedge = sum(mean_t_cover[fg == "sedge"], na.rm=T),
    forb = sum(mean_t_cover[fg == "forb"], na.rm=T),
    fern = sum(mean_t_cover[fg == "fern"], na.rm=T),
    legume = sum(mean_t_cover[fg == "legume"], na.rm=T),
    vine = sum(mean_t_cover[fg == "vine"], na.rm=T),
    shrub = sum(mean_t_cover[fg == "shrub"], na.rm=T),
    tree = sum(mean_t_cover[fg == "tree"], na.rm=T),
    
    grass_s = length(na.omit(taxon [fg == "grass"])),
    sedge_s = length(na.omit(taxon [fg == "sedge"])),
    forb_s = length(na.omit(taxon [fg == "forb"])),
    fern_s = length(na.omit(taxon [fg == "fern"])),
    legume_s = length(na.omit(taxon [fg == "legume"])),
    vine_s = length(na.omit(taxon [fg == "vine"])),
    shrub_s = length(na.omit(taxon [fg == "shrub"])),
    tree_s = length(na.omit(taxon [fg == "tree"])),
    
    prop_grass = grass / knowns,
    prop_sedge = sedge / knowns,
    prop_forb = forb / knowns,
    prop_fern = fern / knowns,
    prop_legume = legume / knowns,
    prop_vine = vine / knowns,
    prop_shrub = shrub / knowns,
    prop_tree = tree / knowns,
    sum_fg = prop_grass+prop_sedge+prop_forb+prop_fern+prop_legume+prop_vine+prop_shrub+prop_tree,
    
    mean_c = mean(c_score,na.rm=T), # this is native mean_c; exotics have NA for c_score
    cwmc = (sum(weighted_c[!is.na(weighted_c)])/native_s),
    fqi = mean_c * sqrt(native_s),
    weighted_fqi = cwmc * sqrt(native_s),
    
    ruderal = sum(mean_t_cover[conservatism == "ruderal"], na.rm=T),
    matrix = sum(mean_t_cover[conservatism == "matrix"], na.rm=T),
    conservative = sum(mean_t_cover[conservatism == "conservative"], na.rm=T),
    
    ruderal_s = length(na.omit(taxon[conservatism == "ruderal"])),
    matrix_s = length(na.omit(taxon[conservatism == "matrix"])),
    conservative_s = length(na.omit(taxon[conservatism == "conservative"])),
    
    prop_ruderal = ruderal / natives,
    prop_matrix = matrix / natives,
    prop_conservative = conservative / natives,
    sum_conservatism = prop_ruderal+prop_matrix+prop_conservative,
    
    seeded = sum(mean_t_cover[seeded_3_6 == "yes"]),
    prop_seeded = seeded / knowns,
    seeded_s = length(na.omit(taxon[seeded_3_6 == "yes"])),
  ) %>%
  ungroup()


transects_7_10_summary <- group_by(transects_7_10_merged, year, transect) %>%
  summarise(
    total_s = length(taxon),
    native_s = length(na.omit(taxon[native == "yes"])),
    exotic_s = length(na.omit(taxon[native == "no"])),
    unknown_s = length(taxon[is.na(native)]),    
    known_s = total_s - unknown_s,
    
    total = sum(mean_t_cover),
    natives = sum(mean_t_cover[native == "yes"],na.rm=T),
    exotics = sum(mean_t_cover[native == "no"],na.rm=T),
    unknowns = sum(mean_t_cover[is.na(native)]),
    knowns = total - unknowns,
    # note: total cover = native + exotic + unknown
    
    prop_unknown = unknowns / total,
    
    # unknown species are excluded from all analyses   
    prop_native = natives / knowns,
    prop_exotic = exotics / knowns,
    
    grass = sum(mean_t_cover[fg == "grass"], na.rm=T),
    sedge = sum(mean_t_cover[fg == "sedge"], na.rm=T),
    forb = sum(mean_t_cover[fg == "forb"], na.rm=T),
    fern = sum(mean_t_cover[fg == "fern"], na.rm=T),
    legume = sum(mean_t_cover[fg == "legume"], na.rm=T),
    vine = sum(mean_t_cover[fg == "vine"], na.rm=T),
    shrub = sum(mean_t_cover[fg == "shrub"], na.rm=T),
    tree = sum(mean_t_cover[fg == "tree"], na.rm=T),
    
    grass_s = length(na.omit(taxon [fg == "grass"])),
    sedge_s = length(na.omit(taxon [fg == "sedge"])),
    forb_s = length(na.omit(taxon [fg == "forb"])),
    fern_s = length(na.omit(taxon [fg == "fern"])),
    legume_s = length(na.omit(taxon [fg == "legume"])),
    vine_s = length(na.omit(taxon [fg == "vine"])),
    shrub_s = length(na.omit(taxon [fg == "shrub"])),
    tree_s = length(na.omit(taxon [fg == "tree"])),
    
    prop_grass = grass / knowns,
    prop_sedge = sedge / knowns,
    prop_forb = forb / knowns,
    prop_fern = fern / knowns,
    prop_legume = legume / knowns,
    prop_vine = vine / knowns,
    prop_shrub = shrub / knowns,
    prop_tree = tree / knowns,
    sum_fg = prop_grass+prop_sedge+prop_forb+prop_fern+prop_legume+prop_vine+prop_shrub+prop_tree,
    
    mean_c = mean(c_score,na.rm=T), # this is native mean_c; exotics have NA for c_score
    cwmc = (sum(weighted_c[!is.na(weighted_c)])/native_s),
    fqi = mean_c * sqrt(native_s),
    weighted_fqi = cwmc * sqrt(native_s),
    
    ruderal = sum(mean_t_cover[conservatism == "ruderal"], na.rm=T),
    matrix = sum(mean_t_cover[conservatism == "matrix"], na.rm=T),
    conservative = sum(mean_t_cover[conservatism == "conservative"], na.rm=T),
    
    ruderal_s = length(na.omit(taxon[conservatism == "ruderal"])),
    matrix_s = length(na.omit(taxon[conservatism == "matrix"])),
    conservative_s = length(na.omit(taxon[conservatism == "conservative"])),
    
    prop_ruderal = ruderal / natives,
    prop_matrix = matrix / natives,
    prop_conservative = conservative / natives,
    sum_conservatism = prop_ruderal+prop_matrix+prop_conservative,
    
    seeded = sum(mean_t_cover[seeded_7_10 == "yes"]),
    prop_seeded = seeded / knowns,
    seeded_s = length(na.omit(taxon[seeded_7_10 == "yes"])),
  ) %>%
  ungroup()


# rbind combines data frames with the same columns
# add mean quadrat richness, transect Simpson's diversity, two columns combining ferns with forbs,
# and a column for seed treatment to our summary data frame
all_transects_summary <- rbind(transects_1_2_summary, transects_3_6_summary, transects_7_10_summary) %>%
  merge(y=quadrat_s, by=c("year","transect"), all.x=T) %>%
  mutate(simp = simp_eff_spp_num,  simp_native = simp_eff_spp_num_native,
         j = pileou_j, j_native = pileou_j_native, 
         shannon = shan_eff_spp_num, shannon_native = shan_eff_spp_num_native,
         forb_2 = forb + fern,
         prop_forb_2 = prop_forb + prop_fern,
         treat = if_else(transect == "01" | transect == "03" | transect == "04" | transect == "08"| transect == "10", "seeded", "control" ),
         treat_bin = if_else (treat == "seeded", 1,0),
         year = as.factor(year))

# just shows that the data have unknown taxa removed prior to analysis
unknown_summary <-group_by(all_transects_summary, year) %>%
  summarise(
    max_prop_unknown = max(prop_unknown),
    mean_prop_unknown = mean(prop_unknown)
  ) %>%
  ungroup()
  
summary_2019 <- filter (all_transects_summary, year == "2019")


#### data processing 2.5 - data set for analysis of forb conservatism ####

transects_forb_merged <- merge(x=transect_mean, y=imls_full_raw[,c( "taxon",
                                                                   "native", "duration", "fg", "c_score",
                                                                   "conservatism")], by="taxon", all.x =T) %>%
  merge(y=total_cover, by.x=c("year","transect"), by.y=c("year","transect")) %>%
  mutate(prop = mean_t_cover / sum_cover,
         weighted_c = prop * c_score )

transects_forb_summary <- group_by(transects_forb_merged, year, transect) %>%
  filter(fg=="forb") %>%
  summarise(
    total = sum(mean_t_cover),
    natives = sum(mean_t_cover[native == "yes"],na.rm=T),
    native_s = length(na.omit(taxon[native == "yes"])),
    mean_c = mean(c_score,na.rm=T), #this is native mean_c; exotics have NA for c_score
    cwmc = (sum(weighted_c[!is.na(weighted_c)])/native_s),
    fqi = mean_c * sqrt(native_s),
    weighted_fqi = cwmc * sqrt(native_s),
    ruderal = sum(mean_t_cover[conservatism == "ruderal"], na.rm=T),
    matrix = sum(mean_t_cover[conservatism == "matrix"], na.rm=T),
    conservative = sum(mean_t_cover[conservatism == "conservative"], na.rm=T),
    ruderal_s = length(na.omit(taxon[conservatism == "ruderal"])),
    matrix_s = length(na.omit(taxon[conservatism == "matrix"])),
    conservative_s = length(na.omit(taxon[conservatism == "conservative"])),
    prop_ruderal = ruderal / natives,
    prop_matrix = matrix / natives,
    prop_conservative = conservative / natives,
  ) %>%
  ungroup()


all_transects_forb_summary <- (transects_forb_summary) %>%
  mutate(treat = if_else(transect == "01" | transect == "03" | transect == "04" | transect == "08"| transect == "10", "seeded", "control" ),
         treat_bin = if_else (treat == "seeded", 1,0),
         year = as.factor(year))

summary_forb_2019 <- filter (all_transects_forb_summary, year == "2019")

#### data processing 3 - seed mix data ####

seeded_1_2 <- filter(seeded_spp, seeded_1_2 == "yes") %>%
              select(sp_id, duration, fg, conservatism, c_score,
                    taxon )
transects_1_2_wide <- pivot_wider(transects_1_2, names_from = c("year","transect"), values_from = mean_t_cover)
seeded_1_2_merged <- merge(x=seeded_1_2, y=transects_1_2_wide, by="taxon" , all.x=T) %>% 
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>% #remove NAs for abundance data
  mutate(present_1_2 = 1 , 
         tot_1_2_2017 = `2017_01` + `2017_02`,
         tot_1_2_2018 = `2018_01` + `2018_02`,
         tot_1_2_2019 = `2019_01` + `2019_02`,
         tot_seeded_1_2_2018 = `2018_01`,
         tot_seeded_1_2_2019 = `2019_01`,
         tot_non_seeded_1_2_2018 = `2018_02`,
         tot_non_seeded_1_2_2019 = `2019_02`)
  
seeded_1_2_merged_long <- melt(seeded_1_2_merged, id= c("taxon", "duration", "fg", "conservatism", "c_score", "sp_id")) %>%
  mutate(sample = variable, abundance = value) %>%
  select(-variable, -value)

mix_1_2_summary <- group_by (seeded_1_2_merged_long, sample) %>%
    summarise(
            total_seeded_s = length(abundance [abundance > 0]),
            mean_c = mean(c_score [abundance >0], na.rm=T),
            seeded_ruderal_s = length(na.omit(abundance [conservatism == "ruderal" & abundance > 0])),
            seeded_matrix_s = length(na.omit(abundance [conservatism == "matrix" & abundance > 0])),
            seeded_conservative_s = length(na.omit(abundance [conservatism == "conservative" & abundance > 0])),
            seeded_forb_s = length(na.omit(abundance [fg == "forb" & abundance > 0])),
            seeded_grass_s = length(na.omit(abundance [fg == "grass" & abundance > 0])),
            seeded_sedge_s = length(na.omit(abundance [fg == "sedge" & abundance > 0])),
            seeded_legume_s = length(na.omit(abundance [fg == "legume" & abundance > 0])),
            seeded_an_s = length(na.omit(abundance [duration == "annual" & abundance > 0])),
            seeded_bienn_s = length(na.omit(abundance [duration == "biennial" & abundance > 0])),
            seeded_perenn_s = length(na.omit(abundance [duration == "perennial" & abundance > 0]))
) %>%
  ungroup() %>%
  mutate (mix = "1_2")


seeded_3_6 <- filter(seeded_spp, seeded_3_6 == "yes") %>%
              select(sp_id, duration, fg, conservatism, c_score,
                    taxon )
transects_3_6_wide <- pivot_wider(transects_3_6, names_from = c("year","transect"), values_from = mean_t_cover)
seeded_3_6_merged <- merge(x=seeded_3_6, y=transects_3_6_wide, by="taxon" , all.x=T) %>%
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>% 
  mutate(present_3_6 = 1 , 
         tot_3_6_2017 = `2017_03` + `2017_04` + `2017_05` + `2017_06`,
         tot_3_6_2018 = `2018_03` + `2018_04` + `2018_05` + `2018_06`,
         tot_3_6_2019 = `2019_03` + `2019_04` + `2019_05` + `2019_06`,
         tot_seeded_3_6_2018 = `2018_03` + `2018_04`, # transects 3 & 4 got seed, 5 & 6 did not 
         tot_seeded_3_6_2019 = `2019_03` + `2019_04`,
         tot_non_seeded_3_6_2018 = `2018_05` + `2018_06`,  
         tot_non_seeded_3_6_2019 = `2019_05` + `2019_06`)

seeded_3_6_merged_long <- melt(seeded_3_6_merged, id= c("taxon", "duration", "fg", "conservatism", "c_score", "sp_id")) %>%
  mutate(sample = variable, abundance = value) %>%
  select(-variable, -value)

mix_3_6_summary <- group_by (seeded_3_6_merged_long, sample) %>%
  summarise(
    total_seeded_s = length(abundance [abundance > 0]),
    mean_c = mean(c_score [abundance >0], na.rm=T),
    seeded_ruderal_s = length(na.omit(abundance [conservatism == "ruderal" & abundance > 0])),
    seeded_matrix_s = length(na.omit(abundance [conservatism == "matrix" & abundance > 0])),
    seeded_conservative_s = length(na.omit(abundance [conservatism == "conservative" & abundance > 0])),
    seeded_forb_s = length(na.omit(abundance [fg == "forb" & abundance > 0])),
    seeded_grass_s = length(na.omit(abundance [fg == "grass" & abundance > 0])),
    seeded_sedge_s = length(na.omit(abundance [fg == "sedge" & abundance > 0])),
    seeded_legume_s = length(na.omit(abundance [fg == "legume" & abundance > 0])),
    seeded_an_s = length(na.omit(abundance [duration == "annual" & abundance > 0])),
    seeded_bienn_s = length(na.omit(abundance [duration == "biennial" & abundance > 0])),
    seeded_perenn_s = length(na.omit(abundance [duration == "perennial" & abundance > 0]))
  ) %>%
  ungroup() %>%
  mutate (mix = "3_6")



seeded_7_10 <- filter(seeded_spp, seeded_7_10 == "yes") %>%
              select(sp_id, duration, fg, conservatism, c_score,
                    taxon )
transects_7_10_wide <- pivot_wider(transects_7_10, names_from = c("year","transect"), values_from = mean_t_cover)
seeded_7_10_merged <- merge(x=seeded_7_10, y=transects_7_10_wide, by="taxon" , all.x=T) %>%
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>% #remove NAs for abundance data
  mutate(present_7_10 = 1 , 
         tot_7_10_2017 = `2017_07` + `2017_08` + `2017_09` + `2017_10`,
         tot_7_10_2018 = `2018_07` + `2018_08` + `2018_09` + `2018_10`,
         tot_7_10_2019 = `2019_07` + `2019_08` + `2019_09` + `2019_10`,
         tot_seeded_7_10_2018 = `2018_08` + `2018_10`, # transects 8 & 10 got seed, 7 & 9 did not 
         tot_seeded_7_10_2019 = `2019_08` + `2019_10`,
         tot_non_seeded_7_10_2018 = `2018_07` + `2018_09`,  
         tot_non_seeded_7_10_2019 = `2019_07` + `2019_09`)

seeded_7_10_merged_long <- melt(seeded_7_10_merged, id= c("taxon", "duration", "fg", "conservatism", "c_score", "sp_id")) %>%
  mutate(sample = variable, abundance = value) %>%
  select(-variable, -value)

mix_7_10_summary <- group_by (seeded_7_10_merged_long, sample) %>%
  summarise(
    total_seeded_s = length(abundance [abundance > 0]),
    mean_c = mean(c_score [abundance >0], na.rm=T),
    seeded_ruderal_s = length(na.omit(abundance [conservatism == "ruderal" & abundance > 0])),
    seeded_matrix_s = length(na.omit(abundance [conservatism == "matrix" & abundance > 0])),
    seeded_conservative_s = length(na.omit(abundance [conservatism == "conservative" & abundance > 0])),
    seeded_forb_s = length(na.omit(abundance [fg == "forb" & abundance > 0])),
    seeded_grass_s = length(na.omit(abundance [fg == "grass" & abundance > 0])),
    seeded_sedge_s = length(na.omit(abundance [fg == "sedge" & abundance > 0])),
    seeded_legume_s = length(na.omit(abundance [fg == "legume" & abundance > 0])),
    seeded_an_s = length(na.omit(abundance [duration == "annual" & abundance > 0])),
    seeded_bienn_s = length(na.omit(abundance [duration == "biennial" & abundance > 0])),
    seeded_perenn_s = length(na.omit(abundance [duration == "perennial" & abundance > 0]))
  ) %>%
  ungroup() %>%
  mutate (mix = "7_10")


all_mix_summary <- rbind(mix_1_2_summary, mix_3_6_summary, mix_7_10_summary)
write_csv(all_mix_summary, "G:/Postdoc Projects/Project 1 - IMLS/IMLS Analyses/IMLS_Analysis/all_mix_summary.csv")

#### data processing 4 - seed mix recruitment by quadrat ####

seeded_long <- select(seeded_spp, -family, -duration, -fg, -conservatism, -native, -c_score ) %>% 
  melt(id= c("taxon", "seeded_1_2", "seeded_3_6", "seeded_7_10", "sp_id")) %>%
  mutate(sample = variable, cover = value) %>%
  select(-variable, -value) %>%
  separate(col=sample, into= c('transect','quadrat','year','season'))
# remove the "X" in transect that showed up upon data import
seeded_long$transect <- substring(seeded_long$transect, 2)      


# calculate max cover for each species in each quadrat within year
seeded_year_max <- group_by(seeded_long, year,transect, quadrat, taxon ) %>% 
  summarise(max_cover = max(cover)) %>%
  ungroup()

seeded_year_max_merged <- merge(x=seeded_year_max, 
                                y=seeded_spp[,c("taxon", "seeded_1_2", "seeded_3_6", "seeded_7_10")], 
                                by="taxon", all.x = T) %>%
                          filter(!is.na(max_cover))

# remove non-seeded transects and count how many quadrats contain each sown species
# independently for 2018 and 2019 data
# recruitment <- filter(seeded_year_max_merged, transect != "02" , transect != "05" ,
#                transect != "06" , transect != "07", transect != "09") 

seeded_year_max_1_2 <- select(seeded_year_max_merged, -seeded_3_6, -seeded_7_10) %>%
  filter(transect == "01" | transect == "02", seeded_1_2 =="yes")


seeded_year_max_1_2_summary <- seeded_year_max_1_2 %>%
  group_by(taxon) %>% 
  summarise(
   sampled_quads_seeded_1_2_2017 = length(max_cover [transect == "01" & year=="2017" & max_cover >0]),
   sampled_quads_control_1_2_2017 = length(max_cover [transect == "02" & year=="2017" & max_cover >0]),
   sampled_quads_seeded_1_2_2018 = length(max_cover [transect == "01" & year=="2018" & max_cover >0]),
   sampled_quads_control_1_2_2018 = length(max_cover [transect == "02" & year=="2018" & max_cover >0]),
   sampled_quads_seeded_1_2_2019 = length(max_cover [transect == "01" & year=="2019" & max_cover >0]),
   sampled_quads_control_1_2_2019 = length(max_cover [transect == "02" & year=="2019" & max_cover >0]),
    ) %>%
  ungroup()


seeded_year_max_3_6 <- select(seeded_year_max_merged, -seeded_1_2, -seeded_7_10) %>%
  filter(transect == "03" | transect == "04" | transect == "05" | transect == "06", 
         seeded_3_6 == "yes")

seeded_year_max_3_6_summary <- seeded_year_max_3_6 %>%
  group_by(taxon) %>% 
  summarise(
    sampled_quads_seeded_3_6_2017 = length(max_cover [transect == "03" & year=="2017" & max_cover >0 
                                                      | transect == "04" & year=="2017" & max_cover >0]),
    sampled_quads_control_3_6_2017 = length(max_cover [transect == "05" & year=="2017" & max_cover >0 
                                                      | transect == "06" & year=="2017" & max_cover >0]),
    sampled_quads_seeded_3_6_2018 = length(max_cover [transect == "03" & year=="2018" & max_cover >0 
                                                      | transect == "04" & year=="2018" & max_cover >0]),
    sampled_quads_control_3_6_2018 = length(max_cover [transect == "05" & year=="2018" & max_cover >0 
                                                       | transect == "06" & year=="2018" & max_cover >0]),
    sampled_quads_seeded_3_6_2019 = length(max_cover [transect == "03" & year=="2019" & max_cover >0 
                                                      | transect == "04" & year=="2019" & max_cover >0]),
    sampled_quads_control_3_6_2019 = length(max_cover [transect == "05" & year=="2019" & max_cover >0 
                                                       | transect == "06" & year=="2019" & max_cover >0]),
  ) %>%
  ungroup()



seeded_year_max_7_10 <- select(seeded_year_max_merged, -seeded_1_2, -seeded_3_6) %>%
  filter(transect == "07" | transect == "08" | transect == "09" | transect == "10", seeded_7_10 =="yes")

seeded_year_max_7_10_summary <- seeded_year_max_7_10 %>%
  group_by(taxon) %>% 
  summarise(
    sampled_quads_seeded_7_10_2017 = length(max_cover [transect == "08" & year=="2017" & max_cover >0 
                                                      | transect == "10" & year=="2017" & max_cover >0]),
    sampled_quads_control_7_10_2017 = length(max_cover [transect == "07" & year=="2017" & max_cover >0 
                                                       | transect == "09" & year=="2017" & max_cover >0]),
    sampled_quads_seeded_7_10_2018 = length(max_cover [transect == "08" & year=="2018" & max_cover >0 
                                                      | transect == "10" & year=="2018" & max_cover >0]),
    sampled_quads_control_7_10_2018 = length(max_cover [transect == "07" & year=="2018" & max_cover >0 
                                                       | transect == "09" & year=="2018" & max_cover >0]),
    sampled_quads_seeded_7_10_2019 = length(max_cover [transect == "08" & year=="2019" & max_cover >0 
                                                      | transect == "10" & year=="2019" & max_cover >0]),
    sampled_quads_control_7_10_2019 = length(max_cover [transect == "07" & year=="2019" & max_cover >0 
                                                       | transect == "09" & year=="2019" & max_cover >0]),
  ) %>%
  ungroup()


# create column to sum hits across all samples for each seeded species
seeded_sampled <- seeded_spp %>%
  mutate(
    sum = select(.,X01.A.2017.Spring:X10.E.2019.Fall)%>% rowSums(na.rm=T), 
    present_ever = ifelse(sum>0,1,0)
  )
sum(seeded_sampled[,"sum"] !=0)
sum(seeded_sampled[,"present_ever"] !=0)
# 71 seeded species were detected at any point in time in the sampling


mix_performance <- seeded_spp[,1:11] %>%
  merge(y=seeded_year_max_1_2_summary, by="taxon",all.x = T) %>%
  merge(y=seeded_year_max_3_6_summary, by="taxon",all.x = T) %>%
  merge(y=seeded_year_max_7_10_summary, by="taxon",all.x = T) %>%
  merge(y=imls_seed_mixes, by="taxon",all.x = T) %>%
  merge(y=seeded_sampled[,c("taxon","present_ever")], by="taxon",all.x = T) %>%
  replace(is.na(.),0) %>%
  mutate(total_quads_2017 = sampled_quads_seeded_1_2_2017 + sampled_quads_seeded_3_6_2017 + sampled_quads_seeded_7_10_2017,
         total_quads_2018 = sampled_quads_seeded_1_2_2018 + sampled_quads_seeded_3_6_2018 + sampled_quads_seeded_7_10_2018,
         total_quads_2019 = sampled_quads_seeded_1_2_2019 + sampled_quads_seeded_3_6_2019 + sampled_quads_seeded_7_10_2019,
         total_quads_2017_control = sampled_quads_control_1_2_2017 + sampled_quads_control_3_6_2017 + sampled_quads_control_7_10_2017,
         total_quads_2018_control = sampled_quads_control_1_2_2018 + sampled_quads_control_3_6_2018 + sampled_quads_control_7_10_2018,
         total_quads_2019_control = sampled_quads_control_1_2_2019 + sampled_quads_control_3_6_2019 + sampled_quads_control_7_10_2019,
         present_17_control = ifelse(total_quads_2017_control >0,1,0),
         present_18_control = ifelse(total_quads_2018_control >0,1,0),
         present_19_control = ifelse(total_quads_2019_control >0,1,0),
         recruit_18 = total_quads_2018 / quadrats_seeded,
         recruit_19 = total_quads_2019 / quadrats_seeded,
         present_17_bin = ifelse(total_quads_2017 >0,1,0),
         present_18_bin = ifelse(total_quads_2018 >0,1,0),
         present_19_bin = ifelse(total_quads_2019 >0,1,0),
         present_18_or_19_bin = as.numeric(pmax(present_18_bin, present_19_bin)),
         present_17_18_19 = as.numeric(pmax(present_17_bin, present_18_bin, present_19_bin)),
         present_17_18_19_control = as.numeric(pmax(present_17_control, present_18_control, present_19_control)),
         control_transects_only = ifelse(present_ever > present_17_18_19,1,0),
         sampled_ever_within_unit = as.numeric(pmax(present_17_bin, present_18_bin, present_19_bin, 
                                        present_17_control, present_18_control, present_19_control)),
         mix_abundance =  t_1_2_oz + t_3_6_oz + t_7_10_oz,
         l_mix_abundance = log(mix_abundance),
         seeds_added = mix_abundance * seeds_oz
  ) 
  

# are recruitment rates of speceis different based on their source?
source_effect <- group_by(mix_performance, source)%>%
  summarise(count = length(source),
            recruit = length(source[recruit_19>0])
            )

mix_performance_fresh_seed <- mix_performance %>% filter(source != "stored")

# colSums(mix_performance !=0)
# colSums(mix_performance_fresh_seed !=0)


# species that recruited into transects where they were not sown
surprising_recruits <- mix_performance %>% 
  filter(control_transects_only>0) %>% 
  merge(y=imls_full_raw [,c("taxon", "latin")], by="taxon", all.x=T) %>%
  relocate(latin, .after=taxon)

write_csv(mix_performance, "G:/Postdoc Projects/Project 1 - IMLS/IMLS Analyses/IMLS_Analysis/mix_performance.csv")
write_csv(surprising_recruits, "G:/Postdoc Projects/Project 1 - IMLS/IMLS Analyses/IMLS_Analysis/surprising_recruits.csv")

# 56 seeded species were observed in seeded quadrats in 2019
# 58 seeded species were observed at any point in 2018 or 2019, meaning two recruited in 2018 but didn't persist to 2019
# 58 seeded species were observed at any point in '17-'19 meaning no species were detected prior to seeding that fell out
# 62 seeded species were observed at any point within their sown management unit
# since only 58 species were detected in seeded transects, 
# that means 13 unique seeded species ONLY recruited into non-seeded transects

#### Analysis step 1 - Two-way repeated measures ANOVA ####
# two-way repeated measures ANVOA tests whether a response variable changes
# over time, and/or between groups

# examine summary statistics

# colnames(all_transects_summary)
response_vars <- c("total_s", "native_s", "mean_quad_s","mean_nat_quad_s",
                   "prop_native", "total", "natives", "prop_seeded", "simp",
                   "simp_native", "shannon", "shannon_native", "j", "j_native")
                  
response_stats <- all_transects_summary %>%
  group_by(treat, year) %>%
  get_summary_stats(response_vars, type= "mean_se") %>%
  ungroup() %>%
  mutate(treat = replace(treat, treat == "control","Control")) %>%
  mutate(treat = replace(treat, treat == "seeded","Seeded")) %>%
  mutate(Treatment = treat)

response_stats_wide<- pivot_wider(response_stats, names_from = variable, values_from = c(mean,se))

# basic visualization
# response_bxp <- ggboxplot(
#  all_transects_summary, x = "year", y = response_vars, color = "treat", add = "point")
# response_bxp

# individual plots for each response variable

nat_quad_rich_plot <- ggplot(response_stats_wide, aes(x=year, y=mean_mean_nat_quad_s, color=Treatment)) + 
  geom_line(aes(group=Treatment)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=mean_mean_nat_quad_s-se_mean_nat_quad_s, ymax=mean_mean_nat_quad_s+se_mean_nat_quad_s), width=.3,
                position=position_dodge(0.05)) +
  labs(x="", y = "Richness (Quadrat)") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual(values=c('#999999','chartreuse3'))
print(nat_quad_rich_plot)

nat_trans_rich_plot <- ggplot(response_stats_wide, aes(x=year, y=mean_native_s, color=Treatment)) + 
  geom_line(aes(group=Treatment)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=mean_native_s-se_native_s, ymax=mean_native_s+se_native_s), width=.3,
                position=position_dodge(0.05)) +
  labs(x="", y = "Richness (Transect)") +
  theme_classic() +
  theme(legend.position = "right") +
  scale_color_manual(values=c('#999999','chartreuse3'))
print(nat_trans_rich_plot)

ggsave(file="transect_richness_graphical_abstract.png", nat_trans_rich_plot, units="in", width=4,height=3, scale=1, dpi=300)

nat_diversity_plot <- ggplot(response_stats_wide, aes(x=year, y=mean_simp_native, color=Treatment)) + 
  geom_line(aes(group=Treatment)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=mean_simp_native-se_simp_native, ymax=mean_simp_native+se_simp_native), width=.3,
                position=position_dodge(0.05)) +
  labs(x="", y = "Simpson's Diversity") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual(values=c('#999999','chartreuse3'))
print(nat_diversity_plot)

nat_shannon_diversity_plot <- ggplot(response_stats_wide, aes(x=year, y=mean_shannon_native, color=Treatment)) + 
  geom_line(aes(group=Treatment)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=mean_shannon_native-se_shannon_native, ymax=mean_shannon_native+se_shannon_native), width=.3,
                position=position_dodge(0.05)) +
  labs(x="", y = "Shannon's Diversity") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual(values=c('#999999','chartreuse3'))
print(nat_shannon_diversity_plot)

nat_evenness_plot <- ggplot(response_stats_wide, aes(x=year, y=mean_j_native, color=Treatment)) + 
  geom_line(aes(group=Treatment)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=mean_j_native-se_j_native, ymax=mean_j_native+se_j_native), width=.3,
                position=position_dodge(0.05)) +
  labs(x="", y = "Native Pielou's Evenness") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual(values=c('#999999','chartreuse3'))
print(nat_evenness_plot)

prop_native_plot <- ggplot(response_stats_wide, aes(x=year, y=mean_prop_native, color=Treatment)) + 
  geom_line(aes(group=Treatment)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=mean_prop_native-se_prop_native, ymax=mean_prop_native+se_prop_native), width=.3,
                position=position_dodge(0.05)) +
  scale_y_continuous(limits = c(0,1)) +
  labs(x="", y = "Proportion Native") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual(values=c('#999999','chartreuse3'))
print(prop_native_plot)

prop_seeded_plot <- ggplot(response_stats_wide, aes(x=year, y=mean_prop_seeded, color=Treatment)) + 
  geom_line(aes(group=Treatment)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=mean_prop_seeded-se_prop_seeded, ymax=mean_prop_seeded+se_prop_seeded), width=.3,
                position=position_dodge(0.05)) +
  scale_y_continuous(limits = c(0,1)) +
  labs(x="", y = "Proportion Seeded") +
  theme_classic() +
  theme(legend.position="right") +
  scale_color_manual(values=c('#999999','chartreuse3'))
print(prop_seeded_plot)

nat_cover_plot <- ggplot(response_stats_wide, aes(x=year, y=mean_natives, color=Treatment)) + 
  geom_line(aes(group=Treatment)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=mean_natives-se_natives, ymax=mean_natives+se_natives), width=.3,
                position=position_dodge(0.05)) +
  labs(x="", y = "Native Cover") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual(values=c('#999999','chartreuse3'))
print(nat_cover_plot)

figure_1 <- ggarrange(nat_quad_rich_plot, nat_trans_rich_plot, 
          nat_diversity_plot, prop_native_plot, prop_seeded_plot, nat_cover_plot, 
          labels = c("A", "B", "C", "D", "E", "F"),
          vjust = 1,
          common.legend = T,
          legend = "bottom",
          ncol = 2, nrow = 3)
figure_1
ggsave(file="figure_1.png", figure_1, units="in", width=6,height=6, scale=1.2, dpi=300)


# test for assumption of no outliers 

    # since we have so few data points within each treat*year (n=5) there is one outlier for total species richness, however examining the 
    # box plots, they appear to be reasonable and not erroneously high or low
all_transects_summary %>%
  group_by(treat, year) %>%
  identify_outliers(total_s) %>%
  ungroup()

    # no outliers for native species richness
all_transects_summary %>%
  group_by(treat, year) %>%
  identify_outliers(native_s) %>%
  ungroup()

    # one in average quadrat richness, but seems within reason
all_transects_summary %>%
  group_by(treat, year) %>%
  identify_outliers(mean_quad_s) %>%
  ungroup()

# no outliers
all_transects_summary %>%
  group_by(treat, year) %>%
  identify_outliers(simp_native) %>%
  ungroup()

# 4 outliers for proportion of cover from native species
all_transects_summary %>%
  group_by(treat, year) %>%
  identify_outliers(prop_native) %>%
  ungroup()

# 3 outliers for summed total cover; all seem reasonable
all_transects_summary %>%
  group_by(treat, year) %>%
  identify_outliers(total) %>%
  ungroup()

# one (pre-treatment) outlier for abundance of sown species. not a problem
all_transects_summary %>%
  group_by(treat, year) %>%
  identify_outliers(prop_seeded) %>%
  ungroup()

# test for assumption of normality

    # total species richness was normally distributed at each time point (p > 0.05)
all_transects_summary %>%
  group_by(treat, year) %>%
  shapiro_test(total_s) %>%
  ungroup()

# all response vars are normally distributed in each treatment except for 3 minor violations
# parametric tests are generally sufficiently robust to minor deviations in assumptions of normality
shapiro <- all_transects_summary %>%
  group_by(treat, year) %>%
  shapiro_test(response_vars) %>%
  ungroup() %>% filter(p<0.05)

# quantile-quantile plots show the correlation between data and the normal distribution
# there are no obvious significant deviations from normality
ggqqplot(all_transects_summary, "total_s", ggtheme = theme_bw()) +
  facet_grid(year ~ treat, labeller = "label_both")
ggqqplot(all_transects_summary, "native_s", ggtheme = theme_bw()) +
  facet_grid(year ~ treat, labeller = "label_both")
ggqqplot(all_transects_summary, "mean_quad_s", ggtheme = theme_bw()) +
  facet_grid(year ~ treat, labeller = "label_both")
ggqqplot(all_transects_summary, "simp_native", ggtheme = theme_bw()) +
  facet_grid(year ~ treat, labeller = "label_both")
ggqqplot(all_transects_summary, "prop_native", ggtheme = theme_bw()) +
  facet_grid(year ~ treat, labeller = "label_both")
ggqqplot(all_transects_summary, "total", ggtheme = theme_bw()) +
  facet_grid(year ~ treat, labeller = "label_both")
ggqqplot(all_transects_summary, "prop_seeded", ggtheme = theme_bw()) +
  facet_grid(year ~ treat, labeller = "label_both")

# test for homogeneity of variances
# significant p-value (<0.05) indicates heteroscedasticity

bartlett.test(native_s ~ treat, data= all_transects_summary)
bartlett.test(simp_native ~ treat, data= all_transects_summary)
bartlett.test(prop_native ~ treat, data= all_transects_summary)
bartlett.test(mean_nat_quad_s ~ treat, data= all_transects_summary)
bartlett.test(prop_seeded ~ treat, data= all_transects_summary) # variance higher in sown transects, but this makes sense
bartlett.test(natives ~ treat, data= all_transects_summary)


# compute repeated measures ANOVA (here as a standard two-way Anova)
# this method does not account for repeated measures co-variance structure
# example with total transect-level species richness
total_s.aov <- anova_test(
  data = all_transects_summary, 
  dv = total_s, 
  wid = transect,
  between = treat,
  within =year)
get_anova_table(total_s.aov)

# We can examine whether our data fit a variance-covariance structure that has
# compound symmetry, unstructured, auto-regressive, or
# auto-regressive with heterogeneous variances
# lets see which structure is best for our total richness variable

total_s_rep_format <- groupedData(total_s ~ as.numeric(treat_bin) * as.numeric(year) | transect, data = all_transects_summary )

# compound symmetry 
total_s.cs <- gls(total_s ~ treat_bin * year, data = total_s_rep_format,
                  corr = corCompSymm(form = ~ 1 | transect)) 
summary(total_s.cs)

# unstructured 
total_s.un <- gls(total_s ~ treat_bin * year, data = total_s_rep_format,
                  corr = corCompSymm(form = ~ 1 | transect),
                  weight = varIdent(form = ~ 1 | year))
summary(total_s.un)

# auto-regressive - the most intuitive option for rep measures data 
total_s.ar1 <- gls(total_s ~ treat_bin * year, data = total_s_rep_format,
                  corr = corAR1(form= ~ 1 | transect))
summary(total_s.ar1)

# auto-regressive with heterogeneous variances
total_s.arh1 <- gls(total_s ~ treat_bin * year, data = total_s_rep_format,
                  corr = corAR1(form= ~ 1 | transect),
                  weight = varIdent(form = ~ 1 | year))
summary(total_s.arh1)

# model comparison
anova(total_s.cs, total_s.un)
anova(total_s.cs, total_s.ar1)
anova(total_s.cs, total_s.arh1)
anova(total_s.ar1, total_s.un)
anova(total_s.un, total_s.ar1)

# we will use the AR1 structure

# complete analysis for total species richness
anova(total_s.ar1) # significant treat*year effect (p= 0.004)
summary(total_s.ar1)
emmeans_total_s <- emmeans(total_s.ar1, "treat_bin", by = "year")
pairs(emmeans_total_s)


# create models for the other response variables

# native species richness
native_s_rep_format <- groupedData(native_s ~ as.numeric(treat_bin) * as.numeric(year) | transect, data = all_transects_summary )
native_s.ar1 <- gls(native_s ~ treat_bin * year, data = native_s_rep_format,
                   corr = corAR1(form= ~ 1 | transect))
summary(native_s.ar1)
anova(native_s.ar1)
emmeans_native_s <- emmeans(native_s.ar1, "treat_bin", by = "year")
pairs(emmeans_native_s)


# mean quadrat richness
mean_quad_s_rep_format <- groupedData(mean_quad_s ~ as.numeric(treat_bin) * as.numeric(year) | transect, data = all_transects_summary )
mean_quad_s.ar1 <- gls(mean_quad_s ~ treat_bin * year, data = mean_quad_s_rep_format,
                    corr = corAR1(form= ~ 1 | transect))
summary(mean_quad_s.ar1)
anova(mean_quad_s.ar1)
emmeans_mean_quad_s <- emmeans(mean_quad_s.ar1, "treat_bin", by = "year")
pairs(emmeans_mean_quad_s)


# mean native quadrat richness
mean_nat_quad_s_rep_format <- groupedData(mean_nat_quad_s ~ as.numeric(treat_bin) * as.numeric(year) | transect, data = all_transects_summary )
mean_nat_quad_s.ar1 <- gls(mean_nat_quad_s ~ treat_bin * year, data = mean_nat_quad_s_rep_format,
                       corr = corAR1(form= ~ 1 | transect))
summary(mean_nat_quad_s.ar1)
anova(mean_nat_quad_s.ar1)
emmeans_mean_nat_quad_s <- emmeans(mean_nat_quad_s.ar1, "treat_bin", by = "year")
pairs(emmeans_mean_nat_quad_s)


# native diversity (Simpson's effective species number)
diversity_rep_format <- groupedData(simp_native ~ as.numeric(treat_bin) * as.numeric(year) | transect, data = all_transects_summary )
diversity.ar1 <- gls(simp_native ~ treat_bin * year, data = diversity_rep_format,
                    corr = corAR1(form= ~ 1 | transect))
summary(diversity.ar1)
anova(diversity.ar1) # no significant time*treat
emmeans_diversity <- emmeans(diversity.ar1, "year")
pairs(emmeans_diversity) 
# diversity is higher in 2018 and 2019 than in 2017, but seeding is not significant


# proportion cover from native species
prop_native_rep_format <- groupedData(prop_native ~ as.numeric(treat_bin) * as.numeric(year) | transect, data = all_transects_summary )
prop_native.ar1 <- gls(prop_native ~ treat_bin * year, data = prop_native_rep_format,
                    corr = corAR1(form= ~ 1 | transect))
summary(prop_native.ar1)
anova(prop_native.ar1)
emmeans_prop_native <- emmeans(prop_native.ar1, "treat_bin", by =  "year")
pairs(emmeans_prop_native)

emmeans_prop_native_year <- emmeans(prop_native.ar1,  "year")
pairs(emmeans_prop_native_year)
# proportional abundance of native species is higher in 2018 and 2019 than in 2017, but seeding is not significant



# proportional cover of seeded species
prop_seeded_rep_format <- groupedData(prop_seeded ~ as.numeric(treat_bin) * as.numeric(year) | transect, data = all_transects_summary )
prop_seeded.ar1 <- gls(prop_seeded ~ treat_bin * year, data = prop_seeded_rep_format,
                       corr = corAR1(form= ~ 1 | transect))
summary(prop_seeded.ar1)
anova(prop_seeded.ar1)
emmeans_prop_seeded <- emmeans(prop_seeded.ar1, "treat_bin", by = "year")
pairs(emmeans_prop_seeded)




# total native cover 
total_rep_format <- groupedData(natives ~ as.numeric(treat_bin) * as.numeric(year) | transect, data = all_transects_summary )
total.ar1 <- gls(natives ~ treat_bin * year, data = total_rep_format,
                    corr = corAR1(form= ~ 1 | transect))
summary(total.ar1)
anova(total.ar1)
emmeans_total <- emmeans(total.ar1, "year")
pairs(emmeans_total)
# total cover is higher in 2018 and 2019 than in 2017, but seeding is not significant




#### Analysis step 2 - PerMANOVA on 2019 Composition data ####
  
wf.orig <- filter(year_max_wide, year == "2019") %>%
    select( -num_quadrat) %>%
    mutate(Seeded = c(rep("Yes",5),rep("No",5),rep("Yes",5),rep("Yes",5),rep("No",5),rep("No",5),
                                                          rep("No",5),rep("Yes",5),rep("No",4),rep("Yes",5))) %>%
    replace(is.na(.),0) %>% #replace NA with zeros
    filter(quadrat != "D") #Remove all quadrat D to allow for restricted permutations. 

  
  #Make a separate table for factors
  env = wf.orig[,1:3]
  env$transect = as.factor(env$transect)
  env$Seeded = wf.orig$Seeded
  wf = wf.orig[,4:214]  # was 260 with unknowns included
  
  dimcheckMDS(wf, distance='bray', k=6, trymax = 100, autotransform=F) #appears 4 axes is best (which is somewhat unsatisfactory)
  
  ord = metaMDS(wf, distance="bray", k=4,try = 500, trymax=1000,maxit=1000, autotransform=F,engine="monoMDS") 
  stressplot(ord) #linear fit = 0.81, that's OK...
  
site_scores = as.data.frame(scores(ord, display = "sites")) %>%
    cbind(env)

ord_spp_fit <- envfit(ord, wf, permutations = 999)

spp_scores = as.data.frame(scores(ord_spp_fit, display = "vectors")) %>%
  tibble::rownames_to_column("taxon") %>%
  cbind(pval = ord_spp_fit$vectors$pvals) 

sig_spp_scores <- subset(spp_scores, pval<=0.01) #34 significant species

# find the most abundant species in 2019
max_2019 <- year_max %>%
  filter(year == 2019) %>%
  group_by(taxon) %>%
  summarise(sum_max_cover = sum(max_cover)) %>%
  ungroup() 
top_20_2019 <- top_n(max_2019, 20)

abundant_spp <- merge(spp_scores,top_20_2019,by="taxon", all.y=T) %>%
  merge(y=sampled_spp[c("taxon","fg","conservatism")],by="taxon",all.x = T) %>%
  mutate_if(is.factor, fct_explicit_na, na_level = 'exotic')%>% # need to replace NA in conservatism col
  mutate(taxon = toupper(taxon)) %>%
  dplyr::mutate(
    conservatism = factor(stringr::str_to_title(conservatism)),
    fg = factor(stringr::str_to_title(fg)))

all_spp <- merge(spp_scores,max_2019,by="taxon", all.y=T) %>%
  merge(y=sampled_spp[c("taxon","fg","conservatism")],by="taxon",all.x = T) %>%
  mutate_if(is.factor, fct_explicit_na, na_level = 'exotic') %>%
  mutate(taxon = toupper(taxon))
all_spp$fg <-recode_factor(all_spp$fg, `exotic` = "non-native",`fern`="forb", `vine`="woody", `shrub`="woody")

# which species are significant and among the 20 most abundant?
abundant_sig_spp <- merge(sig_spp_scores,top_20_2019,by="taxon") %>%
  merge(y=sampled_spp[c("taxon","fg","conservatism")],by="taxon",all.x = T) %>%
  mutate_if(is.factor, fct_explicit_na, na_level = 'exotic') # need to replace NA in conservatism col

# old ordination
  ggplot(site_scores, aes(x=NMDS1,y=NMDS2)) + theme_classic() +
    geom_point(aes(color=Seeded),size=6) + 
    theme(axis.title.y=element_text(size=20),axis.title.x=element_text(size=20),axis.text.y=element_text(size=15),axis.text.x=element_text(size=15),
          legend.text=element_text(size=15),legend.title=element_text(size=20))+
    geom_text(aes(x=-0.75,y=1.0,hjustvar=0,vjustvar=1,label="stress = 0.1274"))

# new ordination plot 
transect_ordination <- site_scores %>%
group_by(transect) %>%
  summarise(
    NMDS1 = mean(NMDS1),
    NMDS2 = mean(NMDS2)
  ) %>%
ungroup()

trans_ord <-  site_scores %>%
  group_by(transect,Seeded) %>%
  get_summary_stats(c("NMDS1","NMDS2"), type= "mean_se") %>%
  ungroup()

trans_ord_wide <- pivot_wider(trans_ord,names_from = variable, values_from=c(mean,se)) %>%
  mutate(Treatment = ifelse(Seeded=="Yes", "Seeded" ,"Control"),
         unit =c("North","North","Central","Central","Central","Central","Riparian","Riparian","Riparian","Riparian")
         )

trans_ord_plot <- ggplot(trans_ord_wide, aes(x=mean_NMDS1, y=mean_NMDS2, color=Treatment, shape=unit)) + 
  geom_point(size=3) +
  coord_cartesian(xlim = c(-.9, .9), ylim =c(-.9,.9)) +
  geom_errorbar(aes(ymin=mean_NMDS2-se_NMDS2, ymax=mean_NMDS2+se_NMDS2)) +
  geom_errorbarh(aes(xmin=mean_NMDS1-se_NMDS1, xmax=mean_NMDS1+se_NMDS1)) +
  labs(x="NMDS1", y = "NMDS2", shape="Subunit") +
  theme_classic() +
  theme(legend.position = "right") +
  scale_color_manual(values=c('#999999','chartreuse3'))
trans_ord_plot

ggsave(file="ordination.png", trans_ord_plot, units="in", width=6,height=6, scale=1, dpi=300)
# ggsave(file="ordination_alt.png", trans_ord_plot, units="in", width=3,height=3, scale=1.2, dpi=300)

# basic plot with the species with p<0.01 significance
bi_plot <- ggplot(trans_ord_wide, aes(x=mean_NMDS1, y=mean_NMDS2)) +
  labs(x="NMDS1", y = "NMDS2") +
  theme_classic() +
  coord_cartesian(xlim = c(-.9, .9), ylim =c(-.9,.9)) +
  theme(legend.position = "right")+
  geom_segment(data = sig_spp_scores, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = sig_spp_scores, aes(x=NMDS1, y=NMDS2, label = taxon), cex = 3, direction = "both", segment.size = 0.25)
bi_plot

# 13 speices that are both significant and among the 20 most abundant
bi_plot2 <- ggplot(abundant_sig_spp, aes(x=mean_NMDS1, y=mean_NMDS2, color =fg, linetype=conservatism)) +
  labs(x="NMDS1", y = "NMDS2") +
  theme_classic() +
  coord_cartesian(xlim = c(-.9, .9), ylim =c(-.9,.9)) +
  theme(legend.position = "right")+
  geom_segment(data = abundant_sig_spp, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), lwd=0.3) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = abundant_sig_spp, aes(x=NMDS1, y=NMDS2, label = toupper(taxon)), cex = 3, direction = "both", segment.size = 0.25) +
  scale_linetype_manual(values=c("solid", "longdash","dashed", "dotted"))+
  labs(color = "funtional group")
bi_plot2

# 20 most abundant species in 2019 plus other species as circles
bi_plot3 <- ggplot(trans_ord_wide, aes(x=mean_NMDS1, y=mean_NMDS2, linetype=conservatism)) +
  labs(x="NMDS1", y = "NMDS2") +
  theme_classic() +
  coord_cartesian(xlim = c(-.9, .9), ylim =c(-.9,.9)) +
  theme(legend.position = "right")+
  geom_point(data=all_spp, aes(x=NMDS1,y=NMDS2, color =fg,),alpha =0.5,cex=2.5) +
  # geom_segment(data = abundant_spp, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), lwd=0.3) + #add vector arrows of significant species
  # scale_linetype_manual(values=c("solid", "longdash","dashed", "dotted")) +
  ggrepel::geom_text_repel(data = abundant_spp, aes(x=NMDS1, y=NMDS2, label = taxon), cex =3, direction = "both", segment.size = 0.3) +
  labs(color = "Funtional group")
bi_plot3


bi_plot4 <- ggplot(abundant_spp, aes(x=mean_NMDS1, y=mean_NMDS2, color=fg, linetype=conservatism)) +
  labs(x="NMDS1", y = "NMDS2") +
  theme_classic() +
  coord_cartesian(ylim =c(-.9,.9)) +
  theme(legend.position = "right")+
  geom_point(data=abundant_spp, aes(x=NMDS1,y=NMDS2, shape =conservatism,),alpha =0.5,cex=2.5) +
  #geom_segment(data = abundant_spp, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), lwd=0.3) + #add vector arrows of significant species
  #scale_linetype_manual(values=c("solid", "longdash","dashed", "dotted")) +
  ggrepel::geom_text_repel(
    data = subset(abundant_spp,NMDS1>0), 
    nudge_x = 1 - subset(abundant_spp,NMDS1>0)$NMDS1,
    aes(x=NMDS1, y=NMDS2, label = taxon),
    cex =3,
    min.segment.length = 0, 
    direction = "y",
    hjust=1) +
  ggrepel::geom_text_repel(
    data = subset(abundant_spp,NMDS1<0), 
    nudge_x = -1.1 - subset(abundant_spp,NMDS1<0)$NMDS1,
    aes(x=NMDS1, y=NMDS2, label = taxon),
    cex =3,
    min.segment.length = 0, 
    direction = "y",
    hjust=0) +
  labs(color = "Functional group", shape="Conservatism") +
  scale_x_continuous(
    breaks = c(-1, -0.5, 0, 0.5, 1),
    limits = c(-1.2, 1.2)) + 
  guides(color = guide_legend(override.aes = aes(label = "")))
bi_plot4


bi_plot5 <- ggplot(abundant_spp, aes(x=mean_NMDS1, y=mean_NMDS2, color=fg, shape= factor(conservatism,levels = c("Conservative", "Matrix", "Ruderal","Exotic")))) +
  labs(x="NMDS1", y = "NMDS2") +
  theme_classic() +
  coord_cartesian(xlim = c(-.9, .9), ylim =c(-.9,.9)) +
  theme(legend.position = "right")+
  scale_shape_manual(values = c(16,17,15,8)) +
  scale_color_manual(values = cb_palette_4) +
  geom_point(data=abundant_spp, aes(x=NMDS1,y=NMDS2),cex=2.5) +
  #geom_segment(data = abundant_spp, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), lwd=0.3) + #add vector arrows of significant species
  #scale_linetype_manual(values=c("solid", "longdash","dashed", "dotted")) +
  ggrepel::geom_text_repel(data = abundant_spp, 
                           aes(x=NMDS1, y=NMDS2, 
                               label = toupper(taxon)), 
                           cex = 3, direction = "both") +
  labs(color = "Functional group", shape="Conservatism") +
  guides(color = guide_legend(override.aes = aes(label = "")))
bi_plot5




figure_3 <- ggarrange(trans_ord_plot, bi_plot5, 
                      labels = c("A", "B"))
figure_3
ggsave(file="figure_3.png", figure_3, units="in", width=10,height=4, scale=1, dpi=300)


  
#Set controls for permutations, this applies to envfit, adonis2 (PERMANOVA)
#First, remove all quadrat D to allow for restricted permutations. 
  wf_alt = wf[wf.orig$quadrat != "D",]
  env_alt = env[env$quadrat != "D",]
  
  ctrl = how(nperm=9999,within= Within(type="none",mirror=F), plots = Plots(strata=env_alt$transect,type="free"))
  #To demonstrate what this does, 
  shuffle(wf_alt,ctrl)
  shuffle(wf_alt,ctrl)
  shuffle(wf_alt,ctrl)
  #look at the output. Transects move around as cohesive blocks, which remain in series (no within transect
  #shuffling) because they are arranged in a spatially explicit line - i.e. 1 closest to 2, 2 to 1 and 3...
  
  fit = adonis2(wf_alt ~ Seeded, env_alt, permutations=ctrl)
  print(fit)
  
#### Analysis step 3 - 2019 Conservatism & FQI ####

# bartlett test for homogenetity of variances
# these data violate the assumption of homogeneity
bartlett.test(prop_ruderal ~ treat, data=summary_2019 )
bartlett.test(prop_matrix ~ treat, data=summary_2019 )
bartlett.test(prop_conservative ~ treat, data=summary_2019 )

# levene test for homogenetity  of variances 
# (less sensitive than bartlett to deviation from normality)
leveneTest(prop_ruderal ~ treat, data=summary_2019 ) 
leveneTest(prop_matrix ~ treat, data=summary_2019 )
leveneTest(prop_conservative ~ treat, data=summary_2019 )


# Mann-whitney U-tests
m1<-wilcox.test(prop_ruderal ~ treat, data=summary_2019, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)
print(m1)
m2<-wilcox.test(prop_matrix ~ treat, data=summary_2019, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)
print(m2)
m3<-wilcox.test(prop_conservative ~ treat, data=summary_2019, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)
print(m3)

#Hodges Lehmann Estimator
m1$estimate
m2$estimate
m3$estimate


### analysis for forbs ###
leveneTest(prop_ruderal ~ treat, data=summary_forb_2019 ) 
leveneTest(prop_matrix ~ treat, data=summary_forb_2019 )
leveneTest(prop_conservative ~ treat, data=summary_forb_2019 )
           
wilcox.test(prop_ruderal ~ treat, data=summary_forb_2019, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)
wilcox.test(prop_matrix ~ treat, data=summary_forb_2019, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)
wilcox.test(prop_conservative ~ treat, data=summary_forb_2019, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)



# MANOVA approach
# NOT used in manuscript
# assumptions: no outliers, normality, homogeneity of variances

# conservatism_man <- manova(cbind(conservative, matrix, ruderal) ~ treat, data = summary_2019 )
# summary(conservatism_man)
# summary.aov(conservatism_man)

# test for multivariate normality is significant indicating deviation from assumptions
# mshapiro.test(t(summary_2019[,c("ruderal","matrix","conservative")]))

# 
# mshapiro.test(t(summary_2019[,c("prop_ruderal","prop_matrix","prop_conservative")]))

# prop_conservatism_man <- manova(cbind(prop_conservative, prop_matrix, prop_ruderal) ~ treat, data = summary_2019 )
# summary(prop_conservatism_man)
# summary.aov(prop_conservatism_man) # compare treatments
# summary.lm(prop_conservatism_man) # extract model R-squared?



# plot proportional abundances in seeded and unseeded transects in 2019

conservatism_stats <- summary_2019 %>%
  group_by(treat) %>%
  get_summary_stats(c("prop_ruderal","prop_matrix","prop_conservative"), type= "mean_se") %>%
  ungroup() %>%
  mutate(treat = replace(treat, treat == "control","Control")) %>%
  mutate(treat = replace(treat, treat == "seeded","Seeded")) %>%
  mutate(variable = replace(variable, variable == "prop_conservative","Conservative")) %>%
  mutate(variable = replace(variable, variable == "prop_matrix","Matrix")) %>%
  mutate(variable = replace(variable, variable == "prop_ruderal","Ruderal")) %>%
  mutate(Treatment = treat)

conservatism_plot <- ggplot(conservatism_stats, aes(x=variable, y=mean, fill=Treatment)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  scale_y_continuous(limits = c(0,1), expand = expansion(mult = c(0, .1))) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,
                position=position_dodge(.9)) +
  labs(x="Conservatism Group", y = "Proportional Cover") +
  theme_classic() +
  theme(
    axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")), # add space between axis labels and tick labels
    axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm")),
    legend.position = "c(0.15,0.85)"
  ) +
  scale_fill_manual(values=c('#999999','chartreuse3')) +
  geom_signif(y_position = c(0.8, 0.8), xmin = c(1.75,2.75), xmax = c(2.25,3.25),
              annotation = c("**","**"), tip_length = 0.03                    
  )
conservatism_plot
# ggsave(file="conservatism_groups.png", conservatism_plot, units="in", width=6,height=6, scale=.7, dpi=300)
# ggsave(file="conservatism_groups_alt.png", conservatism_plot, units="in", width=3,height=3, scale=.9, dpi=300)

# MANOVA on conservatism groups just for forbs

# test for multivariate normality is significant indicating deviation from assumptions
# mshapiro.test(t(summary_forb_2019[,c("prop_ruderal","prop_matrix","prop_conservative")]))

# prop_conservatism_forb_man <- manova(cbind(prop_conservative, prop_matrix, prop_ruderal) ~ treat, data = summary_forb_2019 )
# summary(prop_conservatism_forb_man, test='Pillai')
# summary(prop_conservatism_forb_man, test='Wilks')
# summary.aov(prop_conservatism_forb_man) # compare treatments
# summary.lm(prop_conservatism_forb_man) # extract model R-squared?

#### figure for just forbs

forb_conservatism_stats <- summary_forb_2019 %>% 
  group_by(treat) %>%
  get_summary_stats(c("prop_ruderal","prop_matrix","prop_conservative"), type= "mean_se") %>%
  ungroup() %>%
  mutate(treat = replace(treat, treat == "control","Control")) %>%
  mutate(treat = replace(treat, treat == "seeded","Seeded")) %>%
  mutate(variable = replace(variable, variable == "prop_conservative","Conservative")) %>%
  mutate(variable = replace(variable, variable == "prop_matrix","Matrix")) %>%
  mutate(variable = replace(variable, variable == "prop_ruderal","Ruderal")) %>%
  mutate(Treatment = treat)

forb_conservatism_plot <- ggplot(forb_conservatism_stats, aes(x=variable, y=mean, fill=Treatment)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  scale_y_continuous(limits = c(0,1), expand = expansion(mult = c(0, .1))) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,
                position=position_dodge(.9)) +
  labs(x="Conservatism Group", y = "Forb Proportional Cover") +
  theme_classic() +
  theme(
    axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")), # add space between axis labels and tick labels
    axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm")),
    legend.position = c(0.15,0.85)
  ) +
  scale_fill_manual(values=c('#999999','chartreuse3')) +
  geom_signif(y_position = c(0.6, 0.95), xmin = c(1.75,2.75), xmax = c(2.25,3.25),
              annotation = c("**","**"), tip_length = 0.03                    
  )
forb_conservatism_plot
# ggsave(file="forb_conservatism_groups.png", forb_conservatism_plot, units="in", width=6,height=6, scale=0.7, dpi=300)

# combine the two conservatism plots
conservatsim_plots <- ggarrange(conservatism_plot, forb_conservatism_plot, 
                      labels = c("A", "B"), common.legend = T)
conservatsim_plots
# ggsave(file="figure_2.png", conservatsim_plots, units="in", width=10,height=5, scale=1, dpi=300)




# t-tests on FQI between seeded and control transects in 2019

bartlett.test(fqi ~ treat, data=summary_2019) # test for violation homoscedasticity is not significant; good
t.test(fqi ~ treat, data=summary_2019, var.equal=TRUE, conf.level=0.95) # two sample t-test
histogram(~ fqi | treat, data=summary_2019, layout=c(2,1))
boxplot(fqi ~ treat, data = summary_2019, names=c("control","seeded"), ylab="fqi")

bartlett.test(weighted_fqi ~ treat, data=summary_2019) # test for violation homoscedasticity is not significant; good
t.test(weighted_fqi ~ treat, data=summary_2019, var.equal=TRUE, conf.level=0.95) # two sample t-test
histogram(~ weighted_fqi | treat, data=summary_2019, layout=c(2,1))
boxplot(weighted_fqi ~ treat, data = summary_2019, names=c("control","seeded"), ylab="weighted fqi")

bartlett.test(mean_c ~ treat, data=summary_2019) # test for violation homoscedasticity is not significant; good
t.test(mean_c ~ treat, data=summary_2019, var.equal=TRUE, conf.level=0.95) # two sample t-test
histogram(~ mean_c | treat, data=summary_2019, layout=c(2,1))
boxplot(mean_c ~ treat, data = summary_2019, names=c("control","seeded"), ylab="Mean C-Score")

bartlett.test(cwmc ~ treat, data=summary_2019) # test for violation homoscedasticity is not significant; good
t.test(cwmc ~ treat, data=summary_2019, var.equal=TRUE, conf.level=0.95) # two sample t-test
histogram(~ cwmc | treat, data=summary_2019, layout=c(2,1))
boxplot(cwmc ~ treat, data = summary_2019, names=c("control","seeded"), ylab="weighted Mean C-Score")



# FORB t-tests on FQI between seeded and control transects in 2019
bartlett.test(fqi ~ treat, data=summary_forb_2019) # test for violation homoscedasticity is not significant; good
t.test(fqi ~ treat, data=summary_forb_2019, var.equal=TRUE, conf.level=0.95) # two sample t-test
histogram(~ fqi | treat, data=summary_forb_2019, layout=c(2,1))
boxplot(fqi ~ treat, data = summary_forb_2019, names=c("control","seeded"), ylab="fqi")

bartlett.test(weighted_fqi ~ treat, data=summary_forb_2019) # test for violation homoscedasticity is not significant; good
t.test(weighted_fqi ~ treat, data=summary_forb_2019, var.equal=TRUE, conf.level=0.95) # two sample t-test
histogram(~ weighted_fqi | treat, data=summary_forb_2019, layout=c(2,1))
boxplot(weighted_fqi ~ treat, data = summary_forb_2019, names=c("control","seeded"), ylab="weighted fqi")

bartlett.test(mean_c ~ treat, data=summary_forb_2019) # test for violation homoscedasticity is not significant; good
t.test(mean_c ~ treat, data=summary_forb_2019, var.equal=TRUE, conf.level=0.95) # two sample t-test
histogram(~ mean_c | treat, data=summary_forb_2019, layout=c(2,1))
boxplot(mean_c ~ treat, data = summary_forb_2019, names=c("control","seeded"), ylab="Mean C-Score")

bartlett.test(cwmc ~ treat, data=summary_forb_2019) # test for violation homoscedasticity is not significant; good
t.test(cwmc ~ treat, data=summary_forb_2019, var.equal=TRUE, conf.level=0.95) # two sample t-test
histogram(~ cwmc | treat, data=summary_forb_2019, layout=c(2,1))
boxplot(cwmc ~ treat, data = summary_forb_2019, names=c("control","seeded"), ylab="Mean C-Score")




# FQI plots - not in manuscript
fqi_stats <- summary_2019 %>%
  group_by(treat) %>%
  get_summary_stats(c("fqi"), type= "mean_se") %>%
  ungroup() %>%
  mutate(treat = replace(treat, treat == "control","Control")) %>%
  mutate(treat = replace(treat, treat == "seeded","Seeded")) %>%
  mutate(Treatment = treat)

fqi_forb_stats <- summary_forb_2019 %>%
  group_by(treat) %>%
  get_summary_stats(c("fqi"), type= "mean_se") %>%
  ungroup() %>%
  mutate(treat = replace(treat, treat == "control","Control")) %>%
  mutate(treat = replace(treat, treat == "seeded","Seeded")) %>%
  mutate(Treatment = treat)

fqi_plot <- ggplot(fqi_stats, aes(x=treat, y=mean)) + 
  geom_point(size=2) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,
                position=position_dodge(0.05)) +
  labs(x="Treatment", y = "Floristic Quality Index") +
  theme_classic() 
print(fqi_plot)

mean_c_stats <- summary_2019 %>%
  group_by(treat) %>%
  get_summary_stats(c("mean_c"), type= "mean_se") %>%
  ungroup() %>%
  mutate(treat = replace(treat, treat == "control","Control")) %>%
  mutate(treat = replace(treat, treat == "seeded","Seeded")) %>%
  mutate(Treatment = treat)

mean_c_forb_stats <- summary_forb_2019 %>%
  group_by(treat) %>%
  get_summary_stats(c("mean_c"), type= "mean_se") %>%
  ungroup() %>%
  mutate(treat = replace(treat, treat == "control","Control")) %>%
  mutate(treat = replace(treat, treat == "seeded","Seeded")) %>%
  mutate(Treatment = treat)


#### Analysis step 4 - 2019 functional groups ####


# test for multivariate normality is significant indicating deviation from assumptions
mshapiro.test(t(summary_2019[,c("grass", "sedge", "forb_2", "legume", "vine", "shrub", "tree")]))
mshapiro.test(t(summary_2019[,c("prop_grass", "prop_sedge", "prop_forb_2", "prop_legume", "prop_vine", "prop_shrub", "prop_tree")]))


# MANOVA appraoch - not in manuscript
# fg_man <- manova(cbind(grass, sedge, forb_2, legume, vine, shrub, tree) ~ treat, data = summary_2019 )
# summary(fg_man) # MANOVA on raw cover is not significant
# summary.aov(fg_man)

# prop_fg_man <- manova(cbind(prop_grass, prop_sedge, prop_forb_2, prop_legume, prop_vine, prop_shrub, prop_tree) ~ treat, data = summary_2019 )
# summary(prop_fg_man, test = 'Pillai') # MANOVA on proportional cover is significant
# summary(prop_fg_man, test = 'Wilks')
# summary(prop_fg_man, test = 'Hotelling-Lawley')
# summary(prop_fg_man, test = 'Roy')
# summary.aov(prop_fg_man)

# library(effectsize)
# eta_squared(prop_fg_man) 


# functional group analyses with non-parametric approach
leveneTest(prop_forb ~ treat, data=summary_2019 )
leveneTest(prop_grass ~ treat, data=summary_2019 )
leveneTest(prop_legume ~ treat, data=summary_2019 )
leveneTest(prop_sedge ~ treat, data=summary_2019 )
leveneTest(prop_shrub ~ treat, data=summary_2019 )
leveneTest(prop_tree ~ treat, data=summary_2019 )

wilcox.test(prop_forb ~ treat, data=summary_2019, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)
wilcox.test(prop_grass ~ treat, data=summary_2019, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)
wilcox.test(prop_legume ~ treat, data=summary_2019, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)
wilcox.test(prop_sedge ~ treat, data=summary_2019, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)
wilcox.test(prop_shrub ~ treat, data=summary_2019, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)
wilcox.test(prop_tree ~ treat, data=summary_2019, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)
wilcox.test(prop_vine ~ treat, data=summary_2019, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)



# basic stacked bar graph showing differences between seeded and control transects in 2019
data_fg_2019 <- group_by(summary_2019, treat ) %>%
  summarise(
    count = n(),
    mean_grass = mean(grass, na.rm = TRUE),
    mean_sedge = mean(sedge, na.rm = TRUE),
    mean_forb_2 = mean(forb_2, na.rm = TRUE),
    mean_legume = mean(legume, na.rm = TRUE),
    mean_shrub = mean(shrub, na.rm = TRUE),
    mean_tree = mean(tree, na.rm = TRUE)
    ) %>%
  pivot_longer( cols = c("mean_grass", "mean_sedge", "mean_forb_2" ,
                         "mean_legume", "mean_shrub", "mean_tree"), 
                names_to = "fg",
                values_to = "cover")

ggplot(data_fg_2019, aes(x = treat, y = cover , fill = fg)) +
  geom_col(position = "fill")


# functional group figure

fg_vars <- c("prop_grass","prop_sedge","prop_forb_2","prop_legume","prop_vine","prop_shrub","prop_tree")
fg_stats <- summary_2019 %>%
  group_by(treat) %>%
  get_summary_stats(fg_vars, type= "mean_se") %>%
  ungroup() %>%
  mutate(treat = replace(treat, treat == "control","Control")) %>%
  mutate(treat = replace(treat, treat == "seeded","Seeded")) %>%
  mutate(variable = replace(variable, variable == "prop_grass","Grass")) %>%
  mutate(variable = replace(variable, variable == "prop_sedge","Sedge")) %>%
  mutate(variable = replace(variable, variable == "prop_forb_2","Forb")) %>%
  mutate(variable = replace(variable, variable == "prop_legume","Legume")) %>%
  mutate(variable = replace(variable, variable == "prop_vine","Vine")) %>%
  mutate(variable = replace(variable, variable == "prop_shrub","Shrub")) %>%
  mutate(variable = replace(variable, variable == "prop_tree","Tree")) %>%
  mutate(Treatment = treat)

fg_plot <- ggplot(fg_stats, aes(x=variable, y=mean, fill=Treatment)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  scale_y_continuous(limits = c(0,0.7),expand = expansion(mult = c(0, .1))) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,
                position=position_dodge(.9)) +
  labs( x="Functional Group", y = "Proportional Cover") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 30, vjust=0.8),
    axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm")), # add space between axis labels and tick labels
    axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm")),
    legend.position = c(0.8,0.8)
  ) +
  scale_fill_manual(values=c('#999999','chartreuse3')) +
  geom_signif(y_position = c(0.48, 0.15), xmin = c(1.75,3.75), xmax = c( 2.25,4.25),
                        annotation = c( "**","**"), tip_length = 0.03                    
                    )
fg_plot
# ggsave(file="functional_groups.png", fg_plot, units="in", width=6,height=6, scale=0.7, dpi=300)



#### Analysis step 5 - Seed Mix Performance ####


# question 1, is establishment rate predicted by mass of seeds introduced (a measurement of propagule pressure) ?

abundance.mod <- glm(recruit_19 ~ l_mix_abundance, weights = quadrats_seeded, data=mix_performance, family=quasibinomial(link="logit"))
summary(abundance.mod)
anova(abundance.mod, test = "F") # recruitment is strongly positively related to abundance (mass) of seed sown


# the beta-binomial distribution may fit the data better than a quasibinomial
# test.mod_1 <- glmmTMB(recruit_19 ~ l_mix_abundance , weights = quadrats_seeded, data=mix_performance, family='betabinomial')
# summary(test.mod_1)
# drop1(test.mod_1, test="Chisq")
# sr_1 <- simulateResiduals(test.mod_1)
# plot(sr_1)


mix_performance_new <- mix_performance %>%
  filter(seeds_added != 0) # remove "ranhis" for which we were unable to find an estimated seed mass

seeds.mod <- glm(recruit_19 ~ log(seeds_added), weights = quadrats_seeded, data=mix_performance_new, family=quasibinomial(link="logit"))
# seeds_beta.mod <- glmmTMB(recruit_19 ~ log(seeds_added), weights = quadrats_seeded, data=mix_performance_new, family='betabinomial')
summary(seeds.mod)
car:::Anova(seeds.mod,type=2)


# are mass of sewn seeds and seed number highly correlated? - Yes ; r=.78
cor.test(log(mix_performance_new$mix_abundance),log(mix_performance_new$seeds_added), method="pearson")


# question 2a, how does species' establishment differ by conservatism ?
# from now on we use only fresh collected species for analysis

conservatism_recruitment_rates <- mix_performance_fresh_seed %>% 
  group_by(conservatism) %>% 
  summarise(recruit = mean(recruit_19, na.rm=T)) # what are the raw recruitment rates?

conservatism.mod_1 <-glm(recruit_19 ~ conservatism * l_mix_abundance, weights = quadrats_seeded, data=mix_performance_fresh_seed, family=quasibinomial(link="logit"))
conservatism.mod_2 <-glm(recruit_19 ~ conservatism + l_mix_abundance, weights = quadrats_seeded, data=mix_performance_fresh_seed, family=quasibinomial(link="logit"))
conservatism.mod_3 <-glm(recruit_19 ~ conservatism, weights = quadrats_seeded, data=mix_performance_fresh_seed, family=quasibinomial(link="logit"))
conservatism.mod_4 <-glm(recruit_19 ~ l_mix_abundance, weights = quadrats_seeded, data=mix_performance_fresh_seed, family=quasibinomial(link="logit"))
anova(conservatism.mod_1, conservatism.mod_2, test= "Chi") # models are not different, so interaction term is not necessary
anova(conservatism.mod_2, conservatism.mod_3, test= "Chi") # reduced model with conservatism and log abundance as a co-variate is the best fit 
anova(conservatism.mod_2, conservatism.mod_4, test= "Chi")
anova(conservatism.mod_3, test = "F") # without accounting for abundance seeded, conservatism is not significant
anova(conservatism.mod_2, test = "F") # TEST RESULTS: even with abundance in the model, conservatism groups are not significant
summary(conservatism.mod_2) # data are very over-dispersed
plot(conservatism.mod_2)

# using betabinomial
# conservatism_beta.mod_1 <-glmmTMB(recruit_19 ~ conservatism * l_mix_abundance , weights = quadrats_seeded, data=mix_performance, family='betabinomial')
# conservatism_beta.mod_2 <-glmmTMB(recruit_19 ~ conservatism + l_mix_abundance , weights = quadrats_seeded, data=mix_performance, family='betabinomial')
# conservatism_beta.mod_3 <-glmmTMB(recruit_19 ~ conservatism , weights = quadrats_seeded, data=mix_performance, family='betabinomial')
# conservatism_beta.mod_4 <-glmmTMB(recruit_19 ~ l_mix_abundance , weights = quadrats_seeded, data=mix_performance, family='betabinomial')
# anova(conservatism_beta.mod_1, conservatism_beta.mod_2, test= "Chi") # models are different, so interaction term is significant
# drop1(conservatism_beta.mod_1, test="Chisq")


# table(mix_performance$recruit_19 >0, mix_performance$conservatism)
# plot(simulateResiduals(conservatism_beta.mod_1)) # looks alright


# plotting effects of conservatism on recruitment
conservatism_pred <- summary(emmeans(conservatism.mod_2, ~ conservatism , type= "response"))
conservatism_means <- emmeans(conservatism.mod_2, pairwise ~ conservatism , type= "response")
conservatism_means

cons_plot = ggplot(data = conservatism_pred, aes(x=reorder(conservatism,-prob), y=prob)) + 
  theme_classic() +
  geom_point(stat="identity", color="black", size=2) +
  geom_errorbar(aes(ymax=prob+SE, ymin=prob-SE), width=0.2, size=1) + #Set errobars = prob +/- SE
  scale_x_discrete(labels = c('Ruderal','Matrix','Conservative')) + #capitalize the x axis categories
  xlab("Conservatism group") + ylab("Recruitment rate") +	#customize axis labels
  #scale_fill_manual("CCCCC") + #Use this to customize colors. 	
  ylim(0,0.25) +	#set y axis limits
  #geom_signif(y_position = c(0.8, 0.8), xmin = c(.75,2.75), xmax = c(2.25, 3.25),
  #            annotation = c("A", "B"), tip_length = 0.05 , size=1, textsize = 6) +
  geom_text(aes(label = c("A", "A", "A"), y=prob+SE, size=6), vjust = -1) +
  theme(panel.grid = element_blank(), #this removes almost everything
        panel.border = element_blank(), #I forgot what this was for
        text = element_text(size = 18), #increase text size
        legend.text=element_text(size=18), #legend text size 
        legend.title=element_blank(), #Legend title size
        axis.ticks.length = unit(0.25, 'cm'),
        axis.line=element_line(size=1.25), #axis line size
        axis.text.y=element_text(size=14), #y-axis text size
        axis.text.x=element_text(size=14, margin = unit(c(0, 0, 3, 0), "mm")), #x axis text size
        axis.title=element_text(size=16), #axis title size
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm")),
        legend.position="none", 	#move legend to top (if you don't need this, move erase this line, add an extra ) and + to the previous line, and remove last comma)
        plot.title = element_text(hjust=0, face="plain", size=20)) +
  ggtitle("")

cons_plot
# ggsave(file="conservatism.png", cons_plot, units="in", width=9,height=6, scale=0.8, dpi=300)	

# plot interactive effects of sown abundance and conservatism group on recruitment
conserv_pred_frame <- with(mix_performance_fresh_seed,
                expand.grid(l_mix_abundance=seq(min(l_mix_abundance), max(l_mix_abundance), length =120), conservatism = levels(conservatism), quadrats_seeded=(quadrats_seeded)))
conserv_pred_frame$recruit <- predict(conservatism.mod_2, newdata = conserv_pred_frame, type="response" )

abundance_plot <- ggplot(mix_performance_fresh_seed, aes(y= recruit_19, x= log(mix_abundance), color = conservatism, shape = conservatism )) + 
  theme_classic() + 
  geom_point(size=3) +
  scale_color_manual(values=cb_palette_2)+
  geom_line(data=conserv_pred_frame, aes(y=recruit, x=l_mix_abundance, col=conservatism), linetype = "dashed", size=1) +
  ylim(0,1.01) +
  theme(panel.grid = element_blank(), #this removes almost everything
        panel.border = element_blank(), #I forgot what this was for
        text = element_text(size = 18), #increase text size
        legend.text=element_text(size=18), #legend text size 
        axis.ticks.length = unit(0.25, 'cm'),
        axis.line=element_line(size=1.25), #axis line size
        axis.text.y=element_text(size=14), #y-axis text size
        axis.text.x=element_text(size=14,margin = unit(c(0, 0, 3, 0), "mm")), #x-axis text size
        axis.title=element_text(size=16), #axis title size
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm")),
        legend.position="right")+ 	#move legend to top 
  coord_cartesian(xlim = c(-2.5, 7)) +
  labs(x="Log abundance sown", y = "Recruitment rate", col="Conservatism", shape="Conservatism")
abundance_plot
# ggsave(file="abundance_plot.png", abundance_plot, units="in",scale=0.8, width=10,height=8, dpi=300)	

# There was not a significant interaction between conservatism and sown abundance 
# same plot but with predictions for all data, not based on conservatism.

conserv_pred_frame_2 <- with(mix_performance_fresh_seed,
              expand.grid(l_mix_abundance=seq(min(l_mix_abundance), max(l_mix_abundance), length =120), conservatism = NA, quadrats_seeded=(quadrats_seeded)))
conserv_pred_frame_2$recruit <- predict(conservatism.mod_4, newdata = conserv_pred_frame_2, type="response" )


abundance_plot_2 <- ggplot(mix_performance_fresh_seed, aes(y= recruit_19, x= log(mix_abundance), color = conservatism, shape = conservatism )) + 
  theme_classic() + 
  geom_point(size=3) +
  geom_line(data=conserv_pred_frame_2, aes(y=recruit, x=l_mix_abundance), col="black", linetype = "dashed", size=1) +
  ylim(0,1.01) +
  scale_color_manual(values=cb_palette_3,labels=c("Conservative","Matrix", "Ruderal")) +
  scale_shape(labels=c("Conservative","Matrix", "Ruderal")) +
  theme(panel.grid = element_blank(), #this removes almost everything
        panel.border = element_blank(), #I forgot what this was for
        text = element_text(size = 18), #increase text size
        legend.text=element_text(size=18), #legend text size 
        axis.ticks.length = unit(0.25, 'cm'),
        axis.line=element_line(size=1.25), #axis line size
        axis.text.y=element_text(size=14), #y-axis text size
        axis.text.x=element_text(size=14,margin = unit(c(0, 0, 3, 0), "mm")), #x-axis text size
        axis.title=element_text(size=16), #axis title size
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm")),
        legend.position="right")+ 	#move legend to top 
  coord_cartesian(xlim = c(-2.5, 7)) +
  labs(x="Log abundance sown", y = "Recruitment rate", col="Conservatism", shape="Conservatism")
abundance_plot_2
# ggsave(file="abundance_plot_2.png", abundance_plot_2, units="in",scale=0.8, width=10,height=8, dpi=300)	


# question 2b, how does FORB establishment differ by conservatism ?

seeded_forbs<- mix_performance_fresh_seed %>% filter (fg == "forb" )
# 116 of the 169 seeded species are forbs
# 79 of the 121 fresh seeded species are forbs

conservatism_forb.mod_1 <-glm(recruit_19 ~ conservatism * log(mix_abundance) , weights = quadrats_seeded, data=seeded_forbs, family=quasibinomial(link="logit"))
conservatism_forb.mod_2 <-glm(recruit_19 ~ conservatism + log(mix_abundance) , weights = quadrats_seeded, data=seeded_forbs, family=quasibinomial(link="logit"))
conservatism_forb.mod_3 <-glm(recruit_19 ~ conservatism, weights = quadrats_seeded, data=seeded_forbs, family=quasibinomial(link="logit"))
conservatism_forb.mod_4 <-glm(recruit_19 ~ l_mix_abundance, weights = quadrats_seeded, data=seeded_forbs, family=quasibinomial(link="logit"))
anova(conservatism_forb.mod_1, conservatism_forb.mod_2, test= "Chi") # interaction is not significant
anova(conservatism_forb.mod_2, conservatism_forb.mod_3, test= "Chi") # reduced model with conservatism and log abundance as a co-variate is the best fit 
anova(conservatism_forb.mod_3, test = "F") # without accounting for abundance seeded, conservatism is not significant
anova(conservatism_forb.mod_2, test = "F") # TEST RESULTS: conservatism and seeded abundance both predict forb recruitment
summary(conservatism_forb.mod_2) # data are very over-dispersed
plot(conservatism_forb.mod_2)

# using beta-binomial
# conservatism_forb_beta.mod_1 <-glmmTMB(recruit_19 ~ conservatism * log(mix_abundance) , weights = quadrats_seeded, data=seeded_forbs, family=betabinomial)
# conservatism_forb_beta.mod_2 <-glmmTMB(recruit_19 ~ conservatism + log(mix_abundance) , weights = quadrats_seeded, data=seeded_forbs, family=betabinomial)
# conservatism_forb_beta.mod_3 <-glmmTMB(recruit_19 ~ conservatism, weights = quadrats_seeded, data=seeded_forbs, family=betabinomial)
# anova(conservatism_forb_beta.mod_1, conservatism_forb_beta.mod_2, test= "Chi") # interaction is almost significant p=.07
# anova(conservatism_forb_beta.mod_2, conservatism_forb_beta.mod_3, test= "Chi")

# table(seeded_forbs$recruit_19 >0, seeded_forbs$conservatism)
# plot(simulateResiduals(conservatism_forb_beta.mod_2)) # looks alright


# plotting effects of conservatism on recruitment
conservatism_forb_pred <- summary(emmeans(conservatism_forb.mod_2, ~ conservatism , type= "response"))
conservatism_forb_means <- emmeans(conservatism_forb.mod_2,pairwise ~ conservatism , type= "response")
conservatism_forb_means

forb_cons_plot = ggplot(data = conservatism_forb_pred, aes(x=reorder(conservatism,-prob), y=prob)) +
  theme_classic() +
  geom_point(stat="identity", color="black",size=2) +
  geom_errorbar(aes(ymax=prob+SE, ymin=prob-SE), width=0.2, size=1) + #Set errobars = prob +/- SE
  scale_x_discrete(labels = c('Ruderal','Matrix','Conservative')) + #capitalize the x axis categories
  xlab("Conservatism group") + ylab("Forb recruitment rate") +	#customize axis labels
  #scale_fill_manual("CCCCC") + #Use this to customize colors. 	
  ylim(0,0.5) +	#set y axis limits
  #geom_signif(y_position = c(0.4, 0.3, 0.2), xmin = c(.8,1.8,2.8), xmax = c(1.2,2.2, 3.2),
  #            annotation = c("A", "AB", "B"), tip_length = 0.05 , size=1, textsize = 6) +
  geom_text(aes(label = c("B", "B", "A"), y=prob+SE, size=6), vjust = -1) +
  theme(panel.grid = element_blank(), #this removes almost everything
        panel.border = element_blank(), #I forgot what this was for
        text = element_text(size = 18), #increase text size
        legend.text=element_text(size=18), #legend text size 
        legend.title=element_blank(), #Legend title size
        axis.ticks.length = unit(0.25, 'cm'),
        axis.line=element_line(size=1.25), #axis line size
        axis.text.y=element_text(size=14), #y-axis text size
        axis.text.x=element_text(size=14, margin = unit(c(0, 0, 3, 0), "mm")), #x axis text size
        axis.title=element_text(size=16), #axis title size
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm")),
        legend.position="none", 	#move legend to top (if you don't need this, move erase this line, add an extra ) and + to the previous line, and remove last comma)
        plot.title = element_text(hjust=0, face="plain", size=20))
forb_cons_plot
# ggsave(file="forb_conservatism.png", forb_cons_plot, units="in", width=9,height=6,scale=0.8,  dpi=300)	


## generate prediction frame

forb_conserv_pframe <- with(seeded_forbs,
               expand.grid(l_mix_abundance=seq(min(l_mix_abundance), max(l_mix_abundance), length =120),
                           conservatism=NA, quadrats_seeded=(quadrats_seeded)))

forb_conserv_pframe$recruit <- predict(conservatism_forb.mod_4,newdata =forb_conserv_pframe, type="response" )

ggplot(seeded_forbs, aes(x=l_mix_abundance,y=recruit_19 , col = conservatism)) +
  geom_point() +
  geom_line(data=forb_conserv_pframe, aes(y=recruit))  

forb_abundance_plot <- ggplot(seeded_forbs, aes(y= recruit_19, x= l_mix_abundance, color = conservatism, shape = conservatism )) + 
  theme_classic() + 
  geom_point(size=3) +
  scale_color_manual(values =cb_palette_3, labels=c("Conservative","Matrix", "Ruderal")) +
  scale_shape(labels=c("Conservative","Matrix", "Ruderal")) +
  geom_line(data=forb_conserv_pframe, aes(y=recruit, x=l_mix_abundance),col="black", linetype = "dashed", size=1) +
  ylim(0,1.01) +
  theme(panel.grid = element_blank(), #this removes almost everything
        panel.border = element_blank(), #I forgot what this was for
        text = element_text(size = 18), #increase text size
        legend.text=element_text(size=18), #legend text size 
        axis.ticks.length = unit(0.25, 'cm'),
        axis.line=element_line(size=1.25), #axis line size
        axis.text.y=element_text(size=14), #y-axis text size
        axis.text.x=element_text(size=14,margin = unit(c(0, 0, 3, 0), "mm")), #x-axis text size
        axis.title=element_text(size=16), #axis title size
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm")),
        legend.position="right")+ 	#move legend to top 
  coord_cartesian(xlim = c(-2.5, 7)) +
  labs(x="Log abundance sown", y = "Forb recruitment rate", col="Conservatism", shape="Conservatism" )
forb_abundance_plot
# ggsave(file="forb_abundance_plot.png", forb_abundance_plot, units="in",scale=0.8, width=10,height=8, dpi=300)
# ggsave(file="forb_abundance_plot_alt.png", forb_abundance_plot, units="in",scale=1.5, width=4,height=3, dpi=300)

# question 3, what proportion of each functional group sewn established ?

table(mix_performance_fresh_seed$recruit_19 >0, mix_performance_fresh_seed$fg)
mix_performance_fresh_seed$fg <- factor(mix_performance_fresh_seed$fg)

fg_performance <- mix_performance_fresh_seed %>%
  select (taxon, fg, mix_abundance, recruit_19, quadrats_seeded) %>%
  filter (fg != "shrub" ) %>%
  mutate(l_mix_abundance=log(mix_abundance),
         centered = scale(l_mix_abundance, scale=F))
fg_performance$fg <- factor(fg_performance$fg)
table(fg_performance$fg)
fg_performance <- droplevels(fg_performance)

sown_max <- group_by(fg_performance,fg) %>%
  summarise (
    max_l_abundance = max(l_mix_abundance),
    max_abundance = max(mix_abundance)
  ) %>%
  ungroup()

sedges <- fg_performance %>% filter (fg=="sedge") %>% mutate (fg_center_l_abundance = scale(l_mix_abundance, center = T, scale=F))
grasses <- fg_performance %>% filter (fg=="grass") %>% mutate (fg_center_l_abundance = scale(l_mix_abundance, center = T, scale=F))
legumes <- fg_performance %>% filter (fg=="legume") %>% mutate (fg_center_l_abundance = scale(l_mix_abundance, center = T, scale=F))
forbs <- fg_performance %>% filter (fg=="forb") %>% mutate (fg_center_l_abundance = scale(l_mix_abundance, center = T, scale=F))
combined <- rbind(sedges,grasses,legumes,forbs)

fg_performance_merged <- fg_performance %>%
  merge(y=sown_max,by="fg") %>%
  mutate(scaled_l_abundance = l_mix_abundance/max_l_abundance,
         scaled_abundance = mix_abundance/max_abundance) %>%
  merge(y=combined[,c("taxon","fg_center_l_abundance")],by="taxon")

rm(sedges,grasses,legumes,forbs,combined)

fg_performance_merged$fg_center_l_abundance = c(fg_performance_merged$fg_center_l_abundance)
fg.mod_1 <-glm(recruit_19 ~ fg * fg_center_l_abundance , weights = quadrats_seeded, data=fg_performance_merged, family=quasibinomial(link="logit"))
fg.mod_2 <-glm(recruit_19 ~ fg + fg_center_l_abundance , weights = quadrats_seeded, data=fg_performance_merged, family=quasibinomial(link="logit"))
fg.mod_3 <-glm(recruit_19 ~ fg, weights = quadrats_seeded, data=fg_performance_merged, family=quasibinomial(link="logit"))
anova(fg.mod_1, fg.mod_2, test= "Chi") # interaction is significant
anova(fg.mod_1, test = "F") 
anova(fg.mod_3, test = "F") # functional groups are significantly different on their own
summary(fg.mod_1) # data are very over-dispersed
plot(fg.mod_1)

fg_pred_1 <- summary(emmeans(fg.mod_1, ~ fg , type= "response"))
fg_pred_2 <- summary(emmeans(fg.mod_2, ~ fg , type= "response"))
fg_pred_3 <- summary(emmeans(fg.mod_3, ~ fg , type= "response"))
plot(fg_pred_1)
plot(fg_pred_2)
plot(fg_pred_3)

fg_means <- emmeans(fg.mod_1, pairwise ~ fg , type= "response")
fg_means


# analysis with beta-binomial distribution

# fg_beta.mod_1 <-glmmTMB(recruit_19 ~ fg * fg_center_l_abundance , weights = quadrats_seeded, data=fg_performance_merged, family=betabinomial)
# fg_beta.mod_2 <-glmmTMB(recruit_19 ~ fg + fg_center_l_abundance , weights = quadrats_seeded, data=fg_performance_merged, family=betabinomial)
# fg_beta.mod_3 <-glmmTMB(recruit_19 ~ fg , weights = quadrats_seeded, data=fg_performance_merged, family=betabinomial)
# anova(fg_beta.mod_1, fg_beta.mod_2, test= "Chi") # interaction is significant
# drop1(fg_beta.mod_1, test="Chisq")

# summary(fg_beta.mod_1)
# plot(simulateResiduals(fg_beta.mod_1))

# fg_pred_beta_1 <- summary(emmeans(fg_beta.mod_1, ~ fg , type= "response"))
# fg_pred_beta_2 <- summary(emmeans(fg_beta.mod_2, ~ fg , type= "response"))
# fg_pred_beta_3 <- summary(emmeans(fg_beta.mod_3, ~ fg , type= "response"))
# plot(fg_pred_beta_1)
# plot(fg_pred_beta_2)
# plot(fg_pred_beta_3)

fg_recruitment_rates <- fg_performance %>% 
  group_by(fg) %>% 
  summarise(recruit = mean(recruit_19,na.rm=T)) # what are the raw recruitment rates?

# fg_means_beta <- emmeans(fg_beta.mod_1, pairwise ~ fg , type= "response")
# fg_means_beta


# visualization of effect of functional group on likelihood of recruitment is dependent on where in the abundance sown axis you sample from
# trial1<-emmeans(fg_beta.mod_1, ~fg, type="response", at = list(scaled_abundance = -0.25))
# trial2<-emmeans(fg_beta.mod_1, ~fg, type="response", at = list(scaled_abundance = 0))
# trial3<-emmeans(fg_beta.mod_1, ~fg, type="response", at = list(scaled_abundance = 0.25))
# trial4<-emmeans(fg_beta.mod_1, ~fg, type="response", at = list(scaled_abundance = 0.5))
# trial5<-emmeans(fg_beta.mod_1, ~fg, type="response", at = list(scaled_abundance = 0.75))
# trial6<-emmeans(fg_beta.mod_1, ~fg, type="response", at = list(scaled_abundance = 1))
# plot(summary(trial1))
# plot(summary(trial2))
# plot(summary(trial3))
# plot(summary(trial4))
# plot(summary(trial5))
# plot(summary(trial6))


# instead we can use emtrends to compare the slopes of functional groups
fg_estimates <- emtrends(fg.mod_1, pairwise ~ fg, var = "fg_center_l_abundance")
fg_estimates
emmip(fg.mod_1, fg ~ fg_center_l_abundance, cov.reduce = range)

fg_figure_data <- as.data.frame(fg_estimates$emtrends) %>%
  mutate (fg = str_to_title(fg))

fg_slope_plot <- ggplot(fg_figure_data, aes(x=fg_center_l_abundance.trend, y=fg)) + 
  geom_line() +
  geom_point(size=2) +
  geom_errorbarh(aes(xmin=asymp.LCL, xmax=asymp.UCL), height=.2,
                position=position_dodge(0.05)) +
  labs(x="Slope Estimate", y = "Functional Group") +
  theme_classic() + 
  theme(panel.grid = element_blank(), #this removes almost everything
        panel.border = element_blank(), #I forgot what this was for
        text = element_text(size = 18), #increase text size
        legend.text=element_text(size=18), #legend text size 
        axis.ticks.length = unit(0.25, 'cm'),
        axis.line=element_line(size=1.25), #axis line size
        axis.text.y=element_text(size=14), #y-axis text size
        axis.text.x=element_text(size=14,margin = unit(c(0, 0, 3, 0), "mm")), #x-axis text size
        axis.title=element_text(size=16), #axis title size
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm")),
        legend.position="right" ) + 	#move legend to top 
  geom_vline(xintercept = 0, linetype = "dotted")
print(fg_slope_plot)


figure_s2 <- ggplot(fg_performance_merged, aes(x=l_mix_abundance, y=fg)) + 
  geom_boxplot() +
  geom_jitter(size = 2, alpha=.3, height = .1) +
  labs(x="log sown abundance", y = "Functional Group") +
  theme_classic() + 
  theme(panel.grid = element_blank(), #this removes almost everything
        panel.border = element_blank(), #I forgot what this was for
        text = element_text(size = 18), #increase text size
        legend.text=element_text(size=18), #legend text size 
        axis.ticks.length = unit(0.25, 'cm'),
        axis.line=element_line(size=1.25), #axis line size
        axis.text.y=element_text(size=14), #y-axis text size
        axis.text.x=element_text(size=14,margin = unit(c(0, 0, 3, 0), "mm")), #x-axis text size
        axis.title=element_text(size=16), #axis title size
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm")),
        legend.position="right" )  	#move legend to top 
print(figure_s2)
ggsave(file="figure_s2.png", figure_s1, units="in",scale=0.8, width=10,height=8, dpi=300)


# plot the effect of fg on recruitment accounting for abundance
fg_pred_plot = ggplot(data = fg_pred_1, aes(x=reorder(fg,-prob), y=prob)) + 
  theme_classic() +
  geom_bar(stat="identity", color="black") +
  geom_errorbar(aes(ymax=prob+SE, ymin=prob-SE), width=0.2, size=1) + #Set errobars = prob +/- SE
  xlab("Funtional group") + ylab("Recruitment rate") +	#customize axis labels
  #scale_fill_manual("CCCCC") + #Use this to customize colors. 	
  ylim(0,1.01) +	#set y axis limits
  theme(panel.grid = element_blank(), #this removes almost everything
        panel.border = element_blank(), #I forgot what this was for
        text = element_text(size = 18), #increase text size
        legend.text=element_text(size=18), #legend text size 
        legend.title=element_blank(), #Legend title size
        axis.ticks.length = unit(0.25, 'cm'),
        axis.line=element_line(size=1.25), #axis line size
        axis.text.y=element_text(size=14), #y-axis text size
        axis.text.x=element_text(size=14,margin = unit(c(0, 0, 3, 0), "mm")), #x axis text size
        axis.title=element_text(size=16)) #axis title size
fg_pred_plot


# plot the interactive effect of fg and sown abundance

fg_pframe <- with(fg_performance_merged,
                  expand.grid(fg_center_l_abundance=seq(min(fg_center_l_abundance), max(fg_center_l_abundance), length =120),
                              fg=levels(fg), quadrats_seeded=(unique(quadrats_seeded))))
fg_pframe$recruit <- predict(fg.mod_1, newdata =fg_pframe, type="response" )



fg_abundance_plot <- ggplot(fg_performance_merged, aes(y= recruit_19, x= fg_center_l_abundance, color = fg, shape= fg)) + 
  theme_classic() + 
  scale_color_manual(values =cb_palette_2, labels=c("Forb","Grass","Legume","Sedge")) +
  geom_point(size=3) +
  scale_shape_manual(values = c(16,17,15,18),labels=c("Forb","Grass","Legume","Sedge")) +
  geom_line(data=fg_pframe, aes(y=recruit, x=fg_center_l_abundance), linetype = "dashed", size=1) +
  ylim(0,1.01) +
  theme(panel.grid = element_blank(), #this removes almost everything
        panel.border = element_blank(), #I forgot what this was for
        text = element_text(size = 18), #increase text size
        legend.text=element_text(size=18), #legend text size 
        axis.ticks.length = unit(0.25, 'cm'),
        axis.line=element_line(size=1.25), #axis line size
        axis.text.y=element_text(size=14), #y-axis text size
        axis.text.x=element_text(size=14,margin = unit(c(0, 0, 3, 0), "mm")), #x-axis text size
        axis.title=element_text(size=16), #axis title size
        axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm")),
        legend.position="right")+ 	#move legend to top 
  coord_cartesian(xlim = c(-4, 5)) +
  labs(x="Centered Log abundance sown", y = "Recruitment rate", color = "Functional group", shape = "Functional group")
fg_abundance_plot


figure_4 <- ggarrange(cons_plot, abundance_plot_2 , forb_cons_plot, forb_abundance_plot, fg_slope_plot, fg_abundance_plot,
                      labels = c("A", "B", "C", "D", "E", "F"),
                      vjust = c(1.5,1.5,0,0,1,1),
                      ncol = 2, nrow = 3, widths = c(1,2))
figure_4
ggsave(file="figure_4.png", figure_4, units="in", width=8, height=6, scale=1.6, dpi=300)


#### extra - Analysis step 6 - effect of grasses on forb recruitment ####

# get data on max cover by species in 2019
quads_19 <- year_max %>%
  filter(year == "2019")

# separate data into 3 groups based on management unit / seed mixes used
quads_t_1 <- filter(quads_19, transect == "01" )
nrow(quads_t_1) # 154
quads_t_3_4 <- filter(quads_19, transect == "03" | transect == "04")
nrow(quads_t_3_4) # 301
quads_t_8_10 <- filter(quads_19, transect == "08" | transect == "10")
nrow(quads_t_8_10) # 252

# merge species traits into our 3 data frames corresponding to management areas
quads_t_1_merged <- merge(x=quads_t_1, y=imls_full_raw[,c( "taxon",
                          "seeded_1_2", "duration", "fg", "c_score", "conservatism")], by="taxon", all.x =T)
quads_t_3_4_merged <- merge(x=quads_t_3_4, y=imls_full_raw[,c( "taxon",
                          "seeded_3_6", "duration", "fg", "c_score", "conservatism")], by="taxon", all.x =T)
quads_t_8_10_merged <- merge(x=quads_t_8_10, y=imls_full_raw[,c( "taxon",
                          "seeded_7_10", "duration", "fg", "c_score", "conservatism")], by="taxon", all.x =T)

# calculate summary statistics by quadrat for each of 3 management areas
quads_t_1_summary <- group_by(quads_t_1_merged, transect, quadrat) %>%
  summarise(
    total_s = length(taxon),
    total = sum(max_cover),
    grass = sum(max_cover [fg == "grass"], na.rm=T),
    forb = sum(max_cover [fg == "forb"], na.rm=T),
    seeded_forb = sum(max_cover [fg == "forb" & seeded_1_2 == "yes"], na.rm=T),
    grass_s = length(na.omit(taxon [fg == "grass"])),
    forb_s = length(na.omit(taxon [fg == "forb"])),
    seeded_forb_s_sampled = length(na.omit(taxon [fg == "forb" & seeded_1_2 == "yes"])),
    seeded = sum(max_cover [seeded_1_2 == "yes"]),
    prop_seeded = seeded / total,
    seeded_s_sampled = length(na.omit(taxon[seeded_1_2 == "yes"])),
  ) %>%
  ungroup()

quads_t_3_4_summary <- group_by(quads_t_3_4_merged, transect, quadrat) %>%
  summarise(
    total_s = length(taxon),
    total = sum(max_cover),
    grass = sum(max_cover [fg == "grass"], na.rm=T),
    forb = sum(max_cover [fg == "forb"], na.rm=T),
    seeded_forb = sum(max_cover [fg == "forb" & seeded_3_6 == "yes"], na.rm=T),
    grass_s = length(na.omit(taxon [fg == "grass"])),
    forb_s = length(na.omit(taxon [fg == "forb"])),
    seeded_forb_s_sampled = length(na.omit(taxon [fg == "forb" & seeded_3_6 == "yes"])),
    seeded = sum(max_cover [seeded_3_6 == "yes"]),
    prop_seeded = seeded / total,
    seeded_s_sampled = length(na.omit(taxon[seeded_3_6 == "yes"])),
  ) %>%
  ungroup()

quads_t_8_10_summary <- group_by(quads_t_8_10_merged, transect, quadrat) %>%
  summarise(
    total_s = length(taxon),
    total = sum(max_cover),
    grass = sum(max_cover [fg == "grass"], na.rm=T),
    forb = sum(max_cover [fg == "forb"], na.rm=T),
    seeded_forb = sum(max_cover [fg == "forb" & seeded_7_10 == "yes"], na.rm=T),
    grass_s = length(na.omit(taxon [fg == "grass"])),
    forb_s = length(na.omit(taxon [fg == "forb"])),
    seeded_forb_s_sampled = length(na.omit(taxon [fg == "forb" & seeded_7_10 == "yes"])),
    seeded = sum(max_cover [seeded_7_10 == "yes"]),
    prop_seeded = seeded / total,
    seeded_s_sampled = length(na.omit(taxon[seeded_7_10 == "yes"])),
  ) %>%
  ungroup()

# combine data from the three management areas
seeded_quads_2019_summary <- rbind(quads_t_1_summary, quads_t_3_4_summary, quads_t_8_10_summary)%>%
  mutate(seeded_forb_s = case_when(transect=="01" ~ 66, 
                                   transect=="03" | transect == "04" ~53, 
                                   transect=="08" | transect == "10" ~49),
         forb_recruit = seeded_forb_s_sampled / seeded_forb_s)

plot(y=seeded_quads_2019_summary$forb_recruit, x=seeded_quads_2019_summary$grass)
plot(y=seeded_quads_2019_summary$seeded_forb_s_sampled, x=seeded_quads_2019_summary$grass)
# correlation tests


#### extra - species summary ####


species_summary <- merge(x=sampled_spp[,1:11], y=transect_mean, by= "taxon" , all.x = T)
sp_summary_wide <- pivot_wider(species_summary, names_from = c(year,transect) , values_from = mean_t_cover )
sp_summary_wide[, 12:41][is.na(sp_summary_wide[, 12:41])] <- 0
sp_summary_wide <- sp_summary_wide %>% 
  mutate (sum_17 = `2017_01` + `2017_02` + `2017_03`+ `2017_04`+ `2017_05`+ `2017_06`+ `2017_07`+ `2017_08`+ `2017_09`+ `2017_10` ,
          sum_18 = `2018_01` + `2018_02` + `2018_03`+ `2018_04`+ `2018_05`+ `2018_06`+ `2018_07`+ `2018_08`+ `2018_09`+ `2018_10` ,
          sum_19 = `2019_01` + `2019_02` + `2019_03`+ `2019_04`+ `2019_05`+ `2019_06`+ `2019_07`+ `2019_08`+ `2019_09`+ `2019_10`) %>%
  select (-c(12:41))
write_csv(sp_summary_wide, "species_summary.csv")



#### figure 2 conservatism and fg ratios ####
# combine the two conservatism plots

legend <- cowplot::get_legend(forb_conservatism_plot)
gg_legend <- ggplotify::as.ggplot(legend)
gg_legend

figure_2_new <- ggarrange(conservatism_plot, forb_conservatism_plot, fg_plot,gg_legend, 
                      labels = c("A", "B", "C"),
                      ncol = 2, nrow = 2,
                      common.legend = T,
                      legend = "none")
figure_2_new
ggsave(file="figure_2.png", figure_2_new, units="in", width=5,height=5, scale=1.2, dpi=300)




