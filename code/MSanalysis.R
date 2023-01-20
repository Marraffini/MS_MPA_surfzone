#' @title Analyses of surf zone seine data for manuscript "Responses of Surf Zone Fish Communities to Marine Protected Areas" 
#' @description Script imports data from DataOne, summarizes it, and preforms analysis (ANOVAs, nMDS, PERMANOVA) to examine the role of MPA status and region on surf zone fish abundance, richness, and biomass observed by seines and BRUVs. These act as examples of analyses presented in the manuscript. 
#' @author Michelle 


#Libraries
library("readxl")
library("tidyverse")
library("dplyr")
library("ggpubr")
library("vegan")
library("ggrepel")

# set seed so analyses are reproducible 
set.seed(123)

# data can be found on DataOne
# Jenifer Dugan and Michelle Marraffini. 2022. California coast, ecosystem surveys of sandy beaches and surf zones August 2019-February 2020. California Ocean Protection Council Data Repository. doi:10.25494/P6989J.
# seine observations
fish_seines <- read.csv("fishseine_dataone.csv", row.names=1) 

# site information
MPA_names <- read.csv("MPABeachFishMPAnames.csv") 

# add meta data
seines <- left_join(fish_seines, MPA_names)

seines$month_day <- paste0(seines$month_num, "_", seines$day)
seines$month_day_haul_no <- paste0(seines$month_day, "_", seines$haul_number)

# Convert biomass in kg 
seines$fish_weight_individual <- seines$fish_weight_individual/1000
seines$fish_weight <-seines$fish_weight/1000

# Series of code to generate a data frame of response metrics: CPUE, RPUE, BPUE

# generates a df of the number of hauls pulled per sampling day at each site
site_samplingdate_no_haul <- seines %>%
  group_by(region, mpa_status, site_name, site_pair, year, month_day, month) %>%
  dplyr::summarise( no_hauls = length(unique(month_day_haul_no))) 

# counts abundance, richness, and biomass for each haul 
catch_summary_by_site_haul <- seines %>%
  filter(species_code !="NOSP") %>% # remove observations that did not catch any fish
  group_by(region, mpa_status, site_name, site_pair, haul_number, year, month_day, month) %>%
  dplyr::summarise(no_fish = sum(count), 
                   sp_richness = length(unique(species_code)),
                   bio = sum(fish_weight, na.rm = TRUE)) # how much is caught in each haul

# join the two previous data frames to get final metrics 
catch_summary_by_site_sampling <- 
  site_samplingdate_no_haul[,c("region", "mpa_status", "site_name", "site_pair", "no_hauls", "year", "month_day", "month")] %>%
  left_join(catch_summary_by_site_haul) %>%
  group_by(region, mpa_status, site_name, site_pair, year, month_day, month) %>% 
  dplyr::summarise( no_hauls = no_hauls, 
                    # abundance
                    abundance = sum(no_fish), # add up all the hauls on a given sampling date
                    meanfish_perhaul = abundance/no_hauls, # all sites have different haul numbers
                    # richness
                    species_richness = sum(sp_richness),
                    meanrichness_perhaul=species_richness/no_hauls,
                    # biomass
                    biomass = sum(bio),
                    meanbiomass_perhaul= biomass/no_hauls
  ) %>% distinct()

catch_summary_by_site_sampling[is.na(catch_summary_by_site_sampling)] <- 0 # fill in the sampling dates at sites where no fish were caught

# Create a column to identify sampling event
catch_summary_by_site_sampling <- catch_summary_by_site_sampling %>%
  group_by(region, site_pair, site_name, year) %>%
  mutate(sampling_event = match(month_day, unique(month_day))) 


# Analysis ----------------------------------------------------------------
# in these ANOVAs replicate is the sampling event (n=3 per site for a total of 154)

# Richness 
# example of the model fitting and checking 
fit_rich <- (aov(log1p(meanrichness_perhaul)~mpa_status*region, data=catch_summary_by_site_sampling))

plot(fit_rich)  # examine assumptions on residuals
summary(fit_rich) # standard F test
# adjusted sums of squares
car::Anova(fit_rich, type=2)
TukeyHSD(fit_rich)


# elasmobranchs -----------------------------------------------------------
# repeat the above organization of data and analysis for a specific group. 
# this is an example from the manuscript, could repeat for an individual species or differnt group by changing the filter line
catch_elasmo_haul <- seines %>%
  filter(class == "Elasmobranchii") %>% 
  group_by(region, mpa_status, site_name,CA_MPA_Name_Short, site_pair, haul_number,year, month_day, month) %>%
  dplyr::summarise(no_fish = sum(count), 
                   sp_richness = length(unique(species_code)),
                   trophic_rich = length(unique(broadtrophic)),
                   bio = sum(fish_weight, na.rm = TRUE)) 

elasmo_catch_summary_sampling_date <-
  site_samplingdate_no_haul[,c("region", "mpa_status", "site_name", "site_pair", "no_hauls", "year", "month_day", "month")] %>% # uses the general site and haul information data 
  left_join(catch_elasmo_haul) %>%
  group_by(region, mpa_status, site_name, site_pair, year, month_day, month) %>%
  dplyr::summarise( no_hauls = no_hauls, # all sites have different haul numbers
                    abundance = sum(no_fish), # total abundance per site
                    meanfish_perhaul = abundance/no_hauls,
                    # richness
                    species_richness = sum(sp_richness),
                    meanrichness_perhaul=species_richness/no_hauls,
                    # biomass
                    biomass = sum(bio),
                    meanbiomass_perhaul= biomass/no_hauls
  ) %>%
  distinct()

elasmo_catch_summary_sampling_date[is.na(elasmo_catch_summary_sampling_date)] <- 0

# Elasmobranchs don't actively occur in the surf zone past Drake's MPA and were most abundant in the South coast. In the manuscript we only analyze that south region
elasmo_catch_summary_sampling_dateS <- elasmo_catch_summary_sampling_date %>% filter(region=="South")

# ANOVA
fit_elasmo_bio <- aov(log1p(meanbiomass_perhaul) ~ mpa_status, data = elasmo_catch_summary_sampling_dateS)

plot(fit_elasmo_bio)
car::Anova(fit_elasmo_bio, type=2)

# does not meet assumptions of normality
kruskal.test(meanbiomass_perhaul~ mpa_status, data=elasmo_catch_summary_sampling_dateS)


# Figure ------------------------------------------------------------------
# example of figures in the manuscript

# function to generate error bars
seFunc <- function(x){
  se <- sd(x, na.rm =T)/sqrt(sum(!is.na(x)))
  lims <- c(mean(x) - se, mean(x) + se)
  names(lims) <- c('ymin', 'ymax')
  return(lims)
}

# plot theme
MPA_theme_xaxislabels <- 
  theme(axis.title=element_text(size=14, color="black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(vjust=0.95, size = 12),
        axis.text.x = element_text(size=12, color="black", vjust = 0.95, angle = 45,  hjust = 0.95),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.ticks = element_line(color="black"),
        axis.line = element_line(color="black"),
        panel.background=element_blank(),
        panel.grid.major=element_blank(),
        legend.position="right",
        legend.text = element_text(size=12),
        strip.background=element_blank(),
        strip.text=element_text(size=1, color = "white", vjust=1, face = "bold")
  )


catch_summary_by_site_sampling$region <- factor(catch_summary_by_site_sampling$region, levels = c("North", "Central", "South")) #set levels for figure

richness_state <- ggplot(catch_summary_by_site_sampling, aes(mpa_status, meanrichness_perhaul, fill = mpa_status))+
  scale_y_continuous(expand = c(0, 0)) +
  stat_summary(geom = 'bar', fun=mean, color="black", position=position_dodge(.90))+
  stat_summary(geom = 'errorbar', fun.data=seFunc, width=0, position=position_dodge(.90))+
  ylab("Mean Richness (Species/seine)")+
  scale_fill_manual(name = " ",
                    values=c("#990000", "#3399FF"))+
  MPA_theme_xaxislabels

richness_region <- ggplot(catch_summary_by_site_sampling, aes(region, meanrichness_perhaul, fill = mpa_status))+
  scale_y_continuous(expand = c(0, 0)) +
  stat_summary(geom = 'bar', fun=mean, color="black", position=position_dodge(.90))+
  stat_summary(geom = 'errorbar', fun.data=seFunc, width=0, position=position_dodge(.90))+
  ylab("Mean Richness (Species/seine)")+
  scale_fill_manual(name = " ",
                    values=c("#990000", "#3399FF"))+
  MPA_theme_xaxislabels

# combine into a multipanel figure
figure_rich_seines <- ggpubr::ggarrange(richness_state, richness_region,  # list of plots
                                        labels = NULL, # labels
                                        common.legend = T, # COMMON LEGEND
                                        legend = "right", # legend position
                                        align = "hv", # Align them both, horizontal and vertical
                                        nrow = 1)  # number of rows

figure_rich_seines


# Supplemental LRR analysis -----------------------------------------------
#BPUE
pair_lrr_biomass <- catch_summary_by_site_sampling[,c("region", "mpa_status","year" , "sampling_event", "site_pair", "meanbiomass_perhaul")] %>%
  pivot_wider(
    names_from = mpa_status, values_from = meanbiomass_perhaul, values_fill = 0) 
pair_lrr_biomass <-pair_lrr_biomass%>% filter(MPA>0 & Reference>0) # removes 10 sampling events

pair_lrr_biomass$LRR <- log(pair_lrr_biomass$MPA/pair_lrr_biomass$Reference)

# add names of MPAs
MPA_names_only <- MPA_names %>%filter(mpa_status=="MPA")
pair_lrr_biomass <- merge(pair_lrr_biomass, MPA_names_only)

# t tests -----------------------------------------------------------------
# if p-value is greater than the significance level 0.05 the distribution of the differences (d) are not significantly different from normal distribution. 
hist(pair_lrr_biomass$LRR, breaks=15) #normal 
shapiro.test(pair_lrr_biomass$LRR)# normal

t.test(pair_lrr_biomass$LRR, alternative="two.sided") # p>a not sig different


# Community analysis ------------------------------------------------------
# In the manuscript we combine seine and bruv data into one analysis

# load in BRUV data
fish_bruvs <- read.csv("surf_zone_bruv_data.csv") 

bruvs <- left_join(fish_bruvs, MPA_names, by=c("site_code", "site_type", "site_name", "region", "site_pair", "mpa_status"))


# Seine -------------------------------------------------------------------
seine_community <- seines %>%
  filter(common_name !="no species", common_name !="unspecified") %>% # remove hauls without fish
  group_by(region, site_name, mpa_status, site_pair, common_name) %>%
  summarise(abundance = sum(count)) 

# BRUV --------------------------------------------------------------------
bruv_community <- bruvs %>%
  filter(class == "Actinopteri" | class == "Elasmobranchii") %>% # remove invertebrates
  filter(common_name !="no species", common_name !="unspecified") %>%
  group_by(region, site_name, mpa_status, site_pair, common_name) %>%
  summarise(abundance = sum(max_n, na.rm = TRUE)) 


# combine datasets --------------------------------------------------------
meta_data <- c("site_name",  "region", "mpa_status", "survey", "site_pair")

seine_comm_wide <- pivot_wider(seine_community, names_from = common_name, values_from = abundance, values_fill = 0)

bruv_comm_wide <- pivot_wider(bruv_community, names_from = common_name, values_from = abundance, values_fill =0)

seine_comm_wide$survey <- "Seine"
bruv_comm_wide$survey <- "BRUV"
#combine both surveys 
community_data_allfishsurveys <- full_join(seine_comm_wide,bruv_comm_wide, values_fill=0)

# relabel unknowns see manuscript for details
community_data_allfishsurveys$surfperch_spp <- paste0(community_data_allfishsurveys$`redtail surfperch`+community_data_allfishsurveys$`silver surfperch`+ community_data_allfishsurveys$`unknown surfperch`)

community_data_allfishsurveys <- community_data_allfishsurveys[,!colnames(community_data_allfishsurveys) %in% c("redtail surfperch", "silver surfperch", "unknown surfperch")]

colnames(community_data_allfishsurveys)[colnames(community_data_allfishsurveys)=="unknown rockfish"] <- "Sabastes spp"

colnames(community_data_allfishsurveys)[colnames(community_data_allfishsurveys)=="unknown croaker"] <- "Sciaenidae"


community_data_allfishsurveys[,!colnames(community_data_allfishsurveys) %in% meta_data] <- lapply(community_data_allfishsurveys[,!colnames(community_data_allfishsurveys) %in% meta_data], as.numeric) # make sure all abundances are numeric. NAs are introduced by coercion but we then replace them with 0s

community_data_allfishsurveys[is.na(community_data_allfishsurveys)] <- 0 # set any NAs to 0

community_data_allfishsurveys <- community_data_allfishsurveys[rowSums(community_data_allfishsurveys[,!colnames(community_data_allfishsurveys) %in% meta_data])>5 , ] # removes 3 observations

# Do we have any sites with no fish?
zero <- community_data_allfishsurveys[rowSums(community_data_allfishsurveys[,!colnames(community_data_allfishsurveys) %in% meta_data]) == 0, ] # No

# Calculate relative abundances  -----------------------------------------------
# use relative abundances to set allow comparisons across surveys and different numbers of samples
catch_values_both <- community_data_allfishsurveys[,!colnames(community_data_allfishsurveys) %in% meta_data] # make a matrix of only species and their abundances

# use this code instead to scale with SD 
# scale to st dev units: scale() without centering uses root-mean-square not SD
out <- catch_values_both %>% 
  map_df(~ scale(.x, center = FALSE, scale = sd(.x, na.rm = TRUE)))

#_____ compute NMDS scores for each species ________#
fish_NMDS <- metaMDS(out) # bray distance
stressplot(fish_NMDS)
fish_NMDS 

# PERMANOVA ---------------------------------------------------------------
fish_dist <- vegdist(out, method='bray')
perm <- how(nperm = 9999)
fish_div <- adonis2(fish_dist ~ survey*mpa_status*region, data=community_data_allfishsurveys, permutations = perm, method="bray")
fish_div 

# nMDS figure -------------------------------------------------------------
# gather data for figure
plot_df <- scores(fish_NMDS, display = "sites") %>% # extract NMDS 1 and NMDS 2 coordinates
  as.data.frame() %>% 
  rownames_to_column("sit") %>%
  cbind(community_data_allfishsurveys[, colnames(community_data_allfishsurveys) %in% meta_data]) # Add factors to new NMDS score dataframe for plotting

#__________ compute vectors for species to add to nMDS plot __________#
fish_fit <- envfit(fish_NMDS, community_data_allfishsurveys, perm = 999)
fish_fit # vector loading for each species - NOTE - r2 values are used to scale the NMDS1 and NMDS2 values!

# extract p-values for each species
fish_pvals <- fish_fit$vectors$pvals %>% 
  as.data.frame() %>% 
  rownames_to_column("species") %>% 
  dplyr::rename("pvals" = ".")

# extract coordinates for species, only keep species with p-val < 0.01 because too many at p=0.05
fit_spp <- fish_fit %>% 
  scores(., display = "vectors") %>% 
  as.data.frame() %>% 
  rownames_to_column("species") %>% 
  full_join(., fish_pvals, by = "species") %>% 
  filter(pvals < 0.01)


# Make figures ------------------------------------------------------------
# _________ plot NMDS scores using ggplot ___________ #
MPA_nMDS <- ggplot(plot_df, aes(NMDS1, NMDS2)) +
  geom_point(aes(colour = mpa_status, shape = survey), size = 5, alpha = 0.8) +
  coord_fixed() +
  stat_ellipse(aes(colour = mpa_status), linetype = 2, size = 1) +
  scale_color_manual(name = "MPA",
                     values=c("#990000", "#3399FF"))+
  labs(caption="PERMANOVA: MPA effect: p = 0.80,  stress = 0.19",
       shape="Survey")+
  geom_segment(data = fit_spp,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour="black") +
  geom_text_repel(data = fit_spp, aes(x = NMDS1, y = NMDS2,
                                        label = species), size = 4, min.segment.length = 1, point.padding = 0.2, max.overlaps = 20) +
    theme(axis.title=element_text(size=12, color="black"),
          axis.text = element_text(size = 12, color = "black"),
          axis.ticks = element_line(color="black"),
          axis.line = element_line(color="black"),
          panel.background=element_blank(),
          panel.grid.major=element_blank(),
          legend.key = element_blank(),
          legend.position="right",
          legend.text = element_text(size=12),
          strip.background=element_blank(),
          strip.text=element_text(size=12, color = "white", vjust=1))

region_nMDS <- ggplot(plot_df, aes(NMDS1, NMDS2)) +
  geom_point(aes(colour = region, shape = mpa_status), size = 5, alpha = 0.8) +
  coord_fixed() +
  stat_ellipse(aes(colour = region), linetype = 2, size = 1) +
  scale_color_manual(name = "Region",
                     values=c("#438830", "#745886", "#0094B2"))+
  labs(caption="PERMANOVA: Region effect: p = 0.0001")+
  geom_segment(data = fit_spp,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour="black") +
  geom_text_repel(data = fit_spp, aes(x = NMDS1, y = NMDS2,
                                      label = species), size = 4, min.segment.length = 1, point.padding = 0.2, max.overlaps = 20) +
  labs(shape="MPA")+
  theme(axis.title=element_text(size=12, color="black"),
        axis.text = element_text(size = 12, color = "black"),
        axis.ticks = element_line(color="black"),
        axis.line = element_line(color="black"),
        panel.background=element_blank(),
        panel.grid.major=element_blank(),
        legend.key = element_blank(),
        legend.position="right",
        legend.text = element_text(size=12),
        strip.background=element_blank(),
        strip.text=element_text(size=12, color = "white", vjust=1))


survey_nMDS <- ggplot(plot_df, aes(NMDS1, NMDS2)) +
  geom_point(aes(colour = survey, shape = region), size = 5, alpha = 0.8) +
  coord_fixed() +
  stat_ellipse(aes(colour = survey), linetype = 2, size = 1) +
  scale_color_manual(name = "Survey",
                     values=c("#0066FF", "#FF6600"))+
  labs(caption="PERMANOVA: Survey effect: p = 0.0001, Survey x Region effect: p = 0.0001")+
  geom_segment(data = fit_spp,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour="black") +
  geom_text_repel(data = fit_spp, aes(x = NMDS1, y = NMDS2, 
                                      label = species), size = 4, min.segment.length = 1, point.padding = 0.2,  max.overlaps = 20) +
  labs(shape="Region")+
  theme(axis.title=element_text(size=12, color="black"),
        axis.text = element_text(size = 12, color = "black"),
        axis.ticks = element_line(color="black"),
        axis.line = element_line(color="black"),
        panel.background=element_blank(),
        panel.grid.major=element_blank(),
        legend.key = element_blank(),
        legend.position="right",
        legend.text = element_text(size=12),
        strip.background=element_blank(),
        strip.text=element_text(size=12, color = "white", vjust=1))


