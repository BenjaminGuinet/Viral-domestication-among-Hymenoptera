#BEGIN


# All needed functions

detachAllPackages <- function() {

  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")

  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]

  package.list <- setdiff(package.list,basic.packages)

  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)

}

library(treeio)

#The Hymenoptera tree
tree = read.newick("/Users/bguinet/Desktop/Papier_scientifique/Revbayes_anestral_states_estimates/data/TimeTree_mcmc_10_test_with_outgroup.nwk")
expand.grid.unique <- function(x, y, include.equals=FALSE)
{
  x <- unique(x)

  y <- unique(y)

  g <- function(i)
  {
    z <- setdiff(y, x[seq_len(i-include.equals)])

    if(length(z)) cbind(x[i], z, deparse.level=0)
  }

  do.call(rbind, lapply(seq_along(x), g))
}

library(stringi)
ggname <- function (prefix, grob) {
  grob$name <- grobName(grob, prefix)
  grob
}

geom_label2 <- function(mapping = NULL, data = NULL,
                        stat = "identity", position = "identity",
                        ...,
                        parse = FALSE,
                        nudge_x = 0,
                        nudge_y = 0,
                        label.padding = unit(0.25, "lines"),
                        label.r = unit(0.15, "lines"),
                        label.size = 0.25,
                        na.rm = FALSE,
                        show.legend = NA,
                        inherit.aes = TRUE) {
  if (!missing(nudge_x) || !missing(nudge_y)) {
    if (!missing(position)) {
      stop("Specify either `position` or `nudge_x`/`nudge_y`", call. = FALSE)
    }

    position <- position_nudge(nudge_x, nudge_y)
  }

  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomLabel2,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      parse = parse,
      label.padding = label.padding,
      label.r = label.r,
      label.size = label.size,
      na.rm = na.rm,
      ...
    )
  )
}

GeomLabel2 <- ggproto("GeomLabel2", Geom,
                      required_aes = c("x", "y", "label"),

                      default_aes = aes(
                        colour = "black", fill = "white", size = 3.88, angle = 0,
                        hjust = 0.5, vjust = 0.5, alpha = NA, family = "", fontface = 1,
                        lineheight = 1.2
                      ),

                      draw_panel = function(self, data, panel_params, coord, parse = FALSE,
                                            na.rm = FALSE,
                                            label.padding = unit(0.25, "lines"),
                                            label.r = unit(0.15, "lines"),
                                            label.size = 0.25) {
                        lab <- data$label
                        if (parse) {
                          lab <- parse(text = as.character(lab))
                        }

                        data <- coord$transform(data, panel_params)
                        if (is.character(data$vjust)) {
                          data$vjust <- compute_just(data$vjust, data$y)
                        }
                        if (is.character(data$hjust)) {
                          data$hjust <- compute_just(data$hjust, data$x)
                        }

                        grobs <- lapply(1:nrow(data), function(i) {
                          row <- data[i, , drop = FALSE]
                          labelGrob2(lab[i],
                                     x = unit(row$x, "native"),
                                     y = unit(row$y, "native"),
                                     just = "center",
                                     padding = label.padding,
                                     r = label.r,
                                     text.gp = gpar(
                                       col = row$colour,
                                       fontsize = row$size * .pt,
                                       fontfamily = row$family,
                                       fontface = row$fontface,
                                       lineheight = row$lineheight
                                     ),
                                     rect.gp = gpar(
                                       col = row$colour,
                                       fill = alpha(row$fill, row$alpha),
                                       lwd = label.size * .pt
                                     )
                          )
                        })
                        class(grobs) <- "gList"

                        ggname("geom_label", grobTree(children = grobs))
                      },

                      draw_key = draw_key_label
)

labelGrob2 <- function(label, x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                       just = "center", padding = unit(0.10, "lines"), r = unit(0.1, "snpc"),
                       default.units = "npc", name = NULL,
                       text.gp = gpar(), rect.gp = gpar(fill = "white"), vp = NULL) {

  stopifnot(length(label) == 1)

  if (!is.unit(x))
    x <- unit(x, default.units)
  if (!is.unit(y))
    y <- unit(y, default.units)

  gTree(label = label, x = x, y = y, just = just, padding = padding, r = r,
        name = name, text.gp = text.gp, rect.gp = rect.gp, vp = vp, cl = "labelgrob2")
}

makeContent.labelgrob2 <- function(x) {
  hj <- resolveHJust(x$just, NULL)
  vj <- resolveVJust(x$just, NULL)

  t <- textGrob(
    x$label,
    x$x + 1 * (0.55 - hj) * unit(5, "mm"),
    x$y + 2 * (0.55 - vj) * x$padding,
    just = "center",
    gp = x$text.gp,
    name = "text"
  )

  r <- roundrectGrob(x$x, x$y, default.units = "native",
                     width =  1.5 * unit(max(stri_width(x$x)) + 1, "mm"),
                     height = grobHeight(t) + 2 * x$padding,
                     just = c(hj, vj),
                     r = x$r,
                     gp = x$rect.gp,
                     name = "box"
  )

  setChildren(x, gList(r, t))
}










detachAllPackages()
library(data.table)
library(stringr)
library(randomcoloR)
library(dplyr)
library(phytools)
library(phylogram)
library(tibble)
library(fpc)
library(gplots)
library(ape)
#library(geiger)
#library(nlme)
library(rlist)
#library(tidyverse)
library(ggtree)


#
All_information_df<- read.csv("/Users/bguinet/Desktop/Papier_scientifique/Viralprot_vs_Viral_loci_result_all_match_and_cluster_taxid_ontology_2LCA_HMMER_Augustus_Repeat_Env_pseudo_HSP_Score_Filtred_Mono_TPM_dNdS_synteny_ORFs5.m8",sep=";",header=T,na.strings=c("","NA"))
#subdsDNA_tab2 <- All_information_df[c("Clustername","query","Event","Nsp","HSP_name")]

#All_information_df<-All_information_df[!All_information_df$Clustername=="Cluster3513",]

## hispanica
All_information_df$cov_depth_BUSCO[All_information_df$scaffold_name=="SGAY01008337.1"] <- 20.67932
All_information_df$pvalue_cov[All_information_df$scaffold_name=="SGAY01008337.1"]  <- 0.195438708
All_information_df$FDR_pvalue_cov[All_information_df$scaffold_name=="SGAY01008337.1"]<- "True"
All_information_df$cov_depth_BUSCO[All_information_df$scaffold_name=="SGAY01009995.1"] <- 180.65
All_information_df$pvalue_cov[All_information_df$scaffold_name=="SGAY01009995.1"] <- 0.000000000
All_information_df$FDR_pvalue_cov[All_information_df$scaffold_name=="SGAY01009995.1"]<- "False"
All_information_df$cov_depth_BUSCO[All_information_df$scaffold_name=="SGAY01012126.1"] <- 72.89557
All_information_df$pvalue_cov[All_information_df$scaffold_name=="SGAY01012126.1"]<- 0.013016845
All_information_df$FDR_pvalue_cov[All_information_df$scaffold_name=="SGAY01012126.1"]<- "True"

All_information_df$cov_depth_BUSCO[All_information_df$Species_name=="Dinoponera_quadriceps"] <- "No_cov"
All_information_df$cov_depth_BUSCO[All_information_df$Species_name=="Camponotus_floridanus"] <- "No_cov"
All_information_df$cov_depth_candidat[All_information_df$Species_name=="Dinoponera_quadriceps"] <- "No_cov"
All_information_df$cov_depth_candidat[All_information_df$Species_name=="Camponotus_floridanus"] <- "No_cov"
All_information_df$pvalue_cov[All_information_df$Species_name=="Dinoponera_quadriceps"] <- "No_cov"
All_information_df$pvalue_cov[All_information_df$Species_name=="Camponotus_floridanus"] <- "No_cov"

All_information_df$Mean_dNdS[All_information_df$query_bis %in% c("scaffold_78_286199-286358_-__Meteorus_cinctellus","scaffold_59_101366-101546_+__Meteorus_colon_F","scaffold_146_125820-126000_+__Meteorus_colon_M")]<-0.0001
All_information_df$SE_dNdS[All_information_df$query_bis %in% c("scaffold_78_286199-286358_-__Meteorus_cinctellus","scaffold_59_101366-101546_+__Meteorus_colon_F","scaffold_146_125820-126000_+__Meteorus_colon_M")]<-0.02637698
All_information_df$Pvalue_dNdS[All_information_df$query_bis %in% c("scaffold_78_286199-286358_-__Meteorus_cinctellus","scaffold_59_101366-101546_+__Meteorus_colon_F","scaffold_146_125820-126000_+__Meteorus_colon_M")]<-0.003818
All_information_df$FDR_pvalue_dNdS[All_information_df$query_bis %in% c("scaffold_78_286199-286358_-__Meteorus_cinctellus","scaffold_59_101366-101546_+__Meteorus_colon_F","scaffold_146_125820-126000_+__Meteorus_colon_M")]<-1


All_information_df$Scaffold_score[All_information_df$query=="QANI01000034.1:914891-916274(+):Camponotus_floridanus"]<-"B"
All_information_df$Scaffold_score[All_information_df$query=="QANI01000237.1:118618-119563(+):Camponotus_floridanus"]<-"B"
All_information_df$Scaffold_score[All_information_df$query=="KQ468257.1:12848-13502(+):Dinoponera_quadriceps"]<-"B"
All_information_df$Scaffold_score[All_information_df$query=="KQ466978.1_2-HSPs_+__Dinoponera_quadriceps"]<-"B"
All_information_df$Scaffold_score[All_information_df$query=="KQ466978.1_2-HSPs_+__Dinoponera_quadriceps"]<-"B"


All_information_df<-All_information_df %>% mutate(Scaffold_species = coalesce(Scaffold_species, Scaff_name))



All_information_df <- All_information_df[!All_information_df$best_family_per_query=="Mimiviridae",]
#Remove two cluster containgin Ants species ich are false positive
All_information_df <- All_information_df[!All_information_df$Clustername %in% c("Cluster10968"),]




# Add family to the uknown RNA seq

Shi_RNA_classification<- read.csv("/Users/bguinet/Desktop/RNA_classification.csv",sep=";")

Shi_RNA_classification$Virus.Name..or.published.sequences.<- gsub(" (ref)","",Shi_RNA_classification$Virus.Name..or.published.sequences.)
Shi_RNA_classification$Virus.Name..or.published.sequences.<- gsub(" ","_",Shi_RNA_classification$Virus.Name..or.published.sequences.)

All_information_df$best_family_per_query2 <- Shi_RNA_classification$Classification[match(All_information_df$species, Shi_RNA_classification$Virus.Name..or.published.sequences.)]

All_information_df$best_family_per_query<- ifelse(All_information_df$best_family_per_query=="Unknown", All_information_df$best_family_per_query2, All_information_df$best_family_per_query)



All_information_df$genomic_structure[All_information_df$species=='Abisko_virus']<-'ssRNA'
All_information_df$best_family_per_query[All_information_df$species=='Abisko_virus']<-'Abisko-like'
All_information_df$best_family_per_Event[All_information_df$species=='Abisko_virus']<-'Abisko-like'
All_information_df$family[All_information_df$species=='Abisko_virus']<-'Abisko-like'
All_information_df$consensus_family[All_information_df$species=='Abisko_virus']<-'Abisko-like'


All_information_df$best_family_per_Event[All_information_df$best_family_per_query=="Partiti-Picobirna"] <- "Partiti-Picobirna"
All_information_df$family[All_information_df$best_family_per_query=="Partiti-Picobirna"] <- "Partiti-Picobirna"
All_information_df$genomic_structure[All_information_df$best_family_per_query=="Partiti-Picobirna"] <- "ssRNA"

All_information_df$best_family_per_Event[All_information_df$best_family_per_query=="Narna-Levi"] <- "Narna-Levi"
All_information_df$family[All_information_df$best_family_per_query=="Narna-Levi"] <- "Narna-Levi"

All_information_df$best_family_per_Event[All_information_df$best_family_per_query=="Mono-Chu"] <- "Mono-Chu"
All_information_df$family[All_information_df$best_family_per_query=="Mono-Chu"] <- "Mono-Chu"

All_information_df$best_family_per_Event[All_information_df$best_family_per_query=="Luteo-Sobemo"] <- "Luteo-Sobemo"
All_information_df$family[All_information_df$best_family_per_query=="Luteo-Sobemo"] <- "Luteo-Sobemo"

All_information_df$best_family_per_Event[All_information_df$best_family_per_query=="Bunya-Arenao"] <- "Bunya-Arena"
All_information_df$family[All_information_df$best_family_per_query=="Bunya-Arena"] <- "Bunya-Arena"

All_information_df$best_family_per_Event[All_information_df$best_family_per_query=="Picorna-Calici"] <- "Picorna-Calici"
All_information_df$family[All_information_df$best_family_per_query=="Picorna-Calici"] <- "Picorna-Calici"


All_information_df$best_family_per_Event[All_information_df$best_family_per_query=="Tombusviridae"] <- "Tombusviridae"
All_information_df$family[All_information_df$best_family_per_query=="Tombusviridae"] <- "Tombusviridae"
All_information_df$genomic_structure[All_information_df$best_family_per_query=="Tombusviridae"] <- "ssRNA"

All_information_df$genomic_structure[All_information_df$best_family_per_query=="Phenuiviridae"] <- "ssRNA"
All_information_df$best_family_per_Event[All_information_df$best_family_per_query=="Phenuiviridae"] <- "Phenuiviridae"
All_information_df$family[All_information_df$best_family_per_query=="Phenuiviridae"] <- "Phenuiviridae"



All_information_df$best_family_per_Event[All_information_df$All_information_df$best_family_per_query="Partiti-Picobirna"] <- "Partiti-Picobirna"
All_information_df$family[All_information_df$All_information_df$best_family_per_query=="Partiti-Picobirna"] <- "Partiti-Picobirna"
All_information_df$genomic_structure[All_information_df$All_information_df$best_family_per_query=="Partiti-Picobirna"] <- "dsRNA"

All_information_df$genomic_structure[All_information_df$best_family_per_query=="Chuviridae"] <-"ssRNA"
All_information_df$genomic_structure[All_information_df$best_family_per_query=="Cruciviridae"] <-"ssDNA"

All_information_df$best_family_per_query[All_information_df$species=="Dipteran_anphevirus"] <- "Xinmoviridae"
All_information_df$best_family_per_Event[All_information_df$species=="Dipteran_anphevirus"] <- "Xinmoviridae"
All_information_df$family[All_information_df$species=="Dipteran_anphevirus"] <- "Xinmoviridae"

All_information_df$genomic_structure[All_information_df$best_family_per_query=="Xinmoviridae"]<-"ssRNA"
All_information_df$genomic_structure[All_information_df$best_family_per_query=="Chuviridiae"]<-"ssRNA"
All_information_df$genomic_structure[All_information_df$best_family_per_query=="Lispiviridae"]<-"ssRNA"
All_information_df$genomic_structure[All_information_df$best_family_per_query=="Mono-Chu"]<-"ssRNA"
All_information_df$genomic_structure[All_information_df$best_family_per_query=="Narna-Levi"]<-"ssRNA"
All_information_df$genomic_structure[All_information_df$best_family_per_query=="Luteo-Sobemo"]<-"ssRNA"

All_information_df$genomic_structure[All_information_df$best_family_per_query=="Partiti-Picobirna"]<-"dsRNA"

All_information_df$genomic_structure[All_information_df$best_family_per_query=="Metaviridae"]<-"dsDNA"


All_information_df$genomic_structure[All_information_df$species=="Loreto_virus"]<-"ssRNA"
All_information_df$family[All_information_df$species=="Loreto_virus"]<-"Negevirus-like"
All_information_df$best_family_per_Event[All_information_df$species=="Loreto_virus"]<-"Negevirus-like"
All_information_df$best_family_per_query[All_information_df$species=="Loreto_virus"]<-"Negevirus-like"
All_information_df$genomic_structure[All_information_df$species=="Piura_virus"]<-"ssRNA"
All_information_df$family[All_information_df$species=="Piura_virus"]<-"Negevirus-like"
All_information_df$best_family_per_Event[All_information_df$species=="Piura_virus"]<-"Negevirus-like"
All_information_df$best_family_per_query[All_information_df$species=="Piura_virus"]<-"Negevirus-like"

All_information_df$best_family_per_Event[All_information_df$best_family_per_query =="Hepe-Virga"] <-"Hepe-Virga"
All_information_df$genomic_structure[All_information_df$best_family_per_query =="Hepe-Virga"] <-"ssRNA"

All_information_df$best_family_per_Event[All_information_df$species=="Shahe_heteroptera_virus_3"]<-"Bunya-Arena"
All_information_df$best_family_per_query[All_information_df$species=="Shahe_heteroptera_virus_3"]<-"Bunya-Arena"
All_information_df$family[All_information_df$species=="Shahe_heteroptera_virus_3"]<-"Bunya-Arena"
All_information_df$genomic_structure[All_information_df$species=="Shahe_heteroptera_virus_3"]<-"ssRNA"

All_information_df$family<- ifelse(All_information_df$family=="Unknown", All_information_df$best_family_per_query2, All_information_df$family)
All_information_df$consensus_family<- ifelse(All_information_df$consensus_family=="Unknown", All_information_df$best_family_per_query2, All_information_df$consensus_family)
# recover <-


All_information_df$genomic_structure[All_information_df$species=="Shuangao_insect_virus_7"]<-"ssRNA"
All_information_df$family[All_information_df$species=="Shuangao_insect_virus_7"]<-"Unknown"
All_information_df$best_family_per_query[All_information_df$species=="Shuangao_insect_virus_7"]<-"Unknown"

All_information_df$genomic_structure[All_information_df$species=="Wuhan_flea_virus" ]<-"ssRNA"
All_information_df$family[All_information_df$species=="Wuhan_flea_virus" ]<-"Unknown"
All_information_df$best_family_per_query[All_information_df$species=="Wuhan_flea_virus"]<-"Unknown"

All_information_df$genomic_structure[All_information_df$species=="Wuhan_aphid_virus_1"]<-"ssRNA"
All_information_df$family[All_information_df$species=="Wuhan_aphid_virus_1"]<-"Unknown"
All_information_df$best_family_per_query[All_information_df$species=="Wuhan_aphid_virus_1"]<-"Unknown"

All_information_df$genomic_structure[All_information_df$species=="Circulifer_tenellus_virus_1"]<-"dsRNA"
All_information_df$family[All_information_df$species=="Circulifer_tenellus_virus_1"]<-"Unknown"
All_information_df$best_family_per_query[All_information_df$species=="Circulifer_tenellus_virus_1"]<-"Unknown"

All_information_df$genomic_structure[All_information_df$species=="Persimmon_latent_virus"]<-"dsRNA"
All_information_df$family[All_information_df$species=="Persimmon_latent_virus"]<-"Unknown"
All_information_df$best_family_per_query[All_information_df$species=="Persimmon_latent_virus"]<-"Unknown"

All_information_df$genomic_structure[All_information_df$species=="Shuangao_lacewing_virus_2"]<-"ssRNA"
All_information_df$family[All_information_df$species=="Shuangao_lacewing_virus_2"]<-"Unknown"
All_information_df$best_family_per_query[All_information_df$species=="Shuangao_lacewing_virus_2"]<-"Unknown"

All_information_df$genomic_structure[All_information_df$species=="Shayang_fly_virus_4"]<-"ssRNA"
All_information_df$family[All_information_df$species=="Shayang_fly_virus_4"]<-"Unknown"
All_information_df$best_family_per_query[All_information_df$species=="Shayang_fly_virus_4"]<-"Unknown"

All_information_df$genomic_structure[All_information_df$species=="Spissistilus_festinus_virus_1"]<-"dsRNA"
All_information_df$family[All_information_df$species=="Spissistilus_festinus_virus_1"]<-"Unknown"
All_information_df$best_family_per_query[All_information_df$species=="Spissistilus_festinus_virus_1"]<-"Unknown"

All_information_df$genomic_structure[All_information_df$species=="Spissistilus_festinus_virus_1"]<-"dsRNA"
All_information_df$family[All_information_df$species=="Spissistilus_festinus_virus_1"]<-"Unknown"
All_information_df$best_family_per_query[All_information_df$species=="Spissistilus_festinus_virus_1"]<-"Unknown"

All_information_df$genomic_structure[is.na(All_information_df$genomic_structure)]<-"Unknown"



All_information_df$best_family_per_Event[All_information_df$family == "Unknown" & All_information_df$genomic_structure =="ssRNA"]<-"Unknown"
All_information_df$best_family_per_query[All_information_df$family == "Unknown" & All_information_df$genomic_structure =="ssRNA"]<-"Unknown"
All_information_df$family[All_information_df$family == "Unknown" & All_information_df$genomic_structure =="ssRNA"]<-"Unknown"
All_information_df$best_family_per_Event[All_information_df$family == "Unknown" & All_information_df$genomic_structure =="dsRNA"]<-"Unknown"
All_information_df$best_family_per_query[All_information_df$family == "Unknown" & All_information_df$genomic_structure =="dsRNA"]<-"Unknown"
All_information_df$family[All_information_df$family == "Unknown" & All_information_df$genomic_structure =="dsRNA"]<-"Unknown"
# Add consensus_genomic_structure

All_information_df$consensus_genomic_structure<-NULL
All_information_df <-All_information_df %>%
  filter(genomic_structure != "Unknown") %>%
  group_by(Clustername,Species_name) %>%
  arrange(evalue, desc(bits)) %>%
  dplyr::slice(1) %>%
  select(Clustername, Species_name, consensus_genomic_structure = genomic_structure) %>%
  right_join(All_information_df, by = c("Clustername", "Species_name")) %>%
  relocate(consensus_genomic_structure, .after = genomic_structure)




#All_information_df <- All_information_df[!All_information_df$query %in% c("QANH01000200.1:31242-31851(-):Harpegnathos_saltator","QANH01000009.1:319477-320014(-):Harpegnathos_saltator","QANH01000368.1:67727-68264(+):Harpegnathos_saltator","NJRQ01003872.1:1811-3299(-):Aphaenogaster_ashmeadi"),]

#If HSP, keep HSP fussioned names and remove the the HSPs
#All_information_df <-All_information_df %>%
#  filter(!grepl('-HSPs', HSP_name))

All_information_df$X<-NULL


All_information_df$consensus_family<-NULL
All_information_df <-All_information_df %>%
  filter(family != "Unknown") %>%
  group_by(Clustername,Species_name) %>%
  arrange(evalue, desc(bits)) %>%
  dplyr::slice(1) %>%
  select(Clustername, Species_name, consensus_family = family) %>%
  right_join(All_information_df, by = c("Clustername", "Species_name")) %>%
  relocate(consensus_family, .after = family)


All_information_df$consensus_family[is.na(All_information_df$consensus_family)]<-"Unknown"



#Generate consensus family wich is simply for each species within clusters the hist family with the lowest evalue
All_information_df $family <-  gsub('LbFV-like family', 'LbFV_like', All_information_df$family)

#All_information_df<-All_information_df[order(All_information_df$evalue, decreasing = F),]
#All_information_df_bis<-All_information_df[!duplicated(All_information_df[c("Clustername","Species_name")]),]
#All_information_df_bis<-select(All_information_df_bis,"Clustername",'Species_name','family')

#colnames(All_information_df_bis)<- c("Clustername",'Species_name','consensus_family')

#All_information_df<-merge(All_information_df, All_information_df_bis, by=c("Clustername",'Species_name'))

All_information_df$consensus_family <-  gsub('LbFV-like family', 'LbFV_like', All_information_df $consensus_family)
All_information_df$consensus_family <-  gsub('unknown', 'Unknown', All_information_df $consensus_family)

#Remove clusters without any phylogenies (and then no monophyletic assignment)

#Monop_tab<-read.table("/Users/bguinet/Desktop/these/Monophyletic_tab_all.tab",sep=";",h=T)
#list_cluster_to_keep<-unique(Monop_tab$Clustername)

#All_information_df<-dplyr::filter(All_information_df, Clustername %in% list_cluster_to_keep)


All_information_df$cov_depth_candidat<-gsub('No_cov', '-1', All_information_df$cov_depth_candidat)
library(data.table)
All_information_df$cov_depth_candidat<-as.numeric(All_information_df$cov_depth_candidat)
setDT(All_information_df)[cov_depth_candidat < 20, query := sub("Platygaster_orseoliae", "Unknown_species", query)]
setDT(All_information_df)[cov_depth_candidat < 20, query_bis := sub("Platygaster_orseoliae", "Unknown_species", query_bis)]
All_information_df$cov_depth_candidat<-gsub('-1', 'No_cov', All_information_df$cov_depth_candidat)

All_information_df_save<-All_information_df
#Change dN/dS > 10 to 999

All_information_df$Mean_dNdS[All_information_df$Pvalue_dNdS=="No_value"] <- NA
All_information_df$SE_dNdS[All_information_df$Pvalue_dNdS=="No_value"] <- NA
All_information_df$Pvalue_dNdS[All_information_df$Pvalue_dNdS=="No_value"] <- NA

All_information_df$SE_dNdS<-as.numeric(as.character(All_information_df$SE_dNdS))
All_information_df$Mean_dNdS<-as.numeric(as.character(All_information_df$Mean_dNdS))

All_information_df$Mean_dNdS[All_information_df$Mean_dNdS > 10]<-999

#Fill manually virus famillies where NaN values are present for :

All_information_df$genomic_structure<-as.character(All_information_df$genomic_structure)
All_information_df$family<-as.character(All_information_df$family)


All_information_df_consensus<-All_information_df %>% mutate(consensus_family = coalesce(consensus_family, family_2LCA))



#write.table(All_information_df_consensus,"/Users/bguinet/Desktop/these/All_information_df_consensus.txt",sep=";",row.names=FALSE)
#All_information_df_consensus<-read.table("/Users/bguinet/Desktop/these/All_information_df_consensus.txt",sep=";",h=T)


#Change Events mistakes number for HSPs

#Reasigne consensus_families level for controls Fopius etc
All_information_df_consensus$consensus_family<-as.character(All_information_df_consensus$consensus_family)
All_information_df_consensus$family<-as.character(All_information_df_consensus$family)

#All_information_df_consensus <-All_information_df_consensus %>%
#  group_by(Clustername) %>% mutate(consensus_family= case_when(any(grepl('Burke',target) ) ~ 'Nudiviridae',TRUE ~ consensus_family ))

#All_information_df_consensus <-All_information_df_consensus %>%
#  group_by(Clustername) %>% mutate(consensus_family= case_when(any(grepl('cotesia',target) ) ~ 'Nudiviridae',TRUE ~ consensus_family ))

#All_information_df_consensus <-All_information_df_consensus %>%
#  group_by(Clustername) %>% mutate(consensus_family= case_when(any(grepl('micro',target) ) ~ 'Nudiviridae',TRUE ~ consensus_family ))

#All_information_df_consensus <-All_information_df_consensus %>%
#  group_by(Clustername) %>% mutate(consensus_family= case_when(any(grepl('Pichon',target) ) ~ 'Nudiviridae',TRUE ~ consensus_family ))

#All_information_df_consensus <-All_information_df_consensus %>% mutate(family= case_when(grepl('Burke',target)  ~ 'Nudiviridae',TRUE ~ family ))

#All_information_df_consensus <-All_information_df_consensus %>% mutate(family= case_when(grepl('cotesia',target)  ~ 'Nudiviridae',TRUE ~ family ))

#All_information_df_consensus <-All_information_df_consensus %>%mutate(family= case_when(grepl('micro',target)  ~ 'Nudiviridae',TRUE ~ family ))

All_information_df_consensus <-All_information_df_consensus %>% mutate(consensus_genomic_structure= case_when(grepl('Lispiviridae',consensus_family)  ~ 'ssRNA',TRUE ~ consensus_genomic_structure ))
All_information_df_consensus <-All_information_df_consensus %>% mutate(consensus_genomic_structure= case_when(grepl('Chuviridae',consensus_family)  ~ 'ssRNA',TRUE ~ consensus_genomic_structure ))


All_information_df_consensus<-All_information_df_consensus[order(All_information_df_consensus$evalue, decreasing = F),]

library(data.table)
Resricted_df<-All_information_df_consensus
Resricted_df <- setDT(Resricted_df)
Resricted_df_heatmap <- setDT(Resricted_df)
Resricted_df_heatmap$New_query_bis2<-Resricted_df_heatmap$query
#Delete the scaffolds names in the query names
Resricted_df_heatmap[,scaffold_name := str_extract(query,"[^:]+")]
Resricted_df_heatmap[,query := str_extract(query,"(?<=[0-9]\\([+-]\\):)[A-z ]+")]
#Resricted_df_heatmap[,query := gsub(".*__","",query)]
#Replace scaffolds with a low coverage by the Mistery species name

#All_information_df_consensus<-NULL
#All_information_df<-NULL

Resricted_df_heatmap$consensus_genomic_structure<- as.character(Resricted_df_heatmap$consensus_genomic_structure)


#Subset the dataframe per genomic structure :
List_dsDNA_genomic_structure=c("dsDNA","dsDNA-RT")
List_ssDNA_genomic_structure=c( "ssDNA","ssDNA(+/-)","ssDNA(-)","ssDNA(+)")
List_dsRNA_genomic_structure=c("dsRNA")
List_ssRNA_genomic_structure=c("ssRNA-RT" ,"ssRNA(-)" ,"ssRNA(+)" ,"ssRNA(+/-)","ssRNA")
List_Unknown_genomic_structure=c("Unknown","NA")
List_UnknownRNA_genomic_structure=c("Unknown_RNA")
#combined_list<-c(List_dsDNA_genomic_structure,List_ssDNA_genomic_structure,List_ssRNA_genomic_structure,List_ssDNA_genomic_structure)


#List_DNA_genomic_structure=c("dsDNA","dsDNA-RT","ssDNA","ssDNA(+/-)")
#List_RNA_genomic_structure=c("ssRNA-RT" ,"ssRNA(-)" ,"ssRNA(+)" ,"ssRNA(+/-)","dsRNA","dsRNA","Unknown_RNA")

#putatitive_list<-c("dsDNA","dsDNA-RT","ssDNA","ssDNA(+/-)","ssRNA-RT" ,"ssRNA(-)" ,"ssRNA(+)" ,"ssRNA(+/-)","dsRNA","Unknown_RNA","ssRNA","ssDNA","ssDNA(+/-)","ssDNA(-)","ssDNA(+)","Unknown","NA")


#putatitive_list<-c("dsDNA","dsDNA-RT","ssDNA","ssDNA(+/-)","ssDNA(-)","ssDNA(+)","dsRNA","ssRNA-RT" ,"ssRNA(-)" ,"ssRNA(+)" ,"ssRNA(+/-)","ssRNA")

#Resricted_df_heatmap<-subset(Resricted_df_heatmap, consensus_genomic_structure %in% List_dsDNA_genomic_structure)

#Resricted_df_heatmap<-subset(Resricted_df_heatmap, consensus_genomic_structure %in% c("dsRNA","ssNA"))
#Resricted_df_heatmap<-subset(Resricted_df_heatmap, !consensus_genomic_structure %in% putatitive_list)
#Resricted_df_heatmap<-subset(Resricted_df_heatmap, consensus_genomic_structure %in% List_Unknown_genomic_structure)

#Resricted_df_heatmap<-subset(Resricted_df_heatmap, (!consensus_genomic_structure %in% combined_list))



#t<-subset(Resricted_df_heatmap, consensus_genomic_structure %in% c("Unknown_RNA"))
#t<-select(Resricted_df_heatmap,'Clustername','query','species','genomic_structure','consensus_family','evalue','target','no_rank','family','order','consensus_genomic_structure','Nb_repeat','count_eucaryote','superkingdom_2LCA','family','FDR_pvalue_cov','FDR_pvalue_dNdS','Nb_repeat')


#Replace NA by zero pr present/absence of BUSCO
Resricted_df_heatmap[["Number_busco_loci"]][is.na(Resricted_df_heatmap[["Number_busco_loci"]])] <- 0
Resricted_df_heatmap[["count_eucaryote"]][is.na(Resricted_df_heatmap[["count_eucaryote"]])] <- 0


#As NA if covdepth_BUSCO < 5

Resricted_df_heatmap$superkingdom_2LCA[is.na(Resricted_df_heatmap$superkingdom_2LCA)]<-"unknown"


Resricted_df_heatmap$consensus_family[Resricted_df_heatmap$consensus_family=="Unknown"]<-NA


Resricted_df_heatmap<-Resricted_df_heatmap %>%
  mutate(query = coalesce(query,Species_name))

#Resricted_df_heatmap<-Resricted_df_heatmap[!(Resricted_df_heatmap$query=="nan"),]

length(unique(Resricted_df_heatmap$Clustername))

test_df <-Resricted_df_heatmap

test_df2<-test_df[test_df$Scaffold_score %in% c("A","B","C",'D'),]

#Addd genomic structures

test_df <-test_df %>% mutate(consensus_genomic_structure= case_when(grepl('Chuviridae',consensus_family)  ~ 'ssRNA',TRUE ~ consensus_genomic_structure ))
test_df <-test_df %>% mutate(consensus_genomic_structure= case_when(grepl('Lispiviridae',consensus_family)  ~ 'ssRNA',TRUE ~ consensus_genomic_structure ))
test_df <-test_df %>% mutate(consensus_genomic_structure= case_when(grepl('Xinmoviridae',consensus_family)  ~ 'ssRNA',TRUE ~ consensus_genomic_structure ))



#Manually removed clusters
test_df<-test_df[ ! test_df$Clustername %in% c('Cluster18904', 'Cluster24058', 'Cluster14435',
                                               'Cluster18974', 'Cluster15629', 'Cluster13425', 'Cluster12627',
                                               'Cluster25888', 'Cluster1638', 'Cluster24202', 'Cluster348',
                                               'Cluster24072', 'Cluster5044', 'Cluster23899',
                                               'Cluster26052', 'Cluster6294', 'Cluster7451', 'Cluster22889',
                                               'Cluster17424', 'Cluster3482', 'Cluster11858', 'Cluster3712',
                                               'Cluster8698', 'Cluster15850', 'Cluster27673',
                                               'Cluster15568', 'Cluster3966', 'Cluster11808', 'Cluster13439',
                                               'Cluster7069', 'Cluster4987', 'Cluster1861', 'Cluster127',
                                               'Cluster12543', 'Cluster17317', 'Cluster24467', 'Cluster26440',
                                               'Cluster9418', 'Cluster3354', 'Cluster9138', 'Cluster22922',
                                               'Cluster11771', 'Cluster16155', 'Cluster24397', 'Cluster21227',
                                               'Cluster683', 'Cluster2811', 'Cluster24070', 'Cluster3493','Cluster11916','Cluster4185','Cluster17599','Cluster12361','Cluster17920','Cluster15679', 'Cluster4023',
                                               'Cluster6802','Cluster14292','Cluster4912','Cluster25961','Cluster4383','Cluster25923','Cluster24961','Cluster95','Cluster11478','Cluster15329','Cluster13235','Cluster8661',
                                               'Cluster10736','Cluster14602','Cluster19433')]



TE_table<-read.csv("/Users/bguinet/Desktop/ALL_TE_table.txt",sep=";")
TE_table<-TE_table[TE_table$evalue<=0.0000000001,]
TE_table<-TE_table[TE_table$evalue<=0.0000000001,]

dim(TE_table[TE_table$evalue<=0.0000000001,])/dim(TE_table)
dim(TE_table[TE_table$evalue>=0.0000000001,])/dim(TE_table)

test_df$consensus_family [is.na(test_df$consensus_family)] <- "Unknown"


remove_freeliving_clusters <- test_df %>%
  group_by(Clustername) %>%
  filter(all(grepl('X|F', Scaffold_score)))

length(unique(test_df$query_bis))
test_df <- setDT(test_df )

library(dplyr)
#Deal wit TPM data
TPM_table<-tidyr::pivot_wider(test_df, names_from = Clustername, values_from = TPM_all,
                              values_fn = list(TPM_all = ~any(. > 1000)), values_fill = FALSE,id_cols = query)

test_df$FDR_pvalue_cov_median[is.na(test_df$FDR_pvalue_cov_median)] <- 'No_cov'
test_df$Nb_repeat[is.na(test_df$Nb_repeat)] <- 0


dim(Env_table[rowSums(is.na(Env_table)) != ncol(Env_table), ])


#t<- select(test_df,'Clustername','query','consensus_genomic_structure','superkingdom_2LCA','evalue','target','Nsp','pident','alnlen','family','FDR_pvalue_cov','FDR_pvalue_dNdS','Mean_dNdS','SE_dNdS','Nb_repeat','cov_depth_candidat','cov_depth_BUSCO')
#Merge the two info where uppercase mean that the scaffold contain several candidat loci
Env_table<-as.data.frame(test_df)

#Fill Nas consensus families as the family of other mates

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

Env_table<-Env_table %>%
  mutate(consensus_family = as.character(consensus_family)) %>%
  group_by(Clustername, Event) %>%
  mutate(consensus_family = replace(consensus_family, is.na(consensus_species_family)|consensus_family %in% "unknown",
                                    Mode(consensus_family[consensus_family != "unknown"])))



#In order to correct Clusters where there is no phylogenies and then no Events (we add them events numbers)
Env_table$Event<-as.numeric(Env_table$Event)

Env_table$consensus_genomic_structure[Env_table$consensus_genomic_structure=="ssDNA"& Env_table$genomic_structure=="dsDNA"] <-'dsDNA'

Env_table$consensus_family[Env_table$consensus_family=="LbFV-like family"]<-"LbFV_like"
Env_table$best_family_per_query[Env_table$best_family_per_query=="LbFV-like family"]<-"LbFV_like"
Env_table$best_family_per_Event[Env_table$best_family_per_Event=="LbFV-like family"]<-"LbFV_like"
#print some informations
#Number of total paralogs endogenized
length(unique(Env_table[Env_table$Scaffold_score %in% c("A","B","C",'D'),]$query_bis))
length(unique(Env_table[Env_table$Scaffold_score %in% c("A"),]$query_bis))
length(unique(Env_table[Env_table$Scaffold_score %in% c("B"),]$query_bis))
length(unique(Env_table[Env_table$Scaffold_score %in% c("C"),]$query_bis))
length(unique(Env_table[Env_table$Scaffold_score %in% c("D"),]$query_bis))
length(unique(Env_table[Env_table$Scaffold_score %in% c("E"),]$query_bis))
length(unique(Env_table[Env_table$Scaffold_score %in% c("F"),]$query_bis))
length(unique(Env_table[Env_table$Scaffold_score %in% c("X"),]$query_bis))

length(unique(Env_table[Env_table$Scaffold_score %in% c("A","B","C",'D') & Env_table$Scaffold_length_x > 10000  ,]$query_bis))
length(unique(Env_table[Env_table$Scaffold_score %in% c("A")& Env_table$Scaffold_length_x > 10000 ,]$query_bis))
length(unique(Env_table[Env_table$Scaffold_score %in% c("B")& Env_table$Scaffold_length_x > 10000 ,]$query_bis))
length(unique(Env_table[Env_table$Scaffold_score %in% c("C")& Env_table$Scaffold_length_x > 10000 ,]$query_bis))
length(unique(Env_table[Env_table$Scaffold_score %in% c("D")& Env_table$Scaffold_length_x > 10000 ,]$query_bis))
length(unique(Env_table[Env_table$Scaffold_score %in% c("E")& Env_table$Scaffold_length_x > 10000 ,]$query_bis))
length(unique(Env_table[Env_table$Scaffold_score %in% c("F")& Env_table$Scaffold_length_x > 10000 ,]$query_bis))
length(unique(Env_table[Env_table$Scaffold_score %in% c("X")& Env_table$Scaffold_length_x > 10000 ,]$query_bis))

length(unique(Env_table[Env_table$consensus_genomic_structure %in% c("ssRNA","dsRNA") ,]$query_bis))

length(unique(Env_table[Env_table$Scaffold_score %in% c("X","F")& Env_table$consensus_genomic_structure %in% c("ssRNA","dsRNA") ,]$query_bis))

length(unique(Env_table[Env_table$Scaffold_score %in% c("X","F") & Env_table$consensus_genomic_structure %in% c("ssRNA","dsRNA") ,]$Scaffold_name))

subEnv_table <- Env_table[!Env_table$cov_depth_BUSCO=="No_cov",]
sub_RNA_exo<-subEnv_table[subEnv_table$Scaffold_score %in% c("X","F") & subEnv_table$consensus_genomic_structure %in% c("ssRNA","dsRNA") ,]

sub_RNA_exo$cov_depth_BUSCO <- as.numeric(sub_RNA_exo$cov_depth_BUSCO)
sub_RNA_exo$cov_depth_candidat <- as.numeric(sub_RNA_exo$cov_depth_candidat)
sub_RNA_exo$cov_diff <- sub_RNA_exo$cov_depth_candidat - sub_RNA_exo$cov_depth_BUSCO

sub_RNA_exo<-sub_RNA_exo[! duplicated(sub_RNA_exo$Scaff_name),]
ggplot(sub_RNA_exo, aes(x=cov_diff)) +
  geom_histogram(color="black", fill="white") + xlim(-500,1000)

sub_DNA_exo<-subEnv_table[subEnv_table$Scaffold_score %in% c("X","F") & subEnv_table$consensus_genomic_structure %in% c("ssDNA","dsDNA") ,]

sub_DNA_exo$cov_depth_BUSCO <- as.numeric(sub_DNA_exo$cov_depth_BUSCO)
sub_DNA_exo$cov_depth_candidat <- as.numeric(sub_DNA_exo$cov_depth_candidat)
sub_DNA_exo$cov_diff <- sub_DNA_exo$cov_depth_candidat - sub_DNA_exo$cov_depth_BUSCO

sub_DNA_exo<-sub_DNA_exo[! duplicated(sub_DNA_exo$Scaff_name),]


ggplot(sub_DNA_exo, aes(x=cov_diff)) +
  geom_histogram(color="black", fill="white") + xlim(-500,1000)


length(unique(Env_table[Env_table$Scaffold_score %in% c("X","F") & Env_table$consensus_genomic_structure %in% c("ssRNA","dsRNA") & Env_table$count_repeat>=1 ,]$Scaffold_name))

length(unique(Env_table[Env_table$Scaffold_score %in% c("X","F") & Env_table$consensus_genomic_structure %in% c("ssDNA","dsDNA") ,]$Scaffold_name))

length(unique(Env_table[Env_table$Scaffold_score %in% c("X","F") & Env_table$consensus_genomic_structure %in% c("ssDNA","dsDNA") & Env_table$count_repeat>=1 ,]$Scaffold_name))


length(unique(Env_table[Env_table$Scaffold_score %in% c("A","B","C","D")& Env_table$consensus_genomic_structure %in% c("ssRNA","dsRNA") ,]$query_bis))

length(unique(Env_table[Env_table$Scaffold_score %in% c("A","B","C","D") & Env_table$consensus_genomic_structure %in% c("ssRNA","dsRNA") ,]$Scaffold_name))

length(unique(Env_table[Env_table$Scaffold_score %in% c("A","B","C","D") & Env_table$consensus_genomic_structure %in% c("ssRNA","dsRNA") & Env_table$count_repeat>=1 ,]$Scaffold_name))

length(unique(Env_table[Env_table$Scaffold_score %in% c("A","B","C","D") & Env_table$consensus_genomic_structure %in% c("ssDNA","dsDNA") ,]$Scaffold_name))

length(unique(Env_table[Env_table$Scaffold_score %in% c("A","B","C","D") & Env_table$consensus_genomic_structure %in% c("ssDNA","dsDNA") & Env_table$count_repeat>=1 ,]$Scaffold_name))

#& Env_table$Scaffold_length_x > 10000

#Remove scaffold noted as D presenting coverage lower than the BUSCO

Env_table$cov_depth_BUSCO[Env_table$cov_depth_BUSCO=="No_cov"]<- -1
Env_table$cov_depth_BUSCO <- as.numeric(Env_table$cov_depth_BUSCO)


Env_table$cov_depth_candidat[Env_table$cov_depth_candidat=="No_cov"]<- -1
Env_table$cov_depth_candidat <- as.numeric(Env_table$cov_depth_candidat)

Env_table   <-Env_table[!(Env_table$Scaffold_score=="D" & Env_table$cov_depth_BUSCO > Env_table$cov_depth_candidat),]

#Depending on condition
Env_table<-Env_table[Env_table$Scaffold_score %in% c("A","B","C","D"),]
#Env_table<-Env_table[Env_table$Scaffold_score %in% c("A"),]
#Env_table<-Env_table[Env_table$Scaffold_score %in% c("X","F","E"),]

#Count number of domesticated EVEs

write.table(Env_table,"/Users/bguinet/Desktop/Papier_scientifique/Viralprot_vs_Viral_loci_result_all_match_and_cluster_taxid_ontology_2LCA_HMMER_Augustus_Repeat_Env_pseudo_HSP_Score_Filtred_Mono_TPM_dNdS_synteny_ORFs_R_Free_and_EVEs.m8",sep=";")


#Extract sequence from HSP sequences
library("seqinr")
fastafile<- read.fasta(file = "/Users/bguinet/Desktop/Papier_scientifique/Candidate_loci_filtred.aa", seqtype = "AA",as.string = TRUE, set.attributes = FALSE)

for (i in Env_table$HSP_name[!is.na(Env_table$HSP_name)]){
  print(i)
  Env_table$sequence[Env_table$HSP_name==i]<- as.character(fastafile[c(which(names(fastafile) %in% i))][[1]])
}

for (i in Env_table$New_query_bis2[is.na(Env_table$sequence)]){
  if (grepl("HSP",i)){
    print(i)
  }else if(grepl("Unknown",i)){
    a<-i
    i<-gsub("Unknown_species","Platygaster_orseoliae",i)
    Env_table$sequence[Env_table$New_query_bis2==a]<- as.character(fastafile[c(which(names(fastafile) %in% i))][[1]])
  }else{
    Env_table$sequence[Env_table$New_query_bis2==i]<- as.character(fastafile[c(which(names(fastafile) %in% i))][[1]])
  }
}

#Correct some event
Env_table$Event[Env_table$Clustername=="Cluster15744" & Env_table$Species_name =="Synergus_japonicus"]<-3

Env_table$lifecycle1[Env_table$Species_name=="Megastigmus_dorsalis"]<-"ectoparasitoide"
Env_table$lifecycle1[Env_table$Species_name=="Megastigmus_stigmatizans"]<-"ectoparasitoide"
Env_table$lifecycle1[Env_table$Species_name=="Cecidostiba_semifascia"]<-"ectoparasitoide"
Env_table$lifecycle1[Env_table$Species_name=="Cecidostiba_fungosa"]<-"ectoparasitoide"
Env_table$lifecycle1[Env_table$Species_name=="Eretmocerus_eremicus"]<-"endoparasitoide"
Env_table$lifecycle1[Env_table$Species_name=="Ormyrus_pomaceus"]<-"ectoparasitoide"
Env_table$lifecycle1[Env_table$Species_name=="Ormyrus_nitidulus"]<-"ectoparasitoide"

Env_table<-Env_table[!(Env_table$consensus_family=="Unknown" & Env_table$genomic_structure=="Unknown" & is.na(Env_table$consensus_genomic_structure)) ,]

#We do that here and not only in python since we need to include only species that passed the filter
Env_table <-Env_table %>%
  group_by(Clustername, Event) %>%
  mutate(Events_species = sprintf('[%s]', toString(query)))






detachAllPackages()

library(tidyverse)
library(pastecs)
library(ggplot2)
library(gggenes)
library(data.table)
library(stringr)

dsDNA_tab <-Env_table #read.table("/Users/bguinet/Desktop/these/Table_all_good_scaffols.txt",sep=";",h=T)


###### TO change before submittion !!!!!!!#######

# Psyttalia_lounsbury >  Psyttalia lounsburyi
#suprfam dans fig > Trichopria species: Proctotrupoidea > Diaprioidea

dsDNA_tab$SE_dNdS[dsDNA_tab$SE_dNdS<0 & is.na(dsDNA_tab$SE_dNdS)==F]<- 0
dsDNA_tab$consensus_family <- dsDNA_tab$best_family_per_Event
dsDNA_tab$consensus_family[dsDNA_tab$consensus_family=="LbFV-like family"] <- "LbFV_like"
dsDNA_tab$consensus_genomic_structure[dsDNA_tab$consensus_family=="Unknown" & dsDNA_tab$family=="Unknown" & is.na(dsDNA_tab$genomic_structure)]<-"Unknown"
dsDNA_tab$consensus_genomic_structure[dsDNA_tab$genomic_structure=="Unknown"]<-"Unknown"
dsDNA_tab$best_family_per_query[dsDNA_tab$consensus_family=="Unknown"]<-"Unknown"

dsDNA_tab$Event[dsDNA_tab$query_bis=="scaffold_1340_7196-8042_+__Meteorus_colon_M"] <- 46
dsDNA_tab$Event[dsDNA_tab$query_bis=="NW_014327380.1_1_2-HSPs_+__Copidosoma_floridanum"] <- 27
dsDNA_tab$Event[dsDNA_tab$query_bis=="AOFN02000345.1_3_2-HSPs_+__Athalia_rosae"] <- 11
dsDNA_tab$Event[dsDNA_tab$query_bis=="NJRP01000333.1_2-HSPs_+__Aphaenogaster_floridana"]<- 4
dsDNA_tab$Event[dsDNA_tab$query_bis=="SGBU01000028.1_2-HSPs_+__Nylanderia_fulva"]<- 3
dsDNA_tab$Event[dsDNA_tab$query_bis=="scaffold_4188_5283-7005_+__Ganaspis_ganaspis"] <- 220
dsDNA_tab$Event[dsDNA_tab$query_bis=="scaffold_4188_7008-10152_+__Ganaspis_ganaspis"]<-221
dsDNA_tab$Event[dsDNA_tab$query_bis=="NJRL01003002.1_15722-17663_-__Aphaenogaster_rudis"]<-222
dsDNA_tab$Event[dsDNA_tab$query_bis=="NJRL01003002.1_12567-15621_-__Aphaenogaster_rudis"]<-223
dsDNA_tab$Event[dsDNA_tab$query_bis=="KQ435020.1_194985-196764_-__Dufourea_novaeangliae"]<-47


dsDNA_tab$Event[dsDNA_tab$New_query_bis2=="scaffold_983:16235-16619(-):Platygaster_orseoliae"]<-3


#remove duplicates within events
dsDNA_tab<-dsDNA_tab%>%
  group_by(Clustername,Event)%>%
  distinct(query, .keep_all = TRUE)

#Count number of EVEs

dsDNA_tab$ORF_perc_target <- (((dsDNA_tab$len_ORF/3)*100)/dsDNA_tab$tlen)

nrow(dsDNA_tab %>% filter(ORF_perc>1 & Mean_dNdS <= 1 & FDR_pvalue_dNdS==1 & SE_dNdS < 1 | TPM_all>1000 & Mean_dNdS <= 1 & ORF_perc>1))
nrow(dsDNA_tab %>% filter(ORF_perc>1 & Mean_dNdS <= 1 & FDR_pvalue_dNdS==1 & SE_dNdS < 1 &ORF_perc_target >= 80 | TPM_all>1000 & Mean_dNdS <= 1 & ORF_perc>1 & ORF_perc_target >=80 ))



#dsDNA_tab<-dsDNA_tab[dsDNA_tab$bits >=50,]
write.table(Env_table,"/Users/bguinet/Desktop/Papier_scientifique/Check_candidates.txt",sep=";")


names(dsDNA_tab)[names(dsDNA_tab)=="query"] <- "New_query_bis"
dsDNA_tab$consensus_genomic_structure[is.na(dsDNA_tab$consensus_genomic_structure)] <- "Unknown"


#If we want or not to remove controls

#Remove_controls<-"yes"
Remove_controls<-"no"

if (Remove_controls=="yes"){
  Control_comparaison_table<-read.table("/Users/bguinet/Desktop/Papier_scientifique/Domesticated_vs_Candidate_viral_segment_ALL.tab",sep=";",header=T)
  dsDNA_tab<-dsDNA_tab[!dsDNA_tab$New_query_bis2 %in% Control_comparaison_table$query,]
  dsDNA_tab <-dsDNA_tab%>%
    group_by(New_query_bis,consensus_family) %>%
    filter(!any(grepl('Venturia',New_query_bis2) & grepl('Pichon',target)  ))%>%
    filter(!any(grepl('Cotesia',New_query_bis2) & grepl('Bezier',target)  ))%>%
    filter(!any(grepl('Leptopilina',New_query_bis2) & grepl('Varaldi',target)  ))%>%
    filter(!any(grepl('Microplitis',New_query_bis2) & grepl('Bezier',target)  ))%>%
    filter(!any(grepl('Microplitis',New_query_bis2) & grepl('Pichon',target)  ))%>%
    filter(!any(grepl('Cotesia',New_query_bis2) & grepl('Pichon',target)  ))%>%
    filter(!any(grepl('Venturia',New_query_bis2) & grepl('Bezie',target)  ))%>%
    filter(!any(grepl('Venturia',New_query_bis2) & grepl('Burke',target)  ))%>%
    filter(!any(grepl('Fopius',New_query_bis2) & grepl('Bezier',target)  ))%>%
    filter(!any(grepl('Fopius',New_query_bis2) & grepl('Burke',target)  ))
}


t<-dsDNA_tab[!dsDNA_tab$New_query_bis2 %in% Control_comparaison_table$query,]
dsDNA_tab<-dsDNA_tab[!duplicated(dsDNA_tab[ , c("New_query_bis","Clustername","Event")]),]


#If we want or not to only count multiple EVEs
Only_mulptiple_EVE<-"yes"
Only_mulptiple_EVE<-"no"


#nrow(Env_table[Env_table$ORF_perc>=90 & !grepl("\\["",Env_table$sequence)  & Env_table$pseudogenized ==0 & Env_table$Mean_dNdS <= 1 & Env_table$FDR_pvalue_dNdS==1 & Env_table$SE_dNdS < 1 ,])


Restricted_candidat<- dsDNA_tab %>%
  group_by(Clustername)%>%
  #filter(any(Clustername %in% List_cluster_to_analyse_Platygaster_LbFV))  %>%
  #filter(stringr::str_detect(New_query_bis, species_name))%>%
  #filter(New_query_bis %in% c("Eurytoma_brunniventris","Eurytoma_adleriae",'Cotesia_vestalis','Microplitis_demolitor'))%>%
  mutate(query=New_query_bis2)%>%
  #filter(stringr::str_detect(consensus_family, 'LbFV|Hytrosa'))%>%
  #filter(stringr::str_detect(consensus_family, 'Nudi|Baculo'))%>%
  filter(!duplicated (New_query_bis2))%>%
  separate(New_query_bis2, c("Scaff_name","coordinates","query1"), ":") %>%
  separate(coordinates, c("coordinates","strand"), "[(]") %>%
  mutate(strand = gsub(")", "", strand)) %>%
  separate(coordinates, c("start","end"), "-")%>%
  #filter(any(pvalue_cov>= 0.05) & any(evalue< 0.0000005)) %>%
  select(Scaff_name,consensus_family,Scaffold_score,New_query_bis,query,Clustername,Species_name,qstart,qend,target,evalue,consensus_genomic_structure,family,Domain_description,Gene.ontology..GO.,Gene.names,Protein.names,cov_depth_candidat,cov_depth_BUSCO,FDR_pvalue_cov,count_repeat,count_eucaryote,Event,Boot,Nsp,FDR_pvalue_dNdS,Pvalue_dNdS,Mean_dNdS,SE_dNdS,TPM_all,pseudogenized)

Restricted_candidat$consensus_family[is.na(Restricted_candidat$consensus_family)] <- "Unknown"
#Change the Unknow family names by adding a number for each different target Unknown since they are probably not the same virus
setDT(Restricted_candidat)[consensus_family == "Unknown", consensus_family := paste0(consensus_family, "_", .GRP), by=target]

#Now we will produce several tables (EVEs, dEVES, EVEs_event, dEVEs_event)


look_up_group <- function(one_group, lookup_list) {
  matched_list <- map(lookup_list, function(x) { intersect(x, one_group) } )
  index <- which(unlist(map(matched_list, function(x) { length(x) > 0 })))
  sort(unique(unlist(lookup_list[index])))
}


#Change the Event number by adding the clustername to it
Restricted_candidat$Event <- paste0(Restricted_candidat$Event,"_",Restricted_candidat$Clustername)


#For events EVEs shared only
#This script allows to group together putatives events based on the fact that EVEs are in the same scaffold or by regrouping EVEs coming from the same virus family donnor

####################################################
#Try to capture number of EVEs within shared events
####################################################

table_noduplicate_shared_EVEs <- Restricted_candidat%>%
  group_by(Clustername,Event)%>%
  distinct(New_query_bis, .keep_all = TRUE)%>%
  mutate(nrows=n())%>%
  filter(any(nrows >=2))

table_noduplicate_shared_EVEs<-table_noduplicate_shared_EVEs[!duplicated(table_noduplicate_shared_EVEs[ , c("New_query_bis","Clustername","Event")]),]
table_noduplicate_shared_EVEs <-select(as.data.frame(table_noduplicate_shared_EVEs),'query','Clustername','New_query_bis','consensus_family','consensus_genomic_structure','Species_name','Event','Scaff_name')

#Deal with the fact that duplicated EVEs can be into different Event, the idea is just to check the boostrap value that separates two putative duplicates (g.e EVEs within same cluster and same specie)
table_noduplicate_shared_EVEs<-table_noduplicate_shared_EVEs[!duplicated(table_noduplicate_shared_EVEs[ , c("New_query_bis","Clustername","Event")]),]

#Sort two columns
table_noduplicate_shared_EVEs <- table_noduplicate_shared_EVEs[order(table_noduplicate_shared_EVEs$New_query_bis),]
table_noduplicate_shared_EVEs<-table_noduplicate_shared_EVEs[!duplicated(list_df[ , c("Clustername","Event")]),]

library(dplyr)
library(tidyr)

list_all_query1 <- table_noduplicate_shared_EVEs$query

table_noduplicate_shared_EVEs  <- table_noduplicate_shared_EVEs %>%
  group_by(Clustername, Event)  %>%
  summarize(across(c(New_query_bis, query, consensus_family,consensus_genomic_structure,Species_name,Scaff_name), paste0, collapse = ","), .groups = "drop") %>%
  filter(!duplicated( cbind(New_query_bis,Clustername)))  %>%
  separate_rows(New_query_bis, query, consensus_family,consensus_genomic_structure,Species_name,Scaff_name, sep = ",", convert = TRUE)

list_all_query2<- table_noduplicate_shared_EVEs$query

list_query_remove_because_redundancy<-list_all_query1[!(list_all_query1 %in% list_all_query2)]



list_df<-expand.grid(table_noduplicate_shared_EVEs$Clustername,table_noduplicate_shared_EVEs$Species_name)
colnames(list_df)<-c("Clustername","Species_name")
list_df<-list_df[!duplicated(list_df[ , c("Clustername","Species_name")]),]
list_df$present<-"yes"

list_df<-merge(select(table_noduplicate_shared_EVEs,"Clustername","Species_name"),list_df,by=c("Clustername","Species_name"),all.x = TRUE)

table_noduplicate_shared_EVEs$BootValue<-"NA"
library(ape)
library(phytools)
library(phylobase)
list_df<-list_df[!duplicated(list_df[ , c("Clustername","Species_name")]),]


#Correct Bootstrapvalues
for(i in unique(list_df$Clustername)) {
  sub_table_noduplicate_shared_EVEs<-table_noduplicate_shared_EVEs[with(table_noduplicate_shared_EVEs, Clustername == i ), ]
  sub_table_noduplicate_shared_EVEs<-sub_table_noduplicate_shared_EVEs[!is.na(sub_table_noduplicate_shared_EVEs$Clustername), ]
  sub_table_noduplicate_shared_EVEs$query<-gsub(":","_",sub_table_noduplicate_shared_EVEs$query)
  sub_table_noduplicate_shared_EVEs$query<-gsub("\\(","_",sub_table_noduplicate_shared_EVEs$query)
  sub_table_noduplicate_shared_EVEs$query<-gsub("\\)","_",sub_table_noduplicate_shared_EVEs$query)
  sub_table_noduplicate_shared_EVEs$query<-gsub(" ","",sub_table_noduplicate_shared_EVEs$query)
  for (a in unique(sub_table_noduplicate_shared_EVEs$Event)){
    sub_sub_table_noduplicate_shared_EVEs<-sub_table_noduplicate_shared_EVEs[sub_table_noduplicate_shared_EVEs$Event == a, ]
    if (length(dim(sub_sub_table_noduplicate_shared_EVEs)) >=1 ){
      Event_number<-a
      if (file.exists(paste0("/Users/bguinet/Desktop/these/Cluster_phylogeny_filtred/",unique(sub_table_noduplicate_shared_EVEs$Clustername),"_AA.dna.treefile"))) {
        print(i)
        PhyloTree=read.tree(paste0("/Users/bguinet/Desktop/these/Cluster_phylogeny_filtred/",unique(sub_table_noduplicate_shared_EVEs$Clustername),"_AA.dna.treefile"))
        PhyloTree_labels<-PhyloTree$tip.label
        tax_to_remove<-PhyloTree_labels[grepl("HSP",PhyloTree_labels)]
        tax_to_remove<-gsub("_[0-9]-HSPs.*", "",tax_to_remove)
        list_candidate<-sub_sub_table_noduplicate_shared_EVEs$query
        tryCatch({
          PhyloTree <- as(PhyloTree, "phylo4")
        }, error=function(e){})
        tryCatch({
          mm<- MRCA(PhyloTree,list_candidate)
        }, error=function(e){})
        bootstrap_value<-names(mm)
        print(list_candidate)
        print(bootstrap_value)
        table_noduplicate_shared_EVEs$BootValue[table_noduplicate_shared_EVEs$Event==Event_number] <- bootstrap_value
        if (length(PhyloTree_labels) %in% c(3,4)){
          table_noduplicate_shared_EVEs$BootValue[table_noduplicate_shared_EVEs$Event==Event_number] <- -2
        }
      }
    }
  }
}


for(i in 1:nrow(list_df)) {
  sub_table_noduplicate_shared_EVEs<-table_noduplicate_shared_EVEs[with(table_noduplicate_shared_EVEs, Clustername == list_df$Clustername[[i]] & Species_name == list_df$Species_name[[i]]), ]
  if (length(dim(sub_table_noduplicate_shared_EVEs)) >=1 ){
    sub_table_noduplicate_shared_EVEs<-sub_table_noduplicate_shared_EVEs[!is.na(sub_table_noduplicate_shared_EVEs$Clustername), ]
    sub_table_noduplicate_shared_EVEs$query<-gsub(":","_",sub_table_noduplicate_shared_EVEs$query)
    sub_table_noduplicate_shared_EVEs$query<-gsub("\\(","_",sub_table_noduplicate_shared_EVEs$query)
    sub_table_noduplicate_shared_EVEs$query<-gsub("\\)","_",sub_table_noduplicate_shared_EVEs$query)
    sub_table_noduplicate_shared_EVEs$query<-gsub(" ","",sub_table_noduplicate_shared_EVEs$query)
    if (file.exists(paste0("/Users/bguinet/Desktop/these/Cluster_phylogeny_filtred/",unique(sub_table_noduplicate_shared_EVEs$Clustername),"_AA.dna.treefile"))) {
      print(i)
      PhyloTree=read.tree(paste0("/Users/bguinet/Desktop/these/Cluster_phylogeny_filtred/",unique(sub_table_noduplicate_shared_EVEs$Clustername),"_AA.dna.treefile"))
      PhyloTree_labels<-PhyloTree$tip.label
      tax_to_remove<-PhyloTree_labels[grepl("HSP",PhyloTree_labels)]
      tax_to_remove<-gsub("_[0-9]-HSPs.*", "",tax_to_remove)
      list_candidate<-sub_table_noduplicate_shared_EVEs$query
      if (length(tax_to_remove )> 1) {
        list_candidate<-list_candidate[!grepl(paste(tax_to_remove, collapse = "|"), list_candidate)]
      }
      if (length(list_candidate)>1){
        PhyloTree<-midpoint.root(PhyloTree)
        tryCatch({
          PhyloTree <- as(PhyloTree, "phylo4")
        }, error=function(e){})
        mm<- MRCA(PhyloTree,list_candidate)
        bootstrap_value<-names(mm)
        #print(list_candidate)
        print(bootstrap_value)
        print(list_df$Species_name[[i]])
        print(list_df$Clustername[[i]])
        table_noduplicate_shared_EVEs <-table_noduplicate_shared_EVEs%>%
          group_by(Clustername,Species_name)%>%
          mutate(BootValue = case_when(any(grepl(list_df$Clustername[[i]],Clustername) & Species_name %in% list_df$Species_name[[i]]  ) ~ bootstrap_value,TRUE ~ BootValue))
        print("Cluster7293" %in% table_noduplicate_shared_EVEs$Clustername)
      }
    }
  }
}

#Remove duplicated EVEs when the bootstrap support is below 80
table_noduplicate_shared_EVEs$BootValue[table_noduplicate_shared_EVEs$BootValue=="Root"] <-0
table_noduplicate_shared_EVEs$BootValue[table_noduplicate_shared_EVEs$BootValue=="NA"] <- -1
table_noduplicate_shared_EVEs$BootValue <- as.numeric(table_noduplicate_shared_EVEs$BootValue)
#table_noduplicate_shared_EVEs  <-table_noduplicate_shared_EVEs[duplicated(table_noduplicate_shared_EVEs) | table_noduplicate_shared_EVEs$BootValue < 80,]

Shared_events_splitted <- table_noduplicate_shared_EVEs[table_noduplicate_shared_EVEs$BootValue < 80 & table_noduplicate_shared_EVEs$BootValue > 0,]
#Chang the Event number for those one in order to add them within the alone Event part
for (event in unique(Shared_events_splitted$Event)){
  sub_event_nb <- 1
  subShared_events_splitted <- Shared_events_splitted[Shared_events_splitted$Event==event,]
  for (i in unique(subShared_events_splitted$query)){
    Restricted_candidat$Event[Restricted_candidat$query==i] <- paste0(event,"_bis_",sub_event_nb)
    sub_event_nb= sub_event_nb +1
  }
}


table_noduplicate_shared_EVEs  <- table_noduplicate_shared_EVEs[table_noduplicate_shared_EVEs$BootValue > 80 | table_noduplicate_shared_EVEs$BootValue  == -2 ,]


#Add lifecyle
anoleData<-read.csv("/Users/bguinet/Documents/All_informations_table_fam_busco_life_style.csv",sep=";",header=T)
anoleData$lifecycle1[anoleData$Species_name=="Megastigmus_dorsalis"]<-"ectoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Megastigmus_stigmatizans"]<-"ectoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Cecidostiba_semifascia"]<-"ectoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Cecidostiba_fungosa"]<-"ectoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Eretmocerus_eremicus"]<-"endoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Ormyrus_pomaceus"]<-"ectoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Ormyrus_nitidulus"]<-"ectoparasitoide"
anoleData<-select(anoleData,'Species_name','lifecycle1')
colnames(anoleData)<-c('Species_name','lifecycle1')
table_noduplicate_shared_EVEs <-merge(table_noduplicate_shared_EVEs ,anoleData,by=c("Species_name"))



#Now we will map the EVEs event within the phylogeny, bu for that we need to create a dataframe where we add each Event number for each node number
col.names = c("Node_number", "Event",'Cluster','consensus_family','Clustername','species','Nb_EVEs','Nb_EVEs_dsDNA','Nb_EVEs_ssDNA','Nb_EVEs_ssRNA','Nb_EVEs_dsRNA','Nb_EVEs_Unclassified','Nb_freeliving_EVEs','Nb_ecto_EVEs','Nb_endo_EVEs','lifecycle1')
EVEs_shared_event_noduplicate_df<- read.table(text = "",
                                              col.names = col.names)

EVEs_shared_event_noduplicate_df$lifecycle1<- as.character(EVEs_shared_event_noduplicate_df$lifecycle1)
EVEs_shared_event_noduplicate_df$consensus_family<- as.character(EVEs_shared_event_noduplicate_df$consensus_family)
EVEs_shared_event_noduplicate_df$Clustername<- as.character(EVEs_shared_event_noduplicate_df$Clustername)
EVEs_shared_event_noduplicate_df$species<- as.character(EVEs_shared_event_noduplicate_df$species)
EVEs_shared_event_noduplicate_df$Event<- as.character(EVEs_shared_event_noduplicate_df$Event)

detachAllPackages()

library(tidyverse)
library(pastecs)
library(ggplot2)
library(gggenes)
#library(gsubfn)
library(data.table)
library(stringr)



table_noduplicate_shared_EVEs <-table_noduplicate_shared_EVEs %>%
  group_by(Clustername, Event) %>%
  mutate(Events_species = sprintf('[%s]', toString(Species_name)))

table_noduplicate_shared_EVEs$Events_species2<-NA
for(i in 1:nrow(table_noduplicate_shared_EVEs)) {
  row <- table_noduplicate_shared_EVEs$Events_species[[i]]
  row<-gsub("\\[","",row)
  row<-gsub("\\]","",row)
  row<-unlist(strsplit(as.character(row),","))
  #row<-list(row)
  row<-gsub(" ","",row)
  row<-unique(row)
  row=row<-sort(row)
  #print(row)
  table_noduplicate_shared_EVEs$Events_species2[[i]]<-paste(row, collapse=',')
}



table_shared_EVEs_event<-table_noduplicate_shared_EVEs %>%
  group_by(Clustername,Event)%>%
  distinct(New_query_bis, .keep_all = TRUE)%>%
  mutate(nrows=n())%>%
  filter(any(nrows >=2))%>%
  arrange(Scaff_name,New_query_bis,consensus_family,consensus_genomic_structure,Events_species2) %>%
  # create a Group_2 which is combination of all Group for each family
  group_by(consensus_family,New_query_bis,Events_species2) %>%
  mutate(Group_2 = list(Scaff_name)) %>%
  ungroup() %>%
  # Create Group_3 which is the full combined Group for all intersect Group
  mutate(Group_3 = map(.[["Group_2"]], function(x) { look_up_group(one_group = x, lookup_list = .[["Group_2"]]) })) %>%
  # Combine all Group_3 into a Group_final
  mutate(Group_final = unlist(map(Group_3, function(x) { paste (x, collapse = ",")} ))) %>%
  mutate(Species =New_query_bis )%>%
  # Finally put them all together.
  select(Species,New_query_bis,Group_final, consensus_family,Event,Clustername,consensus_genomic_structure,Events_species2) %>%
  group_by(Group_final) %>%
  summarize(family = paste(consensus_family, collapse = ","),species=paste(unique(Species), collapse = ","),Events_species2=paste(unique(Events_species2), collapse = ","),Events=paste(unique(Event), collapse = ","),Clusters=paste(unique(Clustername), collapse = ","),consensus_genomic_structure=paste(unique(consensus_genomic_structure), collapse = ","), .groups = "drop")


setDT(table_shared_EVEs_event)

#Here we take into account the other species within an event, if for instance two virus families are within the same scaffold in SP1, then I also merge them for SP2.
col.names = c("family","species","Events","Clusters","consensus_genomic_structure")
table_shared_EVEs_event2<- read.table(text = "",
                                      col.names = col.names)

table_shared_EVEs_event2$family<-as.character(table_shared_EVEs_event2$family)
table_shared_EVEs_event2$species<-as.character(table_shared_EVEs_event2$species)
table_shared_EVEs_event2$Events<-as.character(table_shared_EVEs_event2$Events)
table_shared_EVEs_event2$Clusters<-as.character(table_shared_EVEs_event2$Clusters)
table_shared_EVEs_event2$Events_species2<-as.character(table_shared_EVEs_event2$Events_species2)
table_shared_EVEs_event2$consensus_genomic_structure<-as.character(table_shared_EVEs_event2$consensus_genomic_structure)
#Subset the dataframe per genomic structure :
List_dsDNA_genomic_structure=c("dsDNA","dsDNA-RT")
List_ssDNA_genomic_structure=c( "ssDNA","ssDNA(+/-)","ssDNA(-)","ssDNA(+)")
List_dsRNA_genomic_structure=c("dsRNA")
List_ssRNA_genomic_structure=c("ssRNA-RT" ,"ssRNA(-)" ,"ssRNA(+)" ,"ssRNA(+/-)","ssRNA")
List_Unknown_genomic_structure=c("Unknown","NA")
List_UnknownRNA_genomic_structure=c("Unknown_RNA")

all_lists<-list(List_dsDNA_genomic_structure,List_ssDNA_genomic_structure,List_dsRNA_genomic_structure,List_ssRNA_genomic_structure,List_Unknown_genomic_structure,List_UnknownRNA_genomic_structure)
#The idea here is to look if we wan gather in the same event within species  EVEs coming from the same event
for(i in 1:length(all_lists)) {
  table_shared_EVEs_event_sub<-subset(table_shared_EVEs_event, consensus_genomic_structure %in% all_lists[[i]])
  #Resricted_df_heatmap<-subset(Resricted_df_heatmap, !consensus_genomic_structure %in% putatitive_list)
  #Resrict
  if(dim(table_shared_EVEs_event_sub)[1] >=1){
    print(table_shared_EVEs_event_sub)
    setDT(table_shared_EVEs_event_sub)
    table_shared_EVEs_event_sub[, family := as.character(family)]
    table_shared_EVEs_event_sub[, Events := as.character(Events)]
    table_shared_EVEs_event_sub[, Clusters := as.character(Clusters)]
    table_shared_EVEs_event_sub[, Events_species2 := as.character(Events_species2)]
    table_shared_EVEs_event_sub[, consensus_genomic_structure:= as.character(consensus_genomic_structure)]
    table_shared_EVEs_event_sub[, n := row.names(.SD) ]

    dt3 <- merge(table_shared_EVEs_event_sub,
                 table_shared_EVEs_event_sub,
                 by = "species",allow.cartesian=TRUE)[n.x >= n.y]

    dt3[, testfamily := length(intersect(unlist(strsplit(family.x, ",")), unlist(strsplit(family.y, ",")))) > 0, by = 1:nrow(dt3)]
    dt3[, testEventss := length(intersect(unlist(strsplit(Events.x, ",")), unlist(strsplit(Events.y, ",")))) > 0, by = 1:nrow(dt3)]
    dt3[, testEvents_species2 := length(intersect(unlist(strsplit(Events_species2.x, ",")), unlist(strsplit(Events_species2.y, ",")))) > 0, by = 1:nrow(dt3)]
    dt3[, testconsensus_genomic_structure:= length(intersect(unlist(strsplit(consensus_genomic_structure.x, ",")), unlist(strsplit(consensus_genomic_structure.y, ",")))) > 0, by = 1:nrow(dt3)]
    dt3[, testClusters := length(intersect(unlist(strsplit(Clusters.x, ",")), unlist(strsplit(Clusters.y, ",")))) > 0, by = 1:nrow(dt3)]

    dt3[, Events_species2 := paste0(unique(unlist(strsplit(paste(Events_species2.x, Events_species2.y, sep = ","), ","))), collapse = ","), by = 1:nrow(dt3)]
    dt3[, family := paste0(unique(unlist(strsplit(paste(family.x, family.y, sep = ","), ","))), collapse = ","), by = 1:nrow(dt3)]
    dt3[, Events := paste0(unique(unlist(strsplit(paste(Events.x, Events.y, sep = ","), ","))), collapse = ","), by = 1:nrow(dt3)]
    dt3[, consensus_genomic_structure:= paste0(unique(unlist(strsplit(paste(consensus_genomic_structure.x, consensus_genomic_structure.y, sep = ","), ","))), collapse = ","), by = 1:nrow(dt3)]
    dt3[, Clusters := paste0(unique(unlist(strsplit(paste(Clusters.x, Clusters.y, sep = ","), ","))), collapse = ","), by = 1:nrow(dt3)]

    dt3<-dt3[which(dt3$testClusters==TRUE & dt3$testfamily == TRUE & dt3$testClusters == TRUE & dt3$testEventss == TRUE  & dt3$testEvents_species2 & dt3$testconsensus_genomic_structure == TRUE),]

    dtx<-dt3[(testfamily == TRUE | testClusters == TRUE | testEventss == TRUE  |  testEvents_species2 == TRUE | testconsensus_genomic_structure ==TRUE)
             & !n.x %in% dt3[(testfamily == TRUE | testClusters== TRUE | testEventss == TRUE | testEvents_species2 == TRUE | testconsensus_genomic_structure ==TRUE ) & n.x != n.y, n.y]
             & !n.y %in% dt3[(testfamily == TRUE | testClusters == TRUE | testEventss == TRUE | testEvents_species2 ==TRUE| testconsensus_genomic_structure==TRUE) & n.x != n.y, n.x],
             .(species, family, Events, Clusters,consensus_genomic_structure,Events_species2)]

    #table_shared_EVEs_event2<-rbind(table_shared_EVEs_event2,dtx)
    table_shared_EVEs_event2 <-
      full_join(table_shared_EVEs_event2, dtx)
  }
}


#Add lifecyle
anoleData<-read.csv("/Users/bguinet/Documents/All_informations_table_fam_busco_life_style.csv",sep=";",header=T)
anoleData$lifecycle1[anoleData$Species_name=="Megastigmus_dorsalis"]<-"ectoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Megastigmus_stigmatizans"]<-"ectoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Cecidostiba_semifascia"]<-"ectoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Cecidostiba_fungosa"]<-"ectoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Eretmocerus_eremicus"]<-"endoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Ormyrus_pomaceus"]<-"ectoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Ormyrus_nitidulus"]<-"ectoparasitoide"
anoleData<-select(anoleData,'Species_name','lifecycle1')
colnames(anoleData)<-c('species','lifecycle1')
table_shared_EVEs_event2[table_shared_EVEs_event2$species=="Trichopria_sp_",]$species <- "Trichopria_sp_970989"
table_shared_EVEs_event2[table_shared_EVEs_event2$Species=="Trichopria_sp_",]$Species <- "Trichopria_sp_970989"
table_shared_EVEs_event2<-merge(table_shared_EVEs_event2,anoleData,by=c("species"))

#Then the table_shared_EVEs_event2 stores the events candidates for each species

#Now will gather all species from the same event
table_shared_EVEs_event2  <- table_shared_EVEs_event2 %>%
  mutate(across(c(family, Events, Clusters,consensus_genomic_structure,lifecycle1,Events_species2), ~ strsplit(as.character(.), split = ',')))

#Add MRCA number


library(phytools)
table_shared_EVEs_event2$MRCA<-"NA"
for(i in 1:nrow(table_shared_EVEs_event2)) {
  row <- table_shared_EVEs_event2$Events_species2[[i]]
  #row<-gsub("\\[","",row)
  #row<-gsub("\\]","",row)
  #row<-list(row)
  #row<-gsub(" ","",row)
  row<-unlist(strsplit(as.character(row),","))
  if(length(unique(row)) >1){
    #print(row)
    MRCA_numb<-findMRCA(tree,row)
    if(length(MRCA_numb)) {
      #print(MRCA_numb)
      table_shared_EVEs_event2$MRCA[i]<-MRCA_numb
      table_shared_EVEs_event2$Events_species2[[i]]<-paste(row, collapse=',')
    }
  }
}

table_shared_EVEs_event2 <-setDT(table_shared_EVEs_event2)
collapse_rows <- function(i) {
  rows_collapse <- pmap_lgl(table_shared_EVEs_event2, function(family, Events, Clusters,consensus_genomic_structure,MRCA,Events_species2, ...)
    any(table_shared_EVEs_event2$family[[i]] %in% family) & any(table_shared_EVEs_event2$Events[[i]] %in% Events) & any(table_shared_EVEs_event2$Clusters[[i]] %in% Clusters)& any(table_shared_EVEs_event2$consensus_genomic_structure[[i]] %in% consensus_genomic_structure) & any(table_shared_EVEs_event2$MRCA[[i]] %in% MRCA & any(table_shared_EVEs_event2$Events_species2 %in% table_shared_EVEs_event2$Events_species2)))
  table_shared_EVEs_event2 %>%
    filter(rows_collapse) %>%
    mutate(across(everything(), ~ paste(sort(unique(unlist(.))), collapse = ',')))
}


table_shared_EVEs_event<- map_dfr(1:nrow(table_shared_EVEs_event2), collapse_rows) %>% distinct



table_shared_EVEs_event$MRCA<-NULL

colnames(table_shared_EVEs_event)<-c("species","consensus_family","Events","Clusters","consensus_genomic_structure","Events_species2","lifecycle1")


#Now we will map the EVEs event within the phylogeny, bu for that we need to create a dataframe where we add each Event number for each node number
col.names = c("Node_number", "Event",'consensus_genomic_structure','Cluster','consensus_family','Clustername','species','Nb_EVEs','Nb_EVEs_dsDNA','Nb_EVEs_ssDNA','Nb_EVEs_ssRNA','Nb_EVEs_dsRNA','Nb_EVEs_Unclassified','Nb_EVEs_dsDNA_Events','Nb_EVEs_ssDNA_Events','Nb_EVEs_ssRNA_Events','Nb_EVEs_dsRNA_Events','Nb_EVEs_Unclassified_Events','Nb_freeliving_EVEs','Nb_ecto_EVEs','Nb_endo_EVEs','Nb_freeliving_EVEs_Events','Nb_ecto_EVEs_Events','Nb_endo_EVEs_Events','lifecycle1')
EVEs_shared_event_df<- read.table(text = "",
                                  col.names = col.names)

EVEs_shared_event_df$lifecycle1<- as.character(EVEs_shared_event_df$lifecycle1)
EVEs_shared_event_df$consensus_family<- as.character(EVEs_shared_event_df$consensus_family)
EVEs_shared_event_df$Clustername<- as.character(EVEs_shared_event_df$Clustername)
EVEs_shared_event_df$species<- as.character(EVEs_shared_event_df$species)
EVEs_shared_event_df$Event<- as.character(EVEs_shared_event_df$Event)
EVEs_shared_event_df$consensus_genomic_structure<- as.character(EVEs_shared_event_df$consensus_genomic_structure)
EVEs_shared_event_df$Events_species2<- as.character(EVEs_shared_event_df$Events_species2)


library(tidyverse)
library(phytools)
for(i in 1:nrow(table_shared_EVEs_event)) {
  Nb_EVEs<-length(unlist(strsplit(as.character(table_shared_EVEs_event$Clusters[[i]]),",")))
  if(grepl("free",table_shared_EVEs_event$lifecycle1[[i]])){
    Nb_freeliving_EVEs<-length(unlist(strsplit(as.character(table_shared_EVEs_event$Clusters[[i]]),",")))
    Nb_ecto_EVEs<-0
    Nb_endo_EVEs<-0
    Nb_freeliving_EVEs_Events<-1
    Nb_ecto_EVEs_Events<-0
    Nb_endo_EVEs_Events<-0
  }
  if(grepl("endo",table_shared_EVEs_event$lifecycle1[[i]])){
    Nb_freeliving_EVEs<-0
    Nb_ecto_EVEs<-0
    Nb_endo_EVEs<-length(unlist(strsplit(as.character(table_shared_EVEs_event$Clusters[[i]]),",")))
    Nb_freeliving_EVEs_Events<-0
    Nb_ecto_EVEs_Events<-0
    Nb_endo_EVEs_Events<-1
  }
  if(grepl("ecto",table_shared_EVEs_event$lifecycle1[[i]])){
    Nb_freeliving_EVEs<-0
    Nb_ecto_EVEs<-length(unlist(strsplit(as.character(table_shared_EVEs_event$Clusters[[i]]),",")))
    Nb_endo_EVEs<-0
    Nb_freeliving_EVEs_Events<-0
    Nb_ecto_EVEs_Events<-1
    Nb_endo_EVEs_Events<-0
  }
  if(grepl("dsDNA",table_shared_EVEs_event$consensus_genomic_structure[[i]])){
    Nb_EVEs_dsDNA<-length(unlist(strsplit(as.character(table_shared_EVEs_event$Clusters[[i]]),",")))
    Nb_EVEs_ssDNA<-0
    Nb_EVEs_ssRNA<-0
    Nb_EVEs_dsRNA<-0
    Nb_EVEs_Unclassified<-0
    Nb_EVEs_dsDNA_Events<-1
    Nb_EVEs_ssDNA_Events<-0
    Nb_EVEs_ssRNA_Events<-0
    Nb_EVEs_dsRNA_Events<-0
    Nb_EVEs_Unclassified_Events<-0
  }
  if(grepl("ssDNA",table_shared_EVEs_event$consensus_genomic_structure[[i]])){
    Nb_EVEs_ssDNA<-length(unlist(strsplit(as.character(table_shared_EVEs_event$Clusters[[i]]),",")))
    Nb_EVEs_dsDNA<-0
    Nb_EVEs_ssRNA<-0
    Nb_EVEs_dsRNA<-0
    Nb_EVEs_Unclassified<-0
    Nb_EVEs_dsDNA_Events<-0
    Nb_EVEs_ssDNA_Events<-1
    Nb_EVEs_ssRNA_Events<-0
    Nb_EVEs_dsRNA_Events<-0
    Nb_EVEs_Unclassified_Events<-0
  }
  if(grepl("ssRNA",table_shared_EVEs_event$consensus_genomic_structure[[i]])){
    Nb_EVEs_ssRNA<-length(unlist(strsplit(as.character(table_shared_EVEs_event$Clusters[[i]]),",")))
    Nb_EVEs_ssDNA<-0
    Nb_EVEs_dsDNA<-0
    Nb_EVEs_dsRNA<-0
    Nb_EVEs_Unclassified<-0
    Nb_EVEs_dsDNA_Events<-0
    Nb_EVEs_ssDNA_Events<-0
    Nb_EVEs_ssRNA_Events<-1
    Nb_EVEs_dsRNA_Events<-0
    Nb_EVEs_Unclassified_Events<-0
  }
  if(grepl("dsRNA",table_shared_EVEs_event$consensus_genomic_structure[[i]])){
    Nb_EVEs_dsRNA<-length(unlist(strsplit(as.character(table_shared_EVEs_event$Clusters[[i]]),",")))
    Nb_EVEs_ssDNA<-0
    Nb_EVEs_ssRNA<-0
    Nb_EVEs_dsDNA<-0
    Nb_EVEs_Unclassified<-0
    Nb_EVEs_dsDNA_Events<-0
    Nb_EVEs_ssDNA_Events<-0
    Nb_EVEs_ssRNA_Events<-0
    Nb_EVEs_dsRNA_Events<-1
    Nb_EVEs_Unclassified_Events<-0
  }
  if(grepl("Unknown",table_shared_EVEs_event$consensus_genomic_structure[[i]])){
    Nb_EVEs_Unclassified<-length(unlist(strsplit(as.character(table_shared_EVEs_event$Clusters[[i]]),",")))
    Nb_EVEs_ssDNA<-0
    Nb_EVEs_ssRNA<-0
    Nb_EVEs_dsRNA<-0
    Nb_EVEs_dsDNA<-0
    Nb_EVEs_dsDNA_Events<-0
    Nb_EVEs_ssDNA_Events<-0
    Nb_EVEs_ssRNA_Events<-0
    Nb_EVEs_dsRNA_Events<-0
    Nb_EVEs_Unclassified_Events<-1
  }
  lifecycle1<-table_shared_EVEs_event$lifecycle1[[i]]
  Cluster<-table_shared_EVEs_event$Clusters[[i]]
  row <- table_shared_EVEs_event$Events_species2[[i]]
  family<-table_shared_EVEs_event$consensus_family[[i]]
  #row<-gsub("\\[","",row)
  #row<-gsub("\\]","",row)
  #row<-list(row)
  #row<-gsub(" ","",row)
  row<-unlist(strsplit(as.character(row),","))
  row2<-as.character(unique(row))
  row2<-as.character(row2)
  #print(row)
  if(length(unique(row)) >1){
    print(row2)
    print(Cluster)
    MRCA_numb<-findMRCA(tree,row)
    print(MRCA_numb)
    if(length(MRCA_numb)) {
      print(Cluster)
      EVEs_shared_event_df<-EVEs_shared_event_df %>% add_row(Node_number = MRCA_numb, consensus_genomic_structure= table_shared_EVEs_event$consensus_genomic_structure[[i]],Events_species2=table_shared_EVEs_event$Events_species2[[i]],Event =table_shared_EVEs_event$Events[[i]],consensus_family=as.character(family),Clustername=Cluster,species=as.character(row2),Nb_EVEs=Nb_EVEs,Nb_EVEs_dsDNA=Nb_EVEs_dsDNA,Nb_EVEs_ssDNA=Nb_EVEs_ssDNA,Nb_EVEs_ssRNA=Nb_EVEs_ssRNA,Nb_EVEs_dsRNA=Nb_EVEs_dsRNA,Nb_EVEs_Unclassified=Nb_EVEs_Unclassified,Nb_EVEs_dsDNA_Events=Nb_EVEs_dsDNA_Events,Nb_EVEs_ssDNA_Events=Nb_EVEs_ssDNA_Events,Nb_EVEs_ssRNA_Events=Nb_EVEs_ssRNA_Events,Nb_EVEs_dsRNA_Events=Nb_EVEs_dsRNA_Events,Nb_EVEs_Unclassified_Events=Nb_EVEs_Unclassified_Events,Nb_freeliving_EVEs=Nb_freeliving_EVEs,Nb_ecto_EVEs=Nb_ecto_EVEs,Nb_endo_EVEs=Nb_endo_EVEs,Nb_freeliving_EVEs_Events=Nb_freeliving_EVEs_Events,Nb_ecto_EVEs_Events=Nb_ecto_EVEs_Events,Nb_endo_EVEs_Events=Nb_endo_EVEs_Events,lifecycle1=lifecycle1)
      print(MRCA_numb)
    }
  }
  # do stuff with row
}
#########

for(i in 1:nrow(EVEs_shared_event_df)) {
  row <- EVEs_shared_event_df$Events_species2[[i]]
  row<-unlist(strsplit(as.character(row),","))
  row<-unique(row)
  row=row<-sort(row)
  #print(row)
  EVEs_shared_event_df$Events_species2[[i]]<-paste(row, collapse=',')
}



EVEs_shared_event_df<-EVEs_shared_event_df[!duplicated(EVEs_shared_event_df[c("Node_number","Event")]),]


#Count EVEs
#Add Number EVEs
EVEs_shared_event_df$Nb_EVEs <- str_count(EVEs_shared_event_df$Event, ",")+1

#EVENTS according to lifestyles
EVEs_shared_event_df<-EVEs_shared_event_df %>% mutate (
  Nb_endo_EVEs_Events = ifelse(Nb_EVEs >= 1 & lifecycle1 %in% c('endoparasitoide'),1,0 ))
EVEs_shared_event_df<-EVEs_shared_event_df %>% mutate (
  Nb_ecto_EVEs_Events = ifelse(Nb_EVEs >= 1 & lifecycle1 %in% c('ectoparasitoide'),1,0 ))
EVEs_shared_event_df<-EVEs_shared_event_df %>% mutate (
  Nb_freeliving_EVEs_Events = ifelse(Nb_EVEs >= 1 & lifecycle1 %in% c('freeliving'),1,0 ))

#NB EVES according to lifestyles
EVEs_shared_event_df<-EVEs_shared_event_df %>% mutate (
  Nb_ecto_EVEs = ifelse(Nb_EVEs >= 1 & lifecycle1 %in% c('ectoparasitoide'),Nb_EVEs, 0 ))

EVEs_shared_event_df<-EVEs_shared_event_df %>% mutate (
  Nb_freeliving_EVEs = ifelse(Nb_EVEs >= 1 & lifecycle1 %in% c('freeliving'),Nb_EVEs, 0 ))

EVEs_shared_event_df<-EVEs_shared_event_df %>% mutate (
  Nb_endo_EVEs = ifelse(Nb_EVEs >= 1 & lifecycle1 %in% c('endoparasitoide'),Nb_EVEs, 0 ))

#NB EVES according to genomic structure
#dsDNA
EVEs_shared_event_df<-EVEs_shared_event_df %>% mutate (
  Nb_EVEs_dsDNA = ifelse(Nb_EVEs >= 1 & consensus_genomic_structure %in% c('dsDNA'),Nb_EVEs, 0 ))
#ssDNA
EVEs_shared_event_df<-EVEs_shared_event_df %>% mutate (
  Nb_EVEs_ssDNA = ifelse(Nb_EVEs >= 1 & consensus_genomic_structure %in% c('ssDNA'),Nb_EVEs, 0 ))
#dsRNA
EVEs_shared_event_df<-EVEs_shared_event_df %>% mutate (
  Nb_EVEs_dsRNA = ifelse(Nb_EVEs >= 1 & consensus_genomic_structure %in% c('dsRNA'),Nb_EVEs, 0 ))
#ssRNA
EVEs_shared_event_df<-EVEs_shared_event_df %>% mutate (
  Nb_EVEs_ssRNA = ifelse(Nb_EVEs >= 1 & consensus_genomic_structure %in% c('ssRNA'),Nb_EVEs, 0 ))
#Unclassified
EVEs_shared_event_df<-EVEs_shared_event_df %>% mutate (
  Nb_EVEs_Unclassified_Events = ifelse(Nb_EVEs >= 1 & consensus_genomic_structure%in% c('Unknown'),Nb_EVEs, 0 ))

#NB EVES EVENT according to genomic structure
#dsDNA
EVEs_shared_event_df<-EVEs_shared_event_df %>% mutate (
  Nb_EVEs_dsDNA_Events = ifelse(Nb_EVEs >= 1 & consensus_genomic_structure %in% c('dsDNA'),1, 0 ))
#ssDNA
EVEs_shared_event_df<-EVEs_shared_event_df %>% mutate (
  Nb_EVEs_ssDNA_Events  = ifelse(Nb_EVEs >= 1 & consensus_genomic_structure %in% c('ssDNA'),1, 0 ))
#dsRNA
EVEs_shared_event_df<-EVEs_shared_event_df %>% mutate (
  Nb_EVEs_dsRNA_Events  = ifelse(Nb_EVEs >= 1 & consensus_genomic_structure %in% c('dsRNA'),1, 0 ))
#ssRNA
EVEs_shared_event_df<-EVEs_shared_event_df %>% mutate (
  Nb_EVEs_ssRNA_Events  = ifelse(Nb_EVEs >= 1 & consensus_genomic_structure %in% c('ssRNA'),1, 0 ))
#Unclassified
EVEs_shared_event_df<-EVEs_shared_event_df %>% mutate (
  Nb_EVEs_Unclassified_Events = ifelse(Nb_EVEs >= 1 & consensus_genomic_structure%in% c('Unknown'),1, 0 ))


#Count dEVEs
#Add Number dEVEs
library(tidyverse)
Env_table_bis<- Env_table[!duplicated(Env_table$New_query_bis2),]
Env_table_bis<-select(Env_table_bis,"Clustername",'New_query_bis2','Event',"Mean_dNdS","Pvalue_dNdS",'FDR_pvalue_dNdS','pseudogenized','SE_dNdS','TPM_all',"ORF_perc")
Env_table_bis <-Env_table_bis %>%
  group_by(Clustername,New_query_bis2) %>%
  filter(any(FDR_pvalue_dNdS == 1 & as.numeric(as.character(Mean_dNdS))+as.numeric(as.character(SE_dNdS)) < 1 & pseudogenized ==0 & ORF_perc>1 | TPM_all>=1000 & pseudogenized ==0 & ORF_perc>1))


names(Env_table_bis)[names(Env_table_bis)=="Clustername"] <- "Clustername_bis"
names(EVEs_shared_event_df)[names(EVEs_shared_event_df)=="Event"] <- "Event_Groups"

EVEs_shared_event_df$rowname<-NULL

sub<-EVEs_shared_event_df  %>%
  rownames_to_column()  %>%
  mutate(Event_Groups = as.character(Event_Groups)) %>%
  separate_rows(Event_Groups, sep = ",") %>%
  left_join(.,
            Env_table_bis %>%
              unite(col = "Event_Groups", Event, Clustername_bis) %>%
              mutate(count = if_else(Mean_dNdS < 1 & Pvalue_dNdS < 0.05 | TPM_all>=159.559000 & pseudogenized ==0 & ORF_perc>1 ,1L, 0L))) %>%
  distinct(Event_Groups, .keep_all = TRUE) %>%
  group_by(rowname,) %>%
  summarise(Event_Groups = paste(unique(Event_Groups), collapse = ","),
            Nb_dEVEs = sum(count, na.rm=T))

EVEs_shared_event_df$rowname <- rownames(EVEs_shared_event_df)
EVEs_shared_event_df<- merge(x = sub, y =EVEs_shared_event_df, by = "rowname", all = TRUE)
EVEs_shared_event_df$Event_Groups.x <- NULL
EVEs_shared_event_df$Nb_dEVEs[is.na(EVEs_shared_event_df$Nb_dEVEs)] <- 0

names(EVEs_shared_event_df)[names(EVEs_shared_event_df)=="Event_Groups.y"] <- "Event"

#EVENTS according to lifestyles
EVEs_shared_event_df<-EVEs_shared_event_df %>% mutate (
  Nb_endo_dEVEs_Events = ifelse(Nb_dEVEs >= 1 & lifecycle1 %in% c('endoparasitoide'),1,0 ))
EVEs_shared_event_df<-EVEs_shared_event_df %>% mutate (
  Nb_ecto_dEVEs_Events = ifelse(Nb_dEVEs >= 1 & lifecycle1 %in% c('ectoparasitoide'),1,0 ))
EVEs_shared_event_df<-EVEs_shared_event_df %>% mutate (
  Nb_freeliving_dEVEs_Events = ifelse(Nb_dEVEs >= 1 & lifecycle1 %in% c('freeliving'),1,0 ))

#NB EVES according to lifestyles
EVEs_shared_event_df<-EVEs_shared_event_df %>% mutate (
  Nb_ecto_dEVEs = ifelse(Nb_dEVEs >= 1 & lifecycle1 %in% c('ectoparasitoide'),Nb_dEVEs, 0 ))

EVEs_shared_event_df<-EVEs_shared_event_df %>% mutate (
  Nb_freeliving_dEVEs = ifelse(Nb_dEVEs >= 1 & lifecycle1 %in% c('freeliving'),Nb_dEVEs, 0 ))

EVEs_shared_event_df<-EVEs_shared_event_df %>% mutate (
  Nb_endo_dEVEs = ifelse(Nb_dEVEs >= 1 & lifecycle1 %in% c('endoparasitoide'),Nb_dEVEs, 0 ))

#NB EVES according to genomic structure
#dsDNA
EVEs_shared_event_df<-EVEs_shared_event_df %>% mutate (
  Nb_dEVEs_dsDNA = ifelse(Nb_dEVEs >= 1 & consensus_genomic_structure %in% c('dsDNA'),Nb_dEVEs, 0 ))
#ssDNA
EVEs_shared_event_df<-EVEs_shared_event_df %>% mutate (
  Nb_dEVEs_ssDNA = ifelse(Nb_dEVEs >= 1 & consensus_genomic_structure %in% c('ssDNA'),Nb_dEVEs, 0 ))
#dsRNA
EVEs_shared_event_df<-EVEs_shared_event_df %>% mutate (
  Nb_dEVEs_dsRNA = ifelse(Nb_dEVEs >= 1 & consensus_genomic_structure %in% c('dsRNA'),Nb_dEVEs, 0 ))
#ssRNA
EVEs_shared_event_df<-EVEs_shared_event_df %>% mutate (
  Nb_dEVEs_ssRNA = ifelse(Nb_dEVEs >= 1 & consensus_genomic_structure %in% c('ssRNA'),Nb_dEVEs, 0 ))
#Unclassified
EVEs_shared_event_df<-EVEs_shared_event_df %>% mutate (
  Nb_dEVEs_Unclassified_Events = ifelse(Nb_dEVEs >= 1 & consensus_genomic_structure%in% c('Unknown'),Nb_dEVEs, 0 ))

#NB EVES EVENT according to genomic structure
#dsDNA
EVEs_shared_event_df<-EVEs_shared_event_df %>% mutate (
  Nb_dEVEs_dsDNA_Events = ifelse(Nb_dEVEs >= 1 & consensus_genomic_structure %in% c('dsDNA'),1, 0 ))
#ssDNA
EVEs_shared_event_df<-EVEs_shared_event_df %>% mutate (
  Nb_dEVEs_ssDNA_Events  = ifelse(Nb_dEVEs >= 1 & consensus_genomic_structure %in% c('ssDNA'),1, 0 ))
#dsRNA
EVEs_shared_event_df<-EVEs_shared_event_df %>% mutate (
  Nb_dEVEs_dsRNA_Events  = ifelse(Nb_dEVEs >= 1 & consensus_genomic_structure %in% c('dsRNA'),1, 0 ))
#ssRNA
EVEs_shared_event_df<-EVEs_shared_event_df %>% mutate (
  Nb_dEVEs_ssRNA_Events  = ifelse(Nb_dEVEs >= 1 & consensus_genomic_structure %in% c('ssRNA'),1, 0 ))
#Unclassified
EVEs_shared_event_df<-EVEs_shared_event_df %>% mutate (
  Nb_dEVEs_Unclassified_Events = ifelse(Nb_dEVEs >= 1 & consensus_genomic_structure%in% c('Unknown'),1, 0 ))



########

#Save to count in summary :
EVEs_shared_event_df_EVEs_count <- EVEs_shared_event_df


##################

detachAllPackages()

library(tidyverse)
library(pastecs)
library(ggplot2)
library(gggenes)
#library(gsubfn)
library(data.table)
library(stringr)


#For events EVEs alone only
#This script allows to group together putatives events based on the fact that EVEs are in the same scaffold or by regrouping EVEs comming from the same virus family donnor


####################################################
#Try to capture number of EVEs within alone events
####################################################

Restricted_candidat

Restricted_candidat <- Restricted_candidat[order(Restricted_candidat$TPM_all, decreasing = TRUE), ]

table_noduplicate_alone_EVEs <-Restricted_candidat %>%
  group_by(Clustername,query,consensus_genomic_structure)%>%
  distinct(query, .keep_all = TRUE)%>%
  mutate(nrows=n())%>%
  filter(any(nrows <=1))

table_noduplicate_alone_EVEs<-table_noduplicate_alone_EVEs[!duplicated(table_noduplicate_alone_EVEs[ , c("query","Clustername","Event")]),]
table_noduplicate_alone_EVEs <-select(as.data.frame(table_noduplicate_alone_EVEs),'query','Clustername','New_query_bis','consensus_family','consensus_genomic_structure','Species_name','Event','Scaff_name')


#Deal with the fact that duplicated EVEs can be into different Event, the idea is just to check the boostrap value that separates two putative duplicates (g.e EVEs within same cluster and same specie)
table_noduplicate_alone_EVEs<-table_noduplicate_alone_EVEs[!duplicated(table_noduplicate_alone_EVEs[ , c("New_query_bis","Clustername","Event")]),]


list_df<-expand.grid.unique(as.character(table_noduplicate_alone_EVEs$Clustername),table_noduplicate_alone_EVEs$Species_name)
colnames(list_df)<-c("Clustername","Species_name")
list_df<-list_df[!duplicated(list_df[ , c("Clustername","Species_name")]),]
list_df<-as.data.frame(list_df)
list_df$present<-"yes"

list_df<-merge(select(table_noduplicate_alone_EVEs,"Clustername","Species_name"),list_df,by=c("Clustername","Species_name"),all.x = TRUE)

table_noduplicate_alone_EVEs$BootValue<-"NA"
library(ape)
library(phytools)
library(phylobase)
list_df<-list_df[!duplicated(list_df[ , c("Clustername","Species_name")]),]


#Correct Bootstrapvalues
for(i in unique(list_df$Clustername)) {
  sub_table_noduplicate_alone_EVEs<-table_noduplicate_alone_EVEs[with(table_noduplicate_alone_EVEs, Clustername == i ), ]
  sub_table_noduplicate_alone_EVEs<-sub_table_noduplicate_alone_EVEs[!is.na(sub_table_noduplicate_alone_EVEs$Clustername), ]
  sub_table_noduplicate_alone_EVEs$query<-gsub(":","_",sub_table_noduplicate_alone_EVEs$query)
  sub_table_noduplicate_alone_EVEs$query<-gsub("\\(","_",sub_table_noduplicate_alone_EVEs$query)
  sub_table_noduplicate_alone_EVEs$query<-gsub("\\)","_",sub_table_noduplicate_alone_EVEs$query)
  sub_table_noduplicate_alone_EVEs$query<-gsub(" ","",sub_table_noduplicate_alone_EVEs$query)
  for (a in unique(sub_table_noduplicate_alone_EVEs$Event)){
    sub_sub_table_noduplicate_alone_EVEs<-sub_table_noduplicate_alone_EVEs[sub_table_noduplicate_alone_EVEs$Event == a, ]
    if (length(dim(sub_sub_table_noduplicate_alone_EVEs)) >=1 ){
      Event_number<-a
      if (file.exists(paste0("/Users/bguinet/Desktop/these/Cluster_phylogeny_filtred/",unique(sub_table_noduplicate_alone_EVEs$Clustername),"_AA.dna.treefile"))) {
        print(i)
        PhyloTree=read.tree(paste0("/Users/bguinet/Desktop/these/Cluster_phylogeny_filtred/",unique(sub_table_noduplicate_alone_EVEs$Clustername),"_AA.dna.treefile"))
        PhyloTree_labels<-PhyloTree$tip.label
        tax_to_remove<-PhyloTree_labels[grepl("HSP",PhyloTree_labels)]
        tax_to_remove<-gsub("_[0-9]-HSPs.*", "",tax_to_remove)
        list_candidate<-sub_sub_table_noduplicate_alone_EVEs$query
        tryCatch({
          PhyloTree <- as(PhyloTree, "phylo4")
        }, error=function(e){})
        tryCatch({
          mm<- MRCA(PhyloTree,list_candidate)
        }, error=function(e){})
        bootstrap_value<-names(mm)
        print(list_candidate)
        print(bootstrap_value)
        table_noduplicate_alone_EVEs$BootValue[table_noduplicate_alone_EVEs$Event==Event_number] <- bootstrap_value
      }
    }
  }
}
for(i in 1:nrow(list_df)) {
  sub_table_noduplicate_alone_EVEs<-table_noduplicate_alone_EVEs[with(table_noduplicate_alone_EVEs, Clustername == list_df$Clustername[[i]] & Species_name == list_df$Species_name[[i]]), ]
  if (length(dim(sub_table_noduplicate_alone_EVEs)) >=1 ){
    sub_table_noduplicate_alone_EVEs<-sub_table_noduplicate_alone_EVEs[!is.na(sub_table_noduplicate_alone_EVEs$Clustername), ]
    sub_table_noduplicate_alone_EVEs$query<-gsub(":","_",sub_table_noduplicate_alone_EVEs$query)
    sub_table_noduplicate_alone_EVEs$query<-gsub("\\(","_",sub_table_noduplicate_alone_EVEs$query)
    sub_table_noduplicate_alone_EVEs$query<-gsub("\\)","_",sub_table_noduplicate_alone_EVEs$query)
    sub_table_noduplicate_alone_EVEs$query<-gsub(" ","",sub_table_noduplicate_alone_EVEs$query)
    if (file.exists(paste0("/Users/bguinet/Desktop/these/Cluster_phylogeny_filtred/",unique(sub_table_noduplicate_alone_EVEs$Clustername),"_AA.dna.treefile"))) {
      PhyloTree=read.tree(paste0("/Users/bguinet/Desktop/these/Cluster_phylogeny_filtred/",unique(sub_table_noduplicate_alone_EVEs$Clustername),"_AA.dna.treefile"))
      PhyloTree_labels<-PhyloTree$tip.label
      tax_to_remove<-PhyloTree_labels[grepl("HSP",PhyloTree_labels)]
      tax_to_remove<-gsub("_[0-9]-HSPs.*", "",tax_to_remove)
      list_candidate<-sub_table_noduplicate_alone_EVEs$query
      list_candidate<-list_candidate[!grepl(paste(tax_to_remove, collapse = "|"), list_candidate)]
      if (length(list_candidate)>1){
        PhyloTree<-midpoint.root(PhyloTree)
        tryCatch({
          PhyloTree <- as(PhyloTree, "phylo4")
        }, error=function(e){})
        mm<- MRCA(PhyloTree,list_candidate)
        bootstrap_value<-names(mm)
        #print(list_candidate)
        print(bootstrap_value)
        print(list_df$Species_name[[i]])
        print(list_df$Clustername[[i]])
        table_noduplicate_alone_EVEs<-table_noduplicate_alone_EVEs%>%
          group_by(Clustername,Species_name)%>%
          mutate(BootValue = case_when(any(grepl(list_df$Clustername[[i]],Clustername) & Species_name %in% list_df$Species_name[[i]]  ) ~ bootstrap_value,TRUE ~ BootValue))
      }
    }
  }
}

#Remove duplicated EVEs when the bootstrap support is below 80
table_noduplicate_alone_EVEs$BootValue[table_noduplicate_alone_EVEs$BootValue=="Root"] <-0
table_noduplicate_alone_EVEs$BootValue[table_noduplicate_alone_EVEs$BootValue=="NA"] <- -1
table_noduplicate_alone_EVEs$BootValue <- as.numeric(table_noduplicate_alone_EVEs$BootValue)
#table_noduplicate_alone_EVEs<-table_noduplicate_alone_EVEs[ table_noduplicate_alone_EVEs$BootValue < 80,]


#Add lifecycle
anoleData<-read.csv("/Users/bguinet/Documents/All_informations_table_fam_busco_life_style.csv",sep=";",header=T)

#Correct lifestyles
anoleData$lifecycle1[anoleData$Species_name=="Megastigmus_dorsalis"]<-"ectoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Megastigmus_stigmatizans"]<-"ectoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Cecidostiba_semifascia"]<-"ectoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Cecidostiba_fungosa"]<-"ectoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Eretmocerus_eremicus"]<-"endoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Ormyrus_pomaceus"]<-"ectoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Ormyrus_nitidulus"]<-"ectoparasitoide"


anoleData<-select(anoleData,'Species_name','lifecycle1')
colnames(anoleData)<-c('Species_name','lifecycle1')
table_noduplicate_alone_EVEs <-merge(table_noduplicate_alone_EVEs ,anoleData,by=c("Species_name"))



#Now we will map the EVEs event within the phylogeny, but for that we need to create a dataframe where we add each Event number for each node number
col.names = c("Node_number", "Event",'consensus_genomic_structure','Cluster','consensus_family','Clustername','species','Nb_EVEs','Nb_EVEs_dsDNA','Nb_EVEs_ssDNA','Nb_EVEs_ssRNA','Nb_EVEs_dsRNA','Nb_EVEs_Unclassified','Nb_EVEs_dsDNA_Events','Nb_EVEs_ssDNA_Events','Nb_EVEs_ssRNA_Events','Nb_EVEs_dsRNA_Events','Nb_EVEs_Unclassified_Events','Nb_freeliving_EVEs','Nb_ecto_EVEs','Nb_endo_EVEs','Nb_freeliving_EVEs_Events','Nb_ecto_EVEs_Events','Nb_endo_EVEs_Events','lifecycle1')
EVEs_alone_event_noduplicate_df<- read.table(text = "",
                                             col.names = col.names)

EVEs_alone_event_noduplicate_df$lifecycle1<- as.character(EVEs_alone_event_noduplicate_df$lifecycle1)
EVEs_alone_event_noduplicate_df$consensus_family<- as.character(EVEs_alone_event_noduplicate_df$consensus_family)
EVEs_alone_event_noduplicate_df$Clustername<- as.character(EVEs_alone_event_noduplicate_df$Clustername)
EVEs_alone_event_noduplicate_df$species<- as.character(EVEs_alone_event_noduplicate_df$species)
EVEs_alone_event_noduplicate_df$Event<- as.character(EVEs_alone_event_noduplicate_df$Event)
EVEs_alone_event_noduplicate_df$consensus_genomic_structure<- as.character(EVEs_alone_event_noduplicate_df$consensus_genomic_structure)

library(tidyverse)
library(phytools)
for(i in 1:nrow(table_noduplicate_alone_EVEs)) {
  Nb_EVEs<-length(unlist(strsplit(as.character(table_noduplicate_alone_EVEs$Clustername[[i]]),",")))
  if(grepl("free",table_noduplicate_alone_EVEs$lifecycle1[[i]])){
    Nb_freeliving_EVEs<-1
    Nb_ecto_EVEs<-0
    Nb_endo_EVEs<-0
  }
  if(grepl("endo",table_noduplicate_alone_EVEs$lifecycle1[[i]])){
    Nb_freeliving_EVEs<-0
    Nb_ecto_EVEs<-0
    Nb_endo_EVEs<-1
  }
  if(grepl("ecto",table_noduplicate_alone_EVEs$lifecycle1[[i]])){
    Nb_freeliving_EVEs<-0
    Nb_ecto_EVEs<-1
    Nb_endo_EVEs<-0
  }
  if(grepl("dsDNA",table_noduplicate_alone_EVEs$consensus_genomic_structure[[i]])){
    Nb_EVEs_dsDNA<-1
    Nb_EVEs_ssDNA<-0
    Nb_EVEs_ssRNA<-0
    Nb_EVEs_dsRNA<-0
    Nb_EVEs_Unclassified<-0
  }
  if(grepl("ssDNA",table_noduplicate_alone_EVEs$consensus_genomic_structure[[i]])){
    Nb_EVEs_ssDNA<-1
    Nb_EVEs_dsDNA<-0
    Nb_EVEs_ssRNA<-0
    Nb_EVEs_dsRNA<-0
    Nb_EVEs_Unclassified<-0
  }
  if(grepl("ssRNA",table_noduplicate_alone_EVEs$consensus_genomic_structure[[i]])){
    Nb_EVEs_ssRNA<-1
    Nb_EVEs_ssDNA<-0
    Nb_EVEs_dsDNA<-0
    Nb_EVEs_dsRNA<-0
    Nb_EVEs_Unclassified<-0
  }
  if(grepl("dsRNA",table_noduplicate_alone_EVEs$consensus_genomic_structure[[i]])){
    Nb_EVEs_dsRNA<-1
    Nb_EVEs_ssDNA<-0
    Nb_EVEs_ssRNA<-0
    Nb_EVEs_dsDNA<-0
    Nb_EVEs_Unclassified<-0
  }
  if(grepl("Unknown",table_noduplicate_alone_EVEs$consensus_genomic_structure[[i]])){
    Nb_EVEs_Unclassified<-1
    Nb_EVEs_ssDNA<-0
    Nb_EVEs_ssRNA<-0
    Nb_EVEs_dsRNA<-0
    Nb_EVEs_dsDNA<-0
  }
  lifecycle1<-table_noduplicate_alone_EVEs$lifecycle1[[i]]
  Cluster<-table_noduplicate_alone_EVEs$Clustername[[i]]
  row <- table_noduplicate_alone_EVEs$Events_species[[i]]
  family<-table_noduplicate_alone_EVEs$consensus_family[[i]]
  consensus_genomic_structure<-table_noduplicate_alone_EVEs$consensus_genomic_structure[[i]]
  #row<-gsub("\\[","",row)
  #row<-gsub("\\]","",row)
  #row<-list(row)
  #row<-gsub(" ","",row)
  row<-table_noduplicate_alone_EVEs$Species_name[[i]]
  EVEs_alone_event_noduplicate_df<-EVEs_alone_event_noduplicate_df %>% add_row(Event = gsub('_.*','',table_noduplicate_alone_EVEs$Event[[i]]),consensus_genomic_structure=table_noduplicate_alone_EVEs$consensus_genomic_structure[[i]],consensus_family=as.character(family),Clustername=Cluster,species=as.character(row),Nb_EVEs=Nb_EVEs,Nb_EVEs_dsDNA=Nb_EVEs_dsDNA,Nb_EVEs_ssDNA=Nb_EVEs_ssDNA,Nb_EVEs_ssRNA=Nb_EVEs_ssRNA,Nb_EVEs_dsRNA=Nb_EVEs_dsRNA,Nb_EVEs_Unclassified=Nb_EVEs_Unclassified,Nb_freeliving_EVEs=Nb_freeliving_EVEs,Nb_ecto_EVEs=Nb_ecto_EVEs,Nb_endo_EVEs=Nb_endo_EVEs,lifecycle1=lifecycle1)
}

EVEs_alone_event_noduplicate_df<-EVEs_alone_event_noduplicate_df[!duplicated(EVEs_alone_event_noduplicate_df[c("Clustername","Event")]),]

detachAllPackages()

library(tidyverse)
library(pastecs)
library(ggplot2)
library(gggenes)
#library(gsubfn)
library(data.table)
library(stringr)

table_alone_EVEs_event<-table_noduplicate_alone_EVEs %>%
  group_by(Clustername,Event,consensus_genomic_structure)%>%
  distinct(New_query_bis, .keep_all = TRUE)%>%
  mutate(nrows=n())%>%
  filter(any(nrows <=1))%>%
  #select(Clustername,New_query_bis,query,Event,nrows,consensus_family)%>%
  arrange(Scaff_name,New_query_bis,consensus_family,consensus_genomic_structure) %>%
  # create a Group_2 which is combination of all Group for each family
  group_by(consensus_family,New_query_bis,consensus_genomic_structure) %>%
  mutate(Group_2 = list(Scaff_name)) %>%
  ungroup() %>%
  # Create Group_3 which is the full combined Group for all intersect Group
  mutate(Group_3 = map(.[["Group_2"]], function(x) { look_up_group(one_group = x, lookup_list = .[["Group_2"]]) })) %>%
  # Combine all Group_3 into a Group_final
  mutate(Group_final = unlist(map(Group_3, function(x) { paste (x, collapse = ",")} ))) %>%
  mutate(Species =New_query_bis )%>%
  # Finally put them all together.
  select(Species,New_query_bis,Group_final, consensus_family,Event,Clustername,consensus_genomic_structure) %>%
  group_by(Group_final,consensus_genomic_structure,Species) %>%
  summarize(family = paste(consensus_family, collapse = ","),species=paste(unique(Species), collapse = ","),Events=paste(unique(Event), collapse = ","),Clusters=paste(unique(Clustername), collapse = ","),consensus_genomic_structure=paste(unique(consensus_genomic_structure), collapse = ","), .groups = "drop")


#Keep only Event with mutiple EVEs
if (Only_mulptiple_EVE =="yes"){
  table_alone_EVEs_event<-table_alone_EVEs_event %>% filter(grepl(',',Clusters))
  savetable_alone_EVEs_event<-table_alone_EVEs_event
}


#Add lifecycle
anoleData<-read.csv("/Users/bguinet/Documents/All_informations_table_fam_busco_life_style.csv",sep=";",header=T)
anoleData$lifecycle1[anoleData$Species_name=="Megastigmus_dorsalis"]<-"ectoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Megastigmus_stigmatizans"]<-"ectoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Cecidostiba_semifascia"]<-"ectoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Cecidostiba_fungosa"]<-"ectoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Eretmocerus_eremicus"]<-"endoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Ormyrus_pomaceus"]<-"ectoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Ormyrus_nitidulus"]<-"ectoparasitoide"
anoleData<-select(anoleData,'Species_name','lifecycle1')
colnames(anoleData)<-c('species','lifecycle1')


anoleData<-anoleData[!anoleData$species %in% c("Andricus_quercusramuli","Neuroterus_quercusbaccarum","Andricus_inflator",
                                               "Andricus_curvator","Pseudoneuroterus_saliens","Trichopria_drosophilae","Andricus_grossulariae"),]


table_alone_EVEs_event[table_alone_EVEs_event$species=="Trichopria_sp_",]$species <- "Trichopria_sp_970989"
table_alone_EVEs_event[table_alone_EVEs_event$Species=="Trichopria_sp_",]$Species <- "Trichopria_sp_970989"
table_alone_EVEs_event<-merge(table_alone_EVEs_event,anoleData,by=c("species"))

Families_EVEs_alone_event_df<-table(table_alone_EVEs_event$family)

#To remo
table_alone_EVEs_event <-table_alone_EVEs_event %>%
  filter(!(grepl('Cluster69840|Cluster51209', Clusters) & species=='Microplitis_demolitor'))
table_alone_EVEs_event <-table_alone_EVEs_event %>%
  filter(!(grepl('Cluster69840', Clusters) & species=='Cotesia_vestalis'))

#####




#Add Number EVEs
table_alone_EVEs_event$Nb_EVEs <- str_count(table_alone_EVEs_event$Clusters, ",")+1


#EVENTS according to lifestyles
table_alone_EVEs_event<-table_alone_EVEs_event %>% mutate (
  Nb_endo_EVEs_Events = ifelse(Nb_EVEs >= 1 & lifecycle1 %in% c('endoparasitoide'),1,0 ))
table_alone_EVEs_event<-table_alone_EVEs_event %>% mutate (
  Nb_ecto_EVEs_Events = ifelse(Nb_EVEs >= 1 & lifecycle1 %in% c('ectoparasitoide'),1,0 ))
table_alone_EVEs_event<-table_alone_EVEs_event %>% mutate (
  Nb_freeliving_EVEs_Events = ifelse(Nb_EVEs >= 1 & lifecycle1 %in% c('freeliving'),1,0 ))
#NB EVES according to lifestyles
table_alone_EVEs_event %>% mutate (
  Nb_ecto_EVEs = ifelse(Nb_EVEs >= 1 & lifecycle1 %in% c('ectoparasitoide'),Nb_EVEs, 0 ))
table_alone_EVEs_event<-table_alone_EVEs_event %>% mutate (
  Nb_freeliving_EVEs = ifelse(Nb_EVEs >= 1 & lifecycle1 %in% c('freeliving'),Nb_EVEs, 0 ))
table_alone_EVEs_event<-table_alone_EVEs_event %>% mutate (
  Nb_endo_EVEs = ifelse(Nb_EVEs >= 1 & lifecycle1 %in% c('endoparasitoide'),Nb_EVEs, 0 ))

#NB EVEs according to genomic structure
#dsDNA
table_alone_EVEs_event<-table_alone_EVEs_event %>% mutate (
  Nb_EVEs_dsDNA = ifelse(Nb_EVEs >= 1 & consensus_genomic_structure %in% c('dsDNA'),Nb_EVEs, 0 ))
#ssDNA
table_alone_EVEs_event<-table_alone_EVEs_event %>% mutate (
  Nb_EVEs_ssDNA = ifelse(Nb_EVEs >= 1 & consensus_genomic_structure %in% c('ssDNA'),Nb_EVEs, 0 ))
#dsRNA
table_alone_EVEs_event<-table_alone_EVEs_event %>% mutate (
  Nb_EVEs_dsRNA = ifelse(Nb_EVEs >= 1 & consensus_genomic_structure %in% c('dsRNA'),Nb_EVEs, 0 ))
#ssRNA
table_alone_EVEs_event<-table_alone_EVEs_event %>% mutate (
  Nb_EVEs_ssRNA = ifelse(Nb_EVEs >= 1 & consensus_genomic_structure %in% c('ssRNA'),Nb_EVEs, 0 ))
#Unclassified
table_alone_EVEs_event<-table_alone_EVEs_event %>% mutate (
  Nb_EVEs_Unclassified = ifelse(Nb_EVEs >= 1 & consensus_genomic_structure%in% c('Unknown'),Nb_EVEs, 0 ))


#NB EVES EVENT according to genomic structure
#dsDNA
table_alone_EVEs_event<-table_alone_EVEs_event %>% mutate (
  Nb_EVEs_dsDNA_Events = ifelse(Nb_EVEs >= 1 & consensus_genomic_structure %in% c('dsDNA'),1, 0 ))
#ssDNA
table_alone_EVEs_event<-table_alone_EVEs_event %>% mutate (
  Nb_EVEs_ssDNA_Events  = ifelse(Nb_EVEs >= 1 & consensus_genomic_structure %in% c('ssDNA'),1, 0 ))
#dsRNA
table_alone_EVEs_event<-table_alone_EVEs_event %>% mutate (
  Nb_EVEs_dsRNA_Events  = ifelse(Nb_EVEs >= 1 & consensus_genomic_structure %in% c('dsRNA'),1, 0 ))
#ssRNA
table_alone_EVEs_event<-table_alone_EVEs_event %>% mutate (
  Nb_EVEs_ssRNA_Events  = ifelse(Nb_EVEs >= 1 & consensus_genomic_structure %in% c('ssRNA'),1, 0 ))
#Unclassified
table_alone_EVEs_event<-table_alone_EVEs_event %>% mutate (
  Nb_EVEs_Unclassified_Events = ifelse(Nb_EVEs >= 1 & consensus_genomic_structure%in% c('Unknown'),1, 0 ))



#Add dEVES count
library(tidyverse)
Env_table_bis<- Env_table[!duplicated(Env_table$New_query_bis2),]
Env_table_bis<-select(Env_table_bis,"Clustername",'New_query_bis2','Event',"Mean_dNdS","Pvalue_dNdS",'FDR_pvalue_dNdS','pseudogenized','SE_dNdS','TPM_all','ORF_perc')
Env_table_bis <-Env_table_bis %>%
  group_by(Clustername,New_query_bis2) %>%
  filter(any(FDR_pvalue_dNdS == 1 & as.numeric(as.character(Mean_dNdS))+as.numeric(as.character(SE_dNdS)) < 1 & pseudogenized ==0 & ORF_perc>1| TPM_all >=1000 & pseudogenized ==0 & ORF_perc>1))

names(Env_table_bis)[names(Env_table_bis)=="Clustername"] <- "Clustername_bis"
names(table_alone_EVEs_event)[names(table_alone_EVEs_event)=="Events"] <- "Event_Groups"

table_alone_EVEs_event$rowname<-NULL

sub<-table_alone_EVEs_event  %>%
  rownames_to_column()  %>%
  mutate(Event_Groups = as.character(Event_Groups)) %>%
  separate_rows(Event_Groups, sep = ",") %>%
  left_join(.,
            Env_table_bis %>%
              unite(col = "Event_Groups", Event, Clustername_bis) %>%
              mutate(count = if_else(Mean_dNdS < 1 & Pvalue_dNdS < 0.05 | TPM_all>= 1000 & pseudogenized ==0 & ORF_perc>1 ,1L, 0L))) %>%
  distinct(Event_Groups, .keep_all = TRUE) %>%
  group_by(rowname,) %>%
  summarise(Event_Groups = paste(unique(Event_Groups), collapse = ","),
            Nb_dEVEs = sum(count, na.rm=T))

table_alone_EVEs_event$rowname <- rownames(table_alone_EVEs_event)
table_alone_EVEs_event<- merge(x = sub, y =table_alone_EVEs_event, by = "rowname", all = TRUE)
table_alone_EVEs_event$Event_Groups.x <- NULL
table_alone_EVEs_event$Nb_dEVEs[is.na(table_alone_EVEs_event$Nb_dEVEs)] <- 0

names(table_alone_EVEs_event)[names(table_alone_EVEs_event)=="Event_Groups.y"] <- "Event"

#EVENTS according to lifestyles
table_alone_EVEs_event<-table_alone_EVEs_event %>% mutate (
  Nb_endo_dEVEs_Events = ifelse(Nb_dEVEs >= 1 & lifecycle1 %in% c('endoparasitoide'),1,0 ))
table_alone_EVEs_event<-table_alone_EVEs_event %>% mutate (
  Nb_ecto_dEVEs_Events = ifelse(Nb_dEVEs >= 1 & lifecycle1 %in% c('ectoparasitoide'),1,0 ))
table_alone_EVEs_event<-table_alone_EVEs_event %>% mutate (
  Nb_freeliving_dEVEs_Events = ifelse(Nb_dEVEs >= 1 & lifecycle1 %in% c('freeliving'),1,0 ))
#NB EVES according to lifestyles
table_alone_EVEs_event %>% mutate (
  Nb_ecto_dEVEs = ifelse(Nb_dEVEs >= 1 & lifecycle1 %in% c('ectoparasitoide'),Nb_dEVEs, 0 ))
table_alone_EVEs_event<-table_alone_EVEs_event %>% mutate (
  Nb_freeliving_dEVEs = ifelse(Nb_dEVEs >= 1 & lifecycle1 %in% c('freeliving'),Nb_dEVEs, 0 ))
table_alone_EVEs_event<-table_alone_EVEs_event %>% mutate (
  Nb_endo_dEVEs = ifelse(Nb_dEVEs >= 1 & lifecycle1 %in% c('endoparasitoide'),Nb_dEVEs, 0 ))

#NB dEVES according to genomic structure
#dsDNA
table_alone_EVEs_event<-table_alone_EVEs_event %>% mutate (
  Nb_dEVEs_dsDNA = ifelse(Nb_dEVEs >= 1 & consensus_genomic_structure %in% c('dsDNA'),Nb_dEVEs, 0 ))
#ssDNA
table_alone_EVEs_event<-table_alone_EVEs_event %>% mutate (
  Nb_dEVEs_ssDNA = ifelse(Nb_dEVEs >= 1 & consensus_genomic_structure %in% c('ssDNA'),Nb_dEVEs, 0 ))
#dsRNA
table_alone_EVEs_event<-table_alone_EVEs_event %>% mutate (
  Nb_dEVEs_dsRNA = ifelse(Nb_dEVEs >= 1 & consensus_genomic_structure %in% c('dsRNA'),Nb_dEVEs, 0 ))
#ssRNA
table_alone_EVEs_event<-table_alone_EVEs_event %>% mutate (
  Nb_dEVEs_ssRNA = ifelse(Nb_dEVEs >= 1 & consensus_genomic_structure %in% c('ssRNA'),Nb_dEVEs, 0 ))
#Unclassified
table_alone_EVEs_event<-table_alone_EVEs_event %>% mutate (
  Nb_dEVEs_Unclassified = ifelse(Nb_dEVEs >= 1 & consensus_genomic_structure%in% c('Unknown'),Nb_dEVEs, 0 ))


#NB EVES EVENT according to genomic structure
#dsDNA
table_alone_EVEs_event<-table_alone_EVEs_event %>% mutate (
  Nb_dEVEs_dsDNA_Events = ifelse(Nb_dEVEs >= 1 & consensus_genomic_structure %in% c('dsDNA'),1, 0 ))
#ssDNA
table_alone_EVEs_event<-table_alone_EVEs_event %>% mutate (
  Nb_dEVEs_ssDNA_Events  = ifelse(Nb_dEVEs >= 1 & consensus_genomic_structure %in% c('ssDNA'),1, 0 ))
#dsRNA
table_alone_EVEs_event<-table_alone_EVEs_event %>% mutate (
  Nb_dEVEs_dsRNA_Events  = ifelse(Nb_dEVEs >= 1 & consensus_genomic_structure %in% c('dsRNA'),1, 0 ))
#ssRNA
table_alone_EVEs_event<-table_alone_EVEs_event %>% mutate (
  Nb_dEVEs_ssRNA_Events  = ifelse(Nb_dEVEs >= 1 & consensus_genomic_structure %in% c('ssRNA'),1, 0 ))
#Unclassified
table_alone_EVEs_event<-table_alone_EVEs_event %>% mutate (
  Nb_dEVEs_Unclassified_Events = ifelse(Nb_dEVEs >= 1 & consensus_genomic_structure%in% c('Unknown'),1, 0 ))





#table_alone_EVEs_event<-table_alone_EVEs_event[!duplicated(table_alone_EVEs_event[c('Species','Clusters')]),]

EVEs_alone_event_df_EVEs_count<- table_alone_EVEs_event


##################

Event_alone_df <- EVEs_alone_event_df_EVEs_count
Event_shared_df <- EVEs_shared_event_df_EVEs_count


#Correct issue on C vestalis alone
Event_alone_df$Nb_EVEs[ grepl( "3_Cluster24469", Event_alone_df$Event)]<-4

####ADDD the MAXIMUM bayesian names
library(treeio)
library(phytools)
#State_tree<- read.nhx("/Users/bguinet/Desktop/Papier_scientifique/ase_mk.tree")
#State_table<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Table_with_MA_states.txt",h=T,sep=";") #Change Maximum bayesian states

State_tree<- read.nhx("/Users/bguinet/Desktop/Papier_scientifique/Revbayes_anestral_states_estimates/output2/ase_mk.tree")
State_table<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Table_with_MA_states2.txt",h=T,sep=";") #Change Maximum bayesian states


findMRCA(as.phylo(State_table),c('Atta_colombica','Leptopilina_boulardi'))

Event_shared_df$lifecycle1<-'NA'
for (i in 1:nrow(Event_shared_df)){
  MRCA_number<-findMRCA(as.phylo(State_tree),unlist(strsplit(split=",",Event_shared_df$Events_species2[[i]])))
  print(MRCA_number)

  lifestyle<-State_table$anc_state_1[State_table$node==MRCA_number]

  lifestyle<-gsub("B","endoparasitoide",lifestyle)
  lifestyle<-gsub("A","freeliving",lifestyle)
  lifestyle<-gsub("C","ectoparasitoide",lifestyle)
  Event_shared_df$lifecycle1[i]<- lifestyle
}

#Lets gather EVEs and dEVES shared event resuts :


#Nb EVEs event/lifecycles
#dsDNA
#freeliving dsDNA
Nb_EVEs_dsDNA_fl_Event<-sum(Event_alone_df[which(Event_alone_df$lifecycle1=="freeliving"),]$Nb_EVEs_dsDNA_Events,na.rm=T)+sum(Event_shared_df[which(Event_shared_df$lifecycle1=="freeliving"),]$Nb_EVEs_dsDNA_Events,na.rm=T)
#ectoparasitoide dsDNA
Nb_EVEs_dsDNA_ecto_Event<-sum(Event_alone_df[which(Event_alone_df$lifecycle1=="ectoparasitoide"),]$Nb_EVEs_dsDNA_Events,na.rm=T)+sum(Event_shared_df[which(Event_shared_df$lifecycle1=="ectoparasitoide"),]$Nb_EVEs_dsDNA_Events,na.rm=T)
#endoparasitoide dsDNA
Nb_EVEs_dsDNA_endo_Event<-sum(Event_alone_df[which(Event_alone_df$lifecycle1=="endoparasitoide"),]$Nb_EVEs_dsDNA_Events,na.rm=T)+sum(Event_shared_df[which(Event_shared_df$lifecycle1=="endoparasitoide"),]$Nb_EVEs_dsDNA_Events,na.rm=T)
#Unknown dsDNA
Nb_EVEs_dsDNA_Unkn_Event<-sum(Event_alone_df[which(Event_alone_df$lifecycle1=="Unknown"),]$Nb_EVEs_dsDNA_Events,na.rm=T)+sum(Event_shared_df[which(Event_shared_df$lifecycle1=="Unknown"),]$Nb_EVEs_dsDNA_Events,na.rm=T)

#ssDNA
#freeliving ssDNA
Nb_EVEs_ssDNA_fl_Event<-sum(Event_alone_df[which(Event_alone_df$lifecycle1=="freeliving"),]$Nb_EVEs_ssDNA_Events,na.rm=T)+sum(Event_shared_df[which(Event_shared_df$lifecycle1=="freeliving"),]$Nb_EVEs_ssDNA_Events,na.rm=T)
#ectoparasitoide ssDNA
Nb_EVEs_ssDNA_ecto_Event<-sum(Event_alone_df[which(Event_alone_df$lifecycle1=="ectoparasitoide"),]$Nb_EVEs_ssDNA_Events,na.rm=T)+sum(Event_shared_df[which(Event_shared_df$lifecycle1=="ectoparasitoide"),]$Nb_EVEs_ssDNA_Events,na.rm=T)
#endoparasitoide ssDNA
Nb_EVEs_ssDNA_endo_Event<-sum(Event_alone_df[which(Event_alone_df$lifecycle1=="endoparasitoide"),]$Nb_EVEs_ssDNA_Events,na.rm=T)+sum(Event_shared_df[which(Event_shared_df$lifecycle1=="endoparasitoide"),]$Nb_EVEs_ssDNA_Events,na.rm=T)
#Unknown ssDNA
Nb_EVEs_ssDNA_Unkn_Event<-sum(Event_alone_df[which(Event_alone_df$lifecycle1=="Unknown"),]$Nb_EVEs_ssDNA_Events,na.rm=T)+sum(Event_shared_df[which(Event_shared_df$lifecycle1=="Unknown"),]$Nb_EVEs_ssDNA_Events,na.rm=T)

#ssRNA
#freeliving ssRNA
Nb_EVEs_ssRNA_fl_Event<-sum(Event_alone_df[which(Event_alone_df$lifecycle1=="freeliving"),]$Nb_EVEs_ssRNA_Events,na.rm=T)+sum(Event_shared_df[which(Event_shared_df$lifecycle1=="freeliving"),]$Nb_EVEs_ssRNA_Events,na.rm=T)
#ectoparasitoide ssRNA
Nb_EVEs_ssRNA_ecto_Event<-sum(Event_alone_df[which(Event_alone_df$lifecycle1=="ectoparasitoide"),]$Nb_EVEs_ssRNA_Events,na.rm=T)+sum(Event_shared_df[which(Event_shared_df$lifecycle1=="ectoparasitoide"),]$Nb_EVEs_ssRNA_Events,na.rm=T)
#endoparasitoide ssRNA
Nb_EVEs_ssRNA_endo_Event<-sum(Event_alone_df[which(Event_alone_df$lifecycle1=="endoparasitoide"),]$Nb_EVEs_ssRNA_Events,na.rm=T)+sum(Event_shared_df[which(Event_shared_df$lifecycle1=="endoparasitoide"),]$Nb_EVEs_ssRNA_Events,na.rm=T)
#Unknown ssRNA
Nb_EVEs_ssRNA_Unkn_Event<-sum(Event_alone_df[which(Event_alone_df$lifecycle1=="Unknown"),]$Nb_EVEs_ssRNA_Events,na.rm=T)+sum(Event_shared_df[which(Event_shared_df$lifecycle1=="Unknown"),]$Nb_EVEs_ssRNA_Events,na.rm=T)

#dsRNA
#freeliving dsRNA
Nb_EVEs_dsRNA_fl_Event<-sum(Event_alone_df[which(Event_alone_df$lifecycle1=="freeliving"),]$Nb_EVEs_dsRNA_Events,na.rm=T)+sum(Event_shared_df[which(Event_shared_df$lifecycle1=="freeliving"),]$Nb_EVEs_dsRNA_Events,na.rm=T)
#ectoparasitoide dsRNA
Nb_EVEs_dsRNA_ecto_Event<-sum(Event_alone_df[which(Event_alone_df$lifecycle1=="ectoparasitoide"),]$Nb_EVEs_dsRNA_Events,na.rm=T)+sum(Event_shared_df[which(Event_shared_df$lifecycle1=="ectoparasitoide"),]$Nb_EVEs_dsRNA_Events,na.rm=T)
#endoparasitoide dsRNA
Nb_EVEs_dsRNA_endo_Event<-sum(Event_alone_df[which(Event_alone_df$lifecycle1=="endoparasitoide"),]$Nb_EVEs_dsRNA_Events,na.rm=T)+sum(Event_shared_df[which(Event_shared_df$lifecycle1=="endoparasitoide"),]$Nb_EVEs_dsRNA_Events,na.rm=T)
#Unknown dsRNA
Nb_EVEs_dsRNA_Unkn_Event<-sum(Event_alone_df[which(Event_alone_df$lifecycle1=="Unknown"),]$Nb_EVEs_dsRNA_Events,na.rm=T)+sum(Event_shared_df[which(Event_shared_df$lifecycle1=="Unknown"),]$Nb_EVEs_dsRNA_Events,na.rm=T)

#Unclassified
#freeliving Unclassified
Nb_EVEs_Uncla_fl_Event<-sum(Event_alone_df[which(Event_alone_df$lifecycle1=="freeliving"),]$Nb_EVEs_Unclassified_Events,na.rm=T)+sum(Event_shared_df[which(Event_shared_df$lifecycle1=="freeliving"),]$Nb_EVEs_Unclassified_Events,na.rm=T)
#ectoparasitoide Unclassified
Nb_EVEs_Uncla_ecto_Event<-sum(Event_alone_df[which(Event_alone_df$lifecycle1=="ectoparasitoide"),]$Nb_EVEs_Unclassified_Events,na.rm=T)+sum(Event_shared_df[which(Event_shared_df$lifecycle1=="ectoparasitoide"),]$Nb_EVEs_Unclassified_Events,na.rm=T)
#endoparasitoide Unclassified
Nb_EVEs_Uncla_endo_Event<-sum(Event_alone_df[which(Event_alone_df$lifecycle1=="endoparasitoide"),]$Nb_EVEs_Unclassified_Events,na.rm=T)+sum(Event_shared_df[which(Event_shared_df$lifecycle1=="endoparasitoide"),]$Nb_EVEs_Unclassified_Events,na.rm=T)
#Unknown Unclassified
Nb_EVEs_Uncla_Unkn_Event<-sum(Event_alone_df[which(Event_alone_df$lifecycle1=="Unknown"),]$Nb_EVEs_Unclassified_Events,na.rm=T)+sum(Event_shared_df[which(Event_shared_df$lifecycle1=="Unknown"),]$Nb_EVEs_Unclassified_Events,na.rm=T)

#Nb dEVEs event/lifecycles
#dsDNA
#freeliving dsDNA
Nb_dEVEs_dsDNA_fl_Event<-sum(Event_alone_df[which(Event_alone_df$lifecycle1=="freeliving"),]$Nb_dEVEs_dsDNA_Events,na.rm=T)+sum(Event_shared_df[which(Event_shared_df$lifecycle1=="freeliving"),]$Nb_dEVEs_dsDNA_Events,na.rm=T)
#ectoparasitoide dsDNA
Nb_dEVEs_dsDNA_ecto_Event<-sum(Event_alone_df[which(Event_alone_df$lifecycle1=="ectoparasitoide"),]$Nb_dEVEs_dsDNA_Events,na.rm=T)+sum(Event_shared_df[which(Event_shared_df$lifecycle1=="ectoparasitoide"),]$Nb_dEVEs_dsDNA_Events,na.rm=T)
#endoparasitoide dsDNA
Nb_dEVEs_dsDNA_endo_Event<-sum(Event_alone_df[which(Event_alone_df$lifecycle1=="endoparasitoide"),]$Nb_dEVEs_dsDNA_Events,na.rm=T)+sum(Event_shared_df[which(Event_shared_df$lifecycle1=="endoparasitoide"),]$Nb_dEVEs_dsDNA_Events,na.rm=T)
#Unknown dsDNA
Nb_dEVEs_dsDNA_Unkn_Event<-sum(Event_alone_df[which(Event_alone_df$lifecycle1=="Unknown"),]$Nb_dEVEs_dsDNA_Events,na.rm=T)+sum(Event_shared_df[which(Event_shared_df$lifecycle1=="Unknown"),]$Nb_dEVEs_dsDNA_Events,na.rm=T)

#ssDNA
#freeliving ssDNA
Nb_dEVEs_ssDNA_fl_Event<-sum(Event_alone_df[which(Event_alone_df$lifecycle1=="freeliving"),]$Nb_dEVEs_ssDNA_Events,na.rm=T)+sum(Event_shared_df[which(Event_shared_df$lifecycle1=="freeliving"),]$Nb_dEVEs_ssDNA_Events,na.rm=T)
#ectoparasitoide ssDNA
Nb_dEVEs_ssDNA_ecto_Event<-sum(Event_alone_df[which(Event_alone_df$lifecycle1=="ectoparasitoide"),]$Nb_dEVEs_ssDNA_Events,na.rm=T)+sum(Event_shared_df[which(Event_shared_df$lifecycle1=="ectoparasitoide"),]$Nb_dEVEs_ssDNA_Events,na.rm=T)
#endoparasitoide ssDNA
Nb_dEVEs_ssDNA_endo_Event<-sum(Event_alone_df[which(Event_alone_df$lifecycle1=="endoparasitoide"),]$Nb_dEVEs_ssDNA_Events,na.rm=T)+sum(Event_shared_df[which(Event_shared_df$lifecycle1=="endoparasitoide"),]$Nb_dEVEs_ssDNA_Events,na.rm=T)
#Unknown ssDNA
Nb_dEVEs_ssDNA_Unkn_Event<-sum(Event_alone_df[which(Event_alone_df$lifecycle1=="Unknown"),]$Nb_dEVEs_ssDNA_Events,na.rm=T)+sum(Event_shared_df[which(Event_shared_df$lifecycle1=="Unknown"),]$Nb_dEVEs_ssDNA_Events,na.rm=T)

#ssRNA
#freeliving ssRNA
Nb_dEVEs_ssRNA_fl_Event<-sum(Event_alone_df[which(Event_alone_df$lifecycle1=="freeliving"),]$Nb_dEVEs_ssRNA_Events,na.rm=T)+sum(Event_shared_df[which(Event_shared_df$lifecycle1=="freeliving"),]$Nb_dEVEs_ssRNA_Events,na.rm=T)
#ectoparasitoide ssRNA
Nb_dEVEs_ssRNA_ecto_Event<-sum(Event_alone_df[which(Event_alone_df$lifecycle1=="ectoparasitoide"),]$Nb_dEVEs_ssRNA_Events,na.rm=T)+sum(Event_shared_df[which(Event_shared_df$lifecycle1=="ectoparasitoide"),]$Nb_dEVEs_ssRNA_Events,na.rm=T)
#endoparasitoide ssRNA
Nb_dEVEs_ssRNA_endo_Event<-sum(Event_alone_df[which(Event_alone_df$lifecycle1=="endoparasitoide"),]$Nb_dEVEs_ssRNA_Events,na.rm=T)+sum(Event_shared_df[which(Event_shared_df$lifecycle1=="endoparasitoide"),]$Nb_dEVEs_ssRNA_Events,na.rm=T)
#Unknown ssRNA
Nb_dEVEs_ssRNA_Unkn_Event<-sum(Event_alone_df[which(Event_alone_df$lifecycle1=="Unknown"),]$Nb_dEVEs_ssRNA_Events,na.rm=T)+sum(Event_shared_df[which(Event_shared_df$lifecycle1=="Unknown"),]$Nb_dEVEs_ssRNA_Events,na.rm=T)

#dsRNA
#freeliving dsRNA
Nb_dEVEs_dsRNA_fl_Event<-sum(Event_alone_df[which(Event_alone_df$lifecycle1=="freeliving"),]$Nb_dEVEs_dsRNA_Events,na.rm=T)+sum(Event_shared_df[which(Event_shared_df$lifecycle1=="freeliving"),]$Nb_dEVEs_dsRNA_Events,na.rm=T)
#ectoparasitoide dsRNA
Nb_dEVEs_dsRNA_ecto_Event<-sum(Event_alone_df[which(Event_alone_df$lifecycle1=="ectoparasitoide"),]$Nb_dEVEs_dsRNA_Events,na.rm=T)+sum(Event_shared_df[which(Event_shared_df$lifecycle1=="ectoparasitoide"),]$Nb_dEVEs_dsRNA_Events,na.rm=T)
#endoparasitoide dsRNA
Nb_dEVEs_dsRNA_endo_Event<-sum(Event_alone_df[which(Event_alone_df$lifecycle1=="endoparasitoide"),]$Nb_dEVEs_dsRNA_Events,na.rm=T)+sum(Event_shared_df[which(Event_shared_df$lifecycle1=="endoparasitoide"),]$Nb_dEVEs_dsRNA_Events,na.rm=T)
#Unknown dsRNA
Nb_dEVEs_dsRNA_Unkn_Event<-sum(Event_alone_df[which(Event_alone_df$lifecycle1=="Unknown"),]$Nb_dEVEs_dsRNA_Events,na.rm=T)+sum(Event_shared_df[which(Event_shared_df$lifecycle1=="Unknown"),]$Nb_dEVEs_dsRNA_Events,na.rm=T)

#Unclassified
#freeliving Unclassified
Nb_dEVEs_Uncla_fl_Event<-sum(Event_alone_df[which(Event_alone_df$lifecycle1=="freeliving"),]$Nb_dEVEs_Unclassified_Events,na.rm=T)+sum(Event_shared_df[which(Event_shared_df$lifecycle1=="freeliving"),]$Nb_dEVEs_Unclassified_Events,na.rm=T)
#ectoparasitoide Unclassified
Nb_dEVEs_Uncla_ecto_Event<-sum(Event_alone_df[which(Event_alone_df$lifecycle1=="ectoparasitoide"),]$Nb_dEVEs_Unclassified_Events,na.rm=T)+sum(Event_shared_df[which(Event_shared_df$lifecycle1=="ectoparasitoide"),]$Nb_dEVEs_Unclassified_Events,na.rm=T)
#endoparasitoide Unclassified
Nb_dEVEs_Uncla_endo_Event<-sum(Event_alone_df[which(Event_alone_df$lifecycle1=="endoparasitoide"),]$Nb_dEVEs_Unclassified_Events,na.rm=T)+sum(Event_shared_df[which(Event_shared_df$lifecycle1=="endoparasitoide"),]$Nb_dEVEs_Unclassified_Events,na.rm=T)
#Unknown Unclassified
Nb_dEVEs_Uncla_Unkn_Event<-sum(Event_alone_df[which(Event_alone_df$lifecycle1=="Unknown"),]$Nb_dEVEs_Unclassified_Events,na.rm=T)+sum(Event_shared_df[which(Event_shared_df$lifecycle1=="Unknown"),]$Nb_dEVEs_Unclassified_Events,na.rm=T)



###Count EVEs and create the geom_bar plot

#Remove duplicated EVEs within events to not overestimated the number of EVEs
dsDNA_tab_reduced <-dsDNA_tab[!duplicated(dsDNA_tab[ , c("Clustername","Event")]),]

#dsDNA
EVEs_endo_dsDNA<-nrow(dsDNA_tab_reduced[which(dsDNA_tab_reduced$lifecycle1=="endoparasitoide" & dsDNA_tab_reduced$consensus_genomic_structure =="dsDNA"),])
dEVEs_endo_dsDNA<-nrow(dsDNA_tab_reduced[which(dsDNA_tab_reduced$lifecycle1=="endoparasitoide" &  dsDNA_tab_reduced$consensus_genomic_structure =="dsDNA" & dsDNA_tab_reduced$FDR_pvalue_dNdS == 1 & as.numeric(as.character(dsDNA_tab_reduced$Mean_dNdS))&as.numeric(as.character(dsDNA_tab_reduced$SE_dNdS)) < 1 & dsDNA_tab_reduced$pseudogenized ==0 & dsDNA_tab_reduced$ORF_perc>1 | dsDNA_tab_reduced$lifecycle1=="endoparasitoide" &  dsDNA_tab_reduced$consensus_genomic_structure =="dsDNA" & dsDNA_tab_reduced$TPM_all>=1000 & dsDNA_tab_reduced$pseudogenized ==0 & dsDNA_tab_reduced$ORF_perc>1) ,])

EVEs_ecto_dsDNA<-nrow(dsDNA_tab_reduced[which(dsDNA_tab_reduced$lifecycle1=="ectoparasitoide" &  dsDNA_tab_reduced$consensus_genomic_structure =="dsDNA" ),])
dEVEs_ecto_dsDNA<-nrow(dsDNA_tab_reduced[which(dsDNA_tab_reduced$lifecycle1=="ectoparasitoide" &  dsDNA_tab_reduced$consensus_genomic_structure =="dsDNA" & dsDNA_tab_reduced$FDR_pvalue_dNdS == 1 & as.numeric(as.character(dsDNA_tab_reduced$Mean_dNdS))&as.numeric(as.character(dsDNA_tab_reduced$SE_dNdS)) < 1 & dsDNA_tab_reduced$pseudogenized ==0 & dsDNA_tab_reduced$ORF_perc>1 | dsDNA_tab_reduced$lifecycle1=="ectoparasitoide" &  dsDNA_tab_reduced$consensus_genomic_structure =="dsDNA" & dsDNA_tab_reduced$TPM_all>=1000 & dsDNA_tab_reduced$pseudogenized ==0 & dsDNA_tab_reduced$ORF_perc>1) ,])

EVEs_fl_dsDNA<-nrow(dsDNA_tab_reduced[which(dsDNA_tab_reduced$lifecycle1=="freeliving" &  dsDNA_tab_reduced$consensus_genomic_structure =="dsDNA" ),])
dEVEs_fl_dsDNA<-nrow(dsDNA_tab_reduced[which(dsDNA_tab_reduced$lifecycle1=="freeliving" &  dsDNA_tab_reduced$consensus_genomic_structure =="dsDNA" & dsDNA_tab_reduced$FDR_pvalue_dNdS == 1 & as.numeric(as.character(dsDNA_tab_reduced$Mean_dNdS))&as.numeric(as.character(dsDNA_tab_reduced$SE_dNdS)) < 1 & dsDNA_tab_reduced$pseudogenized ==0 & dsDNA_tab_reduced$ORF_perc>1 | dsDNA_tab_reduced$lifecycle1=="freeliving" &  dsDNA_tab_reduced$consensus_genomic_structure =="dsDNA" & dsDNA_tab_reduced$TPM_all>=1000& dsDNA_tab_reduced$pseudogenized ==0 & dsDNA_tab_reduced$ORF_perc>1) ,])


#ssDNA
EVEs_endo_ssDNA<-nrow(dsDNA_tab_reduced[which(dsDNA_tab_reduced$lifecycle1=="endoparasitoide" & dsDNA_tab_reduced$consensus_genomic_structure =="ssDNA"),])
dEVEs_endo_ssDNA<-nrow(dsDNA_tab_reduced[which(dsDNA_tab_reduced$lifecycle1=="endoparasitoide" &  dsDNA_tab_reduced$consensus_genomic_structure =="ssDNA" & dsDNA_tab_reduced$FDR_pvalue_dNdS == 1 & as.numeric(as.character(dsDNA_tab_reduced$Mean_dNdS))&as.numeric(as.character(dsDNA_tab_reduced$SE_dNdS)) < 1 & dsDNA_tab_reduced$pseudogenized ==0 & dsDNA_tab_reduced$ORF_perc>1 | dsDNA_tab_reduced$lifecycle1=="endoparasitoide" &  dsDNA_tab_reduced$consensus_genomic_structure =="ssDNA" & dsDNA_tab_reduced$TPM_all>=1000& dsDNA_tab_reduced$pseudogenized ==0 & dsDNA_tab_reduced$ORF_perc>1) ,])

EVEs_ecto_ssDNA<-nrow(dsDNA_tab_reduced[which(dsDNA_tab_reduced$lifecycle1=="ectoparasitoide" &  dsDNA_tab_reduced$consensus_genomic_structure =="ssDNA" ),])
dEVEs_ecto_ssDNA<-nrow(dsDNA_tab_reduced[which(dsDNA_tab_reduced$lifecycle1=="ectoparasitoide" &  dsDNA_tab_reduced$consensus_genomic_structure =="ssDNA" & dsDNA_tab_reduced$FDR_pvalue_dNdS == 1 & as.numeric(as.character(dsDNA_tab_reduced$Mean_dNdS))&as.numeric(as.character(dsDNA_tab_reduced$SE_dNdS)) < 1 & dsDNA_tab_reduced$pseudogenized ==0 & dsDNA_tab_reduced$ORF_perc>1 | dsDNA_tab_reduced$lifecycle1=="ectoparasitoide" &  dsDNA_tab_reduced$consensus_genomic_structure =="ssDNA" & dsDNA_tab_reduced$TPM_all>=1000& dsDNA_tab_reduced$pseudogenized ==0 & dsDNA_tab_reduced$ORF_perc>1) ,])

EVEs_fl_ssDNA<-nrow(dsDNA_tab_reduced[which(dsDNA_tab_reduced$lifecycle1=="freeliving" &  dsDNA_tab_reduced$consensus_genomic_structure =="ssDNA" ),])
dEVEs_fl_ssDNA<-nrow(dsDNA_tab_reduced[which(dsDNA_tab_reduced$lifecycle1=="freeliving" &  dsDNA_tab_reduced$consensus_genomic_structure =="ssDNA" & dsDNA_tab_reduced$FDR_pvalue_dNdS == 1 & as.numeric(as.character(dsDNA_tab_reduced$Mean_dNdS))&as.numeric(as.character(dsDNA_tab_reduced$SE_dNdS)) < 1 & dsDNA_tab_reduced$pseudogenized ==0 & dsDNA_tab_reduced$ORF_perc>1 | dsDNA_tab_reduced$lifecycle1=="freeliving" &  dsDNA_tab_reduced$consensus_genomic_structure =="ssDNA" & dsDNA_tab_reduced$TPM_all>=1000& dsDNA_tab_reduced$pseudogenized ==0 & dsDNA_tab_reduced$ORF_perc>1) ,])



#dsRNA
EVEs_endo_dsRNA<-nrow(dsDNA_tab_reduced[which(dsDNA_tab_reduced$lifecycle1=="endoparasitoide" & dsDNA_tab_reduced$consensus_genomic_structure =="dsRNA"),])
dEVEs_endo_dsRNA<-nrow(dsDNA_tab_reduced[which(dsDNA_tab_reduced$lifecycle1=="endoparasitoide" &  dsDNA_tab_reduced$consensus_genomic_structure =="dsRNA" & dsDNA_tab_reduced$FDR_pvalue_dNdS == 1 & as.numeric(as.character(dsDNA_tab_reduced$Mean_dNdS))&as.numeric(as.character(dsDNA_tab_reduced$SE_dNdS)) < 1 & dsDNA_tab_reduced$pseudogenized ==0 & dsDNA_tab_reduced$ORF_perc>1 | dsDNA_tab_reduced$lifecycle1=="endoparasitoide" &  dsDNA_tab_reduced$consensus_genomic_structure =="dsRNA" & dsDNA_tab_reduced$TPM_all>=1000& dsDNA_tab_reduced$pseudogenized ==0 & dsDNA_tab_reduced$ORF_perc>1 ) ,])

EVEs_ecto_dsRNA<-nrow(dsDNA_tab_reduced[which(dsDNA_tab_reduced$lifecycle1=="ectoparasitoide" &  dsDNA_tab_reduced$consensus_genomic_structure =="dsRNA" ),])
dEVEs_ecto_dsRNA<-nrow(dsDNA_tab_reduced[which(dsDNA_tab_reduced$lifecycle1=="ectoparasitoide" &  dsDNA_tab_reduced$consensus_genomic_structure =="dsRNA" & dsDNA_tab_reduced$FDR_pvalue_dNdS == 1 & as.numeric(as.character(dsDNA_tab_reduced$Mean_dNdS))&as.numeric(as.character(dsDNA_tab_reduced$SE_dNdS)) < 1 & dsDNA_tab_reduced$pseudogenized ==0 & dsDNA_tab_reduced$ORF_perc>1 | dsDNA_tab_reduced$lifecycle1=="ectoparasitoide" &  dsDNA_tab_reduced$consensus_genomic_structure =="dsRNA" & dsDNA_tab_reduced$TPM_all>=1000& dsDNA_tab_reduced$pseudogenized ==0 & dsDNA_tab_reduced$ORF_perc>1 ) ,])

EVEs_fl_dsRNA<-nrow(dsDNA_tab_reduced[which(dsDNA_tab_reduced$lifecycle1=="freeliving" &  dsDNA_tab_reduced$consensus_genomic_structure =="dsRNA" ),])
dEVEs_fl_dsRNA<-nrow(dsDNA_tab_reduced[which(dsDNA_tab_reduced$lifecycle1=="freeliving" &  dsDNA_tab_reduced$consensus_genomic_structure =="dsRNA" & dsDNA_tab_reduced$FDR_pvalue_dNdS == 1 & as.numeric(as.character(dsDNA_tab_reduced$Mean_dNdS))&as.numeric(as.character(dsDNA_tab_reduced$SE_dNdS)) < 1 & dsDNA_tab_reduced$pseudogenized ==0 & dsDNA_tab_reduced$ORF_perc>1 | dsDNA_tab_reduced$lifecycle1=="freeliving" &  dsDNA_tab_reduced$consensus_genomic_structure =="dsRNA" & dsDNA_tab_reduced$TPM_all>=1000& dsDNA_tab_reduced$pseudogenized ==0 & dsDNA_tab_reduced$ORF_perc>1 ) ,])



#ssRNA
EVEs_endo_ssRNA<-nrow(dsDNA_tab_reduced[which(dsDNA_tab_reduced$lifecycle1=="endoparasitoide" & dsDNA_tab_reduced$consensus_genomic_structure =="ssRNA"),])
dEVEs_endo_ssRNA<-nrow(dsDNA_tab_reduced[which(dsDNA_tab_reduced$lifecycle1=="endoparasitoide" &  dsDNA_tab_reduced$consensus_genomic_structure =="ssRNA" & dsDNA_tab_reduced$FDR_pvalue_dNdS == 1 & as.numeric(as.character(dsDNA_tab_reduced$Mean_dNdS))&as.numeric(as.character(dsDNA_tab_reduced$SE_dNdS)) < 1 & dsDNA_tab_reduced$pseudogenized ==0 & dsDNA_tab_reduced$ORF_perc>1 | dsDNA_tab_reduced$lifecycle1=="endoparasitoide" &  dsDNA_tab_reduced$consensus_genomic_structure =="ssRNA" & dsDNA_tab_reduced$TPM_all>=1000& dsDNA_tab_reduced$pseudogenized ==0& dsDNA_tab_reduced$ORF_perc>1) ,])

EVEs_ecto_ssRNA<-nrow(dsDNA_tab_reduced[which(dsDNA_tab_reduced$lifecycle1=="ectoparasitoide" &  dsDNA_tab_reduced$consensus_genomic_structure =="ssRNA" ),])
dEVEs_ecto_ssRNA<-nrow(dsDNA_tab_reduced[which(dsDNA_tab_reduced$lifecycle1=="ectoparasitoide" &  dsDNA_tab_reduced$consensus_genomic_structure =="ssRNA" & dsDNA_tab_reduced$FDR_pvalue_dNdS == 1 & as.numeric(as.character(dsDNA_tab_reduced$Mean_dNdS))&as.numeric(as.character(dsDNA_tab_reduced$SE_dNdS)) < 1 & dsDNA_tab_reduced$pseudogenized ==0 & dsDNA_tab_reduced$ORF_perc>1 | dsDNA_tab_reduced$lifecycle1=="ectoparasitoide" &  dsDNA_tab_reduced$consensus_genomic_structure =="ssRNA" & dsDNA_tab_reduced$TPM_all>=1000& dsDNA_tab_reduced$pseudogenized ==0 & dsDNA_tab_reduced$ORF_perc>1 ) ,])

EVEs_fl_ssRNA<-nrow(dsDNA_tab_reduced[which(dsDNA_tab_reduced$lifecycle1=="freeliving" &  dsDNA_tab_reduced$consensus_genomic_structure =="ssRNA" ),])
dEVEs_fl_ssRNA<-nrow(dsDNA_tab_reduced[which(dsDNA_tab_reduced$lifecycle1=="freeliving" &  dsDNA_tab_reduced$consensus_genomic_structure =="ssRNA" & dsDNA_tab_reduced$FDR_pvalue_dNdS == 1 & as.numeric(as.character(dsDNA_tab_reduced$Mean_dNdS))&as.numeric(as.character(dsDNA_tab_reduced$SE_dNdS)) < 1 & dsDNA_tab_reduced$pseudogenized ==0 & dsDNA_tab_reduced$ORF_perc>1 | dsDNA_tab_reduced$lifecycle1=="freeliving" &  dsDNA_tab_reduced$consensus_genomic_structure =="ssRNA" & dsDNA_tab_reduced$TPM_all>=1000& dsDNA_tab_reduced$pseudogenized ==0 & dsDNA_tab_reduced$ORF_perc>1 ) ,])


#Unclassified
EVEs_endo_Unclassified<-nrow(dsDNA_tab_reduced[which(dsDNA_tab_reduced$lifecycle1=="endoparasitoide" & dsDNA_tab_reduced$consensus_genomic_structure =="Unknown"),])
dEVEs_endo_Unclassified<-nrow(dsDNA_tab_reduced[which(dsDNA_tab_reduced$lifecycle1=="endoparasitoide" &  dsDNA_tab_reduced$consensus_genomic_structure =="Unknown" & dsDNA_tab_reduced$FDR_pvalue_dNdS == 1 & as.numeric(as.character(dsDNA_tab_reduced$Mean_dNdS))&as.numeric(as.character(dsDNA_tab_reduced$SE_dNdS)) < 1 & dsDNA_tab_reduced$pseudogenized ==0 & dsDNA_tab_reduced$ORF_perc>1  | dsDNA_tab_reduced$lifecycle1=="endoparasitoide" &  dsDNA_tab_reduced$consensus_genomic_structure =="Unknown" & dsDNA_tab_reduced$TPM_all>=1000& dsDNA_tab_reduced$pseudogenized ==0 & dsDNA_tab_reduced$ORF_perc>1 ) ,])

EVEs_ecto_Unclassified<-nrow(dsDNA_tab_reduced[which(dsDNA_tab_reduced$lifecycle1=="ectoparasitoide" &  dsDNA_tab_reduced$consensus_genomic_structure =="Unknown" ),])
dEVEs_ecto_Unclassified<-nrow(dsDNA_tab_reduced[which(dsDNA_tab_reduced$lifecycle1=="ectoparasitoide" &  dsDNA_tab_reduced$consensus_genomic_structure =="Unknown" & dsDNA_tab_reduced$FDR_pvalue_dNdS == 1 & as.numeric(as.character(dsDNA_tab_reduced$Mean_dNdS))&as.numeric(as.character(dsDNA_tab_reduced$SE_dNdS)) < 1 & dsDNA_tab_reduced$pseudogenized ==0 & dsDNA_tab_reduced$ORF_perc>1 | dsDNA_tab_reduced$lifecycle1=="ectoparasitoide" &  dsDNA_tab_reduced$consensus_genomic_structure =="Unknown" & dsDNA_tab_reduced$TPM_all>=1000& dsDNA_tab_reduced$pseudogenized ==0 & dsDNA_tab_reduced$ORF_perc>1 ) ,])

EVEs_fl_Unclassified<-nrow(dsDNA_tab_reduced[which(dsDNA_tab_reduced$lifecycle1=="freeliving" &  dsDNA_tab_reduced$consensus_genomic_structure =="Unknown" ),])
dEVEs_fl_Unclassified<-nrow(dsDNA_tab_reduced[which(dsDNA_tab_reduced$lifecycle1=="freeliving" &  dsDNA_tab_reduced$consensus_genomic_structure =="Unknown" & dsDNA_tab_reduced$FDR_pvalue_dNdS == 1 & as.numeric(as.character(dsDNA_tab_reduced$Mean_dNdS))&as.numeric(as.character(dsDNA_tab_reduced$SE_dNdS)) < 1 & dsDNA_tab_reduced$pseudogenized ==0 & dsDNA_tab_reduced$ORF_perc>1 | dsDNA_tab_reduced$lifecycle1=="freeliving" &  dsDNA_tab_reduced$consensus_genomic_structure =="Unknown" & dsDNA_tab_reduced$TPM_all>=1000& dsDNA_tab_reduced$pseudogenized ==0& dsDNA_tab_reduced$ORF_perc>1 ) ,])




col.names = c("dsDNA", 'ssRNA','dsRNA',"ssDNA",'Category','Type')
Summary_table_count<- read.table(text = "",
                                 col.names = col.names)

Summary_table_count$dsDNA<-as.character(Summary_table_count$dsDNA)
Summary_table_count$ssDNA<-as.character(Summary_table_count$ssDNA)
Summary_table_count$dsRNA<-as.character(Summary_table_count$dsRNA)
Summary_table_count$ssRNA<-as.character(Summary_table_count$ssRNA)
Summary_table_count$Category<-as.character(Summary_table_count$Category)
Summary_table_count$Type<-as.character(Summary_table_count$Type)



Summary_table_count<-Summary_table_count %>%
  summarise(dsDNA = as.character(EVEs_endo_dsDNA),
            ssDNA = as.character(EVEs_endo_ssDNA),
            dsRNA = as.character(EVEs_endo_dsRNA),
            ssRNA = as.character(EVEs_endo_ssRNA),
            Type="Endoparasitoid",
            Category="EVEs")%>%
  bind_rows(Summary_table_count)

Summary_table_count<-Summary_table_count %>%
  summarise(dsDNA = as.character(dEVEs_endo_dsDNA),
            ssDNA = as.character(dEVEs_endo_ssDNA),
            dsRNA = as.character(dEVEs_endo_dsRNA),
            ssRNA = as.character(dEVEs_endo_ssRNA),
            Type="Endoparasitoid",
            Category="dEVEs")%>%
  bind_rows(Summary_table_count)


Summary_table_count<-Summary_table_count %>%
  summarise(dsDNA = as.character(EVEs_ecto_dsDNA),
            ssDNA = as.character(EVEs_ecto_ssDNA),
            dsRNA = as.character(EVEs_ecto_dsRNA),
            ssRNA = as.character(EVEs_ecto_ssRNA),
            Type="Ectoparasitoid",
            Category="EVEs")%>%
  bind_rows(Summary_table_count)

Summary_table_count<-Summary_table_count %>%
  summarise(dsDNA = as.character(dEVEs_ecto_dsDNA),
            ssDNA = as.character(dEVEs_ecto_ssDNA),
            dsRNA = as.character(dEVEs_ecto_dsRNA),
            ssRNA = as.character(dEVEs_ecto_ssRNA),
            Type="Ectoparasitoid",
            Category="dEVEs")%>%
  bind_rows(Summary_table_count)

Summary_table_count<-Summary_table_count %>%
  summarise(dsDNA = as.character(EVEs_fl_dsDNA),
            ssDNA = as.character(EVEs_fl_ssDNA),
            dsRNA = as.character(EVEs_fl_dsRNA),
            ssRNA = as.character(EVEs_fl_ssRNA),
            Type="Free-living",
            Category="EVEs")%>%
  bind_rows(Summary_table_count)

Summary_table_count<-Summary_table_count %>%
  summarise(dsDNA = as.character(dEVEs_fl_dsDNA),
            ssDNA = as.character(dEVEs_fl_ssDNA),
            dsRNA = as.character(dEVEs_fl_dsRNA),
            ssRNA = as.character(dEVEs_fl_ssRNA),
            Type="Free-living",
            Category="dEVEs")%>%
  bind_rows(Summary_table_count)

Summary_table_count<-Summary_table_count %>%
  summarise(dsDNA = as.character(dEVEs_fl_dsDNA),
            ssDNA = as.character(dEVEs_fl_ssDNA),
            dsRNA = as.character(dEVEs_fl_dsRNA),
            ssRNA = as.character(dEVEs_fl_ssRNA),
            Type="Free-living",
            Category="dEVEs")%>%
  bind_rows(Summary_table_count)

tab <- as.data.frame(melt(setDT(Summary_table_count), id.vars = c("Category", "Type"), variable.name = "C"))

tab<-tab [order(tab $Category, decreasing = F),]

tab$Category<- factor(tab$Category, levels = c("EVEs","dEVEs"))
tab$Type<- factor(tab$Type, levels = c("Endoparasitoid","Ectoparasitoid","Free-living"))

tab$Category <- as.factor(tab$Category)
tab$Type <- as.factor(tab$Type)
tab$C <- as.factor(tab$C)
tab$value <- as.integer(tab$value)
tab<-tab[!duplicated(tab), ]
colnames(tab)<-c("Category","Lifecycle","Genomic_structure","Nb_EVEs")

tab$ponderate<-NA
tab$ponderate<-as.character(tab$ponderate)
tab$Lifecycle<-as.character(tab$Lifecycle)
tab<-tab %>% mutate(ponderate= case_when(grepl('Free',Lifecycle) ~ '67',TRUE ~ ponderate ))
tab<-tab %>% mutate(ponderate= case_when(grepl('Endo',Lifecycle) ~ '38',TRUE ~ ponderate ))
tab<-tab %>% mutate(ponderate= case_when(grepl('Ecto',Lifecycle) ~ '19',TRUE ~ ponderate ))
tab$ponderate<-as.numeric(tab$ponderate)

tab$Nb_EVEs2 <- tab$Nb_EVEs/tab$ponderate*100
tab$Nb_EVEs2 <- as.numeric(tab$Nb_EVEs2 )

library(ggrepel)


tab$Nb_EVEs2<-round(tab$Nb_EVEs2,digits=1)
ggplot(tab, aes(x = Lifecycle, y = Nb_EVEs, fill = Category)) +
  facet_wrap( ~ Genomic_structure,nrow=1,strip.position = "bottom")+
  geom_col(data = . %>% filter( Category=="dEVEs" & Lifecycle =="Free-living"), position = position_dodge(width = 0.9), alpha = 1,fill="#467fb3",color="black") +
  geom_col(data = . %>% filter( Category=="EVEs"& Lifecycle =="Free-living"), position = position_dodge(width = 0.9), alpha = 0.4,fill="#467fb3",color="black") +
  geom_col(data = . %>% filter( Category=="dEVEs" & Lifecycle =="Endoparasitoid"), position = position_dodge(width = 0.9), alpha = 1,fill="#06642e",color="black") +
  geom_col(data = . %>% filter( Category=="EVEs"& Lifecycle =="Endoparasitoid"), position = position_dodge(width = 0.9), alpha = 0.4,fill="#06642e",color="black") +
  geom_col(data = . %>% filter( Category=="dEVEs" & Lifecycle =="Ectoparasitoid"), position = position_dodge(width = 0.9), alpha = 1,fill="#d8be03",color="black") +
  geom_col(data = . %>% filter( Category=="EVEs"& Lifecycle =="Ectoparasitoid"), position = position_dodge(width = 0.9), alpha = 0.4,fill="#d8be03",color="black") +
  geom_label_repel(aes(label = Nb_EVEs),position = position_stack(vjust = 0.5),size = 2.75,colour = 'black',fill="white")+
  theme_classic() +
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        panel.spacing.y = unit(2, "lines"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + ylim(0, 1400) +  ylab("Number of viral endogenized elements")


####

#Event plot count
col.names = c("dsDNA", 'ssRNA','dsRNA',"ssDNA",'Category','Type')
Summary_table_count_Event<- read.table(text = "",
                                       col.names = col.names)

Summary_table_count_Event$dsDNA<-as.character(Summary_table_count_Event$dsDNA)
Summary_table_count_Event$ssDNA<-as.character(Summary_table_count_Event$ssDNA)
Summary_table_count_Event$dsRNA<-as.character(Summary_table_count_Event$dsRNA)
Summary_table_count_Event$ssRNA<-as.character(Summary_table_count_Event$ssRNA)
Summary_table_count_Event$Category<-as.character(Summary_table_count_Event$Category)
Summary_table_count_Event$Type<-as.character(Summary_table_count_Event$Type)


Summary_table_count_Event<-Summary_table_count_Event %>%
  summarise(dsDNA = as.character(Nb_EVEs_dsDNA_fl_Event),
            ssDNA = as.character(Nb_EVEs_ssDNA_fl_Event),
            dsRNA = as.character(Nb_EVEs_dsRNA_fl_Event),
            ssRNA = as.character(Nb_EVEs_ssRNA_fl_Event),
            Type="Free-living",
            Category="EVEs")%>%
  bind_rows(Summary_table_count_Event)

Summary_table_count_Event<-Summary_table_count_Event %>%
  summarise(dsDNA = as.character(Nb_dEVEs_dsDNA_fl_Event),
            ssDNA = as.character(Nb_dEVEs_ssDNA_fl_Event),
            dsRNA = as.character(Nb_dEVEs_dsRNA_fl_Event),
            ssRNA = as.character(Nb_dEVEs_ssRNA_fl_Event),
            Type="Free-living",
            Category="dEVEs")%>%
  bind_rows(Summary_table_count_Event)


Summary_table_count_Event<-Summary_table_count_Event %>%
  summarise(dsDNA = as.character(Nb_EVEs_dsDNA_ecto_Event),
            ssDNA = as.character(Nb_EVEs_ssDNA_ecto_Event),
            dsRNA = as.character(Nb_EVEs_dsRNA_ecto_Event),
            ssRNA = as.character(Nb_EVEs_ssRNA_ecto_Event),
            Type="Ectoparasitoid",
            Category="EVEs")%>%
  bind_rows(Summary_table_count_Event)

Summary_table_count_Event<-Summary_table_count_Event %>%
  summarise(dsDNA = as.character(Nb_dEVEs_dsDNA_ecto_Event),
            ssDNA = as.character(Nb_dEVEs_ssDNA_ecto_Event),
            dsRNA = as.character(Nb_dEVEs_dsRNA_ecto_Event),
            ssRNA = as.character(Nb_dEVEs_ssRNA_ecto_Event),
            Type="Ectoparasitoid",
            Category="dEVEs")%>%
  bind_rows(Summary_table_count_Event)



Summary_table_count_Event<-Summary_table_count_Event %>%
  summarise(dsDNA = as.character(Nb_EVEs_dsDNA_endo_Event),
            ssDNA = as.character(Nb_EVEs_ssDNA_endo_Event),
            dsRNA = as.character(Nb_EVEs_dsRNA_endo_Event),
            ssRNA = as.character(Nb_EVEs_ssRNA_endo_Event),
            Type="Endoparasitoid",
            Category="EVEs")%>%
  bind_rows(Summary_table_count_Event)

Summary_table_count_Event<-Summary_table_count_Event %>%
  summarise(dsDNA = as.character(Nb_dEVEs_dsDNA_endo_Event),
            ssDNA = as.character(Nb_dEVEs_ssDNA_endo_Event),
            dsRNA = as.character(Nb_dEVEs_dsRNA_endo_Event),
            ssRNA = as.character(Nb_dEVEs_ssRNA_endo_Event),
            Type="Endoparasitoid",
            Category="dEVEs")%>%
  bind_rows(Summary_table_count_Event)

library(reshape2)
library(data.table)
tab <- as.data.frame(melt(setDT(Summary_table_count_Event), id.vars = c("Category", "Type"), variable.name = "C"))

tab<-tab [order(tab $Category, decreasing = F),]

tab$Category<- factor(tab$Category, levels = c("EVEs","dEVEs"))
tab$Type<- factor(tab$Type, levels = c("Endoparasitoid","Ectoparasitoid","Free-living"))

tab$Category <- as.factor(tab$Category)
tab$Type <- as.factor(tab$Type)
tab$C <- as.factor(tab$C)
tab$value <- as.integer(tab$value)
tab<-tab[!duplicated(tab), ]
colnames(tab)<-c("Category","Lifecycle","Genomic_structure","Nb_EVEs")

tab$sumCat[tab$Genomic_structure =="ssDNA"] <- sum(tab$Nb_EVEs[tab$Genomic_structure =="ssDNA" & tab$Category=="EVEs"])
tab$sumCat[tab$Genomic_structure =="dsDNA"] <- sum(tab$Nb_EVEs[tab$Genomic_structure =="dsDNA" & tab$Category=="EVEs"])
tab$sumCat[tab$Genomic_structure =="ssRNA"] <- sum(tab$Nb_EVEs[tab$Genomic_structure =="ssRNA" & tab$Category=="EVEs"])
tab$sumCat[tab$Genomic_structure =="dsRNA"] <- sum(tab$Nb_EVEs[tab$Genomic_structure =="dsRNA" & tab$Category=="EVEs"])

tab$sumCat2[tab$Genomic_structure =="ssDNA"] <- sum(tab$Nb_EVEs[tab$Genomic_structure =="ssDNA" & tab$Category=="dEVEs"])
tab$sumCat2[tab$Genomic_structure =="dsDNA"] <- sum(tab$Nb_EVEs[tab$Genomic_structure =="dsDNA" & tab$Category=="dEVEs"])
tab$sumCat2[tab$Genomic_structure =="ssRNA"] <- sum(tab$Nb_EVEs[tab$Genomic_structure =="ssRNA" & tab$Category=="dEVEs"])
tab$sumCat2[tab$Genomic_structure =="dsRNA"] <- sum(tab$Nb_EVEs[tab$Genomic_structure =="dsRNA" & tab$Category=="dEVEs"])

tab$Nb_EVEs2 <- tab$Nb_EVEs*100/tab$sumCat
tab$Nb_dEVEs2 <- tab$Nb_EVEs*100/tab$sumCat2
tab$Expected[tab$Lifecycle=="Free-living"]<- (65/124)*100
tab$Expected[tab$Lifecycle=="Endoparasitoid"]<- (37/124)*100
tab$Expected[tab$Lifecycle=="Ectoparasitoid"]<- (24/124)*100

library(ggrepel)


subtab1 <- tab[tab$Category == "EVEs",]
EVEs_Events_representativity <- ggplot(subtab1 , aes(x = Lifecycle, y = Nb_EVEs2, fill = Category)) +
  facet_wrap( ~ Genomic_structure,nrow=1,strip.position = "top")+
  geom_col(data = . %>% filter( Category=="dEVEs" & Lifecycle =="Free-living"), position = position_dodge(width = 0.9), alpha = 1,fill='#467fb3',color="black") +
  geom_col(data = . %>% filter( Category=="EVEs"& Lifecycle =="Free-living"), position = position_dodge(width = 0.9), alpha = 1,fill='#467fb3',color="black") +
  geom_col(data = . %>% filter( Category=="dEVEs" & Lifecycle =="Endoparasitoid"), position = position_dodge(width = 0.9), alpha = 1,fill='#06642e',color="black") +
  geom_col(data = . %>% filter( Category=="EVEs"& Lifecycle =="Endoparasitoid"), position = position_dodge(width = 0.9), alpha = 1,fill='#06642e',color="black") +
  geom_col(data = . %>% filter( Category=="dEVEs" & Lifecycle =="Ectoparasitoid"), position = position_dodge(width = 0.9), alpha = 1,fill='#d8be03',color="black") +
  geom_col(data = . %>% filter( Category=="EVEs"& Lifecycle =="Ectoparasitoid"), position = position_dodge(width = 0.9), alpha = 1,fill='#d8be03',color="black") +
  geom_point(size = 1.5,aes(y=Expected, shape="Exposure",fill="black"),shape=3,stroke = 1)+
  #geom_label_repel(aes(label = Nb_EVEs),position = position_stack(vjust = 0.5),size = 2.75,colour = 'black',fill="white")+
  geom_text(aes(label = Nb_EVEs,y=3),vjust = 0.1, nudge_y = .2,color="white")+
  ylab("EVEs %")+
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line())+
  theme(strip.text.x = element_text(size = 13,face="bold"))+
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        panel.spacing.y = unit(2, "lines"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+ ylim(0, 100)+ theme(legend.position = "none")+
  scale_color_manual(values=c("#d8be03","#06642e","#467fb3"))


subtab2 <- tab[tab$Category == "dEVEs",]
dEVEs_Events_representativity <- ggplot(subtab2 , aes(x = Lifecycle, y = Nb_dEVEs2, fill = Category)) +
  facet_wrap( ~ Genomic_structure,nrow=1,strip.position = "top")+
  geom_col(data = . %>% filter( Category=="dEVEs" & Lifecycle =="Free-living"), position = position_dodge(width = 0.9), alpha = 1,fill='#467fb3',color="black") +
  geom_col(data = . %>% filter( Category=="EVEs"& Lifecycle =="Free-living"), position = position_dodge(width = 0.9), alpha = 0.4,fill='#467fb3',color="black") +
  geom_col(data = . %>% filter( Category=="dEVEs" & Lifecycle =="Endoparasitoid"), position = position_dodge(width = 0.9), alpha = 1,fill='#06642e',color="black") +
  geom_col(data = . %>% filter( Category=="EVEs"& Lifecycle =="Endoparasitoid"), position = position_dodge(width = 0.9), alpha = 0.4,fill='#06642e',color="black") +
  geom_col(data = . %>% filter( Category=="dEVEs" & Lifecycle =="Ectoparasitoid"), position = position_dodge(width = 0.9), alpha = 1,fill='#d8be03',color="black") +
  geom_col(data = . %>% filter( Category=="EVEs"& Lifecycle =="Ectoparasitoid"), position = position_dodge(width = 0.9), alpha = 0.4,fill='#d8be03',color="black") +
  geom_point(size = 1.5,aes(y=Expected, shape="Exposure",fill="black"),shape=3,stroke = 1)+
  #geom_label_repel(aes(label = Nb_EVEs),position = position_stack(vjust = 0.5),size = 2.75,colour = 'black',fill="white")+
  geom_text(aes(label = Nb_EVEs,y=2),vjust = 0.1, nudge_y = .2,color="white") + ylim(0, 100) +
  ylab("dEVEs %")+
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line())+
  theme(strip.text.x = element_text(size = 13,face="bold"))+
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        panel.spacing.y = unit(2, "lines"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+ theme(legend.position = "none")+
  scale_color_manual(values=c("#d8be03","#06642e","#467fb3"))

library(ggpubr)
library(cowplot)
plot_grid(EVEs_Events_representativity ,dEVEs_Events_representativity,ncol=2,nrow=1,labels=LETTERS[1:2],label_size = 16)

#####

col.names = c("dsDNA", 'ssRNA','dsRNA',"ssDNA",'Unclassified','Total','Type')
Summary_table<- read.table(text = "",
                           col.names = col.names)

Summary_table$dsDNA<-as.character(Summary_table$dsDNA)
Summary_table$ssDNA<-as.character(Summary_table$ssDNA)
Summary_table$dsRNA<-as.character(Summary_table$dsRNA)
Summary_table$ssRNA<-as.character(Summary_table$ssRNA)
Summary_table$Unclassified<-as.character(Summary_table$Unclassified)
Summary_table$Total<-as.character(Summary_table$Total)
Summary_table$Type<-as.character(Summary_table$Type)


#Events_dEVEs_dsDNA <- sum(Event_alone_df$Nb_dEVEs_dsDNA_Events,na.rm=T)+sum(Event_shared_df$Nb_dEVEs_dsDNA_Events,na.rm=T)
#Events_dEVEs_ssDNA <- sum(Event_alone_df$Nb_dEVEs_ssDNA_Events,na.rm=T)+sum(Event_shared_df$Nb_dEVEs_ssDNA_Events,na.rm=T)
#Events_dEVEs_dsRNA <- sum(Event_alone_df$Nb_dEVEs_dsRNA_Events,na.rm=T)+sum(Event_shared_df$Nb_dEVEs_dsRNA_Events,na.rm=T)
#Events_dEVEs_ssRNA <- sum(Event_alone_df$Nb_dEVEs_ssRNA_Events,na.rm=T)+sum(Event_shared_df$Nb_dEVEs_ssRNA_Events,na.rm=T)
Events_dEVEs_Unclassified <- sum(Event_alone_df$Nb_dEVEs_Unclassified_Events,na.rm=T)+sum(Event_shared_df$Nb_dEVEs_Unclassified_Events,na.rm=T)




Expected_NB<- c( ((sum(Nb_EVE_dsDNA,Nb_EVE_ssDNA,Nb_EVE_dsRNA,Nb_EVE_ssRNA)*dsDNA_perc)/100),
                 ((sum(Nb_EVE_dsDNA,Nb_EVE_ssDNA,Nb_EVE_dsRNA,Nb_EVE_ssRNA)*ssDNA_perc)/100),
                 ((sum(Nb_EVE_dsDNA,Nb_EVE_ssDNA,Nb_EVE_dsRNA,Nb_EVE_ssRNA)*dsRNA_perc)/100),
                 ((sum(Nb_EVE_dsDNA,Nb_EVE_ssDNA,Nb_EVE_dsRNA,Nb_EVE_ssRNA)*ssRNA_perc)/100))
Expected_Events<- c( ((sum(Events_EVEs_dsDNA,Events_EVEs_ssDNA,Events_EVEs_dsRNA,Events_EVEs_ssRNA)*dsDNA_perc)/100),
                     ((sum(Events_EVEs_dsDNA,Events_EVEs_ssDNA,Events_EVEs_dsRNA,Events_EVEs_ssRNA)*ssDNA_perc)/100),
                     ((sum(Events_EVEs_dsDNA,Events_EVEs_ssDNA,Events_EVEs_dsRNA,Events_EVEs_ssRNA)*dsRNA_perc)/100),
                     ((sum(Events_EVEs_dsDNA,Events_EVEs_ssDNA,Events_EVEs_dsRNA,Events_EVEs_ssRNA)*ssRNA_perc)/100))

Summary_table2 <- data.frame(Genomic_structures,Type,value,Expected_NB,Expected_Events)

Bar_plot_NB<-Summary_table2[!grepl("Event",Summary_table2$Type),] %>%
  ggplot( aes(y=value, x=factor(Genomic_structures , levels = c("dsDNA", "ssDNA", "dsRNA", "ssRNA")) ,fill=factor(Type, levels = c("dEVEs","EVEs")),label=Genomic_structures)) +
  geom_bar(stat="identity",position = "identity",color="black",size=0.2)+
  geom_point(size = 1.5,aes(y=Expected_NB, shape="Exposure",fill="black",stroke = 1),shape=3)+
  scale_fill_manual(values = c("EVEs"="#FFD460","dEVEs"="#F07B3F"))+  theme_minimal() +
  theme(
    panel.grid.major.x  = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line())+
  ylab("Number of EVEs") + xlab("Viral genomic structures") +
  theme(axis.text.y = element_text( color="black",size=9))+
  theme(legend.position = "none") +scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line())+
  theme(text = element_text(size=17)) + theme(strip.text.x = element_text(size = 13,face="bold"))+
  theme(axis.text.y = element_text( color="black",size=10))+theme(axis.title.x = element_blank())

#Add dEVEs Event count
Summary_table<-Summary_table %>%
  summarise(dsDNA = paste0(Events_dEVEs_dsDNA,"(",Nb_dEVEs_dsDNA_fl_Event,",",Nb_dEVEs_dsDNA_endo_Event,",",Nb_dEVEs_dsDNA_ecto_Event,",",Nb_dEVEs_dsDNA_Unkn_Event,")"),
            ssDNA = paste0(Events_dEVEs_ssDNA,"(",Nb_dEVEs_ssDNA_fl_Event,",",Nb_dEVEs_ssDNA_endo_Event,",",Nb_dEVEs_ssDNA_ecto_Event,",",Nb_dEVEs_ssDNA_Unkn_Event,")"),
            dsRNA = paste0(Events_dEVEs_dsRNA,"(",Nb_dEVEs_dsRNA_fl_Event,",",Nb_dEVEs_dsRNA_endo_Event,",",Nb_dEVEs_dsRNA_ecto_Event,",",Nb_dEVEs_dsRNA_Unkn_Event,")"),
            ssRNA = paste0(Events_dEVEs_ssRNA,"(",Nb_dEVEs_ssRNA_fl_Event,",",Nb_dEVEs_ssRNA_endo_Event,",",Nb_dEVEs_ssRNA_ecto_Event,",",Nb_dEVEs_ssRNA_Unkn_Event,")"),
            Unclassified = paste0(Events_dEVEs_Unclassified,"(",Nb_dEVEs_Uncla_fl_Event,",",Nb_dEVEs_Uncla_endo_Event,",",Nb_dEVEs_Uncla_ecto_Event,",",Nb_dEVEs_Uncla_Unkn_Event,")"),
            Type="dEVEs_Events")%>%
  bind_rows(Summary_table)

#Add EVEs Event count


Events_EVEs_ssDNA <-length(unique(Env_table2$ID[Env_table2$consensus_genomic_structure=="ssDNA" ]))-1
Events_EVEs_dsDNA <-length(unique(Env_table2$ID[Env_table2$consensus_genomic_structure=="dsDNA" ]))-1
Events_EVEs_ssRNA <-length(unique(Env_table2$ID[Env_table2$consensus_genomic_structure=="ssRNA" ]))-1
Events_EVEs_dsRNA <-length(unique(Env_table2$ID[Env_table2$consensus_genomic_structure=="dsRNA" ]))-1


Events_EVEs_Unclassified <- sum(Event_alone_df$Nb_EVEs_Unclassified_Events,na.rm=T)+sum(Event_shared_df$Nb_EVEs_Unclassified_Events,na.rm=T)

Summary_table<-Summary_table %>%
  summarise(dsDNA = paste0(Events_EVEs_dsDNA,"(",Nb_EVEs_dsDNA_fl_Event,",",Nb_EVEs_dsDNA_endo_Event,",",Nb_EVEs_dsDNA_ecto_Event,",",Nb_EVEs_dsDNA_Unkn_Event,")"),
            ssDNA = paste0(Events_EVEs_ssDNA,"(",Nb_EVEs_ssDNA_fl_Event,",",Nb_EVEs_ssDNA_endo_Event,",",Nb_EVEs_ssDNA_ecto_Event,",",Nb_EVEs_ssDNA_Unkn_Event,")"),
            dsRNA = paste0(Events_EVEs_dsRNA,"(",Nb_EVEs_dsRNA_fl_Event,",",Nb_EVEs_dsRNA_endo_Event,",",Nb_EVEs_dsRNA_ecto_Event,",",Nb_EVEs_dsRNA_Unkn_Event,")"),
            ssRNA = paste0(Events_EVEs_ssRNA,"(",Nb_EVEs_ssRNA_fl_Event,",",Nb_EVEs_ssRNA_endo_Event,",",Nb_EVEs_ssRNA_ecto_Event,",",Nb_EVEs_ssRNA_Unkn_Event,")"),
            Unclassified = paste0(Events_EVEs_Unclassified,"(",Nb_EVEs_Uncla_fl_Event,",",Nb_EVEs_Uncla_endo_Event,",",Nb_EVEs_Uncla_ecto_Event,",",Nb_EVEs_Uncla_Unkn_Event,")"),
            Type="EVEs_Events")%>%
  bind_rows(Summary_table)


##########


#Add dEVEs count


dsDNA_tab_reduced<- dsDNA_tab_reduced[!is.na(dsDNA_tab_reduced$Species_name),]


library(dplyr)


NB_dEVEs_dsDNA <- nrow(filter (dsDNA_tab_reduced, consensus_genomic_structure %in% c("dsDNA") &   Mean_dNdS < 1 & SE_dNdS < 1 & pseudogenized ==0 | consensus_genomic_structure %in% c("dsDNA") & TPM_all > 1000 ))
NB_dEVEs_ssDNA <- nrow(filter (dsDNA_tab_reduced, consensus_genomic_structure %in% c("ssDNA") &   Mean_dNdS < 1 & SE_dNdS < 1 & pseudogenized ==0 | consensus_genomic_structure %in% c("ssDNA") & TPM_all > 1000 ))
NB_dEVEs_dsRNA <- nrow(filter (dsDNA_tab_reduced, consensus_genomic_structure %in% c("dsRNA") &   Mean_dNdS < 1 & SE_dNdS < 1 & pseudogenized ==0 | consensus_genomic_structure %in% c("dsRNA") & TPM_all > 1000 ))
NB_dEVEs_ssRNA <- nrow(filter (dsDNA_tab_reduced, consensus_genomic_structure %in% c("ssRNA") &   Mean_dNdS < 1 & SE_dNdS < 1 & pseudogenized ==0 | consensus_genomic_structure %in% c("ssRNA") & TPM_all > 1000 ))
NB_dEVEs_Unknown<- nrow(filter (dsDNA_tab_reduced, consensus_genomic_structure %in% c("Unknown") &   Mean_dNdS < 1 & SE_dNdS < 1 & pseudogenized ==0 | consensus_genomic_structure %in% c("Unknown") & TPM_all > 1000 ))


Summary_table<-Summary_table %>%
  summarise(dsDNA = as.character(NB_dEVEs_dsDNA),
            ssDNA = as.character(NB_dEVEs_ssDNA),
            dsRNA = as.character(NB_dEVEs_dsRNA),
            ssRNA = as.character(NB_dEVEs_ssRNA),
            Unclassified = as.character(NB_dEVEs_Unknown),
            Type="dEVEs")%>%
  bind_rows(Summary_table)
#Add EVEs count


NB_EVEs_dsDNA <- nrow(filter (dsDNA_tab_reduced, consensus_genomic_structure %in% c("dsDNA") ))
NB_EVEs_ssDNA <- nrow(filter (dsDNA_tab_reduced, consensus_genomic_structure %in% c("ssDNA") ))
NB_EVEs_dsRNA <- nrow(filter (dsDNA_tab_reduced, consensus_genomic_structure %in% c("dsRNA") ))
NB_EVEs_ssRNA <- nrow(filter (dsDNA_tab_reduced, consensus_genomic_structure %in% c("ssRNA") ))
NB_EVEs_Unknown<- nrow(filter (dsDNA_tab_reduced, consensus_genomic_structure %in% c("Unknown")))


Summary_table<-Summary_table %>%
  summarise(dsDNA = as.character(NB_EVEs_dsDNA),
            ssDNA = as.character(NB_EVEs_ssDNA),
            dsRNA = as.character(NB_EVEs_dsRNA),
            ssRNA = as.character(NB_EVEs_ssRNA),
            Unclassified = as.character(NB_EVEs_Unknown),
            Type="EVEs")%>%
  bind_rows(Summary_table)
#Add total
Summary_table$Total<-as.numeric(gsub("\\(.*","",Summary_table$dsDNA))+as.numeric(gsub("\\(.*","",Summary_table$ssDNA))+as.numeric(gsub("\\(.*","",Summary_table$dsRNA))+as.numeric(gsub("\\(.*","",Summary_table$ssRNA))+as.numeric(gsub("\\(.*","",Summary_table$Unclassified))



library(gridExtra)
pdf("/Users/bguinet/Desktop/Papier_scientifique/Table_count_summary_freeliving.pdf", height=11, width=10)
grid.table(Summary_table)
dev.off()

library(tidyverse)
library(ggrepel)
library(forcats)
library(dplyr)

###################
###Paper output ###
##################

Event_all_df <- rbind(Event_alone_df[c('Event')],Event_shared_df[c('Event')])
Event_all_df$ID <- seq.int(nrow(Event_all_df))
library(dplyr)
library(tidyr) # unnest, separate

Env_table2<-Env_table #[c("New_query_bis2","Clustername","target","Event","Boot","Species_name","pident","alnlen","evalue","bits","Scaffold_length_x","best_family_per_query","family","species","genomic_structure","Scaffold_score","Scaffold_score2","GC_content_scaffold","pvalue_gc","cov_depth_BUSCO","cov_depth_candidat","pvalue_cov","count_repeat","count_eucaryote","Number_viral_loci","Number_busco_loci","Gene.names","Protein.names","Gene.ontology..GO.","consensus_Protein_names","consensus_Gene_names","consensus_Domain_description","Gene.ontology..biological.process.","Gene.ontology..molecular.function.","Gene.ontology..cellular.component.","Gene.ontology.IDs","Pfam_acc","Hmmscan_evalue","Domain_description","Synteny_Windows","Synteny_Nb_HSPs","Synteny_Tot_alnlen","NSPtot","Nviraltot","Mean_dNdS","Pvalue_dNdS","SE_dNdS","FDR_pvalue_dNdS","pseudogenized","ORF_perc","len_ORF","start_ORF","end_ORF","perc_tlen_alnlen","sequence","Events_species","lifecycle1")]
length(unique(Env_table2[Env_table2$best_family_per_query =="Parvoviridae",]$query))

Env_table2$Event <- as.integer(Env_table2$Event)

Event_all_df<-Event_all_df %>%
  mutate(Event = strsplit(Event, "[ ,]+")) %>%
  unnest(Event) %>%
  separate(Event, into = c("Event", "Clustername")) %>%
  mutate(Event = as.integer(Event))

Env_table2<-merge(x = Event_all_df, y = Env_table2, by = c("Clustername","Event"), all = TRUE)

Env_table2<-Env_table2[!Env_table2$New_query_bis2 %in% list_query_remove_because_redundancy,]

Env_table2$Event <- as.character(Env_table2$Event)

Env_table2 <-Env_table2[c("New_query_bis2","consensus_Protein_names","consensus_Domain_description","Gene.ontology..GO.")]%>%
  #filter(!is.na(Protein.names)) %>%
  pivot_longer(cols = -New_query_bis2) %>%
  filter(!is.na(value) & value != "Uncharacterized protein") %>%
  group_by(New_query_bis2) %>%
  summarise(Consensus_function= first(value)) %>%
  ungroup() %>%
  right_join(Env_table2) %>%
  rowwise() %>%
  mutate(Consensus_function= replace(Consensus_function, is.na(Consensus_function), consensus_Protein_names)) %>%
  relocate(Consensus_function, .after = Gene.ontology..GO.) %>%
  arrange(New_query_bis2)

Env_table2 <-Env_table2 %>%
  group_by(Clustername,Event) %>%
  slice(1) %>%
  select(Clustername, Event, Consensus_function_Event = Consensus_function) %>%
  right_join(Env_table2, by = c("Clustername", "Event")) %>%
  relocate(Consensus_function_Event, .after = Consensus_function)

Env_table2<-Env_table2[!is.na(Env_table2$ID),]


#Set consensus family per ID Events
Env_table2 <- Env_table2%>%
  group_by(ID) %>%
  add_count(best_family_per_query) %>%
  top_n(1, n) %>%
  distinct(consensus_family_ID = best_family_per_query) %>%
  right_join(Env_table2)

#Change Porseoliae to Unknown
Env_table2$Species_name[grepl("Unknown_species",Env_table2$New_query_bis2)]<-"Unknown_sp"

#Change Campopleginae
Env_table2$family[grepl("Pichon",Env_table2$target) & Env_table2$Species_name =="Campopleginae" ] <- "IVSPERs"
Env_table2$best_family_per_query[grepl("Pichon",Env_table2$target) & Env_table2$Species_name =="Campopleginae" ] <- "IVSPERs"


Env_table2<-Env_table2[!duplicated(Env_table2[c("New_query_bis2","Event")]),]


Env_table3<-Env_table2[c("New_query_bis2","Clustername","target","Event","ID","Boot","Species_name","pident","alnlen","evalue","bits","Scaffold_length_x","best_family_per_query","family","consensus_family","species","genomic_structure","consensus_genomic_structure","Scaffold_score","GC_content_scaffold","pvalue_gc","cov_depth_BUSCO","cov_depth_candidat","pvalue_cov","count_repeat","count_eucaryote","Number_viral_loci","Number_busco_loci","Gene.names","Protein.names","Gene.ontology..GO.","consensus_Protein_names","consensus_Gene_names","consensus_Domain_description","Consensus_function","Gene.ontology..biological.process.","Gene.ontology..molecular.function.","Gene.ontology..cellular.component.","Gene.ontology.IDs","Pfam_acc","Hmmscan_evalue","Domain_description","Synteny_Windows","Synteny_Nb_HSPs","Synteny_Tot_alnlen","NSPtot","Nviraltot","Mean_dNdS","Pvalue_dNdS","SE_dNdS","FDR_pvalue_dNdS","TPM_all","pseudogenized","ORF_perc","len_ORF","start_ORF","end_ORF","perc_tlen_alnlen","sequence","Events_species","lifecycle1")]

Env_table3$consensus_genomic_structure[Env_table3$consensus_family=="Unknown" & Env_table3$family=="Unknown" & is.na(Env_table3$genomic_structure)]<-"Unclassified_RNA"

write.table(Env_table3,"/Users/bguinet/Desktop/Papier_scientifique/Env_table3.txt",sep=";")


#Add Sequence when there is not ORFs
library(Biostrings)
Fasta_sequences<- readAAStringSet("/Users/bguinet/Desktop/Papier_scientifique/Candidate_loci_filtred.aa")

for (seq in unique(Env_table3$New_query_bis2)){
  if (is.na(Env_table3$sequence[Env_table3$New_query_bis2==seq]) ==TRUE){
    Env_table3$sequence[Env_table3$New_query_bis2==seq] <- as.character(Fasta_sequences[[seq]])

  }
}

write.table(Env_table3,"/Users/bguinet/Desktop/Papier_scientifique/Env_table2.txt",sep=",",row.names = F)

##### Create control heatmap results table

Control_comparaison_table<-read.table("/Users/bguinet/Desktop/Papier_scientifique/Domesticated_vs_Candidate_viral_segment_ALL.tab",sep=";",header=T)

Control_comparaison_table$Control_Names<-gsub("GbNV gp","GbNV_gp",Control_comparaison_table$Control_Names)
#Add informations from Big table
Control_comparaison_table<- merge(Control_comparaison_table, Env_table3[c("Clustername","New_query_bis2","Scaffold_score","Pvalue_dNdS","FDR_pvalue_dNdS","Mean_dNdS","pseudogenized","TPM_all","SE_dNdS")],  by.y="New_query_bis2", by.x = "query",outer=T)

#If the standard error of the dN/dS value is to high
Control_comparaison_table$SE_dNdS_plus<- Control_comparaison_table$SE_dNdS+Control_comparaison_table$Mean_dNdS

#Add control to the big table

Env_table3<-merge(x=Env_table3, y =Control_comparaison_table[c("query","Control_Sp_Names")], by.x="New_query_bis2",by.y="query",all.x=T)

Control_comparaison_table <-Control_comparaison_table %>%
  group_by(Clustername, query)  %>%
  mutate(dNdS = case_when(
    any(Pvalue_dNdS>=0.05) ~ "2Purifying++",
    any(is.na(Pvalue_dNdS))~ "6None",
    any(FDR_pvalue_dNdS==1 & as.numeric(as.character(Mean_dNdS))+as.numeric(as.character(SE_dNdS)) < 0.98 & pseudogenized ==0 ) ~ "1*Purifying++",
    #any(FDR_pvalue_dNdS=="True" & between(Mean_dNdS, 0.50,0.98)) ~ "*Purifying++",
    any(FDR_pvalue_dNdS==0 & as.numeric(as.character(Mean_dNdS))+as.numeric(as.character(SE_dNdS)) < 0.98 & pseudogenized ==0 ) ~ "2Purifying++",
    #any(FDR_pvalue_dNdS=="False" & between(Mean_dNdS, 0.50,0.98)) ~ "Purifying++",
    any(FDR_pvalue_dNdS==1 & as.numeric(as.character(Mean_dNdS))+as.numeric(as.character(SE_dNdS)) < 10 & pseudogenized ==0) ~ "3*Diversifying++",
    any(Mean_dNdS==999) ~ "5Incomputable",
    any(SE_dNdS_plus >0.80) ~ "5Incomputable",
    any(FDR_pvalue_dNdS==0 & as.numeric(as.character(Mean_dNdS))+as.numeric(as.character(SE_dNdS)) < 10  & pseudogenized ==0) ~ "3Diversifying++",
    TRUE ~ NA_character_)
  ) %>%
  group_by(Clustername, query)


#Remove paralogues
Control_comparaison_table<-Control_comparaison_table[order(Control_comparaison_table$dNdS),]
Control_comparaison_table$dNdS <- gsub("1\\*Purifying\\++","\\*Purifying\\++",Control_comparaison_table$dNdS)
Control_comparaison_table$dNdS <- gsub("2Purifying\\++","Purifying\\++",Control_comparaison_table$dNdS)
Control_comparaison_table$dNdS <- gsub("3\\*Diversifying\\++","\\*Diversifying\\++",Control_comparaison_table$dNdS)
Control_comparaison_table$dNdS <- gsub("5Incomputable","\\*Incomputable",Control_comparaison_table$dNdS)
Control_comparaison_table$dNdS <- gsub("6None","None",Control_comparaison_table$dNdS)



#### Nudi heatmap

Control_comparaison_table_nudi<-Control_comparaison_table [!duplicated(Control_comparaison_table[c('Control_Sp_Names','Control_Names')]),]

#Keep only if present in proteins
df_order_nudi<-read.csv("/Users/bguinet/Desktop/these/Tab_order_candidate_name_nudiviridae.csv",header=T,stringsAsFactors = FALSE)
colnames(df_order_nudi) <- c("Name")
Control_comparaison_table_nudi<-Control_comparaison_table_nudi[Control_comparaison_table_nudi$Control_Names %in% df_order_nudi$Name,]
#Format

Control_comparaison_table_nudi<- Control_comparaison_table_nudi[c("Control_Sp_Names","Control_Names","dNdS","Scaffold_score","TPM_all")]
Control_comparaison_table_nudi$Groups<- Control_comparaison_table_nudi$Control_Sp_Names
Control_comparaison_table_nudi$Groups_f<- Control_comparaison_table_nudi$Control_Sp_Names
colnames(Control_comparaison_table_nudi)<-c("Species_name","Name","dNdS","Genomic_env_level2", "TPM_all","Groups","Groups_f")

Control_comparaison_table_nudi$Species_name <- ifelse(Control_comparaison_table_nudi$Species_name == "Farisanus", "F.arisanus*", as.character(Control_comparaison_table_nudi$Species_name))
Control_comparaison_table_nudi$Species_name <- ifelse(Control_comparaison_table_nudi$Species_name == "Vcanescens", "V.canescens*", as.character(Control_comparaison_table_nudi$Species_name))
Control_comparaison_table_nudi$Species_name <- ifelse(Control_comparaison_table_nudi$Species_name == "Mdemolitor", "M.demolitor*", as.character(Control_comparaison_table_nudi$Species_name))

library(tidyverse)
eg <- expand.grid(Species_name = unique(Control_comparaison_table_nudi$Species_name), Name = unique( df_order_nudi$Name), stringsAsFactors = FALSE)
Control_comparaison_table_nudi<-merge(Control_comparaison_table_nudi, eg, by = c("Species_name", "Name"), all = TRUE)

Controls_check_nudi_nudi<-read.table("/Users/bguinet/Desktop/these/Controls_only_Nudiviridae.txt",sep="\t",header=T,stringsAsFactors = FALSE)

Controls_check_nudi_nudi$Genomic_env_level2<-""
Controls_check_nudi_nudi$TPM_all<-FALSE
Controls_check_nudi_nudi$Groups<-"Free_living_viruses"
Controls_check_nudi_nudi$Groups_f<-"Free_living_viruses"
colnames(Controls_check_nudi_nudi) <- c("Name", "Species_name","dNdS","Genomic_env_level2", "TPM_all","Groups","Groups_f")
Controls_check_nudi_nudi$dNdS <- as.character(Controls_check_nudi_nudi$dNdS)

control_check_nudi_and_comparaison <- rbind(Control_comparaison_table_nudi,Controls_check_nudi_nudi)

control_check_nudi_and_comparaison<-control_check_nudi_and_comparaison[!control_check_nudi_and_comparaison$Name=='HvAV3g_orf95-like',]

control_check_nudi_and_comparaison$Species_name <- ifelse(control_check_nudi_and_comparaison$Species_name == "Cotesia_Chelonus", "C.vestalis*", as.character(control_check_nudi_and_comparaison$Species_name))
control_check_nudi_and_comparaison$Species_name <- ifelse(control_check_nudi_and_comparaison$Species_name == "MdBV", "M.demolitor", as.character(control_check_nudi_and_comparaison$Species_name))
control_check_nudi_and_comparaison$Species_name <- ifelse(control_check_nudi_and_comparaison$Species_name == "Vcanescens", "V.canescens", as.character(control_check_nudi_and_comparaison$Species_name))
control_check_nudi_and_comparaison$Species_name <- ifelse(control_check_nudi_and_comparaison$Species_name == "Farisanus", "F.arisanus", as.character(control_check_nudi_and_comparaison$Species_name))
control_check_nudi_and_comparaison$Species_name <- ifelse(control_check_nudi_and_comparaison$Species_name == "CcBV", "C.congregata", as.character(control_check_nudi_and_comparaison$Species_name))
control_check_nudi_and_comparaison$Species_name <- ifelse(control_check_nudi_and_comparaison$Species_name == "CiBV", "C.inanitus", as.character(control_check_nudi_and_comparaison$Species_name))
control_check_nudi_and_comparaison$Species_name <- ifelse(control_check_nudi_and_comparaison$Species_name == "Alphanudivirus_core", "alpha-Nv core", as.character(control_check_nudi_and_comparaison$Species_name))
control_check_nudi_and_comparaison$Species_name <- ifelse(control_check_nudi_and_comparaison$Species_name == "Betanudivirus_core", "beta-Nv core", as.character(control_check_nudi_and_comparaison$Species_name))
control_check_nudi_and_comparaison$Species_name <- ifelse(control_check_nudi_and_comparaison$Species_name == "Baculovirus_core", "Bv core", as.character(control_check_nudi_and_comparaison$Species_name))
control_check_nudi_and_comparaison$Species_name <- ifelse(control_check_nudi_and_comparaison$Species_name == "Nudivirus_core", "Nv core", as.character(control_check_nudi_and_comparaison$Species_name))

control_check_nudi_and_comparaison$Species_name<-factor(control_check_nudi_and_comparaison$Species_name , levels = c("Bv core"   , "Nv core"   ,  "alpha-Nv core","beta-Nv core","V.canescens","V.canescens*","F.arisanus","F.arisanus*","C.inanitus","C.congregata","C.vestalis","C.vestalis*","M.demolitor","M.demolitor*") )

control_check_nudi_and_comparaison$Name<-factor(control_check_nudi_and_comparaison$Name , levels = df_order_nudi$Name )

control_check_nudi_and_comparaison$Genomic_env_level2[control_check_nudi_and_comparaison$Genomic_env_level2=="NONE"]<-"N"
control_check_nudi_and_comparaison$dNdS[is.na(control_check_nudi_and_comparaison$dNdS)] <- 0
#deal with TPM
control_check_nudi_and_comparaison$TPM_all<- as.numeric(control_check_nudi_and_comparaison$TPM_all)
control_check_nudi_and_comparaison$TPM_all2<-NA
control_check_nudi_and_comparaison$TPM_all2[ control_check_nudi_and_comparaison$TPM_all >= 1000 ] <- "TRUE"
control_check_nudi_and_comparaison$TPM_all2[control_check_nudi_and_comparaison$TPM_all2 != "TRUE" ] <- "FALSE"
control_check_nudi_and_comparaison$TPM_all2[is.na(control_check_nudi_and_comparaison$TPM_all2)] <- "FALSE"

#Create group to separate them
control_check_nudi_and_comparaison$tile_group <- NA
control_check_nudi_and_comparaison$tile_group[ grepl("V.canescens",control_check_nudi_and_comparaison$Species_name)] <- "Venturia"
control_check_nudi_and_comparaison$tile_group[ grepl("F.arisanus",control_check_nudi_and_comparaison$Species_name)] <- "Fopius"
control_check_nudi_and_comparaison$tile_group[ grepl("C.",control_check_nudi_and_comparaison$Species_name)] <- "Cotesia"
control_check_nudi_and_comparaison$tile_group[ grepl("M.dem",control_check_nudi_and_comparaison$Species_name)] <- "Microplitis"
control_check_nudi_and_comparaison$tile_group[ grepl("core",control_check_nudi_and_comparaison$Species_name)] <- "Free_living_viruses"

cols <- c("endogenized"="dark green","1" = "dark green","0"="white","None"="grey","*Purifying++"="red","*Purifying+"="orange","Purifying++"="#3a3a3a","Purifying+"="#3a3a3a","Incomputable"="grey","*Diversifying++"="#0174DF")

control_check_nudi_and_comparaison$tile_group= factor(control_check_nudi_and_comparaison$tile_group, levels=c('Free_living_viruses','Venturia','Fopius','Cotesia','Microplitis'))
control_check_nudi_and_comparaison<-control_check_nudi_and_comparaison[!is.na(control_check_nudi_and_comparaison$tile_group),]
control_check_nudi_and_comparaison[is.na(control_check_nudi_and_comparaison)]<-0



#ODVE66 Farisanus <-
#  HvAV3g_orf95-like <- to many euk
#
p_nudi <- ggplot(control_check_nudi_and_comparaison,aes(y=Name,x=Species_name,fill=dNdS)) + theme_minimal()+
  scale_fill_manual(values = cols)+
  geom_tile(aes(fill = dNdS), color = "grey20",size=0.8)+
  facet_grid(~tile_group, scale="free") +
  theme(axis.text.x = element_text(angle=30, vjust = 0, hjust = 0))+
  theme(strip.text.x = element_blank())+scale_x_discrete(position = "top")+
  #theme(axis.text.x = element_text(size=18),axis.text.y = element_text(size=30))+
  geom_text(aes(label = Genomic_env_level2,fontface=2),size=11,col="white")  +
  theme(strip.text.x = element_blank()) +geom_tile(data=control_check_nudi_and_comparaison, mapping=aes(colour=factor(TPM_all2, c(TRUE, FALSE)),size=factor(TPM_all2, c(TRUE, FALSE))),size=2, alpha=0)+
  theme(legend.position="none")+scale_colour_manual("TPM_all2", values=c("black",adjustcolor("red", alpha.f = 0))) +
  scale_size_manual("TPM_all2", values=c(5, 0)) +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle=30, vjust = 0, hjust = -0.3),
        axis.text.y = element_text(size = 27)) + scale_x_discrete(position = "top")


### LbFV heatmap



#### LbFV heatmap

#Keep only if present in proteins
df_order_LbFV<-read.csv("/Users/bguinet/Desktop/these/Tab_order_candidate_name_LbFV.csv",header=T,stringsAsFactors = FALSE)

Control_comparaison_table_LbFV<-Control_comparaison_table
Control_comparaison_table_LbFV$Control_Names[ grepl("ORF78",Control_comparaison_table_LbFV$Control_Names)]<-"ORF78 (lef-9)"
Control_comparaison_table_LbFV$Control_Names[ grepl("ORF87",Control_comparaison_table_LbFV$Control_Names)]<-"ORF87"
Control_comparaison_table_LbFV$Control_Names[ grepl("ORF94",Control_comparaison_table_LbFV$Control_Names)]<-"ORF94"
Control_comparaison_table_LbFV$Control_Names[ grepl("ORF92",Control_comparaison_table_LbFV$Control_Names)]<-"ORF92 (lef-3)"

colnames(df_order_LbFV) <- c("Name")
Control_comparaison_table_LbFV<-Control_comparaison_table_LbFV[Control_comparaison_table_LbFV$Control_Names %in% df_order_LbFV$Name,]


#Format
Control_comparaison_table_LbFV$Control_Sp_Names <- gsub(".*:","",Control_comparaison_table_LbFV$query)
Control_comparaison_table_LbFV<- Control_comparaison_table_LbFV[c("Control_Sp_Names","Control_Names","dNdS","Scaffold_score","TPM_all")]
Control_comparaison_table_LbFV$Groups<- Control_comparaison_table_LbFV$Control_Sp_Names
Control_comparaison_table_LbFV$Groups_f<- Control_comparaison_table_LbFV$Control_Sp_Names
colnames(Control_comparaison_table_LbFV)<-c("Species_name","Name","dNdS","Genomic_env_level2", "TPM_all","Groups","Groups_f")

Control_comparaison_table_LbFV$Species_name <- ifelse(Control_comparaison_table_LbFV$Species_name == "Leptopilina_boulardi", "L.boulardi*", as.character(Control_comparaison_table_LbFV$Species_name))
Control_comparaison_table_LbFV$Species_name <- ifelse(Control_comparaison_table_LbFV$Species_name == "Leptopilina_heterotoma", "L.heterotoma*", as.character(Control_comparaison_table_LbFV$Species_name))
Control_comparaison_table_LbFV$Species_name <- ifelse(Control_comparaison_table_LbFV$Species_name == "Leptopilina_clavipes", "L.clavipes*", as.character(Control_comparaison_table_LbFV$Species_name))

library(tidyverse)
eg <- expand.grid(Species_name = unique(Control_comparaison_table_LbFV$Species_name), Name = unique( df_order_LbFV$Name), stringsAsFactors = FALSE)
Control_comparaison_table_LbFV<-merge(Control_comparaison_table_LbFV, eg, by = c("Species_name", "Name"), all = TRUE)

Controls_check_LbFV<-read.table("/Users/bguinet/Desktop/these/Controls_only_LbFV.txt",sep=";",header=T)

Controls_check_LbFV$Genomic_env_level2<-""
Controls_check_LbFV$TPM_all<-FALSE
Controls_check_LbFV$Groups<-"Free_living_viruses"
Controls_check_LbFV$Groups_f<-"Free_living_viruses"
colnames(Controls_check_LbFV) <- c("Name", "Species_name","dNdS","Genomic_env_level2", "TPM_all","Groups","Groups_f")
Controls_check_LbFV$dNdS <- as.character(Controls_check_LbFV$dNdS)

control_check_LbFV_and_comparaison <- rbind(Control_comparaison_table_LbFV,Controls_check_LbFV)

control_check_LbFV_and_comparaison$Genomic_env_level2[control_check_LbFV_and_comparaison$Genomic_env_level2=="NONE"]<-"N"
control_check_LbFV_and_comparaison$dNdS[is.na(control_check_LbFV_and_comparaison$dNdS)] <- 0
#deal with TPM
control_check_LbFV_and_comparaison$TPM_all<- as.numeric(control_check_LbFV_and_comparaison$TPM_all)
control_check_LbFV_and_comparaison$TPM_all2<-NA
control_check_LbFV_and_comparaison$TPM_all2[ control_check_LbFV_and_comparaison$TPM_all >= 1000 ] <- "TRUE"
control_check_LbFV_and_comparaison$TPM_all2[control_check_LbFV_and_comparaison$TPM_all2 != "TRUE" ] <- "FALSE"
control_check_LbFV_and_comparaison$TPM_all2[is.na(control_check_LbFV_and_comparaison$TPM_all2)] <- "FALSE"

#Create group to separate them
control_check_LbFV_and_comparaison$tile_group <- NA
control_check_LbFV_and_comparaison$tile_group[ grepl("L.boulardi",control_check_LbFV_and_comparaison$Species_name)] <- "boulardi"
control_check_LbFV_and_comparaison$tile_group[ grepl("L.clavipes",control_check_LbFV_and_comparaison$Species_name)] <- "clavipes"
control_check_LbFV_and_comparaison$tile_group[ grepl("L.heterotoma",control_check_LbFV_and_comparaison$Species_name)] <- "heterotoma"


cols <- c("endogenized"="dark green","1" = "dark green","0"="white","None"="grey","*Purifying++"="red","*Purifying+"="orange","Purifying++"="#3a3a3a","Purifying+"="#3a3a3a","Incomputable"="grey","*Diversifying++"="#0174DF")

control_check_LbFV_and_comparaison$tile_group= factor(control_check_LbFV_and_comparaison$tile_group, levels=c('boulardi','clavipes','heterotoma'))
control_check_LbFV_and_comparaison<-control_check_LbFV_and_comparaison[!is.na(control_check_LbFV_and_comparaison$tile_group),]
control_check_LbFV_and_comparaison[is.na(control_check_LbFV_and_comparaison)]<-0


control_check_LbFV_and_comparaison$Name<-factor(control_check_LbFV_and_comparaison$Name , levels = df_order_LbFV$Name )
control_check_LbFV_and_comparaison<-control_check_LbFV_and_comparaison[! duplicated(control_check_LbFV_and_comparaison),]


control_check_LbFV_and_comparaison<-control_check_LbFV_and_comparaison[!duplicated(control_check_LbFV_and_comparaison[c('Species_name','Name')]),]


#ODVE66 Farisanus <-
#  HvAV3g_orf95-like <- to many euk
#
p_LbFV <- ggplot(control_check_LbFV_and_comparaison,aes(y=Name,x=Species_name,fill=dNdS)) + theme_minimal()+
  scale_fill_manual(values = cols)+
  geom_tile(aes(fill = dNdS), color = "grey20",size=0.8)+
  facet_grid(~tile_group, scale="free") +
  theme(axis.text.x = element_text(angle=30, vjust = 0, hjust = 0))+
  theme(strip.text.x = element_blank())+scale_x_discrete(position = "top")+
  #theme(axis.text.x = element_text(size=18),axis.text.y = element_text(size=30))+
  geom_text(aes(label = Genomic_env_level2,fontface=2),size=11,col="white")  +
  theme(strip.text.x = element_blank()) +geom_tile(data=control_check_LbFV_and_comparaison, mapping=aes(colour=factor(TPM_all2, c(TRUE, FALSE)),size=factor(TPM_all2, c(TRUE, FALSE))),size=2, alpha=0)+
  theme(legend.position="none")+scale_colour_manual("TPM_all2", values=c("black",adjustcolor("red", alpha.f = 0))) +
  scale_size_manual("TPM_all2", values=c(5, 0)) +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle=30, vjust = 0, hjust = -0.3),
        axis.text.y = element_text(size = 27)) + scale_x_discrete(position = "top")


## add both heatmap together

p_nudi_lbfv=cowplot::plot_grid(p_nudi, p_LbFV,nrow=2,align = "v",axis = "lr", rel_heights=c(0.8,0.2))


# Create output control dataframe

Ntot_EVEs_Vcanescens= nrow(control_check_nudi_and_comparaison[control_check_nudi_and_comparaison$Species_name=="V.canescens*" & control_check_nudi_and_comparaison$Genomic_env_level2 %in% c("A","B","C","D"),])
Ntot_EVEs_Farisanus= nrow(control_check_nudi_and_comparaison[control_check_nudi_and_comparaison$Species_name=="F.arisanus*" & control_check_nudi_and_comparaison$Genomic_env_level2 %in% c("A","B","C","D"),])
Ntot_EVEs_Cvestalis=  nrow(control_check_nudi_and_comparaison[control_check_nudi_and_comparaison$Species_name=="C.vestalis*" & control_check_nudi_and_comparaison$Genomic_env_level2 %in% c("A","B","C","D"),])
Ntot_EVEs_Mdemolitor=  nrow(control_check_nudi_and_comparaison[control_check_nudi_and_comparaison$Species_name=="M.demolitor*" & control_check_nudi_and_comparaison$Genomic_env_level2 %in% c("A","B","C","D"),])
Ntot_EVEs_Lboulardi= nrow(control_check_LbFV_and_comparaison[control_check_LbFV_and_comparaison$Species_name=="L.boulardi*" & control_check_LbFV_and_comparaison$Genomic_env_level2 %in% c("A","B","C","D"),])
Ntot_EVEs_Lclavipes= nrow(control_check_LbFV_and_comparaison[control_check_LbFV_and_comparaison$Species_name=="L.clavipes*" & control_check_LbFV_and_comparaison$Genomic_env_level2 %in% c("A","B","C","D"),])
Ntot_EVEs_Lheterotoma= nrow(control_check_LbFV_and_comparaison[control_check_LbFV_and_comparaison$Species_name=="L.heterotoma*" & control_check_LbFV_and_comparaison$Genomic_env_level2 %in% c("A","B","C","D"),])
sum_Ntot_EVEs <- sum(Ntot_EVEs_Vcanescens,Ntot_EVEs_Farisanus,Ntot_EVEs_Cvestalis, Ntot_EVEs_Mdemolitor,Ntot_EVEs_Lboulardi, Ntot_EVEs_Lclavipes,Ntot_EVEs_Lheterotoma)


Ntot_dEVEs_Vcanescens= nrow(subset(control_check_nudi_and_comparaison , Species_name=="V.canescens*" & Genomic_env_level2 %in% c("A","B","C","D") & dNdS %in% c("*Purifying++","*Purifying+") | Species_name=="V.canescens*" & Genomic_env_level2 %in% c("A","B","C","D") & TPM_all2 %in% c("TRUE")))
Ntot_dEVEs_Cvestalis= nrow(subset(control_check_nudi_and_comparaison , Species_name=="C.vestalis*" & Genomic_env_level2 %in% c("A","B","C","D") & dNdS %in% c("*Purifying++","*Purifying+") | Species_name=="C.vestalis*" & Genomic_env_level2 %in% c("A","B","C","D") & TPM_all2 %in% c("TRUE")))
Ntot_dEVEs_Vcanescens= nrow(subset(control_check_nudi_and_comparaison , Species_name=="M.demolitor*" & Genomic_env_level2 %in% c("A","B","C","D") & dNdS %in% c("*Purifying++","*Purifying+") | Species_name=="M.demolitor*" & Genomic_env_level2 %in% c("A","B","C","D") & TPM_all2 %in% c("TRUE")))
Ntot_dEVEs_Mdemolitor= nrow(subset(control_check_nudi_and_comparaison , Species_name=="V.canescens*" & Genomic_env_level2 %in% c("A","B","C","D") & dNdS %in% c("*Purifying++","*Purifying+") | Species_name=="V.canescens*" & Genomic_env_level2 %in% c("A","B","C","D") & TPM_all2 %in% c("TRUE")))
Ntot_dEVEs_Farisanus= nrow(subset(control_check_nudi_and_comparaison , Species_name=="F.arisanus*" & Genomic_env_level2 %in% c("A","B","C","D") & dNdS %in% c("*Purifying++","*Purifying+") | Species_name=="F.arisanus*" & Genomic_env_level2 %in% c("A","B","C","D") & TPM_all2 %in% c("TRUE")))
Ntot_dEVEs_Lboulardi= nrow(subset(control_check_LbFV_and_comparaison , Species_name=="L.boulardi*" & Genomic_env_level2 %in% c("A","B","C","D") & dNdS %in% c("*Purifying++","*Purifying+") | Species_name=="L.boulardi*" & Genomic_env_level2 %in% c("A","B","C","D") & TPM_all2 %in% c("TRUE")))
Ntot_dEVEs_Lclavipes= nrow(subset(control_check_LbFV_and_comparaison , Species_name=="L.clavipes*" & Genomic_env_level2 %in% c("A","B","C","D") & dNdS %in% c("*Purifying++","*Purifying+") | Species_name=="L.clavipes*" & Genomic_env_level2 %in% c("A","B","C","D") & TPM_all2 %in% c("TRUE")))
Ntot_dEVEs_Lheterotoma= nrow(subset(control_check_LbFV_and_comparaison , Species_name=="L.heterotoma*" & Genomic_env_level2 %in% c("A","B","C","D") & dNdS %in% c("*Purifying++","*Purifying+") | Species_name=="L.heterotoma*" & Genomic_env_level2 %in% c("A","B","C","D") & TPM_all2 %in% c("TRUE")))

sum_Ntot_dEVEs <- sum(Ntot_dEVEs_Vcanescens,Ntot_dEVEs_Farisanus,Ntot_dEVEs_Cvestalis, Ntot_dEVEs_Mdemolitor,Ntot_dEVEs_Lboulardi, Ntot_dEVEs_Lclavipes,Ntot_dEVEs_Lheterotoma)


statistics_table_control=data.frame (Type= c("V.canescens", 'F.arisanus','C.vestalis',"M.demolitor","L.boulardi","L.clavipes","L.heterotoma","Nb_expected"),
                                     Raw_EVEs= c(paste0(" ",Ntot_EVEs_Vcanescens, "/40"),paste0(" ",Ntot_EVEs_Farisanus, "/47"),paste0(" ",Ntot_EVEs_Cvestalis, "/21"),paste0(" ",Ntot_EVEs_Mdemolitor, "/25"),paste0(" ",Ntot_EVEs_Lboulardi, "/13"),paste0(" ",Ntot_EVEs_Lclavipes, "/13"),paste0(" ",Ntot_EVEs_Lheterotoma, "/13"), sum_Ntot_EVEs/sum(40,47,21,25,13,13,13)*100),
                                     Scaffold_A= c(paste0(" ",nrow(control_check_nudi_and_comparaison[control_check_nudi_and_comparaison$Species_name=="V.canescens*" & control_check_nudi_and_comparaison$Genomic_env_level2 %in% c("A"),]), "/",Ntot_EVEs_Vcanescens),
                                                   paste0(" ",nrow(control_check_nudi_and_comparaison[control_check_nudi_and_comparaison$Species_name=="F.arisanus*" & control_check_nudi_and_comparaison$Genomic_env_level2 %in% c("A"),]), "/",Ntot_EVEs_Farisanus),
                                                   paste0(" ",nrow(control_check_nudi_and_comparaison[control_check_nudi_and_comparaison$Species_name=="C.vestalis*" & control_check_nudi_and_comparaison$Genomic_env_level2 %in% c("A"),]), "/",Ntot_EVEs_Cvestalis),
                                                   paste0(" ",nrow(control_check_nudi_and_comparaison[control_check_nudi_and_comparaison$Species_name=="M.demolitor*" & control_check_nudi_and_comparaison$Genomic_env_level2 %in% c("A"),]), "/",Ntot_EVEs_Mdemolitor),
                                                   paste0(" ",nrow(control_check_LbFV_and_comparaison[control_check_LbFV_and_comparaison$Species_name=="L.boulardi*" & control_check_LbFV_and_comparaison$Genomic_env_level2 %in% c("A"),]), "/",Ntot_EVEs_Lboulardi),
                                                   paste0(" ",nrow(control_check_LbFV_and_comparaison[control_check_LbFV_and_comparaison$Species_name=="L.clavipes*" & control_check_LbFV_and_comparaison$Genomic_env_level2 %in% c("A"),]), "/",Ntot_EVEs_Lclavipes),
                                                   paste0(" ",nrow(control_check_LbFV_and_comparaison[control_check_LbFV_and_comparaison$Species_name=="L.heterotoma*" & control_check_LbFV_and_comparaison$Genomic_env_level2 %in% c("A"),]), "/",Ntot_EVEs_Lheterotoma),
                                                   sum_Ntot_EVEs),
                                     Scaffold_B= c(paste0(" ",nrow(control_check_nudi_and_comparaison[control_check_nudi_and_comparaison$Species_name=="V.canescens*" & control_check_nudi_and_comparaison$Genomic_env_level2 %in% c("B"),]), "/",Ntot_EVEs_Vcanescens),
                                                   paste0(" ",nrow(control_check_nudi_and_comparaison[control_check_nudi_and_comparaison$Species_name=="F.arisanus*" & control_check_nudi_and_comparaison$Genomic_env_level2 %in% c("B"),]), "/",Ntot_EVEs_Farisanus),
                                                   paste0(" ",nrow(control_check_nudi_and_comparaison[control_check_nudi_and_comparaison$Species_name=="C.vestalis*" & control_check_nudi_and_comparaison$Genomic_env_level2 %in% c("B"),]), "/",Ntot_EVEs_Cvestalis),
                                                   paste0(" ",nrow(control_check_nudi_and_comparaison[control_check_nudi_and_comparaison$Species_name=="M.demolitor*" & control_check_nudi_and_comparaison$Genomic_env_level2 %in% c("B"),]), "/",Ntot_EVEs_Mdemolitor),
                                                   paste0(" ",nrow(control_check_LbFV_and_comparaison[control_check_LbFV_and_comparaison$Species_name=="L.boulardi*" & control_check_LbFV_and_comparaison$Genomic_env_level2 %in% c("B"),]), "/",Ntot_EVEs_Lboulardi),
                                                   paste0(" ",nrow(control_check_LbFV_and_comparaison[control_check_LbFV_and_comparaison$Species_name=="L.clavipes*" & control_check_LbFV_and_comparaison$Genomic_env_level2 %in% c("B"),]), "/",Ntot_EVEs_Lclavipes),
                                                   paste0(" ",nrow(control_check_LbFV_and_comparaison[control_check_LbFV_and_comparaison$Species_name=="L.heterotoma*" & control_check_LbFV_and_comparaison$Genomic_env_level2 %in% c("B"),]), "/",Ntot_EVEs_Lheterotoma),
                                                   sum_Ntot_EVEs),
                                     Scaffold_B= c(paste0(" ",nrow(control_check_nudi_and_comparaison[control_check_nudi_and_comparaison$Species_name=="V.canescens*" & control_check_nudi_and_comparaison$Genomic_env_level2 %in% c("B"),]), "/",Ntot_EVEs_Vcanescens),
                                                   paste0(" ",nrow(control_check_nudi_and_comparaison[control_check_nudi_and_comparaison$Species_name=="F.arisanus*" & control_check_nudi_and_comparaison$Genomic_env_level2 %in% c("B"),]), "/",Ntot_EVEs_Farisanus),
                                                   paste0(" ",nrow(control_check_nudi_and_comparaison[control_check_nudi_and_comparaison$Species_name=="C.vestalis*" & control_check_nudi_and_comparaison$Genomic_env_level2 %in% c("B"),]), "/",Ntot_EVEs_Cvestalis),
                                                   paste0(" ",nrow(control_check_nudi_and_comparaison[control_check_nudi_and_comparaison$Species_name=="M.demolitor*" & control_check_nudi_and_comparaison$Genomic_env_level2 %in% c("B"),]), "/",Ntot_EVEs_Mdemolitor),
                                                   paste0(" ",nrow(control_check_LbFV_and_comparaison[control_check_LbFV_and_comparaison$Species_name=="L.boulardi*" & control_check_LbFV_and_comparaison$Genomic_env_level2 %in% c("B"),]), "/",Ntot_EVEs_Lboulardi),
                                                   paste0(" ",nrow(control_check_LbFV_and_comparaison[control_check_LbFV_and_comparaison$Species_name=="L.clavipes*" & control_check_LbFV_and_comparaison$Genomic_env_level2 %in% c("B"),]), "/",Ntot_EVEs_Lclavipes),
                                                   paste0(" ",nrow(control_check_LbFV_and_comparaison[control_check_LbFV_and_comparaison$Species_name=="L.heterotoma*" & control_check_LbFV_and_comparaison$Genomic_env_level2 %in% c("B"),]), "/",Ntot_EVEs_Lheterotoma),
                                                   sum_Ntot_EVEs),
                                     Scaffold_C= c(paste0(" ",nrow(control_check_nudi_and_comparaison[control_check_nudi_and_comparaison$Species_name=="V.canescens*" & control_check_nudi_and_comparaison$Genomic_env_level2 %in% c("C"),]), "/",Ntot_EVEs_Vcanescens),
                                                   paste0(" ",nrow(control_check_nudi_and_comparaison[control_check_nudi_and_comparaison$Species_name=="F.arisanus*" & control_check_nudi_and_comparaison$Genomic_env_level2 %in% c("C"),]), "/",Ntot_EVEs_Farisanus),
                                                   paste0(" ",nrow(control_check_nudi_and_comparaison[control_check_nudi_and_comparaison$Species_name=="C.vestalis*" & control_check_nudi_and_comparaison$Genomic_env_level2 %in% c("C"),]), "/",Ntot_EVEs_Cvestalis),
                                                   paste0(" ",nrow(control_check_nudi_and_comparaison[control_check_nudi_and_comparaison$Species_name=="M.demolitor*" & control_check_nudi_and_comparaison$Genomic_env_level2 %in% c("C"),]), "/",Ntot_EVEs_Mdemolitor),
                                                   paste0(" ",nrow(control_check_LbFV_and_comparaison[control_check_LbFV_and_comparaison$Species_name=="L.boulardi*" & control_check_LbFV_and_comparaison$Genomic_env_level2 %in% c("C"),]), "/",Ntot_EVEs_Lboulardi),
                                                   paste0(" ",nrow(control_check_LbFV_and_comparaison[control_check_LbFV_and_comparaison$Species_name=="L.clavipes*" & control_check_LbFV_and_comparaison$Genomic_env_level2 %in% c("C"),]), "/",Ntot_EVEs_Lclavipes),
                                                   paste0(" ",nrow(control_check_LbFV_and_comparaison[control_check_LbFV_and_comparaison$Species_name=="L.heterotoma*" & control_check_LbFV_and_comparaison$Genomic_env_level2 %in% c("C"),]), "/",Ntot_EVEs_Lheterotoma),
                                                   sum_Ntot_EVEs),
                                     Scaffold_D= c(paste0(" ",nrow(control_check_nudi_and_comparaison[control_check_nudi_and_comparaison$Species_name=="V.canescens*" & control_check_nudi_and_comparaison$Genomic_env_level2 %in% c("D"),]), "/",Ntot_EVEs_Vcanescens),
                                                   paste0(" ",nrow(control_check_nudi_and_comparaison[control_check_nudi_and_comparaison$Species_name=="F.arisanus*" & control_check_nudi_and_comparaison$Genomic_env_level2 %in% c("D"),]), "/",Ntot_EVEs_Farisanus),
                                                   paste0(" ",nrow(control_check_nudi_and_comparaison[control_check_nudi_and_comparaison$Species_name=="C.vestalis*" & control_check_nudi_and_comparaison$Genomic_env_level2 %in% c("D"),]), "/",Ntot_EVEs_Cvestalis),
                                                   paste0(" ",nrow(control_check_nudi_and_comparaison[control_check_nudi_and_comparaison$Species_name=="M.demolitor*" & control_check_nudi_and_comparaison$Genomic_env_level2 %in% c("D"),]), "/",Ntot_EVEs_Mdemolitor),
                                                   paste0(" ",nrow(control_check_LbFV_and_comparaison[control_check_LbFV_and_comparaison$Species_name=="L.boulardi*" & control_check_LbFV_and_comparaison$Genomic_env_level2 %in% c("D"),]), "/",Ntot_EVEs_Lboulardi),
                                                   paste0(" ",nrow(control_check_LbFV_and_comparaison[control_check_LbFV_and_comparaison$Species_name=="L.clavipes*" & control_check_LbFV_and_comparaison$Genomic_env_level2 %in% c("D"),]), "/",Ntot_EVEs_Lclavipes),
                                                   paste0(" ",nrow(control_check_LbFV_and_comparaison[control_check_LbFV_and_comparaison$Species_name=="L.heterotoma*" & control_check_LbFV_and_comparaison$Genomic_env_level2 %in% c("D"),]), "/",Ntot_EVEs_Lheterotoma),
                                                   sum_Ntot_EVEs),
                                     Scaffold_E_F_X= c(paste0(" ",nrow(control_check_nudi_and_comparaison[control_check_nudi_and_comparaison$Species_name=="V.canescens*" & control_check_nudi_and_comparaison$Genomic_env_level2 %in% c("E,F,X"),]), "/",Ntot_EVEs_Vcanescens),
                                                       paste0(" ",nrow(control_check_nudi_and_comparaison[control_check_nudi_and_comparaison$Species_name=="F.arisanus*" & control_check_nudi_and_comparaison$Genomic_env_level2 %in% c("E,F,X"),]), "/",Ntot_EVEs_Farisanus),
                                                       paste0(" ",nrow(control_check_nudi_and_comparaison[control_check_nudi_and_comparaison$Species_name=="C.vestalis*" & control_check_nudi_and_comparaison$Genomic_env_level2 %in% c("E,F,X"),]), "/",Ntot_EVEs_Mdemolitor),
                                                       paste0(" ",nrow(control_check_nudi_and_comparaison[control_check_nudi_and_comparaison$Species_name=="M.demolitor*" & control_check_nudi_and_comparaison$Genomic_env_level2 %in% c("E,F,X"),]), "/",Ntot_EVEs_Lboulardi),
                                                       paste0(" ",nrow(control_check_LbFV_and_comparaison[control_check_LbFV_and_comparaison$Species_name=="L.boulardi*" & control_check_LbFV_and_comparaison$Genomic_env_level2 %in% c("E,F,X"),]), "/",Ntot_EVEs_Lclavipes),
                                                       paste0(" ",nrow(control_check_LbFV_and_comparaison[control_check_LbFV_and_comparaison$Species_name=="L.clavipes*" & control_check_LbFV_and_comparaison$Genomic_env_level2 %in% c("E,F,X"),]), "/",Ntot_EVEs_Lclavipes),
                                                       paste0(" ",nrow(control_check_LbFV_and_comparaison[control_check_LbFV_and_comparaison$Species_name=="L.heterotoma*" & control_check_LbFV_and_comparaison$Genomic_env_level2 %in% c("E,F,X"),]), "/",Ntot_EVEs_Lheterotoma),
                                                       sum_Ntot_EVEs),
                                     Raw_dEVEs = c(paste0(" ",Ntot_dEVEs_Vcanescens, "/",Ntot_EVEs_Vcanescens),paste0(" ",Ntot_dEVEs_Farisanus, "/",Ntot_EVEs_Farisanus),paste0(" ",Ntot_dEVEs_Cvestalis, "/",Ntot_EVEs_Cvestalis),paste0(" ",Ntot_dEVEs_Mdemolitor, "/",Ntot_EVEs_Mdemolitor),paste0(" ",Ntot_dEVEs_Lboulardi, "/",Ntot_EVEs_Lboulardi),paste0(" ",Ntot_dEVEs_Lclavipes, "/",Ntot_EVEs_Lclavipes),paste0(" ",Ntot_dEVEs_Lheterotoma, "/",Ntot_dEVEs_Lheterotoma), paste0(" ",sum_Ntot_dEVEs,"/",sum_Ntot_EVEs)),
                                     by_dNdS =   c(paste0(" ",nrow(control_check_nudi_and_comparaison[control_check_nudi_and_comparaison$Species_name=="V.canescens*" & control_check_nudi_and_comparaison$Genomic_env_level2 %in% c("A","B","C","D") & control_check_nudi_and_comparaison$dNdS %in% c("*Purifying++","*Purifying+") ,]), "/",Ntot_EVEs_Vcanescens),
                                                   paste0(" ",nrow(control_check_nudi_and_comparaison[control_check_nudi_and_comparaison$Species_name=="F.arisanus*" & control_check_nudi_and_comparaison$Genomic_env_level2 %in% c("A","B","C","D")& control_check_nudi_and_comparaison$dNdS %in% c("*Purifying++","*Purifying+") ,]), "/",Ntot_EVEs_Farisanus),
                                                   paste0(" ",nrow(control_check_nudi_and_comparaison[control_check_nudi_and_comparaison$Species_name=="C.vestalis*" & control_check_nudi_and_comparaison$Genomic_env_level2 %in% c("A","B","C","D")& control_check_nudi_and_comparaison$dNdS %in% c("*Purifying++","*Purifying+") ,]), "/",Ntot_EVEs_Cvestalis),
                                                   paste0(" ",nrow(control_check_nudi_and_comparaison[control_check_nudi_and_comparaison$Species_name=="M.demolitor*" & control_check_nudi_and_comparaison$Genomic_env_level2 %in% c("A","B","C","D") & control_check_nudi_and_comparaison$dNdS %in% c("*Purifying++","*Purifying+"),]), "/",Ntot_EVEs_Mdemolitor),
                                                   paste0(" ",nrow(control_check_LbFV_and_comparaison[control_check_LbFV_and_comparaison$Species_name=="L.boulardi*" & control_check_LbFV_and_comparaison$Genomic_env_level2 %in% c("A","B","C","D") & control_check_LbFV_and_comparaison$dNdS %in% c("*Purifying++","*Purifying+"),] ), "/",Ntot_EVEs_Lboulardi),
                                                   paste0(" ",nrow(control_check_LbFV_and_comparaison[control_check_LbFV_and_comparaison$Species_name=="L.clavipes*" & control_check_LbFV_and_comparaison$Genomic_env_level2 %in% c("A","B","C","D") & control_check_LbFV_and_comparaison$dNdS %in% c("*Purifying++","*Purifying+"),] ), "/",Ntot_EVEs_Lclavipes),
                                                   paste0(" ",nrow(control_check_LbFV_and_comparaison[control_check_LbFV_and_comparaison$Species_name=="L.heterotoma*" & control_check_LbFV_and_comparaison$Genomic_env_level2 %in% c("A","B","C","D") & control_check_LbFV_and_comparaison$dNdS %in% c("*Purifying++","*Purifying+"),]), "/",Ntot_EVEs_Lheterotoma),
                                                   sum_Ntot_EVEs),
                                     by_TPM =   c(paste0(" ",nrow(control_check_nudi_and_comparaison[control_check_nudi_and_comparaison$Species_name=="V.canescens*" & control_check_nudi_and_comparaison$Genomic_env_level2 %in% c("A","B","C","D") & control_check_nudi_and_comparaison$TPM_all2 %in% c("TRUE") ,]), "/",Ntot_EVEs_Vcanescens),
                                                  paste0(" ",nrow(control_check_nudi_and_comparaison[control_check_nudi_and_comparaison$Species_name=="F.arisanus*" & control_check_nudi_and_comparaison$Genomic_env_level2 %in% c("A","B","C","D")& control_check_nudi_and_comparaison$TPM_all2 %in% c("TRUE") ,]), "/",Ntot_EVEs_Farisanus),
                                                  paste0(" ",nrow(control_check_nudi_and_comparaison[control_check_nudi_and_comparaison$Species_name=="C.vestalis*" & control_check_nudi_and_comparaison$Genomic_env_level2 %in% c("A","B","C","D")& control_check_nudi_and_comparaison$TPM_all2 %in% c("TRUE") ,]), "/",Ntot_EVEs_Cvestalis),
                                                  paste0(" ",nrow(control_check_nudi_and_comparaison[control_check_nudi_and_comparaison$Species_name=="M.demolitor*" & control_check_nudi_and_comparaison$Genomic_env_level2 %in% c("A","B","C","D") & control_check_nudi_and_comparaison$TPM_all2 %in% c("TRUE"),]), "/",Ntot_EVEs_Mdemolitor),
                                                  paste0(" ",nrow(control_check_LbFV_and_comparaison[control_check_LbFV_and_comparaison$Species_name=="L.boulardi*" & control_check_LbFV_and_comparaison$Genomic_env_level2 %in% c("A","B","C","D") & control_check_LbFV_and_comparaison$TPM_all2 %in% c("TRUE"),] ), "/",Ntot_EVEs_Lboulardi),
                                                  paste0(" ",nrow(control_check_LbFV_and_comparaison[control_check_LbFV_and_comparaison$Species_name=="L.clavipes*" & control_check_LbFV_and_comparaison$Genomic_env_level2 %in% c("A","B","C","D") & control_check_LbFV_and_comparaison$TPM_all2 %in% c("TRUE"),]) , "/",Ntot_EVEs_Lclavipes),
                                                  paste0(" ",nrow(control_check_LbFV_and_comparaison[control_check_LbFV_and_comparaison$Species_name=="L.heterotoma*" & control_check_LbFV_and_comparaison$Genomic_env_level2 %in% c("A","B","C","D") & control_check_LbFV_and_comparaison$TPM_all %in% c("TRUE") ,]), "/",Ntot_EVEs_Lheterotoma),
                                                  sum_Ntot_EVEs),
                                     by_TPM100 =   c(paste0(" ",nrow(control_check_nudi_and_comparaison[control_check_nudi_and_comparaison$Species_name=="V.canescens*" & control_check_nudi_and_comparaison$Genomic_env_level2 %in% c("A","B","C","D") & control_check_nudi_and_comparaison$TPM_all >100 ,]), "/",Ntot_EVEs_Vcanescens),
                                                     paste0(" ",nrow(control_check_nudi_and_comparaison[control_check_nudi_and_comparaison$Species_name=="F.arisanus*" & control_check_nudi_and_comparaison$Genomic_env_level2 %in% c("A","B","C","D")& control_check_nudi_and_comparaison$TPM_all >100 ,]), "/",Ntot_EVEs_Farisanus),
                                                     paste0(" ",nrow(control_check_nudi_and_comparaison[control_check_nudi_and_comparaison$Species_name=="C.vestalis*" & control_check_nudi_and_comparaison$Genomic_env_level2 %in% c("A","B","C","D")& control_check_nudi_and_comparaison$TPM_all >100 ,]), "/",Ntot_EVEs_Cvestalis),
                                                     paste0(" ",nrow(control_check_nudi_and_comparaison[control_check_nudi_and_comparaison$Species_name=="M.demolitor*" & control_check_nudi_and_comparaison$Genomic_env_level2 %in% c("A","B","C","D") & control_check_nudi_and_comparaison$TPM_all >100,]), "/",Ntot_EVEs_Mdemolitor),
                                                     paste0(" ",nrow(control_check_LbFV_and_comparaison[control_check_LbFV_and_comparaison$Species_name=="L.boulardi*" & control_check_LbFV_and_comparaison$Genomic_env_level2 %in% c("A","B","C","D") & control_check_LbFV_and_comparaison$TPM_all >100 ,] ), "/",Ntot_EVEs_Lboulardi),
                                                     paste0(" ",nrow(control_check_LbFV_and_comparaison[control_check_LbFV_and_comparaison$Species_name=="L.clavipes*" & control_check_LbFV_and_comparaison$Genomic_env_level2 %in% c("A","B","C","D") & control_check_LbFV_and_comparaison$TPM_all >100 ,] ), "/",Ntot_EVEs_Lclavipes),
                                                     paste0(" ",nrow(control_check_LbFV_and_comparaison[control_check_LbFV_and_comparaison$Species_name=="L.heterotoma*" & control_check_LbFV_and_comparaison$Genomic_env_level2 %in% c("A","B","C","D") & control_check_LbFV_and_comparaison$TPM_all >100 ,]), "/",Ntot_EVEs_Lheterotoma),
                                                     sum_Ntot_EVEs),
                                     Nb_Events =   c(paste0(" ",length(unique(Env_table3$ID[Env_table3$Species_name=="Venturia_canescens" & Env_table3$Control_Sp_Names=="Vcanescens"]))-1 ,"/1"),
                                                     paste0(" ",length(unique(Env_table3$ID[Env_table3$Species_name=="Fopius_arisanus" & Env_table3$Control_Sp_Names=="Farisanus"]))-1 ,"/1"),
                                                     paste0(" ",length(setdiff(unique(Env_table3$ID[Env_table3$Species_name=="Cotesia_vestalis" & Env_table3$Control_Sp_Names=="Cotesia_Chelonus"]),unique(Env_table3$ID[Env_table3$Species_name=="Microplitis_demolitor" & Env_table3$Control_Sp_Names=="Mdemolitor"]))) ,"/1"),
                                                     paste0(" ",length(setdiff(unique(Env_table3$ID[Env_table3$Species_name=="Cotesia_vestalis" & Env_table3$Control_Sp_Names=="Cotesia_Chelonus"]),unique(Env_table3$ID[Env_table3$Species_name=="Microplitis_demolitor" & Env_table3$Control_Sp_Names=="Mdemolitor"]))) ,"/1"),
                                                     paste0(" ",length(unique(Env_table3$ID[ grepl('Leptopilina',Env_table3$Species_name) & Env_table3$Control_Sp_Names=="Leptopilina"]))-1 ,"/1"),
                                                     paste0(" ",length(unique(Env_table3$ID[ grepl('Leptopilina',Env_table3$Species_name) & Env_table3$Control_Sp_Names=="Leptopilina"]))-1 ,"/1"),
                                                     paste0(" ",length(unique(Env_table3$ID[ grepl('Leptopilina',Env_table3$Species_name) & Env_table3$Control_Sp_Names=="Leptopilina"]))-1 ,"/1"),
                                                     sum_Ntot_EVEs))



statistics_table_control<-as.data.frame(t(statistics_table_control),)
colnames(statistics_table_control) <- statistics_table_control[1,]
statistics_table_control<-statistics_table_control[-1,]


#write.table(statistics_table_control,"/Users/bguinet/Desktop/Papier_scientifique/Statistics_control_table.txt",sep=";")

################




# BAR PLOT dEVENT
Events_dEVEs_dsDNA <-length(unique(Env_table2$ID[Env_table2$consensus_genomic_structure=="dsDNA" & Env_table2$FDR_pvalue_dNdS == 1 & as.numeric(as.character(Env_table2$Mean_dNdS))+as.numeric(as.character(Env_table2$SE_dNdS)) < 1 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1 | Env_table2$consensus_genomic_structure=="dsDNA" &Env_table2$TPM_all>=1000 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1]))-1
Events_dEVEs_ssDNA <-length(unique(Env_table2$ID[Env_table2$consensus_genomic_structure=="ssDNA" & Env_table2$FDR_pvalue_dNdS == 1 & as.numeric(as.character(Env_table2$Mean_dNdS))+as.numeric(as.character(Env_table2$SE_dNdS)) < 1 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1| Env_table2$consensus_genomic_structure=="ssDNA" &Env_table2$TPM_all>=1000 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1]))-1
Events_dEVEs_ssRNA <- length(unique(Env_table2$ID[Env_table2$consensus_genomic_structure=="ssRNA" & Env_table2$FDR_pvalue_dNdS == 1 & as.numeric(as.character(Env_table2$Mean_dNdS))+as.numeric(as.character(Env_table2$SE_dNdS)) < 1 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1| Env_table2$consensus_genomic_structure=="ssRNA" &Env_table2$TPM_all>=1000 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1]))-1
Events_dEVEs_dsRNA <- length(unique(Env_table2$ID[Env_table2$consensus_genomic_structure=="dsRNA" & Env_table2$FDR_pvalue_dNdS == 1 & as.numeric(as.character(Env_table2$Mean_dNdS))+as.numeric(as.character(Env_table2$SE_dNdS)) < 1 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1| Env_table2$consensus_genomic_structure=="dsRNA" &Env_table2$TPM_all>=1000 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1]))-1

Events_dEVEs_dsDNA_ENDO <-length(unique(Env_table2$ID[Env_table2$consensus_genomic_structure=="dsDNA" & Env_table2$FDR_pvalue_dNdS == 1 & as.numeric(as.character(Env_table2$Mean_dNdS))+as.numeric(as.character(Env_table2$SE_dNdS)) < 1 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1  & Env_table2$lifecycle1 =="endoparasitoide" | Env_table2$consensus_genomic_structure=="dsDNA" &Env_table2$TPM_all>=1000 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1 & Env_table2$lifecycle1 =="endoparasitoide"]))-1
Events_dEVEs_dsDNA_ECTO <-length(unique(Env_table2$ID[Env_table2$consensus_genomic_structure=="dsDNA" & Env_table2$FDR_pvalue_dNdS == 1 & as.numeric(as.character(Env_table2$Mean_dNdS))+as.numeric(as.character(Env_table2$SE_dNdS)) < 1 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1  & Env_table2$lifecycle1 =="ectoparasitoide" | Env_table2$consensus_genomic_structure=="dsDNA" &Env_table2$TPM_all>=1000 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1 & Env_table2$lifecycle1 =="ectoparasitoide"]))-1
Events_dEVEs_dsDNA_FREE <-length(unique(Env_table2$ID[Env_table2$consensus_genomic_structure=="dsDNA" & Env_table2$FDR_pvalue_dNdS == 1 & as.numeric(as.character(Env_table2$Mean_dNdS))+as.numeric(as.character(Env_table2$SE_dNdS)) < 1 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1  & Env_table2$lifecycle1 =="freeliving" | Env_table2$consensus_genomic_structure=="dsDNA" &Env_table2$TPM_all>=1000 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1 & Env_table2$lifecycle1 =="freeliving"]))-1
sum_Events_dEVEs_dsDNA<- Events_dEVEs_dsDNA_ENDO+Events_dEVEs_dsDNA_ECTO+Events_dEVEs_dsDNA_FREE

Events_dEVEs_ssDNA_ENDO <-length(unique(Env_table2$ID[Env_table2$consensus_genomic_structure=="ssDNA" & Env_table2$FDR_pvalue_dNdS == 1 & as.numeric(as.character(Env_table2$Mean_dNdS))+as.numeric(as.character(Env_table2$SE_dNdS)) < 1 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1  & Env_table2$lifecycle1 =="endoparasitoide" | Env_table2$consensus_genomic_structure=="ssDNA" &Env_table2$TPM_all>=1000 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1 & Env_table2$lifecycle1 =="endoparasitoide"]))-1
Events_dEVEs_ssDNA_ECTO <-length(unique(Env_table2$ID[Env_table2$consensus_genomic_structure=="ssDNA" & Env_table2$FDR_pvalue_dNdS == 1 & as.numeric(as.character(Env_table2$Mean_dNdS))+as.numeric(as.character(Env_table2$SE_dNdS)) < 1 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1  & Env_table2$lifecycle1 =="ectoparasitoide" | Env_table2$consensus_genomic_structure=="ssDNA" &Env_table2$TPM_all>=1000 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1 & Env_table2$lifecycle1 =="ectoparasitoide"]))-1
Events_dEVEs_ssDNA_FREE <-length(unique(Env_table2$ID[Env_table2$consensus_genomic_structure=="ssDNA" & Env_table2$FDR_pvalue_dNdS == 1 & as.numeric(as.character(Env_table2$Mean_dNdS))+as.numeric(as.character(Env_table2$SE_dNdS)) < 1 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1  & Env_table2$lifecycle1 =="freeliving" | Env_table2$consensus_genomic_structure=="ssDNA" &Env_table2$TPM_all>=1000 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1 & Env_table2$lifecycle1 =="freeliving"]))-1
sum_Events_dEVEs_ssDNA<- Events_dEVEs_ssDNA_ENDO+Events_dEVEs_ssDNA_ECTO+Events_dEVEs_ssDNA_FREE

Events_dEVEs_ssRNA_ENDO <-length(unique(Env_table2$ID[Env_table2$consensus_genomic_structure=="ssRNA" & Env_table2$FDR_pvalue_dNdS == 1 & as.numeric(as.character(Env_table2$Mean_dNdS))+as.numeric(as.character(Env_table2$SE_dNdS)) < 1 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1  & Env_table2$lifecycle1 =="endoparasitoide" | Env_table2$consensus_genomic_structure=="ssRNA" &Env_table2$TPM_all>=1000 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1 & Env_table2$lifecycle1 =="endoparasitoide"]))-1
Events_dEVEs_ssRNA_ECTO <-length(unique(Env_table2$ID[Env_table2$consensus_genomic_structure=="ssRNA" & Env_table2$FDR_pvalue_dNdS == 1 & as.numeric(as.character(Env_table2$Mean_dNdS))+as.numeric(as.character(Env_table2$SE_dNdS)) < 1 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1  & Env_table2$lifecycle1 =="ectoparasitoide" | Env_table2$consensus_genomic_structure=="ssRNA" &Env_table2$TPM_all>=1000 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1 & Env_table2$lifecycle1 =="ectoparasitoide"]))-1
Events_dEVEs_ssRNA_FREE <-length(unique(Env_table2$ID[Env_table2$consensus_genomic_structure=="ssRNA" & Env_table2$FDR_pvalue_dNdS == 1 & as.numeric(as.character(Env_table2$Mean_dNdS))+as.numeric(as.character(Env_table2$SE_dNdS)) < 1 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1  & Env_table2$lifecycle1 =="freeliving" | Env_table2$consensus_genomic_structure=="ssRNA" &Env_table2$TPM_all>=1000 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1 & Env_table2$lifecycle1 =="freeliving"]))-1
sum_Events_dEVEs_ssRNA<- Events_dEVEs_ssRNA_ENDO+Events_dEVEs_ssRNA_ECTO+Events_dEVEs_ssRNA_FREE

Events_dEVEs_dsRNA_ENDO <-length(unique(Env_table2$ID[Env_table2$consensus_genomic_structure=="dsRNA" & Env_table2$FDR_pvalue_dNdS == 1 & as.numeric(as.character(Env_table2$Mean_dNdS))+as.numeric(as.character(Env_table2$SE_dNdS)) < 1 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1  & Env_table2$lifecycle1 =="endoparasitoide" | Env_table2$consensus_genomic_structure=="dsRNA" &Env_table2$TPM_all>=1000 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1 & Env_table2$lifecycle1 =="endoparasitoide"]))-1
Events_dEVEs_dsRNA_ECTO <-length(unique(Env_table2$ID[Env_table2$consensus_genomic_structure=="dsRNA" & Env_table2$FDR_pvalue_dNdS == 1 & as.numeric(as.character(Env_table2$Mean_dNdS))+as.numeric(as.character(Env_table2$SE_dNdS)) < 1 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1  & Env_table2$lifecycle1 =="ectoparasitoide" | Env_table2$consensus_genomic_structure=="dsRNA" &Env_table2$TPM_all>=1000 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1 & Env_table2$lifecycle1 =="ectoparasitoide"]))-1
Events_dEVEs_dsRNA_FREE <-length(unique(Env_table2$ID[Env_table2$consensus_genomic_structure=="dsRNA" & Env_table2$FDR_pvalue_dNdS == 1 & as.numeric(as.character(Env_table2$Mean_dNdS))+as.numeric(as.character(Env_table2$SE_dNdS)) < 1 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1  & Env_table2$lifecycle1 =="freeliving" | Env_table2$consensus_genomic_structure=="dsRNA" &Env_table2$TPM_all>=1000 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1 & Env_table2$lifecycle1 =="freeliving"]))-1
sum_Events_dEVEs_dsRNA<- Events_dEVEs_dsRNA_ENDO+Events_dEVEs_dsRNA_ECTO+Events_dEVEs_dsRNA_FREE

dEvent_tab<- data.frame(
  Category = c("Events","Events","Events","Events","Events","Events","dEvents","dEvents","dEvents","dEvents","dEvents","dEvents"),
  Lifecycle=c("Endoparasitoid","Ectoparasitoid","Free-living","Endoparasitoid","Ectoparasitoid","Free-living","Endoparasitoid","Ectoparasitoid","Free-living","Endoparasitoid","Ectoparasitoid","Free-living"),
  Genomic_structure = c("dsDNA","dsDNA","dsDNA","ssDNA","ssDNA","ssDNA","ssRNA","ssRNA","ssRNA","dsRNA","dsRNA","dsRNA"),
  Nb_events = c(Events_dEVEs_dsDNA_ENDO,Events_dEVEs_dsDNA_ECTO,Events_dEVEs_dsDNA_FREE,Events_dEVEs_ssDNA_ENDO,Events_dEVEs_ssDNA_ECTO,Events_dEVEs_ssDNA_FREE,Events_dEVEs_ssRNA_ENDO,Events_dEVEs_ssRNA_ECTO,Events_dEVEs_ssRNA_FREE,Events_dEVEs_dsRNA_ENDO,Events_dEVEs_dsRNA_ECTO,Events_dEVEs_dsRNA_FREE),
  Sumcat = c(sum_Events_dEVEs_dsDNA,sum_Events_dEVEs_dsDNA,sum_Events_dEVEs_dsDNA,sum_Events_dEVEs_ssDNA,sum_Events_dEVEs_ssDNA,sum_Events_dEVEs_ssDNA,sum_Events_dEVEs_ssRNA,sum_Events_dEVEs_ssRNA,sum_Events_dEVEs_ssRNA,sum_Events_dEVEs_dsRNA,sum_Events_dEVEs_dsRNA,sum_Events_dEVEs_dsRNA)
)

dEvent_tab$Expected<-NA
dEvent_tab$Expected[dEvent_tab$Lifecycle=="Endoparasitoid"]<- (37/124)*100
dEvent_tab$Expected[dEvent_tab$Lifecycle=="Ectoparasitoid"]<- (24/124)*100
dEvent_tab$Expected[dEvent_tab$Lifecycle=="Free-living"]<- (65/124)*100
dEvent_tab$Nb_events_corrected <- dEvent_tab$Nb_events*100/dEvent_tab$Sumcat

library(ggrepel)
dEvent_tab$Genomic_structure = factor(dEvent_tab$Genomic_structure, levels=c('dsDNA','ssDNA','dsRNA','ssRNA'))


dEvents_representativity <- ggplot(dEvent_tab , aes(x = Lifecycle, y = Nb_events_corrected, fill =Lifecycle)) +
  facet_wrap( ~ Genomic_structure,nrow=1,strip.position = "top")+
  geom_bar(stat="identity",color="black")+
  geom_col(data = . %>% filter( Category=="dEVEs" & Lifecycle =="Free-living"), position = position_dodge(width = 0.9), alpha = 1,fill='#467fb3',color="black") +
  geom_col(data = . %>% filter( Category=="EVEs"& Lifecycle =="Free-living"), position = position_dodge(width = 0.9), alpha = 1,fill='#467fb3',color="black") +
  geom_col(data = . %>% filter( Category=="dEVEs" & Lifecycle =="Endoparasitoid"), position = position_dodge(width = 0.9), alpha = 1,fill='#06642e',color="black") +
  geom_col(data = . %>% filter( Category=="EVEs"& Lifecycle =="Endoparasitoid"), position = position_dodge(width = 0.9), alpha = 1,fill='#06642e',color="black") +
  geom_col(data = . %>% filter( Category=="dEVEs" & Lifecycle =="Ectoparasitoid"), position = position_dodge(width = 0.9), alpha = 1,fill='#d8be03',color="black") +
  geom_col(data = . %>% filter( Category=="EVEs"& Lifecycle =="Ectoparasitoid"), position = position_dodge(width = 0.9), alpha = 1,fill='#d8be03',color="black") +
  geom_point(size = 1.5,aes(y=Expected, shape="Exposure",fill="black"),shape=3,stroke = 1)+
  geom_text(aes(label = Nb_events,y=3),vjust = 0.1, nudge_y = .2,color="white")+
  ylab("EVEs %")+
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line())+
  theme(strip.text.x = element_text(size = 13,face="bold"))+
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        panel.spacing.y = unit(2, "lines"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+ ylim(0, 100)+ theme(legend.position = "none")+
  scale_fill_manual(values=c("#d8be03","#d8be03","#06642e","#467fb3"))

# BAR PLOT dEVENT
Events_EVEs_dsDNA <-length(unique(Env_table2$ID[Env_table2$consensus_genomic_structure=="dsDNA"]))-1
Events_EVEs_ssDNA <-length(unique(Env_table2$ID[Env_table2$consensus_genomic_structure=="ssDNA"]))-1
Events_EVEs_ssRNA <- length(unique(Env_table2$ID[Env_table2$consensus_genomic_structure=="ssRNA"]))-1
Events_EVEs_dsRNA <- length(unique(Env_table2$ID[Env_table2$consensus_genomic_structure=="dsRNA"]))-1


Events_EVEs_dsDNA_ENDO <-length(unique(Env_table2$ID[Env_table2$consensus_genomic_structure=="dsDNA" & Env_table2$lifecycle1 =="endoparasitoide"]))-1
Events_EVEs_dsDNA_ECTO <-length(unique(Env_table2$ID[Env_table2$consensus_genomic_structure=="dsDNA" & Env_table2$lifecycle1 =="ectoparasitoide"]))-1
Events_EVEs_dsDNA_FREE <-length(unique(Env_table2$ID[Env_table2$consensus_genomic_structure=="dsDNA" & Env_table2$lifecycle1 =="freeliving"]))-1
sum_Events_EVEs_dsDNA<- Events_EVEs_dsDNA_ENDO+Events_EVEs_dsDNA_ECTO+Events_EVEs_dsDNA_FREE

Events_EVEs_ssDNA_ENDO <-length(unique(Env_table2$ID[Env_table2$consensus_genomic_structure=="ssDNA" & Env_table2$lifecycle1 =="endoparasitoide"]))-1
Events_EVEs_ssDNA_ECTO <-length(unique(Env_table2$ID[Env_table2$consensus_genomic_structure=="ssDNA" & Env_table2$lifecycle1 =="ectoparasitoide"]))-1
Events_EVEs_ssDNA_FREE <-length(unique(Env_table2$ID[Env_table2$consensus_genomic_structure=="ssDNA" & Env_table2$lifecycle1 =="freeliving"]))-1
sum_Events_EVEs_ssDNA<- Events_EVEs_ssDNA_ENDO+Events_EVEs_ssDNA_ECTO+Events_EVEs_ssDNA_FREE

Events_EVEs_ssRNA_ENDO <-length(unique(Env_table2$ID[Env_table2$consensus_genomic_structure=="ssRNA" & Env_table2$lifecycle1 =="endoparasitoide"]))-1
Events_EVEs_ssRNA_ECTO <-length(unique(Env_table2$ID[Env_table2$consensus_genomic_structure=="ssRNA" & Env_table2$lifecycle1 =="ectoparasitoide"]))-1
Events_EVEs_ssRNA_FREE <-length(unique(Env_table2$ID[Env_table2$consensus_genomic_structure=="ssRNA" & Env_table2$lifecycle1 =="freeliving"]))-1
sum_Events_EVEs_ssRNA<- Events_EVEs_ssRNA_ENDO+Events_EVEs_ssRNA_ECTO+Events_EVEs_ssRNA_FREE

Events_EVEs_dsRNA_ENDO <-length(unique(Env_table2$ID[Env_table2$consensus_genomic_structure=="dsRNA" & Env_table2$lifecycle1 =="endoparasitoide"]))-1
Events_EVEs_dsRNA_ECTO <-length(unique(Env_table2$ID[Env_table2$consensus_genomic_structure=="dsRNA" & Env_table2$lifecycle1 =="ectoparasitoide"]))-1
Events_EVEs_dsRNA_FREE <-length(unique(Env_table2$ID[Env_table2$consensus_genomic_structure=="dsRNA" & Env_table2$lifecycle1 =="freeliving"]))-1
sum_Events_EVEs_dsRNA<- Events_EVEs_dsRNA_ENDO+Events_EVEs_dsRNA_ECTO+Events_EVEs_dsRNA_FREE

Event_tab<- data.frame(
  Category = c("Events","Events","Events","Events","Events","Events","dEvents","dEvents","dEvents","dEvents","dEvents","dEvents"),
  Lifecycle=c("Endoparasitoid","Ectoparasitoid","Free-living","Endoparasitoid","Ectoparasitoid","Free-living","Endoparasitoid","Ectoparasitoid","Free-living","Endoparasitoid","Ectoparasitoid","Free-living"),
  Genomic_structure = c("dsDNA","dsDNA","dsDNA","ssDNA","ssDNA","ssDNA","ssRNA","ssRNA","ssRNA","dsRNA","dsRNA","dsRNA"),
  Nb_events = c(Events_EVEs_dsDNA_ENDO,Events_EVEs_dsDNA_ECTO,Events_EVEs_dsDNA_FREE,Events_EVEs_ssDNA_ENDO,Events_EVEs_ssDNA_ECTO,Events_EVEs_ssDNA_FREE,Events_EVEs_ssRNA_ENDO,Events_EVEs_ssRNA_ECTO,Events_EVEs_ssRNA_FREE,Events_EVEs_dsRNA_ENDO,Events_EVEs_dsRNA_ECTO,Events_EVEs_dsRNA_FREE),
  Sumcat = c(sum_Events_EVEs_dsDNA,sum_Events_EVEs_dsDNA,sum_Events_EVEs_dsDNA,sum_Events_EVEs_ssDNA,sum_Events_EVEs_ssDNA,sum_Events_EVEs_ssDNA,sum_Events_EVEs_ssRNA,sum_Events_EVEs_ssRNA,sum_Events_EVEs_ssRNA,sum_Events_EVEs_dsRNA,sum_Events_EVEs_dsRNA,sum_Events_EVEs_dsRNA)
)

Event_tab$Expected<-NA
Event_tab$Expected[Event_tab$Lifecycle=="Endoparasitoid"]<- (37/124)*100
Event_tab$Expected[Event_tab$Lifecycle=="Ectoparasitoid"]<- (24/124)*100
Event_tab$Expected[Event_tab$Lifecycle=="Free-living"]<- (65/124)*100
Event_tab$Nb_events_corrected <- Event_tab$Nb_events*100/Event_tab$Sumcat

library(ggrepel)

Event_tab$Genomic_structure = factor(Event_tab$Genomic_structure, levels=c('dsDNA','ssDNA','dsRNA','ssRNA'))

Events_representativity <- ggplot(Event_tab , aes(x = Lifecycle, y = Nb_events_corrected, fill =Lifecycle)) +
  facet_wrap( ~ Genomic_structure,nrow=1,strip.position = "top")+
  geom_bar(stat="identity",color="black")+
  geom_col(data = . %>% filter( Category=="EVEs" & Lifecycle =="Free-living"), position = position_dodge(width = 0.9), alpha = 1,fill='#467fb3',color="black") +
  geom_col(data = . %>% filter( Category=="EVEs"& Lifecycle =="Free-living"), position = position_dodge(width = 0.9), alpha = 1,fill='#467fb3',color="black") +
  geom_col(data = . %>% filter( Category=="EVEs" & Lifecycle =="Endoparasitoid"), position = position_dodge(width = 0.9), alpha = 1,fill='#06642e',color="black") +
  geom_col(data = . %>% filter( Category=="EVEs"& Lifecycle =="Endoparasitoid"), position = position_dodge(width = 0.9), alpha = 1,fill='#06642e',color="black") +
  geom_col(data = . %>% filter( Category=="EVEs" & Lifecycle =="Ectoparasitoid"), position = position_dodge(width = 0.9), alpha = 1,fill='#d8be03',color="black") +
  geom_col(data = . %>% filter( Category=="EVEs"& Lifecycle =="Ectoparasitoid"), position = position_dodge(width = 0.9), alpha = 1,fill='#d8be03',color="black") +
  geom_point(size = 1.5,aes(y=Expected, shape="Exposure",fill="black"),shape=3,stroke = 1)+
  geom_text(aes(label = Nb_events,y=3),vjust = 0.1, nudge_y = .2,color="white")+
  ylab("EVEs %")+
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line())+
  theme(strip.text.x = element_text(size = 13,face="bold"))+
  theme(strip.placement = "outside",
        strip.background = element_blank(),
        panel.spacing.y = unit(2, "lines"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+ ylim(0, 100)+ theme(legend.position = "none")+
  scale_fill_manual(values=c("#d8be03","#d8be03","#06642e","#467fb3"))


library(ggpubr)
library(cowplot)
plot_grid(Events_representativity ,dEvents_representativity,ncol=2,nrow=1,labels=LETTERS[1:2],label_size = 16)





#Nb dEVEs total with complet ORF :
length(unique(Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),][Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),]$ORF_perc>0 & Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),]$FDR_pvalue_dNdS == 1 & as.numeric(as.character(Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),]$Mean_dNdS))+as.numeric(as.character(Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),]$SE_dNdS)) < 1 & Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),]$pseudogenized ==0 | Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),]$ORF_perc>0 &  Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),]$TPM_all>=1000 & Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),]$pseudogenized ==0 ,]$New_query_bis2))
length(unique(Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),][Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),]$ORF_perc>0 & Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),]$FDR_pvalue_dNdS == 1 & as.numeric(as.character(Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),]$Mean_dNdS))+as.numeric(as.character(Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),]$SE_dNdS)) < 1 & Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),]$pseudogenized ==0 | Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),]$ORF_perc>1 &  Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),]$TPM_all>=1000 & Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),]$pseudogenized ==0 ,]$New_query_bis2))


#Count number shared Events with synteny
Env_table2<-Env_table2 %>%
  group_by(Clustername,Event) %>%
  mutate(Nb_species_within_events = n_distinct(Species_name))

Count_synteny_tab<-Env_table2[Env_table2$Nb_species_within_events>1,]
#Sort by PASS2
Count_synteny_tab<-Count_synteny_tab[order(Count_synteny_tab$Synteny_PASS2),]
#Remove redondancy
Count_synteny_tab<-Count_synteny_tab[!duplicated(Count_synteny_tab[c("Clustername", "Event")]), ]
#Count total number of shared Events
length(unique(Count_synteny_tab[Count_synteny_tab$Synteny_PASS2=="YES",]$ID))-1

#Count only one EVE/CLusters and species
#Count only one EVE / Clusters and events
Env_table2bis <- Env_table2[!duplicated(Env_table2[c('Clustername','Species_name')]),]
Env_table2bis <- Env_table2bis[!duplicated(Env_table2bis[c('Clustername','Event')]),]
Env_table2bis <- Env_table2bis[order(Env_table2bis$TPM_all, decreasing = TRUE), ]

Env_table2$genomic_structure[is.na(Env_table2$genomic_structure)]<-'Unclassified'
#Get statistics for each genomic structures
Nb_para_EVE_dsDNA_A =length(unique(Env_table2[Env_table2$genomic_structure=="dsDNA" & Env_table2$Scaffold_score=="A",]$New_query_bis2 ))
Nb_para_EVE_dsDNA_B =length(unique(Env_table2[Env_table2$genomic_structure=="dsDNA" & Env_table2$Scaffold_score=="B",]$New_query_bis2 ))
Nb_para_EVE_dsDNA_C =length(unique(Env_table2[Env_table2$genomic_structure=="dsDNA" & Env_table2$Scaffold_score=="C",]$New_query_bis2 ))
Nb_para_EVE_dsDNA_D =length(unique(Env_table2[Env_table2$genomic_structure=="dsDNA" & Env_table2$Scaffold_score=="D",]$New_query_bis2 ))

Nb_para_dEVE_dsDNA =length(unique(Env_table2[Env_table2$FDR_pvalue_dNdS == 1 & Env_table2$genomic_structure=="dsDNA" & as.numeric(as.character(Env_table2$Mean_dNdS))+as.numeric(as.character(Env_table2$SE_dNdS)) < 1 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1| Env_table2$genomic_structure=="dsDNA" &  Env_table2$TPM_all>=1000 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1,]$New_query_bis2 ))

Max_pident_ALL=max(Env_table2$pident,na.rm=T)
Min_pident_ALL=min(Env_table2$pident,na.rm=T)

Mean_pident_dsDNA=mean(Env_table2[Env_table2$genomic_structure=="dsDNA",]$pident,na.rm=T)
SD_pident_dsDNA=sd(Env_table2[Env_table2$genomic_structure=="dsDNA",]$pident,na.rm=T)
Mean_pident_dsDNA<- format(round(Mean_pident_dsDNA, 2), nsmall = 2)
SD_pident_dsDNA<- format(round(SD_pident_dsDNA, 2), nsmall = 2)
Nb_EVE_dsDNA_pident_ge_80=length(unique(Env_table2[Env_table2$genomic_structure=="dsDNA" & Env_table2$pident >=80,]$New_query_bis2))
Nb_EVE_dsDNA_evalue_less_exp20=length(unique(Env_table2[Env_table2$genomic_structure=="dsDNA" & Env_table2$evalue <= 0.00000000000000000001,]$New_query_bis2))
Nb_EVE_dsDNA=length(unique((Env_table2bis[Env_table2bis$genomic_structure=="dsDNA",]$New_query_bis2)))

Nb_EVEs_dsDNA_complet_ORF= length(unique((Env_table2[Env_table2$genomic_structure=="dsDNA" & Env_table2$len_ORF > 36,]$New_query_bis2)))
Nb_EVEs_dsDNA_complet_ORF_ge_80= length(unique((Env_table2[Env_table2$genomic_structure=="dsDNA" & Env_table2$ORF_perc >= 80,]$New_query_bis2)))
Nb_dEVE_dsDNA= length(unique(Env_table2bis[Env_table2bis$genomic_structure=="dsDNA" & Env_table2bis$FDR_pvalue_dNdS == 1 & as.numeric(as.character(Env_table2bis$Mean_dNdS))+as.numeric(as.character(Env_table2bis$SE_dNdS)) < 1 & Env_table2bis$pseudogenized ==0  & Env_table2bis$ORF_perc>1 | Env_table2bis$genomic_structure=="dsDNA" &  Env_table2bis$TPM_all>=1000 & Env_table2bis$pseudogenized ==0 & Env_table2bis$ORF_perc>1,]$New_query_bis2))
Nb_viral_families_dsDNA= length(unique((Env_table2[Env_table2$genomic_structure=="dsDNA" & Env_table2$consensus_family_ID !="Unknown" & ! is.na(Env_table2$consensus_family_ID),]$consensus_family_ID)))
Nb_hymeno_families_dsDNA= length(unique((Env_table2[Env_table2$genomic_structure=="dsDNA",]$Family)))
Nb_Clusters_dsDNA= length(unique((Env_table2[Env_table2$genomic_structure=="dsDNA",]$Clustername)))

##ssDNA
Nb_para_EVE_ssDNA_A =length(unique(Env_table2[Env_table2$genomic_structure=="ssDNA" & Env_table2$Scaffold_score=="A",]$New_query_bis2 ))
Nb_para_EVE_ssDNA_B =length(unique(Env_table2[Env_table2$genomic_structure=="ssDNA" & Env_table2$Scaffold_score=="B",]$New_query_bis2 ))
Nb_para_EVE_ssDNA_C =length(unique(Env_table2[Env_table2$genomic_structure=="ssDNA" & Env_table2$Scaffold_score=="C",]$New_query_bis2 ))
Nb_para_EVE_ssDNA_D =length(unique(Env_table2[Env_table2$genomic_structure=="ssDNA" & Env_table2$Scaffold_score=="D",]$New_query_bis2 ))

Nb_para_dEVE_ssDNA =length(unique(Env_table2[Env_table2$FDR_pvalue_dNdS == 1 & Env_table2$genomic_structure=="ssDNA" &  as.numeric(as.character(Env_table2$Mean_dNdS))+as.numeric(as.character(Env_table2$SE_dNdS)) < 1 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1| Env_table2$genomic_structure=="ssDNA" &  Env_table2$TPM_all>=1000 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1,]$New_query_bis2 ))

Mean_pident_ssDNA=mean(Env_table2[Env_table2$genomic_structure=="ssDNA",]$pident,na.rm=T)
SD_pident_ssDNA=sd(Env_table2[Env_table2$genomic_structure=="ssDNA",]$pident,na.rm=T)
Mean_pident_ssDNA<- format(round(Mean_pident_ssDNA, 2), nsmall = 2)
SD_pident_ssDNA<- format(round(SD_pident_ssDNA, 2), nsmall = 2)
Nb_EVE_ssDNA_pident_ge_80=length(unique(Env_table2[Env_table2$genomic_structure=="ssDNA" & Env_table2$pident >=80,]$New_query_bis2))
Nb_EVE_ssDNA_evalue_less_exp20=length(unique(Env_table2[Env_table2$genomic_structure=="ssDNA" & Env_table2$evalue <= 0.00000000000000000001,]$New_query_bis2))
Nb_EVE_ssDNA=length(unique((Env_table2bis[Env_table2bis$genomic_structure=="ssDNA",]$New_query_bis2)))

Nb_EVEs_ssDNA_complet_ORF= length(unique((Env_table2[Env_table2$genomic_structure=="ssDNA" & Env_table2$len_ORF > 36,]$New_query_bis2)))
Nb_EVEs_ssDNA_complet_ORF_ge_80= length(unique((Env_table2[Env_table2$genomic_structure=="ssDNA" & Env_table2$ORF_perc >= 80,]$New_query_bis2)))
Nb_dEVE_ssDNA= length(unique(Env_table2bis[Env_table2bis$genomic_structure=="ssDNA" & Env_table2bis$FDR_pvalue_dNdS == 1 & as.numeric(as.character(Env_table2bis$Mean_dNdS))+as.numeric(as.character(Env_table2bis$SE_dNdS)) < 1 & Env_table2bis$pseudogenized ==0 & Env_table2bis$ORF_perc>1 | Env_table2bis$genomic_structure=="ssDNA" &  Env_table2bis$TPM_all>=1000 & Env_table2bis$pseudogenized ==0 & Env_table2bis$ORF_perc>1,]$New_query_bis2))
Nb_viral_families_ssDNA= length(unique((Env_table2[Env_table2$genomic_structure=="ssDNA" & Env_table2$consensus_family_ID !="Unknown" & ! is.na(Env_table2$consensus_family_ID),]$consensus_family_ID)))
Nb_hymeno_families_ssDNA= length(unique((Env_table2[Env_table2$genomic_structure=="ssDNA",]$Family)))
Nb_Clusters_ssDNA= length(unique((Env_table2[Env_table2$genomic_structure=="ssDNA",]$Clustername)))

##dsRNA
Nb_para_EVE_dsRNA_A =length(unique(Env_table2[Env_table2$genomic_structure=="dsRNA" & Env_table2$Scaffold_score=="A",]$New_query_bis2 ))
Nb_para_EVE_dsRNA_B =length(unique(Env_table2[Env_table2$genomic_structure=="dsRNA" & Env_table2$Scaffold_score=="B",]$New_query_bis2 ))
Nb_para_EVE_dsRNA_C =length(unique(Env_table2[Env_table2$genomic_structure=="dsRNA" & Env_table2$Scaffold_score=="C",]$New_query_bis2 ))
Nb_para_EVE_dsRNA_D =length(unique(Env_table2[Env_table2$genomic_structure=="dsRNA" & Env_table2$Scaffold_score=="D",]$New_query_bis2 ))

Nb_para_dEVE_dsRNA =length(unique(Env_table2[Env_table2$FDR_pvalue_dNdS == 1 & Env_table2$genomic_structure=="dsRNA" &  as.numeric(as.character(Env_table2$Mean_dNdS))+as.numeric(as.character(Env_table2$SE_dNdS)) < 1 & Env_table2$pseudogenized ==0 | Env_table2$genomic_structure=="dsRNA" &  Env_table2$TPM_all>=1000 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1,]$New_query_bis2 ))


Mean_pident_dsRNA=mean(Env_table2[Env_table2$genomic_structure=="dsRNA",]$pident,na.rm=T)
SD_pident_dsRNA=sd(Env_table2[Env_table2$genomic_structure=="dsRNA",]$pident,na.rm=T)
Mean_pident_dsRNA<- format(round(Mean_pident_dsRNA, 2), nsmall = 2)
SD_pident_dsRNA<- format(round(SD_pident_dsRNA, 2), nsmall = 2)
Nb_EVE_dsRNA_pident_ge_80=length(unique(Env_table2[Env_table2$genomic_structure=="dsRNA" & Env_table2$pident >=80,]$New_query_bis2))
Nb_EVE_dsRNA_evalue_less_exp20=length(unique(Env_table2[Env_table2$genomic_structure=="dsRNA" & Env_table2$evalue <= 0.00000000000000000001,]$New_query_bis2))
Nb_EVE_dsRNA=length(unique((Env_table2bis[Env_table2bis$genomic_structure=="dsRNA",]$New_query_bis2)))

Nb_EVEs_dsRNA_complet_ORF= length(unique((Env_table2[Env_table2$genomic_structure=="dsRNA" & Env_table2$len_ORF > 36,]$New_query_bis2)))
Nb_EVEs_dsRNA_complet_ORF_ge_80= length(unique((Env_table2[Env_table2$genomic_structure=="dsRNA" & Env_table2$ORF_perc >= 80,]$New_query_bis2)))
Nb_dEVE_dsRNA= length(unique(Env_table2bis[Env_table2bis$genomic_structure=="dsRNA" & Env_table2bis$FDR_pvalue_dNdS == 1 & as.numeric(as.character(Env_table2bis$Mean_dNdS))+as.numeric(as.character(Env_table2bis$SE_dNdS)) < 1 & Env_table2bis$pseudogenized ==0 & Env_table2bis$ORF_perc>1 | Env_table2bis$genomic_structure=="dsRNA" &  Env_table2bis$TPM_all>=1000 & Env_table2bis$pseudogenized ==0 & Env_table2bis$ORF_perc>1,]$New_query_bis2))
Nb_viral_families_dsRNA= length(unique((Env_table2[Env_table2$genomic_structure=="dsRNA" & Env_table2$consensus_family_ID !="Unknown" & ! is.na(Env_table2$consensus_family_ID),]$consensus_family_ID)))
Nb_hymeno_families_dsRNA= length(unique((Env_table2[Env_table2$genomic_structure=="dsRNA",]$Family)))
Nb_Clusters_dsRNA= length(unique((Env_table2[Env_table2$genomic_structure=="dsRNA",]$Clustername)))

##ssRNA
Nb_para_EVE_ssRNA_A =length(unique(Env_table2[Env_table2$genomic_structure=="ssRNA" & Env_table2$Scaffold_score=="A",]$New_query_bis2 ))
Nb_para_EVE_ssRNA_B =length(unique(Env_table2[Env_table2$genomic_structure=="ssRNA" & Env_table2$Scaffold_score=="B",]$New_query_bis2 ))
Nb_para_EVE_ssRNA_C =length(unique(Env_table2[Env_table2$genomic_structure=="ssRNA" & Env_table2$Scaffold_score=="C",]$New_query_bis2 ))
Nb_para_EVE_ssRNA_D =length(unique(Env_table2[Env_table2$genomic_structure=="ssRNA" & Env_table2$Scaffold_score=="D",]$New_query_bis2 ))

Nb_para_dEVE_ssRNA =length(unique(Env_table2[Env_table2$FDR_pvalue_dNdS == 1 &Env_table2$genomic_structure=="ssRNA" &  as.numeric(as.character(Env_table2$Mean_dNdS))+as.numeric(as.character(Env_table2$SE_dNdS)) < 1 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1 | Env_table2$genomic_structure=="ssRNA" &  Env_table2$TPM_all>=1000 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1,]$New_query_bis2 ))


Mean_pident_ssRNA=mean(Env_table2[Env_table2$genomic_structure=="ssRNA",]$pident,na.rm=T)
SD_pident_ssRNA=sd(Env_table2[Env_table2$genomic_structure=="ssRNA",]$pident,na.rm=T)
Mean_pident_ssRNA<- format(round(Mean_pident_ssRNA, 2), nsmall = 2)
SD_pident_ssRNA<- format(round(SD_pident_ssRNA, 2), nsmall = 2)
Nb_EVE_ssRNA_pident_ge_80=length(unique(Env_table2[Env_table2$genomic_structure=="ssRNA" & Env_table2$pident >=80,]$New_query_bis2))
Nb_EVE_ssRNA_evalue_less_exp20=length(unique(Env_table2[Env_table2$genomic_structure=="ssRNA" & Env_table2$evalue <= 0.00000000000000000001,]$New_query_bis2))
Nb_EVE_ssRNA=length(unique((Env_table2bis[Env_table2bis$genomic_structure=="ssRNA",]$New_query_bis2)))

Nb_EVEs_ssRNA_complet_ORF= length(unique((Env_table2[Env_table2$genomic_structure=="ssRNA" & Env_table2$len_ORF > 36,]$New_query_bis2)))
Nb_EVEs_ssRNA_complet_ORF_ge_80= length(unique((Env_table2[Env_table2$genomic_structure=="ssRNA" & Env_table2$ORF_perc >= 80,]$New_query_bis2)))
Nb_dEVE_ssRNA= length(unique(Env_table2bis[Env_table2bis$genomic_structure=="ssRNA" & Env_table2bis$FDR_pvalue_dNdS == 1 & as.numeric(as.character(Env_table2bis$Mean_dNdS))+as.numeric(as.character(Env_table2bis$SE_dNdS)) < 1 & Env_table2bis$pseudogenized ==0 &Env_table2bis$ORF_perc>1 | Env_table2bis$genomic_structure=="ssRNA" &  Env_table2bis$TPM_all>=1000 & Env_table2bis$pseudogenized ==0 &Env_table2bis$ORF_perc>1 ,]$New_query_bis2))
Nb_viral_families_ssRNA= length(unique((Env_table2[Env_table2$genomic_structure=="ssRNA" & Env_table2$consensus_family_ID !="Unknown" & ! is.na(Env_table2$consensus_family_ID),]$consensus_family_ID)))
Nb_hymeno_families_ssRNA= length(unique((Env_table2[Env_table2$genomic_structure=="ssRNA",]$Family)))
Nb_Clusters_ssRNA= length(unique((Env_table2[Env_table2$genomic_structure=="ssRNA",]$Clustername)))

##Unclassified

Env_table2bis$genomic_structure[is.na(Env_table2bis$genomic_structure)]<-"Unclassified"

Nb_para_EVE_Unclassified_A =length(unique(Env_table2[Env_table2$genomic_structure=="Unclassified" & Env_table2$Scaffold_score=="A",]$New_query_bis2 ))
Nb_para_EVE_Unclassified_B =length(unique(Env_table2[Env_table2$genomic_structure=="Unclassified" & Env_table2$Scaffold_score=="B",]$New_query_bis2 ))
Nb_para_EVE_Unclassified_C =length(unique(Env_table2[Env_table2$genomic_structure=="Unclassified" & Env_table2$Scaffold_score=="C",]$New_query_bis2 ))
Nb_para_EVE_Unclassified_D =length(unique(Env_table2[Env_table2$genomic_structure=="Unclassified" & Env_table2$Scaffold_score=="D",]$New_query_bis2 ))

Nb_para_dEVE_Unclassified =length(unique(Env_table2[Env_table2$FDR_pvalue_dNdS == 1 & Env_table2$genomic_structure=="Unclassified" &  as.numeric(as.character(Env_table2$Mean_dNdS))+as.numeric(as.character(Env_table2$SE_dNdS)) < 1 & Env_table2$pseudogenized ==0  & Env_table2$ORF_perc>1 | Env_table2$genomic_structure=="Unclassified" &  Env_table2$TPM_all>=1000 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1,]$New_query_bis2 ))


Mean_pident_Unclassified=mean(Env_table2[Env_table2$genomic_structure=="Unclassified",]$pident,na.rm=T)
SD_pident_Unclassified=sd(Env_table2[Env_table2$genomic_structure=="Unclassified",]$pident,na.rm=T)
Mean_pident_Unclassified<- format(round(Mean_pident_Unclassified, 2), nsmall = 2)
SD_pident_Unclassified<- format(round(SD_pident_Unclassified, 2), nsmall = 2)
Nb_EVE_Unclassified_pident_ge_80=length(unique(Env_table2[Env_table2$genomic_structure=="Unclassified" & Env_table2$pident >=80,]$New_query_bis2))
Nb_EVE_Unclassified_evalue_less_exp20=length(unique(Env_table2[Env_table2$genomic_structure=="Unclassified" & Env_table2$evalue <= 0.00000000000000000001,]$New_query_bis2))
Nb_EVE_Unclassified=length(unique((Env_table2bis[Env_table2bis$genomic_structure=="Unclassified",]$New_query_bis2)))
Unknown

Nb_EVEs_Unclassified_complet_ORF= length(unique((Env_table2[Env_table2$genomic_structure=="Unclassified" & Env_table2$len_ORF > 36,]$New_query_bis2)))
Nb_EVEs_Unclassified_complet_ORF_ge_80= length(unique((Env_table2[Env_table2$genomic_structure=="Unclassified" & Env_table2$ORF_perc >= 80,]$New_query_bis2)))
Nb_dEVE_Unclassified= length(unique(Env_table2bis[Env_table2bis$genomic_structure=="Unclassified" & Env_table2bis$FDR_pvalue_dNdS == 1 & as.numeric(as.character(Env_table2bis$Mean_dNdS))+as.numeric(as.character(Env_table2bis$SE_dNdS)) < 1 & Env_table2bis$pseudogenized ==0 & Env_table2bis$ORF_perc>1| Env_table2bis$genomic_structure=="Unclassified" &  Env_table2bis$TPM_all>=1000 & Env_table2bis$pseudogenized ==0 & Env_table2bis$ORF_perc>1,]$New_query_bis2))
Nb_viral_families_Unclassified= length(unique((Env_table2[Env_table2$genomic_structure=="Unclassified" & Env_table2$consensus_family_ID !="Unknown" & ! is.na(Env_table2$consensus_family_ID),]$consensus_family_ID)))
Nb_hymeno_families_Unclassified= length(unique((Env_table2[Env_table2$genomic_structure=="Unclassified",]$Family)))
Nb_Clusters_Unclassified= length(unique((Env_table2[Env_table2$genomic_structure=="Unclassified",]$Clustername)))


statistics_table=data.frame (Type= c("dsDNA", 'ssDNA','dsRNA',"ssRNA","Unclassified RNA","Total"),
                             Raw_EVEs= c(paste0(Nb_para_EVE_dsDNA_A,"|",Nb_para_EVE_dsDNA_B,"|",Nb_para_EVE_dsDNA_C,"|",Nb_para_EVE_dsDNA_D),
                                         paste0(Nb_para_EVE_ssDNA_A,"|",Nb_para_EVE_ssDNA_B,"|",Nb_para_EVE_ssDNA_C,"|",Nb_para_EVE_ssDNA_D),
                                         paste0(Nb_para_EVE_dsRNA_A,"|",Nb_para_EVE_dsRNA_B,"|",Nb_para_EVE_dsRNA_C,"|",Nb_para_EVE_dsRNA_D),
                                         paste0(Nb_para_EVE_ssRNA_A,"|",Nb_para_EVE_ssRNA_B,"|",Nb_para_EVE_ssRNA_C,"|",Nb_para_EVE_ssRNA_D),
                                         paste0(Nb_para_EVE_Unclassified_A,"|",Nb_para_EVE_Unclassified_B,"|",Nb_para_EVE_Unclassified_C,"|",Nb_para_EVE_Unclassified_D),
                                         sum(Nb_para_EVE_dsDNA_A,Nb_para_EVE_dsDNA_B,Nb_para_EVE_dsDNA_C,Nb_para_EVE_dsDNA_D,Nb_para_EVE_ssDNA_A,Nb_para_EVE_ssDNA_B,Nb_para_EVE_ssDNA_C,Nb_para_EVE_ssDNA_D,Nb_para_EVE_ssRNA_A,Nb_para_EVE_ssRNA_B,Nb_para_EVE_ssRNA_C,Nb_para_EVE_ssRNA_D,Nb_para_EVE_dsRNA_A,Nb_para_EVE_dsRNA_B,Nb_para_EVE_dsRNA_C,Nb_para_EVE_dsRNA_D,Nb_para_EVE_Unclassified_A,Nb_para_EVE_Unclassified_B,Nb_para_EVE_Unclassified_C,Nb_para_EVE_Unclassified_D)),
                             Nb_Clusters=c(Nb_Clusters_dsDNA,Nb_Clusters_ssDNA,Nb_Clusters_dsRNA,Nb_Clusters_ssRNA,Nb_Clusters_Unclassified,sum(Nb_Clusters_dsDNA,Nb_Clusters_ssDNA,Nb_Clusters_dsRNA,Nb_Clusters_ssRNA,Nb_Clusters_Unclassified)),
                             Nb_EVEs = c(Nb_EVE_dsDNA,Nb_EVE_ssDNA,Nb_EVE_dsRNA,Nb_EVE_ssRNA,Nb_EVE_Unclassified, sum(Nb_EVE_dsDNA,Nb_EVE_ssDNA,Nb_EVE_dsRNA,Nb_EVE_ssRNA,Nb_EVE_Unclassified)),
                             Mean_pident  = c(paste0(Mean_pident_dsDNA," (",SD_pident_dsDNA,")")  ,paste0(Mean_pident_ssDNA," (",SD_pident_ssDNA,")"),paste0(Mean_pident_dsRNA," (",SD_pident_dsRNA,")"),paste0(Mean_pident_ssRNA," (",SD_pident_ssRNA,")"),paste0(Mean_pident_Unclassified," (",SD_pident_Unclassified,")"), mean(c(as.numeric(Mean_pident_dsDNA),as.numeric(Mean_pident_ssDNA),as.numeric(Mean_pident_dsRNA),as.numeric(Mean_pident_ssRNA),as.numeric(Mean_pident_Unclassified)))),
                             Nb_EVEs_pident_ge_80 =c(Nb_EVE_dsDNA_pident_ge_80,Nb_EVE_ssDNA_pident_ge_80,Nb_EVE_dsRNA_pident_ge_80,Nb_EVE_ssRNA_pident_ge_80,Nb_EVE_Unclassified_pident_ge_80 , sum(Nb_EVE_dsDNA_pident_ge_80,Nb_EVE_ssDNA_pident_ge_80,Nb_EVE_dsRNA_pident_ge_80,Nb_EVE_ssRNA_pident_ge_80,Nb_EVE_Unclassified_pident_ge_80)),
                             Nb_EVEs_evalue_le_exp20 =c(Nb_EVE_dsDNA_evalue_less_exp20,Nb_EVE_ssDNA_evalue_less_exp20,Nb_EVE_dsRNA_evalue_less_exp20,Nb_EVE_ssRNA_evalue_less_exp20,Nb_EVE_Unclassified_evalue_less_exp20, sum(Nb_EVE_dsDNA_evalue_less_exp20,Nb_EVE_ssDNA_evalue_less_exp20,Nb_EVE_dsRNA_evalue_less_exp20,Nb_EVE_ssRNA_evalue_less_exp20,Nb_EVE_Unclassified_evalue_less_exp20)),
                             Nb_Events_EVEs=c(Events_EVEs_dsDNA,Events_EVEs_ssDNA,Events_EVEs_dsRNA,Events_EVEs_ssRNA,Events_EVEs_Unclassified, sum(Events_EVEs_dsDNA,Events_EVEs_ssDNA,Events_EVEs_dsRNA,Events_EVEs_ssRNA,Events_EVEs_Unclassified)),
                             Nb_shared_Events_EVEs=c(sum(Event_shared_df$Nb_EVEs_dsDNA_Events,na.rm=T),sum(Event_shared_df$Nb_EVEs_ssDNA_Events,na.rm=T),sum(Event_shared_df$Nb_EVEs_dsRNA_Events,na.rm=T),sum(Event_shared_df$Nb_EVEs_ssRNA_Events,na.rm=T),sum(Event_shared_df$Nb_EVEs_Unclassified_Events,na.rm=T), sum(sum(Event_shared_df$Nb_EVEs_dsDNA_Events,na.rm=T),sum(Event_shared_df$Nb_EVEs_ssDNA_Events,na.rm=T),sum(Event_shared_df$Nb_EVEs_dsRNA_Events,na.rm=T),sum(Event_shared_df$Nb_EVEs_ssRNA_Events,na.rm=T),sum(Event_shared_df$Nb_EVEs_Unclassified_Events,na.rm=T))),
                             Nb_viral_families= c(Nb_viral_families_dsDNA,Nb_viral_families_ssDNA,Nb_viral_families_dsRNA,Nb_viral_families_ssRNA,Nb_viral_families_Unclassified,sum(Nb_viral_families_dsDNA,Nb_viral_families_ssDNA,Nb_viral_families_dsRNA,Nb_viral_families_ssRNA,Nb_viral_families_Unclassified)),
                             Nb_Hymenoptera_families = c(Nb_hymeno_families_dsDNA,Nb_hymeno_families_ssDNA,Nb_hymeno_families_dsRNA,Nb_hymeno_families_ssRNA,Nb_hymeno_families_Unclassified,sum(Nb_hymeno_families_dsDNA,Nb_hymeno_families_ssDNA,Nb_hymeno_families_dsRNA,Nb_hymeno_families_ssRNA,Nb_hymeno_families_Unclassified)),
                             Raw_dEVEs=c(Nb_para_dEVE_dsDNA,Nb_para_dEVE_ssDNA,Nb_para_dEVE_dsRNA,Nb_para_dEVE_ssRNA,Nb_para_dEVE_Unclassified,sum(Nb_para_dEVE_dsDNA,Nb_para_dEVE_ssDNA,Nb_para_dEVE_dsRNA,Nb_para_dEVE_ssRNA,Nb_para_dEVE_Unclassified)),
                             Nb_dEVEs = c(Nb_dEVE_dsDNA,Nb_dEVE_ssDNA,Nb_dEVE_dsRNA,Nb_dEVE_ssRNA,Nb_dEVE_Unclassified,sum(Nb_dEVE_dsDNA,Nb_dEVE_ssDNA,Nb_dEVE_dsRNA,Nb_dEVE_ssRNA,Nb_dEVE_Unclassified)),
                             Nb_EVEs_with_ORFs = c(Nb_EVEs_dsDNA_complet_ORF,Nb_EVEs_ssDNA_complet_ORF,Nb_EVEs_dsRNA_complet_ORF,Nb_EVEs_ssRNA_complet_ORF,Nb_EVEs_Unclassified_complet_ORF,sum(Nb_EVEs_dsDNA_complet_ORF,Nb_EVEs_ssDNA_complet_ORF,Nb_EVEs_dsRNA_complet_ORF,Nb_EVEs_ssRNA_complet_ORF,Nb_EVEs_Unclassified_complet_ORF)),
                             Nb_EVEs_with_ORFS_ge_80 = c(Nb_EVEs_dsDNA_complet_ORF_ge_80,Nb_EVEs_ssDNA_complet_ORF_ge_80,Nb_EVEs_dsRNA_complet_ORF_ge_80,Nb_EVEs_ssRNA_complet_ORF_ge_80,Nb_EVEs_Unclassified_complet_ORF_ge_80,sum(Nb_EVEs_dsDNA_complet_ORF_ge_80,Nb_EVEs_ssDNA_complet_ORF_ge_80,Nb_EVEs_dsRNA_complet_ORF_ge_80,Nb_EVEs_ssRNA_complet_ORF_ge_80,Nb_EVEs_Unclassified_complet_ORF_ge_80)),
                             Nb_Events_dEVEs=c(Events_dEVEs_dsDNA,Events_dEVEs_ssDNA,Events_dEVEs_dsRNA,Events_dEVEs_ssRNA,Events_dEVEs_Unclassified,sum(Events_dEVEs_dsDNA,Events_dEVEs_ssDNA,Events_dEVEs_dsRNA,Events_dEVEs_ssRNA,Events_dEVEs_Unclassified)))


statistics_table<-as.data.frame(t(statistics_table),)
colnames(statistics_table) <- statistics_table[1,]
statistics_table<-statistics_table[-1,]
write.table(statistics_table,"/Users/bguinet/Desktop/Papier_scientifique/Statistics_endogenization_table.txt",sep=";")



#Get count of events per species
#Median  & Max event per species
median(as.data.frame(table(aggregate(Env_table2, by=list(Env_table2$query,Env_table2$ID), FUN=length)$Group.1))$Freq)
max(as.data.frame(table(aggregate(Env_table2, by=list(Env_table2$query,Env_table2$ID), FUN=length)$Group.1))$Freq)
#Nb species with events
nrow(as.data.frame(table(aggregate(Env_table2, by=list(Env_table2$query,Env_table2$ID), FUN=length)$Group.1)))

#Nb dEVEs with complet ORFs
Env_table2[Env_table2$FDR_pvalue_dNdS == 1 & as.numeric(as.character(Env_table2$Mean_dNdS))+as.numeric(as.character(Env_table2$SE_dNdS)) < 1 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1| Env_table2$TPM_all>=1000 & Env_table2$pseudogenized ==0 & Env_table2$ORF_perc>1,]
length(unique(Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),][Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),]$ORF_perc >=1 & Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),]$FDR_pvalue_dNdS == 1 & as.numeric(as.character(Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),]$Mean_dNdS))+as.numeric(as.character(Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),]$SE_dNdS)) < 1 & Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),]$pseudogenized ==0 & Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),]$ORF_perc>1 |    Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),]$TPM_all>=1000 & Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),]$pseudogenized ==0 & Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),]$ORF_perc>1,]$New_query_bis2))


# Add classification
Env_table2$consensus_family_ID[Env_table2$consensus_family_ID=="Unknown" & Env_table2$consensus_genomic_structure=="ssRNA"]<-"Unclassified_ssRNA"
Env_table2$consensus_family_ID[Env_table2$consensus_family_ID=="Unknown" & Env_table2$consensus_genomic_structure=="dsRNA"]<-"Unclassified_dsRNA"
Env_table2$consensus_family_ID[is.na(Env_table2$consensus_family_ID) & Env_table2$consensus_genomic_structure=="ssDNA"]<-"Unclassified_ssDNA"

Env_table2$consensus_family[is.na(Env_table2$consensus_family) & Env_table2$genomic_structure == "ssDNA"]<-"Unclassified ssDNA"
Env_table2$consensus_family[is.na(Env_table2$consensus_family) & Env_table2$genomic_structure == "ssRNA"]<-"Unclassified ssRNA"
Env_table2$consensus_family[is.na(Env_table2$consensus_family) & Env_table2$genomic_structure == "dsRNA"]<-"Unclassified dsRNA"

Env_table2$consensus_family[Env_table2$consensus_family=="Unknown" & Env_table2$genomic_structure == "ssDNA"]<-"Unclassified ssDNA"
Env_table2$consensus_family[Env_table2$consensus_family=="Unknown" & Env_table2$genomic_structure == "ssRNA"]<-"Unclassified ssRNA"
Env_table2$consensus_family[Env_table2$consensus_family=="Unknown" & Env_table2$genomic_structure == "dsRNA"]<-"Unclassified dsRNA"


Env_table2$consensus_family_ID[is.na(Env_table2$consensus_family_ID) & is.na(Env_table2$consensus_genomic_structure)]<-"Unclassified"
Env_table2$genomic_structure[is.na(Env_table2$genomic_structure)]<-"Unclassified"
Env_table2$genomic_structure[is.na(Env_table2$consensus_genomic_structure)]<-"Unclassified"
Env_table2$consensus_genomic_structure[Env_table2$consensus_family_ID =="Unclassified"]<-"Unclassified"

Env_table2$consensus_genomic_structure[is.na(Env_table2$consensus_family_ID) & is.na(Env_table2$consensus_genomic_structure)]<-"Unclassified"

write.table(Env_table2,"/Users/bguinet/Desktop/Papier_scientifique/Env_table2.txt",sep=";")
write.table(Env_table2,"/Users/bguinet/Desktop/Supplementary_material_Current_Biology/All_EVEs_dEVEs_informations.txt",sep=";")



#Possible bias in EVE detection metrics

ggplot(Env_table2, aes( y=alnlen)) +
  geom_boxplot(notch=TRUE)

nrow(Env_table2[Env_table2$alnlen<=50,])
mean(Env_table2[Env_table2$alnlen<=50,]$tlen,na.rm=TRUE)
mean(Env_table2$alnlen,na.rm=TRUE)
median(Env_table2$alnlen,na.rm=TRUE)

median(Env_table2$alnlen / Env_table2$tlen,na.rm=TRUE)
min(Env_table2$alnlen / Env_table2$tlen,na.rm=TRUE)
sd(Env_table2$alnlen / Env_table2$tlen,na.rm=TRUE)
sd(Env_table2$alnlen,na.rm=TRUE)

#Count the number of different hymenoptera specie within clusters
library(dplyr)
data_count_2 <- Env_table2%>%                              # Applying group_by & summarise
  group_by(Clustername) %>%
  summarise(count = n_distinct(Species_name))



#Possible bias in dEVE detection metrics
length(unique(Env_table2[((Env_table2$Mean_dNdS  +Env_table2$SE_dNdS) < 1) & Env_table2$Pvalue_dNdS>0.05 & Env_table2$TPM_all >1000,]$ID))

length(unique(Env_table2[((Env_table2$Mean_dNdS  +Env_table2$SE_dNdS) > 1) & Env_table2$Pvalue_dNdS<0.05 & Env_table2$TPM_all >1000,]$ID))

#python

# TPM > 1000 et dN/dS < 1

tab1=tab.loc[tab['Mean_dNdS'].lt(0.8) & tab['Pvalue_dNdS'].lt(0.05) & tab['TPM_all'].gt(1000)]
Count_TPMgt1000_and_dNdS_less1=len(tab1.drop_duplicates(['Clustername','Event'],keep= 'first')]

tab2=tab.loc[tab['Mean_dNdS'].lt(0.8) & tab['Pvalue_dNdS'].lt(0.05) & tab['TPM_all'].lt(1000)]
Count_TPMlt1000_and_dNdS_less1=len(tab2.drop_duplicates(['Clustername','Event'],keep= 'first'))

tab3=tab.loc[tab['Mean_dNdS'].gt(0.8) & tab['Pvalue_dNdS'].gt(0.05) & tab['TPM_all'].gt(1000)]
Count_TPMgt1000_and_dNdS_equal1=len(tab3.drop_duplicates(['Clustername','Event'],keep= 'first'))

tab4=tab.loc[tab['Mean_dNdS'].gt(0.8) & tab['Pvalue_dNdS'].gt(0.05) & tab['TPM_all'].lt(1000)]
Count_TPMlt1000_and_dNdS_equal1=len(tab4.drop_duplicates(['Clustername','Event'],keep= 'first'))



#OUtput for Sebastian Lequime

taxo<-read.table("/Users/bguinet/Desktop/Papier_scientifique/Species_sub_sup_families_order.txt", header=T)
taxo<-taxo[taxo$SuperFamily=="Formicoidea",]
#Only focus on ants species
lequimeEnv_table2<- Env_table2[Env_table2$Species_name %in% taxo$Species,]

lequimeEnv_table2<-lequimeEnv_table2[c("ID","Clustername","Event","Species_name","New_query_bis2","start","end","Scaffold_species","target","evalue","bits","alnlen","family","consensus_family","genomic_structure","Protein.names","GC_content_BUSCO","GC_content_scaffold","pvalue_gc","cov_depth_BUSCO","cov_depth_candidat","pvalue_cov","count_repeat","count_eucaryote","Scaffold_score","pseudogenized","TPM_all","Mean_dNdS","Pvalue_dNdS","SE_dNdS","ORF_perc","start_ORF","end_ORF",'sequence')]
write.table(lequimeEnv_table2,file="/Users/bguinet/Desktop/Papier_scientifique/Env_table_for_sebastian.txt",sep=";")

# Count number EVEs Events per species

Nb_Events_per_species<-Env_table2 %>%
  group_by(Species_name) %>%
  summarise(count = n_distinct(ID))

dim(Nb_Events_per_species)
median(Nb_Events_per_species$count)
max(Nb_Events_per_species$count)
min(Nb_Events_per_species$count)


#TE file output
TE_dataframe <- data.frame(matrix(ncol = 19, nrow = 0))
x <- c("query","target","pident","alnlen","mismatch","gapopen","qstart","qend","tstart","tend","evalue","bits",
       "qlen","tlen","tcov","strand","Newqstart" ,"Newqend","qmatchlen")
colnames(TE_dataframe) <- x
for (species in unique(Env_table2$Species_name)){
  if (species =="Unknown_sp"){
    subTE_table<- read.table(paste0("/Users/bguinet/Desktop/Papier_scientifique/Repeat_env_analysis2/BlastX_BUSCO_and_EVEs_on_RepeatPeps_results_Platygaster_orseoliae.tab"),sep=",")
  }else{
    subTE_table<- read.csv(paste0("/Users/bguinet/Desktop/Papier_scientifique/Repeat_env_analysis2/BlastX_BUSCO_and_EVEs_on_RepeatPeps_results_",species,".tab"),sep=";")
  }

  #Only keep scaffold remaining in the data
  subTE_table <- subTE_table[subTE_table$query %in% Env_table2$Scaffold_species,]
  subTE_table<- subTE_table[subTE_table$evalue < 0.0000000001,]
  TE_dataframe<-rbind(TE_dataframe,subTE_table)

}

#Save TE file

TE_dataframe<-TE_dataframe[!(colnames(TE_dataframe) %in% c('Newqstart','Newqend'))]
colnames(TE_dataframe)<-c("Scaffold_species_names","TE_names","pident","alnlen","mismatch","gapopen","Scaffold_start","Scaffold_end","TE_start","TE_end","evalue","bits","Scaffold_length","TE_length","TE_coverage","strand","Scaffold_matchlen")
write.table(TE_dataframe[c("Scaffold_species_names","TE_names","pident","alnlen","mismatch","gapopen","strand","evalue","bits","Scaffold_start","Scaffold_end","TE_start","TE_end","strand","Scaffold_length","TE_length","TE_coverage","Scaffold_matchlen")],file="/Users/bguinet/Desktop/Fichiers_papier_finalises/TE_blastX_table.txt",sep=";")

#Plot Consensus Event functions


#subEnv_table2<-Env_table2[Env_table2$consensus_family=="Nudiviridae",]
subEnv_table2<-Env_table2[!duplicated(Env_table2[c("Clustername","Event")]),]


subEnv_table2 <- subEnv_table2[subEnv_table2$FDR_pvalue_dNdS == 1 & as.numeric(as.character(subEnv_table2$Mean_dNdS))+as.numeric(as.character(subEnv_table2$SE_dNdS)) < 1 & subEnv_table2$pseudogenized ==0 &subEnv_table2$ORF_perc>1 | subEnv_table2$TPM_all>=1000 & subEnv_table2$pseudogenized ==0 &subEnv_table2$ORF_perc>1 ,]

subEnv_table2<- subEnv_table2[ ! is.na(subEnv_table2$Clustername),]

subEnv_table2<- subEnv_table2[c("Consensus_function_Event","consensus_family","consensus_genomic_structure")]

subEnv_table2<-subEnv_table2 %>%
  group_by(Consensus_function_Event, consensus_family) %>%
  add_count(name = "Ndup")

subEnv_table2<-subEnv_table2[!duplicated(subEnv_table2[c("consensus_family","Consensus_function_Event")]),]


Bar_plot_families<-ggplot(subEnv_table2,aes(x = factor(consensus_family),y = Ndup,fill=Consensus_function_Event))+
  facet_grid(consensus_genomic_structure ~ . , scales="free_y",space="free") +
  geom_bar(stat="identity",position = "identity",color="black",size=0.2)+
  theme_bw() +
  theme(axis.text.x=element_text(angle=0,hjust=0.5))+  coord_flip() +theme(
    panel.grid.major.y  = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line())+
  theme(strip.text.y = element_text(size = 7,angle=0,face="bold"))+
  theme(axis.text.y = element_blank())+ ylab("Number of Events") + theme(axis.text.y = element_text( color="black",size=9))+
  theme(legend.position="bottom")+
  scale_y_continuous(breaks=seq(0,75,5)) +
  theme(axis.text.y = element_text(hjust=0.95))


###Output for dsDNA phylogeny
Alone_dsDNA_4_phylo <- Event_alone_df[Event_alone_df$Nb_EVEs >=15 & Event_alone_df$Nb_EVEs < 19,][c("species","Event")] %>%
  separate_rows(Event, sep = ',') %>%
  separate(Event, c('Event', 'Clusters'), sep = '_', convert = TRUE)


Shared_dsDNA_4_phylo <- Event_shared_df[Event_shared_df$Nb_EVEs >=5,][c("species","Event")] %>%
  separate_rows(Event, sep = ',') %>%
  separate(Event, c('Event', 'Clusters'), sep = '_', convert = TRUE)


###########################################################
#Create histogram plots

#dsDNA hitmap
library("tidyverse")


Event_alone_df$Nb_EVEs_dsDNA[ grepl( "3_Cluster24469", Event_alone_df$Event)] <- 4

type1<-append(Event_alone_df$Nb_EVEs_dsDNA , Event_shared_df$Nb_EVEs_dsDNA)
type2<-append(Event_alone_df$Nb_dEVEs_dsDNA, Event_shared_df$Nb_dEVEs_dsDNA)
type2[type2>=1]<-1

list1<- list(EVEs = type1, dEVEs = type2)


grp <- factor(with(list1, fct_collapse(as.character(EVEs),
                                       `>=5` = as.character(EVEs)[EVEs >=5])), levels = c(0:4, ">=5"))
v1 <- table(grp)
v2 <- tapply(list1$dEVEs, grp, FUN = sum)
dsDNA_hist_df<- bind_rows(list(EVEs = stack(v1)[2:1], dEVEs = stack(v2)[2:1]), .id = 'type') %>%
  filter(ind != '0')

dsDNA_hist_df[dsDNA_hist_df==0] <- NA
dsDNA_hist <- ggplot(data = dsDNA_hist_df,
                     aes(x = ind, y = values, fill = factor(type, levels = c("EVEs","dEVEs")))) +
  geom_bar(width = 0.65,
           stat = 'identity',
           position = 'identity',color="black",size=0.2)+
  #geom_label_repel(aes(label = values),
  #position = position_stack(vjust = 0.5),
  #size = 2.75,
  # colour = 'white')
  ggtitle("dsDNA") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line())+
  ylab("Number of Events") + xlab("Number of EVEs within Events")+
  scale_fill_manual(values = c("#FFD460","#F07B3F"))+
  theme(plot.title = element_text( face = "bold"))+ coord_cartesian(ylim=c(0,150))+
  theme(legend.title=element_blank()) +
  theme(
    axis.title.y = element_text(size = 17))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20))+
  theme(axis.title.x=element_blank())

#dsRNA hitmap
type1<-append(Event_alone_df$Nb_EVEs_dsRNA , Event_shared_df$Nb_EVEs_dsRNA)
type2<-append(Event_alone_df$Nb_dEVEs_dsRNA, Event_shared_df$Nb_dEVEs_dsRNA)
type2[type2>=1]<-1

list1<- list(EVEs = type1, dEVEs = type2)
grp <- factor(with(list1, fct_collapse(as.character(EVEs),
                                       `>=5` = as.character(EVEs)[EVEs >=5])), levels = c(0:4, ">=5"))
v1 <- table(grp)
v2 <- tapply(list1$dEVEs, grp, FUN = sum)
dsRNA_hist_df<- bind_rows(list(EVEs = stack(v1)[2:1], dEVEs = stack(v2)[2:1]), .id = 'type') %>%
  filter(ind != '0')

dsRNA_hist_df[dsRNA_hist_df==0] <- NA
dsRNA_hist <- ggplot(data = dsRNA_hist_df,
                     aes(x = ind, y = values, fill = factor(type, levels = c("EVEs","dEVEs")))) +
  geom_bar(width = 0.65,
           stat = 'identity',
           position = 'identity',color="black",size=0.2)+ #geom_label_repel(aes(label = values),
  #position = position_stack(vjust = 0.5),
  #size = 2.75,
  # colour = 'white')
  ggtitle("dsRNA") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line())+
  ylab("Number of Events") + xlab("Number of EVEs within Events")+
  scale_fill_manual(values = c("#FFD460","#F07B3F"))+
  theme(plot.title = element_text( face = "bold"))+coord_cartesian(ylim=c(0,150))+
  theme(legend.title=element_blank())+theme(
    axis.title.y = element_blank())+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20))+
  theme(axis.title.x=element_blank())

#ssDNA hitmap
type1<-append(Event_alone_df$Nb_EVEs_ssDNA , Event_shared_df$Nb_EVEs_ssDNA)
type2<-append(Event_alone_df$Nb_dEVEs_ssDNA, Event_shared_df$Nb_dEVEs_ssDNA)
type2[type2>=1]<-1

list1<- list(EVEs = type1, dEVEs = type2)
grp <- factor(with(list1, fct_collapse(as.character(EVEs),
                                       `>=5` = as.character(EVEs)[EVEs >=5])), levels = c(0:4, ">=5"))
v1 <- table(grp)
v2 <- tapply(list1$dEVEs, grp, FUN = sum)
ssDNA_hist_df<- bind_rows(list(EVEs = stack(v1)[2:1], dEVEs = stack(v2)[2:1]), .id = 'type') %>%
  filter(ind != '0')

ssDNA_hist_df[ssDNA_hist_df==0] <- NA
ssDNA_hist <- ggplot(data = ssDNA_hist_df,
                     aes(x = ind, y = values, fill = factor(type, levels = c("EVEs","dEVEs")))) +
  geom_bar(width = 0.65,
           stat = 'identity',
           position = 'identity',color="black",size=0.2)+ #geom_label_repel(aes(label = values),
  #position = position_stack(vjust = 0.5),
  #size = 2.75,
  # colour = 'white')
  ggtitle("ssDNA") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line())+
  ylab("Number of Events") + xlab("Number of EVEs within Events")+
  scale_fill_manual(values = c("#FFD460","#F07B3F"))+
  theme(plot.title = element_text( face = "bold"))+ coord_cartesian(ylim=c(0,150))+
  theme(legend.title=element_blank())+theme(
    axis.title.y = element_blank())+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20))+
  theme(axis.title.x=element_blank())



#ssRNA hitmap
type1<-append(Event_alone_df$Nb_EVEs_ssRNA , Event_shared_df$Nb_EVEs_ssRNA)
type2<-append(Event_alone_df$Nb_dEVEs_ssRNA, Event_shared_df$Nb_dEVEs_ssRNA)
type2[type2>=1]<-1

list1<- list(EVEs = type1, dEVEs = type2)
grp <- factor(with(list1, fct_collapse(as.character(EVEs),
                                       `>=5` = as.character(EVEs)[EVEs >=5])), levels = c(0:4, ">=5"))
v1 <- table(grp)
v2 <- tapply(list1$dEVEs, grp, FUN = sum)
ssRNA_hist_df<- bind_rows(list(EVEs = stack(v1)[2:1], dEVEs = stack(v2)[2:1]), .id = 'type') %>%
  filter(ind != '0')


ssRNA_hist_df[ssRNA_hist_df==0] <- NA
ssRNA_hist <- ggplot(data = ssRNA_hist_df,
                     aes(x = ind, y = values, fill = factor(type, levels = c("EVEs","dEVEs")))) +
  geom_bar(width = 0.65,
           stat = 'identity',
           position = 'identity',color="black",size=0.2)+ #geom_label_repel(aes(label = values),
  #position = position_stack(vjust = 0.5),
  #size = 2.75,
  # colour = 'white')
  ggtitle("ssRNA") +
  ylab("Number of Events") + xlab("Number of EVEs within Events")+
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line())+
  scale_fill_manual(values = c("#FFD460","#F07B3F"))+
  theme(plot.title = element_text( face = "bold"))+ coord_cartesian(ylim=c(0,150))+
  theme(legend.title=element_blank()) +theme(
    axis.title.y = element_blank())+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 20))+
  theme(axis.title.x=element_blank())



library(ggridges)
library(ggpubr)
library(MASS)
#test

ssRNA_hist_df$structure<- "ssRNA"
dsRNA_hist_df$structure<- "dsRNA"
ssDNA_hist_df$structure<- "ssDNA"
dsDNA_hist_df$structure<- "dsDNA"
All_hist_df <- rbind(ssRNA_hist_df,dsRNA_hist_df,ssDNA_hist_df,dsDNA_hist_df)

All_hist_df$ind=factor(All_hist_df$ind,levels=c("1","2","3","4",">=5"))
All_hist_df$structure=factor(All_hist_df$structure,levels=c("dsDNA","ssDNA","dsRNA","ssRNA"))


# Number of Events with single EVEs
sum(All_hist_df$values[All_hist_df$type=="EVEs" & All_hist_df$ind=="1"],na.rm=T)

# Number of Events with >4 EVEs
sum(All_hist_df$values[All_hist_df$type=="EVEs" & All_hist_df$ind==">=5"],na.rm=T)


ggplot(All_hist_df, aes(x = ind, y = values, fill = factor(type, levels = c("EVEs","dEVEs")))) +
  geom_bar(stat = 'identity', position = 'stack',color="black",size=0.2) +   facet_grid(~ structure) +
  scale_fill_manual(values = c("#FFD460", "#F07B3F")) +
  theme(panel.background = element_rect(fill = 'white'))+
  ylab("Numumber of Events") + xlab("Number of EVEs/dEVEs within Events") + theme(axis.text.y = element_text( color="black",size=10))+
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line())+
  theme(text = element_text(size=17)) + theme(strip.text.x = element_text(size = 13,face="bold"))+
  theme(strip.text.y = element_text(size = 13,face="bold"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


library(grid)
Nb_EVES_distribution_plot<-ggarrange(dsDNA_hist,ssDNA_hist,dsRNA_hist,ssRNA_hist,ncol=4, common.legend = TRUE, legend="right")
Nb_EVES_distribution_plot<-annotate_figure(Nb_EVES_distribution_plot,
                                           bottom = textGrob("Number of Endogenous Viral Elements (EVEs) and domesticated EVEs (dEVEs) within Events", gp = gpar(cex = 1.3)))

dev.print(device = pdf, file = paste0("/Users/bguinet/Desktop/Papier_scientifique/Hist_EVEs_dEVEs_count.pdf"), width = 10,height=10)



#Create barplot with expected number to be found for each genomic structure categoris
Genomic_structures <- c("dsDNA", "ssDNA", "dsRNA", "ssRNA","dsDNA", "ssDNA", "dsRNA", "ssRNA","dsDNA", "ssDNA", "dsRNA", "ssRNA","dsDNA", "ssDNA", "dsRNA", "ssRNA")
Type<-c("EVEs","EVEs","EVEs","EVEs","dEVEs","dEVEs","dEVEs","dEVEs","EVEs_Event","EVEs_Event","EVEs_Event","EVEs_Event","dEVEs_Event","dEVEs_Event","dEVEs_Event","dEVEs_Event")
value <- c(Nb_EVE_dsDNA,Nb_EVE_ssDNA,Nb_EVE_dsRNA,Nb_EVE_ssRNA,
           Nb_dEVE_dsDNA,Nb_dEVE_ssDNA,Nb_dEVE_dsRNA,Nb_dEVE_ssRNA,
           Events_EVEs_dsDNA,Events_EVEs_ssDNA,Events_EVEs_dsRNA,Events_EVEs_ssRNA,
           Events_dEVEs_dsDNA,Events_dEVEs_ssDNA,Events_dEVEs_dsRNA,Events_dEVEs_ssRNA)

ssRNA_perc <- 1844*100/2475
dsRNA_perc <- 401*100/2475
dsDNA_perc <- 155*100/2475
ssDNA_perc <- 75*100/2475

Expected_NB<- c( ((sum(Nb_EVE_dsDNA,Nb_EVE_ssDNA,Nb_EVE_dsRNA,Nb_EVE_ssRNA)*dsDNA_perc)/100),
                 ((sum(Nb_EVE_dsDNA,Nb_EVE_ssDNA,Nb_EVE_dsRNA,Nb_EVE_ssRNA)*ssDNA_perc)/100),
                 ((sum(Nb_EVE_dsDNA,Nb_EVE_ssDNA,Nb_EVE_dsRNA,Nb_EVE_ssRNA)*dsRNA_perc)/100),
                 ((sum(Nb_EVE_dsDNA,Nb_EVE_ssDNA,Nb_EVE_dsRNA,Nb_EVE_ssRNA)*ssRNA_perc)/100))
Expected_Events<- c( ((sum(Events_EVEs_dsDNA,Events_EVEs_ssDNA,Events_EVEs_dsRNA,Events_EVEs_ssRNA)*dsDNA_perc)/100),
                     ((sum(Events_EVEs_dsDNA,Events_EVEs_ssDNA,Events_EVEs_dsRNA,Events_EVEs_ssRNA)*ssDNA_perc)/100),
                     ((sum(Events_EVEs_dsDNA,Events_EVEs_ssDNA,Events_EVEs_dsRNA,Events_EVEs_ssRNA)*dsRNA_perc)/100),
                     ((sum(Events_EVEs_dsDNA,Events_EVEs_ssDNA,Events_EVEs_dsRNA,Events_EVEs_ssRNA)*ssRNA_perc)/100))

Summary_table2 <- data.frame(Genomic_structures,Type,value,Expected_NB,Expected_Events)

Bar_plot_NB<-Summary_table2[!grepl("Event",Summary_table2$Type),] %>%
  ggplot( aes(y=value, x=factor(Genomic_structures , levels = c("dsDNA", "ssDNA", "dsRNA", "ssRNA")) ,fill=factor(Type, levels = c("dEVEs","EVEs")),label=Genomic_structures)) +
  geom_bar(stat="identity",position = "identity",color="black",size=0.2)+
  geom_point(size = 1.5,aes(y=Expected_NB, shape="Exposure",fill="black",stroke = 1),shape=3)+
  scale_fill_manual(values = c("EVEs"="#FFD460","dEVEs"="#F07B3F"))+  theme_minimal() +
  theme(
    panel.grid.major.x  = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line())+
  ylab("Number of EVEs") + xlab("Viral genomic structures") +
  theme(axis.text.y = element_text( color="black",size=9))+
  theme(legend.position = "none") +scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line())+
  theme(text = element_text(size=17)) + theme(strip.text.x = element_text(size = 13,face="bold"))+
  theme(axis.text.y = element_text( color="black",size=10))+theme(axis.title.x = element_blank())

Bar_plot_EVENT<-Summary_table2[grepl("Event",Summary_table2$Type),] %>%
  ggplot( aes(y=value, x=factor(Genomic_structures , levels = c("dsDNA", "ssDNA", "dsRNA", "ssRNA")),fill=factor(Type, levels = c("dEVEs_Event","EVEs_Event")),label=Genomic_structures)) +
  geom_bar(stat="identity",position = "identity",color="black",size=0.2)+
  geom_point(size = 1.5,aes(y=Expected_Events, shape="Exposure",fill="black",stroke = 1),shape=3)+
  scale_fill_manual(values = c("EVEs_Event"="#FFD460","dEVEs_Event"="#F07B3F"))+  theme_minimal() +
  theme(
    panel.grid.major.x  = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line())+
  ylab("Number of Events") + xlab("Viral genomic structures") +
  theme(axis.text.y = element_text( color="black",size=9))+
  theme(legend.position = "none") +scale_y_continuous(breaks = scales::pretty_breaks(n = 20))+
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line())+
  theme(text = element_text(size=17)) + theme(strip.text.x = element_text(size = 13,face="bold"))+
  theme(axis.text.y = element_text( color="black",size=10))+theme(axis.title.x = element_blank()) +
  scale_y_continuous(labels = abs)

library(ggpubr)
Bar_plot_EVENT_NB<-ggarrange(Bar_plot_NB,Bar_plot_EVENT,ncol=1)

####
library(forcats)
library(ggplot2)
library(ggtext)
library(ggnewscale)
#Get consensus Families dEVEs Event count
Env_table2<-Env_table2[order(Env_table2$TPM_all, decreasing = TRUE),]


Families_EVEs_event<- Env_table2[!duplicated(Env_table2[ c("ID")]),]

Families_EVEs_event<-Families_EVEs_event[!Families_EVEs_event$consensus_family_ID=="Unclassified",]
#Families_EVEs_event$consensus_family_ID[is.na(Families_EVEs_event$consensus_family_ID)]<-"Unknown"
#Families_EVEs_event$genomic_structure[is.na(Families_EVEs_event$consensus_family_ID)]<- "Unclassified RNA"


Families_dEVEs_event <- Families_EVEs_event[Families_EVEs_event$FDR_pvalue_dNdS == 1 & as.numeric(as.character(Families_EVEs_event$Mean_dNdS))+as.numeric(as.character(Families_EVEs_event$SE_dNdS)) < 1 & Families_EVEs_event$pseudogenized ==0 & Families_EVEs_event$ORF_perc>1| Families_EVEs_event$TPM_all>=1000 & Families_EVEs_event$pseudogenized ==0 & Families_EVEs_event$ORF_perc>1,]

Families_EVEs_event<-as.data.frame(table(Families_EVEs_event$consensus_family_ID))

Families_dEVEs_event <-as.data.frame(table(Families_dEVEs_event$consensus_family_ID))

#Families_dEVEs_event<-Families_dEVEs_event %>% add_row(Var1 = "Unknown", Freq = 0)

EVEs_dEVEs_table_family<-merge(x = Families_EVEs_event, y = Families_dEVEs_event, by = "Var1", all = TRUE)

EVEs_dEVEs_Familiescount<-reshape2::melt(EVEs_dEVEs_table_family , id=c("Var1"))

colnames(EVEs_dEVEs_Familiescount) <- c("consensus_family_ID",'variable','value')

##Add genomic structure for coloring
sub_fam_tab<-Env_table2[!duplicated(Env_table2$consensus_family_ID),]
sub_fam_tab<-sub_fam_tab[c("consensus_family_ID","genomic_structure")]


EVEs_dEVEs_Familiescount<-merge(sub_fam_tab,EVEs_dEVEs_Familiescount,by="consensus_family_ID")
EVEs_dEVEs_Familiescount$genomic_structure[EVEs_dEVEs_Familiescount$consensus_family_ID=="Unclassified_ssDNA"]<-"ssDNA"
EVEs_dEVEs_Familiescount$value[is.na(EVEs_dEVEs_Familiescount$value)] <- 0
EVEs_dEVEs_Familiescount$variable<-gsub("Freq.x","EVEs",EVEs_dEVEs_Familiescount$variable)
EVEs_dEVEs_Familiescount$variable<-gsub("Freq.y","dEVEs",EVEs_dEVEs_Familiescount$variable)

EVEs_dEVEs_Familiescount$genomic_structure[EVEs_dEVEs_Familiescount$best_family_per_query=="Partiti-Picobirna"]<-"dsRNA"
EVEs_dEVEs_Familiescount$genomic_structure[EVEs_dEVEs_Familiescount$best_family_per_query=="Luteo-Sobemo"]<-"ssRNA"

axis_color <- dplyr::distinct(EVEs_dEVEs_Familiescount, consensus_family_ID, genomic_structure) %>%
  tibble::deframe()

color_list<- c("EVEs"="#FFD460","dEVEs"="#F07B3F")

EVEs_dEVEs_Familiescount$genomic_structure <- factor(EVEs_dEVEs_Familiescount$genomic_structure,      # Reordering group factor levels
                                                               levels = c("dsDNA", "ssDNA", "dsRNA", "ssRNA"))


List_inset_viruses<-c("Polyomaviridae","Lispivirida (-) ","Phasmaviridae","Nimaviridae","Marseilleviridae","IVSPERs","Artoviridae (-) ","Circoviridae","Chuviridae (-) ","Xinmoviridae","Nyamiviridae","Marseivviridae","Chuviridae",
                      "Apis_filamentous-like","Astroviridae","Aspiviridae","Benyviridae","Hepeviridae","Alphatetraviridae","Togaviridae","Bastrovirus","Birnavirdae","Bunyavirales","Arenaviridae","Cystoviridae","Endornaviridae","Flaviviridae","Hypoviridae","Leviviridae","Mononegavirales","Jingchuvirales","Narnaviridae","Botourmiaviridae","Nidovirales","Partitiviridae",
                      "Cruciviridae","Amalgaviridae","Picobirnaviridae","Orthomyxoviridae","Permutotetraviridae","Potyviridae","Genomoviridae",
                      "Picornavirales","Solinviviridaes","Caliciviridae","Marnaviridae","Qinviridae","Solemoviridae","Luteoviridae",
                      "Alvernaviridae","Reoviridae","Tombusviridae","Nodaviridae","Sinaivirus","Carmotetraviridae","Luteoviridae","Totiviridae","Chrysoviridae","Megabirnaviridae","Quadriviridae","Botybirnavirus","Tymovirales","Virgaviridae","Togaviridae","Bromoviridae","Closteroviridae","Idaeovirus","Weivirus","LbFV_like","Nudiviridae","Bidnaviridae","Parvoviridae",
                      "Hytrosaviridae","Ascoviridae","Retroviridae","Nyamiviridae (-) ","Partitiviridae","Birnaviridae","Plasmaviridae","Rhabdoviridae","Dicistroviridae","Iflaviridae","Permutotetraviridae","Flaviviridae","Negevirus","Nodaviridae","Noravirus","Baculoviridae","Iridoviridae","AmFV-like","Poxviridae","Mesoviridae","Bunyaviridae","Mesoniviridae","Sarthroviridae",
                      "Solinviviridae","Peribunyaviridae","Phenuiviridae","Tymoviridae (+) ","Orthomyxoviridae (-) ","Peribunyaviridae (-) ","Phasmaviridae (-) ","Phenuiviridae (-) ","Qinviridae (-) ","Rhabdoviridae (-) ","Tombusviridae (+) ","Xinmoviridae (-) ","Luteo-Sobemo (+) ","Tombus-Noda (+) ","Abisko","Partiti-Picobirna","Mono-Chu (-) ","Narna-Levi (+) ","Bunya-Arena (-) ","Negevirus-like (+) ","Hepe-Virga (+) ","Tymoviridae (-) ","Abisko-like (+) ")


EVEs_dEVEs_Familiescount$consensus_family_ID[EVEs_dEVEs_Familiescount$consensus_family_ID=="Tymoviridae"] <- "Tymoviridae (+) "
EVEs_dEVEs_Familiescount$consensus_family_ID[EVEs_dEVEs_Familiescount$consensus_family_ID=="Mono-Chu"] <- "Mono-Chu (-) "
EVEs_dEVEs_Familiescount$consensus_family_ID[EVEs_dEVEs_Familiescount$consensus_family_ID=="Artoviridae"] <- "Artoviridae (-) "
EVEs_dEVEs_Familiescount$consensus_family_ID[EVEs_dEVEs_Familiescount$consensus_family_ID=="Narna-Levi"] <- "Narna-Levi (+) "
EVEs_dEVEs_Familiescount$consensus_family_ID[EVEs_dEVEs_Familiescount$consensus_family_ID=="Picorna-Calici"] <- "Picorna-Calici (+) "
EVEs_dEVEs_Familiescount$consensus_family_ID[EVEs_dEVEs_Familiescount$consensus_family_ID=="Hepe-Virga"] <- "Hepe-Virga (+) "
EVEs_dEVEs_Familiescount$consensus_family_ID[EVEs_dEVEs_Familiescount$consensus_family_ID=="Tombus-Noda"] <- "Tombus-Noda (+) "
EVEs_dEVEs_Familiescount$consensus_family_ID[EVEs_dEVEs_Familiescount$consensus_family_ID=="Lispiviridae"] <- "Lispivirida (-) "
EVEs_dEVEs_Familiescount$consensus_family_ID[EVEs_dEVEs_Familiescount$consensus_family_ID=="Mypoviridae"] <- "Mypoviridae (-) "
EVEs_dEVEs_Familiescount$consensus_family_ID[EVEs_dEVEs_Familiescount$consensus_family_ID=="Nyamiviridae"] <- "Nyamiviridae (-) "
EVEs_dEVEs_Familiescount$consensus_family_ID[EVEs_dEVEs_Familiescount$consensus_family_ID=="Orthomyxoviridae"] <- "Orthomyxoviridae (-) "
EVEs_dEVEs_Familiescount$consensus_family_ID[EVEs_dEVEs_Familiescount$consensus_family_ID=="Bunya-Arena"] <- "Bunya-Arena (-) "
EVEs_dEVEs_Familiescount$consensus_family_ID[EVEs_dEVEs_Familiescount$consensus_family_ID=="Phasmaviridae"] <- "Phasmaviridae (-) "
EVEs_dEVEs_Familiescount$consensus_family_ID[EVEs_dEVEs_Familiescount$consensus_family_ID=="Phenuiviridae"] <- "Phenuiviridae (-) "
EVEs_dEVEs_Familiescount$consensus_family_ID[EVEs_dEVEs_Familiescount$consensus_family_ID=="Qinviridae"] <- "Qinviridae (-) "
EVEs_dEVEs_Familiescount$consensus_family_ID[EVEs_dEVEs_Familiescount$consensus_family_ID=="Rhabdoviridae"] <- "Rhabdoviridae (-) "
EVEs_dEVEs_Familiescount$consensus_family_ID[EVEs_dEVEs_Familiescount$consensus_family_ID=="Tombusviridae"] <- "Tombusviridae (+) "
EVEs_dEVEs_Familiescount$consensus_family_ID[EVEs_dEVEs_Familiescount$consensus_family_ID=="Xinmoviridae"] <- "Xinmoviridae (-) "
EVEs_dEVEs_Familiescount$consensus_family_ID[EVEs_dEVEs_Familiescount$consensus_family_ID=="Tymoviridae"] <- "Tymoviridae (+) "
EVEs_dEVEs_Familiescount$consensus_family_ID[EVEs_dEVEs_Familiescount$consensus_family_ID=="Bornaviridae"] <- "Bornaviridae (-) "
EVEs_dEVEs_Familiescount$consensus_family_ID[EVEs_dEVEs_Familiescount$consensus_family_ID=="Chuviridae"] <- "Chuviridae (-) "
EVEs_dEVEs_Familiescount$consensus_family_ID[EVEs_dEVEs_Familiescount$consensus_family_ID=="Luteo-Sobemo"] <- "Luteo-Sobemo (+) "
EVEs_dEVEs_Familiescount$consensus_family_ID[EVEs_dEVEs_Familiescount$consensus_family_ID=="Abisko-like"] <- "Abisko-like (+) "
EVEs_dEVEs_Familiescount$consensus_family_ID[EVEs_dEVEs_Familiescount$consensus_family_ID=="Peribunyaviridae"] <- "Peribunyaviridae (-) "
EVEs_dEVEs_Familiescount$consensus_family_ID[EVEs_dEVEs_Familiescount$consensus_family_ID=="Negevirus-like"] <- "Negevirus-like (+) "
EVEs_dEVEs_Familiescount$Insect <- ""
for (i in 1:nrow(EVEs_dEVEs_Familiescount)){
  if (EVEs_dEVEs_Familiescount$consensus_family_ID[i] %in% List_inset_viruses){
    print(EVEs_dEVEs_Familiescount$consensus_family_ID[i])
    EVEs_dEVEs_Familiescount$Insect[i] <- "*"
  }
}


#insect :
#  Chuviridae (Shi, M. et al 2016)
#  Iridoviridae * https://www.scielo.br/j/ne/a/ZJTrQfNCxL5knVf6L7sNzGL/?lang=en
#  Flaviviridae *  (Bradley J. Blitvich1,* and Andrew E. Firth2, 2015)
#  Narnaviridae (Haoming Wu et al 2020)
#  Nudiviridae ICTV
#  Nyamiviridae (Fei Wang  et al 2017 & ICTV )
#  Partitiviridae ( Shi, M. et al 2016 Members of Partitiviridae are known viruses of plants, fungi and protists and the first record from an insect (Drosophila) was only in 2015 [5], with subsequent reports in other invertebrates also using NGS approaches [1,5,37].)
#  Parvoviridae (Densovirinae : ICTV)
#  Orthomyxoviridae :Quaranjaviruses  & Isavirus
#  Peribunyaviridae - ICTV
#  Phenuiviridae - Pontus Öhlund et al 2019
#  Poxviridae (Entomopoxvirinae) (Graziele Pereira Oliveira et al 2017)
#  Qinviridae (Haoming Wu et al 2020)
#  Reoviridae : Viralzone & Haoming Wu et al 2020)
#  Retroviridae (Bryony C. Bonning*)
#  Rhabdoviridae (Bryony C. Bonning; ICTV)
#  Tombusviridae(Andrew J.Bennett et al 2019 ) ?
#  Totiviridae (Kenta Okamoto, et al 2016)
#  Xinmoviridae (Simon Käfer et al 2019)
#  Marseilleviridae (Dehia Sahmi-Bounsiar*)
#  Circoviride (Tristan P. W. Dennis et al 2018)
#  Artoviridae  *  - ICTV
#  Phasmaviridae - Haoming Wu et al 2020)
#  Lispiviridae -  Haoming Wu et al 2020)
#  Cruciviridae  - Rosario,K. et al 2018  - on fire ants
#  Mypoviridae - Thekke-Veetil,T et al 2020)
#  Asfaviridae  - ICTV - on ticks
#  Polyomaviridae  - Rosario,K. et al 2018 (T antigen)

EVEs_dEVEs_Familiescount$consensus_family_ID<- paste0(EVEs_dEVEs_Familiescount$consensus_family_ID,EVEs_dEVEs_Familiescount$Insect)


EVEs_dEVEs_Familiescount$value2 <- EVEs_dEVEs_Familiescount$value
EVEs_dEVEs_Familiescount$value2[EVEs_dEVEs_Familiescount$variable=="dEVEs"]<-0
EVEs_dEVEs_Familiescount$consensus_family_ID <- reorder(EVEs_dEVEs_Familiescount$consensus_family_ID, EVEs_dEVEs_Familiescount$value2)


# Number of Event inferred from virus not known to infect insects

length(unique(Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),][Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),]$consensus_family_ID%in% c("Caulimoviridae" ,"Mypoviridae","Genomoviridae","Herpesviridae"  ,"Papillomaviridae","Phycodnaviridae","Bornaviridae") ,]$ID))

#Get percentage EVE identity of non-insect virus families
mean(Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),][Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),]$consensus_family_ID%in% c("Caulimoviridae" ,"Mypoviridae","Genomoviridae","Herpesviridae"  ,"Papillomaviridae","Phycodnaviridae","Bornaviridae") ,]$pident,na.rm=T)
min(Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),][Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),]$consensus_family_ID %in% c("Caulimoviridae" ,"Mypoviridae","Genomoviridae","Herpesviridae"  ,"Papillomaviridae","Phycodnaviridae","Bornaviridae") ,]$pident,na.rm=T)
max(Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),][Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),]$consensus_family_ID %in% c("Caulimoviridae" ,"Mypoviridae","Genomoviridae","Herpesviridae"  ,"Papillomaviridae","Phycodnaviridae","Bornaviridae") ,]$pident,na.rm=T)

#View(Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),][Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),]$best_family_per_query %in% c("Caulimoviridae","Nimaviridae" ,"Tymoviridae","Genomoviridae","Herpesviridae"  ,"Papillomaviridae","Phycodnaviridae") ,])

median(Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),][Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),]$consensus_family_ID %in% c("Caulimoviridae","Mypoviridae","Genomoviridae","Herpesviridae"  ,"Papillomaviridae","Phycodnaviridae","Bornaviridae"),]$evalue)
min(Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),][Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),]$consensus_family_ID %in% c("Caulimoviridae" ,"Mypoviridae","Genomoviridae","Herpesviridae"  ,"Papillomaviridae","Phycodnaviridae","Bornaviridae"),]$evalue)
max(Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),][Env_table2[!duplicated(Env_table2[c('Clustername','Event')]),]$consensus_family_ID %in% c("Caulimoviridae" ,"Mypoviridae","Genomoviridae","Herpesviridae"  ,"Papillomaviridae","Phycodnaviridae","Bornaviridae"),]$evalue)


#Test different between families

EVEs_dEVEs_Familiescount$genomic_structure[grepl("Tymo",EVEs_dEVEs_Familiescount$consensus_family_ID)] <- "ssRNA"

Bar_plot_families<-ggplot(EVEs_dEVEs_Familiescount,aes(x = factor(consensus_family_ID),y = value,fill=variable))+
  facet_grid(genomic_structure ~ . , scales="free_y",space="free") +
  geom_bar(stat="identity",position = "identity",color="black",size=0.2)+
  geom_bar(data = EVEs_dEVEs_Familiescount[EVEs_dEVEs_Familiescount$variable=="dEVEs",], aes(y = value, x = consensus_family_ID),
           stat = "identity", colour = "black",size=0.2)+
  scale_fill_manual(values = color_list)+
  theme_bw() +
  theme(axis.text.x=element_text(angle=0,hjust=0.5))+  coord_flip() +theme(
    panel.grid.major.y  = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line())+
  theme(strip.text.y = element_text(size = 7,angle=0,face="bold"))+
  theme(axis.text.y = element_blank())+ ylab("Number of Events") + xlab("Viral families") + theme(axis.text.y = element_text( color="black",size=12))+
  theme(legend.position="none")+
  scale_y_continuous(breaks=seq(0,75,5)) +
  theme(axis.text.y = element_text(hjust=0.95))+
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line())+
  theme(text = element_text(size=17)) + theme(strip.text.x = element_text(size = 13,face="bold"))+
  theme(axis.text.y = element_text(face = "italic"))


# Count number events in ssRNA (-) and ssRNA (+)

sum(EVEs_dEVEs_Familiescount$value[grepl("\\(-\\)",EVEs_dEVEs_Familiescount$consensus_family_ID) & EVEs_dEVEs_Familiescount$variable=="EVEs"])
sum(EVEs_dEVEs_Familiescount$value[grepl("\\(\\+\\)",EVEs_dEVEs_Familiescount$consensus_family_ID) & EVEs_dEVEs_Familiescount$variable=="EVEs"])
sum(EVEs_dEVEs_Familiescount$value[EVEs_dEVEs_Familiescount$genomic_structure=="ssRNA" & EVEs_dEVEs_Familiescount$variable=="EVEs"])

#dcast(EVEs_dEVEs_Familiescount, consensus_family , value_var = 'test_result' )
#melt(EVEs_dEVEs_Familiescount, id.vars="consensus_family")



#ggarrange(Bar_plot_EVENT_NB,Bar_plot_families,ncol=2,widths =c(0.5,1))



# Combine all plots
library(cowplot)
ggarrange(
  ggarrange(Bar_plot_EVENT_NB,Bar_plot_families, ncol = 2, widths =c(1,2)),
  Nb_EVES_distribution_plot,nrow = 2, heights  =c(2,1))+
  draw_plot_label(label = c("A", "B","C", "D"), size = 15,
                  x = c(0, 0,0.35, 0), y = c(1, 0.69,1, 0.37))


ggsave("/Users/bguinet/Desktop/Papier_scientifique/Genomic_structure_and_families_Events_distributions2.pdf",scale = 1.4)


#############################################
#Get  distribution of EVEs and dEVEs per species
#############################################

detachAllPackages()
library(ggplot2)
library(tidyverse)
library(dplyr)
Binary_table_alone_EVEs_count<-select(EVEs_alone_event_df_EVEs_count,"species","Clusters","Nb_EVEs")
colnames(Binary_table_alone_EVEs_count)<-c("Events_species","Clustername","Nb_EVEs")

Binary_table_shared_EVEs_count <-  select(EVEs_shared_event_df_EVEs_count,"Events_species2","Clustername","Nb_EVEs")
colnames(Binary_table_shared_EVEs_count)<- c("Events_species","Clustername","Nb_EVEs")

EVEs_matrix_count <-bind_rows(select(Binary_table_shared_EVEs_count,"Events_species","Clustername","Nb_EVEs"),select(Binary_table_alone_EVEs_count,"Events_species","Clustername","Nb_EVEs"))

EVEs_matrix_count<-data.frame(EVEs_matrix_count[rep(seq_len(dim(EVEs_matrix_count)[1]), EVEs_matrix_count$Nb_EVEs), c(1,2,3), drop = FALSE], row.names=NULL)
EVEs_matrix_count$EVE_numbers <- paste0("EVE_",rownames(EVEs_matrix_count))


Species_name_df<-select(anoleData,'species')
colnames(Species_name_df)<-'Events_species'
Binary_EVEs_matrix_count <-select(as.data.frame(EVEs_matrix_count),"Events_species","EVE_numbers") %>% separate_rows(Events_species,sep=",") %>%
  mutate(Val = 1) %>% type.convert(as.is = TRUE)

Binary_EVEs_matrix_count<- as.data.frame(table(Binary_EVEs_matrix_count$Events_species))

Binary_EVEs_matrix_count<-merge(anoleData,Binary_EVEs_matrix_count,by.y="Var1",by.x="species")

library(forcats)
Binary_EVEs_plot_count<-Binary_EVEs_matrix_count %>%
  mutate(species = fct_reorder(species, -desc(Freq))) %>%
  ggplot( aes(x=species, y=Freq, fill=lifecycle1)) +
  geom_bar(stat="identity") +
  coord_flip()+
  theme_bw()+ theme(axis.text.y = element_text(color="black",
                                               size=5))+
  labs(x ="Species names", y = "Nb endogenous viral elements (EVEs)")


#dEVEs matrix count
Binary_table_alone_dEVEs_count<-select(EVEs_alone_event_df_EVEs_count[EVEs_alone_event_df_EVEs_count$Nb_dEVEs>=1,],"species","Clusters","Nb_dEVEs")
colnames(Binary_table_alone_dEVEs_count)<-c("Events_species","Clustername","Nb_dEVEs")

Binary_table_shared_dEVEs_count <-  select(EVEs_shared_event_df_EVEs_count[EVEs_shared_event_df_EVEs_count$Nb_dEVEs>=1,],"Events_species2","Clustername","Nb_dEVEs")
colnames(Binary_table_shared_dEVEs_count)<- c("Events_species","Clustername","Nb_dEVEs")

dEVEs_matrix_count <-bind_rows(select(Binary_table_shared_dEVEs_count,"Events_species","Clustername","Nb_dEVEs"),select(Binary_table_alone_dEVEs_count,"Events_species","Clustername","Nb_dEVEs"))
dEVEs_matrix_count<-data.frame(dEVEs_matrix_count[rep(seq_len(dim(dEVEs_matrix_count)[1]), dEVEs_matrix_count$Nb_dEVEs), c(1,2,3), drop = FALSE], row.names=NULL)
dEVEs_matrix_count$dEVE_numbers <- paste0("dEVE_",rownames(dEVEs_matrix_count))

Species_name_df<-select(anoleData,'species')
colnames(Species_name_df)<-'Events_species'

Binary_dEVEs_matrix_count <-select(as.data.frame(dEVEs_matrix_count),"Events_species","dEVE_numbers") %>% separate_rows(Events_species,sep=",") %>%
  mutate(Val = 1) %>% type.convert(as.is = TRUE)


Binary_dEVEs_matrix_count<- table(Binary_dEVEs_matrix_count$Events_species)
Binary_dEVEs_matrix_count<-merge(anoleData,Binary_dEVEs_matrix_count,by.y="Var1",by.x="species")


Binary_dEVEs_plot_count<-Binary_dEVEs_matrix_count %>%
  mutate(species = fct_reorder(species, -desc(Freq))) %>%
  ggplot( aes(x=species, y=Freq, fill=lifecycle1)) +
  geom_bar(stat="identity") +
  coord_flip()+
  theme_bw()+ theme(axis.text.y = element_text(color="black",
                                               size=7))+
  labs(x ="Species names", y = "Nb domesticated EVEs (dEVEs)")

library(ggpubr)
ggarrange(Binary_EVEs_plot_count,Binary_dEVEs_plot_count,nrow=1, common.legend = TRUE, legend="bottom")

ggsave("/Users/bguinet/Desktop/Papier_scientifique/EVEs_dEVEs_count_species_distributions.pdf",scale = 1.15)



#####Distribution of function around multiple EVEs acquired for each lifestyle


saveEVEs_shared_event_df<-select(saveEVEs_shared_event_df,'species','Clustername')
saveEVEs_shared_event_df<-merge(x = select(anoleData,"species","lifecycle1"), y = saveEVEs_shared_event_df, by = "species", all = TRUE)
saveEVEs_shared_event_df<-saveEVEs_shared_event_df[!is.na(saveEVEs_shared_event_df$Clustername),]
colnames(saveEVEs_shared_event_df)<-c("species","lifecycle1","Clusternames")

savetable_alone_EVEs_event <- select(savetable_alone_EVEs_event,"species",'Clusters')
savetable_alone_EVEs_event <-merge(x = select(anoleData,"species","lifecycle1"), y = savetable_alone_EVEs_event , by = "species", all = TRUE)
saveEVEs_alone_event_df<-savetable_alone_EVEs_event [!is.na(savetable_alone_EVEs_event$Clusters),]
colnames(saveEVEs_alone_event_df)<-c("species","lifecycle1","Clusternames")

saveEVEs_event_df<-rbind(saveEVEs_shared_event_df,saveEVEs_alone_event_df)

list_cluster_name_freeliving<-c()
list_cluster_name_endoparasitoide<-c()
for(i in 1:nrow(saveEVEs_event_df)) {
  row <- saveEVEs_event_df[i,]
  if (row$lifecycle1 == "freeliving"){
    list_cluster_name_freeliving<-c(list_cluster_name_freeliving,strsplit(row$Clusternames,","))
    # do stuff with row
  }
  if(row$lifecycle1=="endoparasitoide"){
    list_cluster_name_endoparasitoide<-c(list_cluster_name_endoparasitoide,strsplit(row$Clusternames,","))
  }
}

paste0("Mean number of EVEs implicated in Events within endoparasitoide sp:",sum(sapply(list_cluster_name_endoparasitoide,length))/length(list_cluster_name_endoparasitoide))
paste0("Mean number of EVEs implicated in Events within freeliving sp:",sum(sapply(list_cluster_name_freeliving,length))/length(list_cluster_name_freeliving))

for (i in list_cluster_name_endoparasitoide){
  print(i)
}
sub_tab_endoparasitoide<- dsDNA_tab%>%
  filter(Clustername %in% unlist(list_cluster_name_endoparasitoide, recursive = FALSE) )

#We remove NA to not include them into the analysis
sub_tab_endoparasitoide<-sub_tab_endoparasitoide[!is.na(sub_tab_endoparasitoide$Domain_description),]
#Sorte by evalue
sub_tab_endoparasitoide<-sub_tab_endoparasitoide[order(sub_tab_endoparasitoide$evalue, decreasing = FALSE),]
sub_tab_endoparasitoide<-sub_tab_endoparasitoide[ !duplicated(sub_tab_endoparasitoide$Clustername), ]


for (i in list_cluster_name_freeliving){
  print(i)
}
sub_tab_freeliving<- dsDNA_tab%>%
  filter(Clustername %in% unlist(list_cluster_name_freeliving, recursive = FALSE) )

#We remove NA to not include them into the analysis
sub_tab_freeliving<-sub_tab_freeliving[!is.na(sub_tab_freeliving$Domain_description),]
#Sorte by evalue
sub_tab_freeliving<-sub_tab_freeliving[order(sub_tab_freeliving$evalue, decreasing = FALSE),]
sub_tab_freeliving<-sub_tab_freeliving[ !duplicated(sub_tab_freeliving$Clustername), ]


savedEVEs_shared_event_df
savetable_alone_dEVEs_event
savetable_alone_EVEs_event



library(Hmisc)
library(knitr)
latex(Summary_table, file="")

kable(Summary_table, "latex")


####################################################################
#Add it into the phylogeny

Event_shared_df<-Event_shared_df%>% mutate (Nb_EVEs_Events = ifelse(Nb_EVEs >= 1 ,1,0 ))
Event_shared_df$Nb_dEVEs_Events<- Event_shared_df$Nb_dEVEs_dsDNA_Events + Event_shared_df$Nb_dEVEs_ssDNA_Events +  Event_shared_df$Nb_dEVEs_ssRNA_Events +  Event_shared_df$Nb_dEVEs_dsRNA_Events +  Event_shared_df$Nb_dEVEs_Unclassified_Events

Event_shared_tree_df<- select(Event_shared_df,"Node_number","Nb_EVEs","Nb_dEVEs",'Nb_EVEs_Events',"Nb_dEVEs_Events") %>%
  group_by(Node_number) %>%
  summarise(Nb_EVEs_shared = sum(Nb_EVEs),Nb_dEVEs_shared = sum(Nb_dEVEs),Nb_EVEs_shared_Events = sum(Nb_EVEs_Events), Nb_dEVEs_shared_Events = sum(Nb_dEVEs_Events) )


# Create piechart table

library(dplyr)
library(tidyverse)
Pie_chart_tab<-Event_shared_df



Pie_chart_tab$consensus_family[Pie_chart_tab$consensus_family=="LbFV_like,Nudiviridae"]<-"Nudiviridae"


Pie_chart_tab<-Pie_chart_tab%>% select (Node_number,consensus_family,Nb_EVEs,Nb_dEVEs) %>%
  group_by(Node_number,consensus_family) %>%
  mutate(Frequency = sum(Nb_EVEs))

detachAllPackages()
library(dplyr)
library(ggtree)
library(ggplot2)
library(tidyverse)
library(dplyr)

Pie_chart_tab_save<-Pie_chart_tab [c("Node_number","consensus_family","Nb_EVEs")]%>%
  group_by(Node_number) %>%
  summarise(across(Nb_EVEs, sum))

library(tidyverse)


library(tidyr)

nodepie2 <- function(data, cols, color, alpha=1) {
  if (! "Node_number" %in% colnames(data)) {
    stop("data should have a column 'node'...")
  }
  type <- value <- NULL
  if (missingArg(color)) {
    color <- NA
  }
  ldf<-split(data , f = data$Node_number )
  #print(ldf)
  lapply(ldf, function(df) ggpie2(df, y=~value,fill=~type, color, alpha))
}

ggpie2 <- function(data, y, fill, color, alpha=1) {

  ldf1 <- data %>%
    mutate(nonNb_dEVEs = Nb_EVEs - Nb_dEVEs) %>%
    pivot_longer(c(Nb_dEVEs, nonNb_dEVEs), names_to = "Groups1", values_to = "Nb_dEVEs") %>%
    mutate(Node_number = interaction(Node_number, Groups1, lex.order = TRUE))

  Nb_fam<- length(unique(data$consensus_family))

  pal <- scales::hue_pal()(Nb_fam)
  names(pal) <- unique(data$consensus_family)
  pal <- c(pal, c(Nb_dEVEs = "black", nonNb_dEVEs = "transparent"))

  pi<-ggplot() +
    geom_col(data =data , aes(x = 1, y = Nb_EVEs, fill = factor(consensus_family)), color = "white",stroke=.1) +
    geom_col(data = ldf1, aes(x = 1.5, y = Nb_dEVEs, fill = Groups1, group = Node_number), width = .2) +
    scale_fill_manual(values=color)+
    coord_polar(theta = "y") +
    theme_void() +
    theme(legend.position="none")
  return(pi)
}


Color_palette<-c("Phycodnaviridae"="#6c919e","Nudiviridae"="#004D80","Iridoviridae"="#67BCEB",'LbFV_like'="#0276BB","Poxviridae"="#00FFFF","Apis_filamentous-like"="#485F92","IVSPERs"="#B5DCE5","Baculoviridae"="#4773B8","Herpesviridae"="#3BABA8","Nimaviridae"="#3BABA8","Papillomaviridae"="#3BABA8","Ascoviridae"="#3BABA8","Asfaviridae"="#3BABA8","Caulimoviridae"="#3BABA8","Tymoviridae"="#3BABA8","Phycodnaviridae"="#3BABA8",
                 "Parvoviridae"="#E52421","Circoviridae"="#F4A09D","Genomoviridae"="#971D14","Partiti-Picobirna"="#734b00",
                 "Artoviridae"="#5BB888","Chuviridae"="#E8E100","Nyamiviridae"="#77FF46","Orthomyxoviridae"="#E9EEBF","Xinmoviridae"="#4F9436","Lispiviridae"="#BFDFCD","Phenuiviridae"="#2C5933","Bornaviridae"="#C0CE4A","Virgaviridae"="#ADF2CE","Bornaviridae"="#ADF2CE","Tombusviridae"="#ADF2CE","Qinviridae"="#ADF2CE","Peribunyaviridae"="#ADF2CE","Mypoviridae"="#ADF2CE","Phasmaviridae"="#ADF2CE","Narna-Levi"="#ADF2CE","Mono-Chu"="#ADF2CE","Hepe-Virga"="#ADF2CE","Abisko-like"="#ADF2CE","Tombusviridae"="#ADF2CE","Luteo-Sobemo"="#ADF2CE","Negevirus-like"="#ADF2CE",
                 "Totiviridae"="#FC990B","Reoviridae"="#C98A4B","Partitiviridae"="#F4D0AB",
                 "Unknown"="black","Unclassified ssDNA"="#d6d6d6","Unclassified ssRNA"="#9e9e9e","Unclassified dsRNA"="#525251","Nb_dEVEs"="black","nonNb_dEVEs"="transparent" )

pies<-nodepie2(Pie_chart_tab,color=Color_palette)


library(ggplot2)
library(dplyr)
library(tidyr)


Event_alone_df$species[Event_alone_df$species=="Unknown_species"]<-"Unknown_sp"


Event_alone_df<-Event_alone_df%>% mutate (Nb_EVEs_Events = ifelse(Nb_EVEs >= 1 ,1,0 ))
Event_alone_df$Nb_dEVEs_Events<- Event_alone_df$Nb_dEVEs_dsDNA_Events + Event_alone_df$Nb_dEVEs_ssDNA_Events +  Event_alone_df$Nb_dEVEs_ssRNA_Events +  Event_alone_df$Nb_dEVEs_dsRNA_Events +  Event_alone_df$Nb_dEVEs_Unclassified_Events

Event_alone_tree_df<- select(Event_alone_df,"species","Nb_EVEs","Nb_dEVEs","Nb_EVEs_Events","Nb_dEVEs_Events") %>%
  group_by(species) %>%
  summarise(Nb_EVEs_alone = sum(Nb_EVEs),Nb_dEVEs_alone = sum(Nb_dEVEs),Nb_EVEs_alone_Events = sum(Nb_EVEs_Events),Nb_dEVEs_alone_Events = sum(Nb_dEVEs_Events) )




#Event shared
detachAllPackages()
library(ape)
library(treeio)
library(dplyr)
library(ape)
library(ggtree)
library(tibble)
library(RevGadgets)

#Load datation MAP tree and ad the lengths
#tree = read.newick("/Users/bguinet/Desktop/TimeTree_mcmc_3.nwk")
#treeNW = as.tibble(read.newick("/Users/bguinet/Desktop/Papier_scientifique/Revbayes_anestral_states_estimates/data/TimeTree_mcmc_10_test_with_outgroup.nwk"))
treeNW = as.tibble(read.newick("/Users/bguinet/Desktop/Papier_scientifique/alignment.fasta.approx.dating_20000_moremoves_noMC3_Exp_CompoundMove_BIS_map.nwk"))
treeNW$label<-gsub("Unknwon_species","Unknown_sp",treeNW$label)


library(ape)
library(treeio)
#treeNX1 <- readTrees(paths = "/Users/bguinet/Desktop/Papier_scientifique/Revbayes_anestral_states_estimates/data/TimeTree_mcmc_10_test_with_outgroup2.tree")
#treeNX1 <- readTrees(paths = "/Users/bguinet/Desktop/Papier_scientifique/alignment.fasta.approx.dating_20000_moremoves_noMC3_Exp_CompoundMove_BIS_map.tree")
treeNX1 <- readTrees(paths = "/Users/bguinet/Desktop/Papier_scientifique/AnalysesTER/alignment.fasta.approx.dating_20000_moremoves_noMC3_Exp_CompoundMove_TER_map.tree")

#treeNX1 <- readTrees(paths = "/Users/bguinet/Desktop/Papier_scientifique/alignment.fasta.approx.dating_20000.tree")

treeNX<-plotFBDTree(tree = treeNX1,
                    timeline = T,
                    time_bars = T,
                    geo_units = "epochs",
                    tip_labels_italics = T,
                    tip_labels_remove_underscore = T,
                    tip_labels_size = 2,
                    tip_age_bars = T,
                    node_age_bars = T,
                    age_bars_colored_by = "posterior",
                    label_sampled_ancs = TRUE) +
  theme(legend.position='none')+
  geom_text(aes(label=index), hjust=-.3)



tree<-merge(subset(treeNX$data, select=-c(label,branch.length,parent)), treeNW,by="node")
#tree<-drop.tip(as.phylo(tree), c("Tribolium_castaneum", "Anoplophora_glabripennis"), root.edge = 0)

#tree<-as.tibble(tree)
#colnames(tree)<- c("parent","node","branch.lengthMA","label" )


tree_bis<-tree


#tree_bis<-drop.tip(as.phylo(tree_bis), "Unknown_sp")

tree2<-merge(as_tibble(tree_bis),Event_shared_tree_df, by.x="node", by.y="Node_number",all = TRUE)
#tree2$node_depth<-node.depth(tree, method = 2)
#tree2$node_depth<-node.depth.edgelength(tree)

Event_alone_tree_df$species <- gsub(" ","_",Event_alone_tree_df$species)
tree3<-merge(x=tree2,y=Event_alone_tree_df, by.x="label", by.y="species",all = TRUE)

tree3<-as.data.frame(tree3)
#Replace zero by NA

tree3$Nb_dEVEs_alone<- as.numeric(tree3$Nb_dEVEs_alone)
tree3$Nb_dEVEs_alone[tree3$Nb_dEVEs_alone==0]<-NA

tree3$Nb_EVEs_alone<- as.numeric(tree3$Nb_EVEs_alone)
tree3$Nb_EVEs_alone[tree3$Nb_EVEs_alone==0]<-NA


tree3$Nb_dEVEs_shared<- as.numeric(tree3$Nb_dEVEs_shared)
tree3$Nb_dEVEs_shared[tree3$Nb_dEVEs_shared==0]<-NA

tree3$Nb_EVEs_shared<- as.numeric(tree3$Nb_EVEs_shared)
tree3$Nb_EVEs_shared[tree3$Nb_EVEs_shared==0]<-NA

tree3$Nb_EVEs_shared_Events<- as.numeric(tree3$Nb_EVEs_shared_Events)
tree3$Nb_EVEs_shared_Events[tree3$Nb_EVEs_shared_Events==0]<-NA

tree3$Nb_EVEs_alone_Events<- as.numeric(tree3$Nb_EVEs_alone_Events)
tree3$Nb_EVEs_alone_Events[tree3$Nb_EVEs_alone_Events==0]<-NA


tree4<-left_join(fortify(tree_bis, ladderize = FALSE),fortify(select(tree3,-c("label","parent"))) , by = c("node" = "node") )
names(tree4)[names(tree4)=="parent.x"] <- "parent"

tree6 <- tree4
#Transforme data
tree6$Nb_EVEs_shared[is.na(tree6$Nb_EVEs_shared)]<-0
tree6$Nb_EVEs_alone[is.na(tree6$Nb_EVEs_alone)]<-0
tree6$Nb_dEVEs_shared[is.na(tree6$Nb_dEVEs_shared)]<-0
tree6$Nb_dEVEs_alone[is.na(tree6$Nb_dEVEs_alone)]<-0

tree6$Nb_EVEs_shared_Events[is.na(tree6$Nb_EVEs_shared_Events)]<-0
tree6$Nb_EVEs_alone_Events[is.na(tree6$Nb_EVEs_alone_Events)]<-0
tree6$Nb_dEVEs_shared_Events[is.na(tree6$Nb_dEVEs_shared_Events)]<-0
tree6$Nb_dEVEs_alone_Events[is.na(tree6$Nb_dEVEs_alone_Events)]<-0

tree6$Nb_EVEs_Events<- tree6$Nb_EVEs_shared_Events + tree6$Nb_EVEs_alone_Events
tree6$Nb_dEVEs_Events<- tree6$Nb_dEVEs_shared_Events + tree6$Nb_dEVEs_alone_Events
tree6$Nb_EVEs2<- tree6$Nb_EVEs_shared + tree6$Nb_EVEs_alone
tree6$Nb_dEVEs2<- tree6$Nb_dEVEs_shared + tree6$Nb_dEVEs_alone

tree6<-tree6[c("node","parent", "branch.length.x","label","index.x","Nb_EVEs_Events","Nb_dEVEs_Events","Nb_EVEs2","Nb_dEVEs2","x.y")]



#write.table(as.data.frame(tree6),"/Users/bguinet/Desktop/Papier_scientifique/All_table_for_GLM2/Matrix_EVEs_dEVEs_dsDNA_A",sep=";")


#ssRNA_A <- NbEVEs2 = 91
#ssRNA_A <- NbdEVEs2 = 30
#ssRNA_A= tree6$Nb_EVEs_Events = 89

#ssDNA_A <- NbEVEs2 = 26
#ssDNA_A <- NbdEVEs2 = 39
#ssDNA_A= tree6$Nb_EVEs_Events = 18

#tree4<-tree4[c("node","parent.x.x", "branch.length.x.x","index.x","isTip.x","Nb_EVEs_Events","Nb_dEVEs_Events","Nb_EVEs2","Nb_dEVEs2")]

#colnames(tree4)
#Skip to Revbayes analysis
#
library(ggtree)
library(treeio)
library(ggtree)
library(ggstance)
library(ggplot2)
library(ggh4x)
library(deeptime)

#tree5<- tree4[c("node","parent","branch.length.x","label",'branch.x','x.x','Nb_EVEs_shared','Nb_dEVEs_shared','Nb_EVEs_alone_Events','Nb_dEVEs_alone_Events','Nb_EVEs_shared_Events','Nb_dEVEs_shared_Events','Nb_dEVEs_alone','Nb_EVEs_alone')]
#colnames(tree5)<-c("node","parent","branch.length",'label','branch','x','Nb_EVEs_shared','Nb_dEVEs_shared','Nb_EVEs_alone_Events','Nb_dEVEs_alone_Events','Nb_EVEs_shared_Events','Nb_dEVEs_shared_Events','Nb_dEVEs_alone','Nb_EVEs_alone')
tree5<- tree4[c("node","parent","branch.length.x","label",'Nb_EVEs_shared','Nb_dEVEs_shared','Nb_EVEs_alone_Events','Nb_dEVEs_alone_Events','Nb_EVEs_shared_Events','Nb_dEVEs_shared_Events','Nb_dEVEs_alone','Nb_EVEs_alone')]
colnames(tree5)<-c("node","parent","branch.length",'label','Nb_EVEs_shared','Nb_dEVEs_shared','Nb_EVEs_alone_Events','Nb_dEVEs_alone_Events','Nb_EVEs_shared_Events','Nb_dEVEs_shared_Events','Nb_dEVEs_alone','Nb_EVEs_alone')


#Add lifestyles infered by bayesian analysis
#Load the MA scenario tree to get posterior probabilities of node states
#Bigtable<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Table_with_MA_states.txt",sep=";",h=T)
Bigtable<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Table_with_MA_states2.txt",sep=";",h=T)
Bigtable$index <- as.character(Bigtable$index)
Bigtable <- Bigtable [c("label","parent","node","anc_state_1")]
colnames(Bigtable) <- c("label_bigtable","parent","node","anc_state_1")
#Merge informations
#Here there is an issue of node number bewwetn Bigtable and tree5 !!!
tree5<-merge(x=tree5,y=Bigtable,by=c("node"),all.y=T)
tree5$parent <- tree5$parent.x
library(phytools)
findMRCA(as.phylo(tree5),c("Eupelmus_urozonus","Eupelmus_kiefferi"))


findMRCA(as.phylo(tree5),c("Apis_mellifera","Platygaster_orseoliae"))

#Get node number
ggtree(as.treedata(tree5)) + geom_text(aes(label=node), hjust=-.3)


#names(tree5)[names(tree5)=="branch.length.x"] <- "branch.length"
#names(tree5)[names(tree5)=="parent.x"] <- "parent"
#tree5<-drop.tip(as.phylo(tree5), c("Tribolium_castaneum", "Anoplophora_glabripennis"), root.edge = 0)


periods_cust <- periods # get the periods object from deeptime
periods_cust$box_fill <- c("#e5e5e5", "white") # add alternating colors to the data frame

tree5$label <- gsub(" ","_",tree5$label)

library(ape)
#tree5<-drop.tip(as.phylo(tree5), "Unknown_sp")
p<-ggtree(as.treedata(tree5),color="grey40") + # plots tree
  geom_tiplab(offset=25,size=7,align=TRUE, linesize=.5,fontface="italic")


p<-flip(p, 187, 133)

#p<- p+geom_inset(pies2, x="branch", width=c(rep(0.06,length(unique(Pie_chart_tab$node))-1),0.12)  ,hjust=7,vjust=0,height=c(rep(0.06,length(unique(Pie_chart_tab$node))-1),0.12) )+
#  coord_cartesian()

p<- p+
  scale_x_reverse("Age (Ma)") +
  geom_rect(data = periods_cust, aes(xmin = -min_age, xmax = -max_age, ymin = -Inf, ymax = Inf),
            fill = periods_cust$box_fill, inherit.aes = FALSE)+# adds tip labels
  coord_geo(xlim = c(-320,150), ylim = c(0.4,127 + 0.4),
            neg = TRUE, abbrv = FALSE) +
  scale_x_continuous(breaks=seq(-380,0,20), labels=abs(seq(-380,0,20))) +
  theme_tree2()+
  theme(axis.text=element_text(size=18))

p$layers <- rev(p$layers)
# put the rectangles behind the tree
p<-revts(p)

anoleData<-read.csv("/Users/bguinet/Documents/All_informations_table_fam_busco_life_style.csv",sep=";",header=T)
anoleData$lifecycle1[anoleData$Species_name=="Megastigmus_dorsalis"]<-"ectoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Megastigmus_stigmatizans"]<-"ectoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Cecidostiba_semifascia"]<-"ectoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Cecidostiba_fungosa"]<-"ectoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Eretmocerus_eremicus"]<-"endoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Ormyrus_pomaceus"]<-"ectoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Ormyrus_nitidulus"]<-"ectoparasitoide"
anoleData<-select(anoleData,'Species_name','lifecycle1')
colnames(anoleData)<-c('species','lifecycle1')


anoleData<-anoleData[!anoleData$species %in% c("Andricus_quercusramuli","Neuroterus_quercusbaccarum","Andricus_inflator",
                                               "Andricus_curvator","Pseudoneuroterus_saliens","Trichopria_drosophilae","Andricus_grossulariae"),]


colors_associated_dsDNA<-c("freeliving"="#467fb3","ectoparasitoide"="#d8be03","endoparasitoide"="#06642e","red"="white")


rownames(anoleData) <- anoleData$species
p<-gheatmap(p, anoleData[c("lifecycle1")] , hjust=0, offset = 15, width=0.02,colnames=FALSE) +
  scale_fill_manual(values =colors_associated_dsDNA)   + theme(legend.position="none")

library(stringi)
library(ggplot2)

z <- as.character(p$data$Nb_EVEs_alone_Events)
z[is.na(z)]<-"NA"
z2 <- stri_pad_both(z, width = 3)
p$data$Nb_EVEs_alone_Events2 <- z2
p$data$Nb_EVEs_alone_Events2[p$data$Nb_EVEs_alone_Events2=="NA "]<-NA
p$data$Nb_EVEs_alone_Events2[p$data$Nb_EVEs_alone_Events2=="10 "]<-"10"
p$data$Nb_EVEs_alone_Events2[p$data$Nb_EVEs_alone_Events2=="11 "]<-"11"
p$data$Nb_EVEs_alone_Events2[p$data$Nb_EVEs_alone_Events2=="12 "]<-"12"
p$data$Nb_EVEs_alone_Events2[p$data$Nb_EVEs_alone_Events2=="13 "]<-"13"
p$data$Nb_EVEs_alone_Events2[p$data$Nb_EVEs_alone_Events2=="14 "]<-"14"
p$data$Nb_EVEs_alone_Events2[p$data$Nb_EVEs_alone_Events2=="15 "]<-"15"
p$data$Nb_EVEs_alone_Events2[p$data$Nb_EVEs_alone_Events2=="16 "]<-"16"

bb <- c(1, 3,5,10,20,40) # define breaks.
ll <- c("1","3","5","10","20","40")


p+ geom_text(aes(label=node), hjust=-.3)


for (i in 1:nrow(Pie_chart_tab_save)){
  if (Pie_chart_tab_save$Nb_EVEs[i]== 19) {
    Pie_chart_tab_save$Nb_2EVEs[i]<-0.08
  }else if (Pie_chart_tab_save$Nb_EVEs[i]== 15){
    Pie_chart_tab_save$Nb_2EVEs[i]<-0.075
  }else if (Pie_chart_tab_save$Nb_EVEs[i]== 6){
    Pie_chart_tab_save$Nb_2EVEs[i]<-0.065
  }else if (Pie_chart_tab_save$Nb_EVEs[i]== 4){
    Pie_chart_tab_save$Nb_2EVEs[i]<-0.055
  }else if (Pie_chart_tab_save$Nb_EVEs[i]== 2){
    Pie_chart_tab_save$Nb_2EVEs[i]<-0.045
  }else if (Pie_chart_tab_save$Nb_EVEs[i]== 1){
    Pie_chart_tab_save$Nb_2EVEs[i]<-0.035
  }
}


p<-  p +
  theme(axis.line.x = element_blank())+ # hide the x axis so it doesn't show up under the tip labels
  geom_label(aes(label = Nb_EVEs_alone_Events2,x=x), fill = "white",size = 5.5,nudge_x = 9.5 )+
  #geom_inset(pies2, x="branch", width=Pie_chart_tab_save2$Nb_2EVEs  ,hjust=0,vjust=0,height=Pie_chart_tab_save2$Nb_2EVEs )+
  geom_inset(pies, x="branch", width=Pie_chart_tab_save$Nb_2EVEs,hjust=0,vjust=0,height=Pie_chart_tab_save$Nb_2EVEs )+
  coord_cartesian(xlim=c(-320,150))+
  #geom_tippoint(aes(size=Nb_EVEs_alone,fill="black"),color='black')+
  #geom_tippoint(aes(size=Nb_dEVEs_alone,fill="red"),color='red')+ # adds tip labels
  #geom_nodepoint(aes(subset = !is.na(Nb_EVEs_shared), size = Nb_EVEs_shared, x = branch)) + # adds your data as points
  ggnewscale::new_scale_fill()+ #geom_nodepoint(aes(subset = !is.na(Nb_dEVEs_shared), size = Nb_dEVEs_shared, x = branch,fill="red"),color='red')+
  ggrepel::geom_label_repel(aes(label = Nb_EVEs_shared_Events,x=branch), fill = "white",size = 7.5,box.padding = 1.5, max.overlaps = Inf)+
  scale_size_continuous(name = "Prop.",
                        breaks = bb,
                        limits = c(0, 47),
                        labels = ll,
                        range = c(0, 15) )+
  theme_tree2() +# theme(legend.position="none") +
  theme(axis.text=element_text(size=18))# +theme(legend.text=element_text(size=30))

p<-p+ theme(plot.margin = unit(c(1,-27,1,1), "mm")) #Reduce the space between the first histogramm an the labels
#p+ geom_text(aes(label=node))
#You need to run the part below before to generate the table EVEs and dEVEs count

#Binary_EVEs_matrix_count<-as.data.frame(table(Binary_EVEs_matrix_count$Events_species))

library(plyr)
Binary_EVEs_matrix_count_table <- ddply(Binary_EVEs_matrix_count, .(Binary_EVEs_matrix_count$Events_species, Binary_EVEs_matrix_count$consensus_genomic_structure), nrow)
names(Binary_EVEs_matrix_count_table ) <- c("Events_species", "consensus_genomic_structure", "Freq")
#Binary_EVEs_matrix_count_table<-Binary_EVEs_matrix_count
#names(Binary_EVEs_matrix_count_table ) <- c("Events_species", "Freq")
Binary_EVEs_matrix_count_table $Events_species<-gsub("Unknown_species","Unknown_sp",Binary_EVEs_matrix_count_table $Events_species)


Binary_dEVEs_matrix_count_table  <- ddply(Binary_dEVEs_matrix_count, .(Binary_dEVEs_matrix_count$Events_species, Binary_dEVEs_matrix_count$consensus_genomic_structure), nrow)
names(Binary_dEVEs_matrix_count_table ) <- c("Events_species", "consensus_genomic_structure", "Freq")
#Binary_dEVEs_matrix_count_table<-Binary_dEVEs_matrix_count
#names(Binary_dEVEs_matrix_count_table ) <- c("Events_species", "Freq")
Binary_dEVEs_matrix_count_table $Events_species<-gsub("Unknown_species","Unknown_sp",Binary_dEVEs_matrix_count_table $Events_species)

#Get the same for events number
Binary_EVEs_matrix<-as.data.frame(table(Binary_EVEs_matrix_save$Events_species))
Binary_dEVEs_matrix<-as.data.frame(table(Binary_dEVEs_matrix_save$Events_species))

library(tidyverse)
library(phytools)
treeNW = as.tibble(read.newick("/Users/bguinet/Desktop/Papier_scientifique/Revbayes_anestral_states_estimates/data/TimeTree_mcmc_10_test_with_outgroup.nwk"))
treeNW$label<-gsub("Unknwon_species","Unknown_sp",treeNW$label)

Binary_EVEs_matrix_count_table <-as.data.frame(Binary_EVEs_matrix_count_table )
Binary_EVEs_dEVEs_matrix_count_table <-merge(Binary_EVEs_matrix_count_table ,Binary_dEVEs_matrix_count_table ,by="Events_species",all = T)
Labels<-as.data.frame(treeNW$label)
Labels<-Labels[!is.na(Labels),]



Binary_EVEs_matrix_count_table<-merge(Binary_EVEs_matrix_count_table,as.data.frame(Labels),by.x="Events_species",by.y="Labels",all = T)
Binary_EVEs_matrix_count_table[is.na(Binary_EVEs_matrix_count_table)]<-0

Binary_dEVEs_matrix_count_table<-merge(Binary_dEVEs_matrix_count_table,as.data.frame(Labels),by.x="Events_species",by.y="Labels",all = T)
Binary_dEVEs_matrix_count_table[is.na(Binary_dEVEs_matrix_count_table)]<-0


list_names<-p$data[order(p$data$y),decreasing = TRUE]$label
list_names<-rev(list_names[!is.na(list_names)])

Events_species <- c("Anoplophora_glabripennis", "Tribolium_castaneum")
consensus_genomic_structure <- c("dsDNA", "dsDNA")
Freq <- c(0,0)

df_outgroup <- data.frame(Events_species,consensus_genomic_structure,Freq)

Binary_EVEs_matrix_count_table<-rbind(Binary_EVEs_matrix_count_table,df_outgroup)
Binary_dEVEs_matrix_count_table<-rbind(Binary_dEVEs_matrix_count_table,df_outgroup)




#Get statistics

Binary_EVEs_matrix_count_table$Events_species <- factor(Binary_EVEs_matrix_count_table$Events_species, levels = rev(list_names))
Binary_dEVEs_matrix_count_table$Events_species <- factor(Binary_dEVEs_matrix_count_table$Events_species, levels = rev(list_names))

max(aggregate(Binary_EVEs_matrix_count_table$Freq, by=list(Category=Binary_EVEs_matrix_count_table$Events_species), FUN=sum)$x)
median(aggregate(Binary_EVEs_matrix_count_table$Freq, by=list(Category=Binary_EVEs_matrix_count_table$Events_species), FUN=sum)$x)
min(aggregate(Binary_EVEs_matrix_count_table$Freq, by=list(Category=Binary_EVEs_matrix_count_table$Events_species), FUN=sum)$x)

max(aggregate(Binary_dEVEs_matrix_count_table$Freq, by=list(Category=Binary_dEVEs_matrix_count_table$Events_species), FUN=sum)$x)
median(aggregate(Binary_dEVEs_matrix_count_table$Freq, by=list(Category=Binary_dEVEs_matrix_count_table$Events_species), FUN=sum)$x)
min(aggregate(Binary_dEVEs_matrix_count_table$Freq, by=list(Category=Binary_dEVEs_matrix_count_table$Events_species), FUN=sum)$x)


#Extract family EVE distribution for each species of the phylogeny

subEnv_table2 <- Env_table2[!duplicated(Env_table2[c('Clustername','Species_name','Event')]),]

subEnv_table2$Species_name <- factor(subEnv_table2 $Species_name, levels = rev(list_names))

#Remove Harpegnathos

#subEnv_table2<-subEnv_table2[! (subEnv_table2$Species_name=="Harpegnathos_saltator" &  subEnv_table2$family=="Iridoviridae" |  subEnv_table2$Species_name=="Harpegnathos_saltator"  & subEnv_table2$family=="Apis_filamentous-like"),]

Color_palette<-c("Phycodnaviridae"="#6c919e","Nudiviridae"="#004D80","Iridoviridae"="#67BCEB",'LbFV_like'="#0276BB","Poxviridae"="#00FFFF","Apis_filamentous-like"="#485F92","IVSPERs"="#B5DCE5","Baculoviridae"="#4773B8","Herpesviridae"="#3BABA8","Nimaviridae"="#3BABA8","Papillomaviridae"="#3BABA8","Ascoviridae"="#3BABA8","Asfaviridae"="#3BABA8","Caulimoviridae"="#3BABA8","Tymoviridae"="#3BABA8","Phycodnaviridae"="#3BABA8",
                 "Parvoviridae"="#E52421","Circoviridae"="#F4A09D","Genomoviridae"="#971D14","Partiti-Picobirna"="#734b00",
                 "Artoviridae"="#5BB888","Chuviridae"="#E8E100","Nyamiviridae"="#77FF46","Orthomyxoviridae"="#E9EEBF","Xinmoviridae"="#4F9436","Lispiviridae"="#BFDFCD","Phenuiviridae"="#2C5933","Bornaviridae"="#C0CE4A","Virgaviridae"="#ADF2CE","Bornaviridae"="#ADF2CE","Tombusviridae"="#ADF2CE","Qinviridae"="#ADF2CE","Peribunyaviridae"="#ADF2CE","Mypoviridae"="#ADF2CE","Phasmaviridae"="#ADF2CE","Narna-Levi"="#ADF2CE","Mono-Chu"="#ADF2CE","Hepe-Virga"="#ADF2CE","Abisko-like"="#ADF2CE","Tombusviridae"="#ADF2CE","Luteo-Sobemo"="#ADF2CE","Negevirus-like"="#ADF2CE",
                 "Totiviridae"="#FC990B","Reoviridae"="#C98A4B","Partitiviridae"="#F4D0AB",
                 "Unknown"="black","Unclassified ssDNA"="#d6d6d6","Unclassified ssRNA"="#9e9e9e","Unclassified dsRNA"="#525251","Nb_dEVEs"="black","nonNb_dEVEs"="transparent" )


# Only keep one EVE per species per cluster

subEnv_table2  <- subEnv_table2 [!duplicated(subEnv_table2[c("ID","Clustername",'Species_name')]),]

subEnv_table2 <- subEnv_table2 %>%
  arrange(desc(Mean_dNdS)) %>%
  arrange(desc(TPM_all))

subEnv_table2<- subEnv_table2[!subEnv_table2$Clustername=="Cluster3513",]
p9 <- ggplot(as.data.frame(table(subEnv_table2$Species_name,subEnv_table2$consensus_family)), aes(y = Freq, x = Var1))+
  geom_col(aes(fill = Var2), width=.5,color="black",size=0.1)  + coord_flip()+
  scale_fill_manual(values=Color_palette)+
  theme_bw()+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank())+
  theme(axis.text=element_text(size=18))+
  theme(legend.position="none")+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 55),breaks = scales::pretty_breaks(n = 20))+
  ylab("Number of EVEs")+
  theme(axis.title.x = element_text(size=22))


subEnv_table2_dEVEs <-  subEnv_table2[subEnv_table2$FDR_pvalue_dNdS == 1 & as.numeric(as.character(subEnv_table2$Mean_dNdS))+as.numeric(as.character(subEnv_table2$SE_dNdS)) < 1 & subEnv_table2$pseudogenized ==0 & subEnv_table2$ORF_perc>1 |  subEnv_table2$TPM_all>=1000 & subEnv_table2$pseudogenized ==0 & subEnv_table2$ORF_perc>1,]

p10 <- ggplot(as.data.frame(table(subEnv_table2_dEVEs$Species_name ,subEnv_table2_dEVEs$consensus_family)), aes(y = Freq, x = Var1))+
  geom_col(aes(fill = Var2), width=.5,color="black",size=0.1)  + coord_flip()+
  theme_bw()+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position="none")+
  theme(panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank())+
  scale_fill_manual(values=Color_palette)+
  theme(axis.text=element_text(size=18))+
 ylab("Number of dEVEs")+
  theme(axis.title.x = element_text(size=22))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 30),breaks = scales::pretty_breaks(n = 20))

#Add famillies info

anoleData<-read.csv("/Users/bguinet/Documents/All_informations_table_fam_busco_life_style.csv",sep=";",header=T)
anoleData$lifecycle1[anoleData$Species_name=="Megastigmus_dorsalis"]<-"ectoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Megastigmus_stigmatizans"]<-"ectoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Cecidostiba_semifascia"]<-"ectoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Cecidostiba_fungosa"]<-"ectoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Eretmocerus_eremicus"]<-"endoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Ormyrus_pomaceus"]<-"ectoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Ormyrus_nitidulus"]<-"ectoparasitoide"
anoleData$Species_name[anoleData$Species_name=="Unknown_species"]<-"Unknown_sp"
anoleData$SuperFamily[anoleData$SuperFamily=="Chalcidoidea ? "] <- "Unknown"
anoleData$SuperFamily[anoleData$SuperFamily=="formicoidea"] <- "Formicoidea"
names(anoleData)[names(anoleData)=="Species_name"] <- "Var1"

plot_subEnv_table2_EVEs<-as.data.frame(table(subEnv_table2$Species_name ,subEnv_table2$family))
plot_subEnv_table2_EVEs<-dplyr::left_join(plot_subEnv_table2_EVEs,anoleData[,c("SuperFamily","Var1")], by = c(Var1 = "Var1"))
plot_subEnv_table2_EVEs$Var1 <- factor(plot_subEnv_table2_EVEs$Var1, levels = unique(as.data.frame(table(subEnv_table2$Species_name,subEnv_table2$family))$Var1))
plot_subEnv_table2_EVEs$SuperFamily <- factor(plot_subEnv_table2_EVEs$SuperFamily, levels = rev(unique(plot_subEnv_table2_EVEs$SuperFamily)))
plot_subEnv_table2_EVEs$Freq<-0

library(ggplot2)
library(ggtree)
p11 <-  ggplot(plot_subEnv_table2_EVEs, aes(y = Freq, x = Var1))+
  geom_col(aes(fill = Var2), width=.5,color="black",size=0.1)  + coord_flip()+
  theme(axis.text=element_text(size=18))+
  theme(axis.title.x = element_text(size=20))+
  scale_y_continuous(expand = c(0, 0), limits = c(0),breaks = scales::pretty_breaks(n = 1))+
  theme_tree2(panel.grid.major   = element_line(color="grey", size=.2),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank())+
  facet_grid(SuperFamily~., scales = "free", space = "free") +
  theme(strip.text.y = element_text(size = 22,angle = 0,face = "italic"))+
  scale_x_discrete(expand = c(0, 0), position = "top") +
  theme(plot.margin = margin(0, 0, 0, 0),
        axis.line.y.right = element_line(size = 3),
        strip.background = element_blank(),
        strip.placement = "outside")+  theme(legend.position="none") +theme(axis.title.y=element_blank(),
                                                                           axis.text.y=element_blank())+
   theme(panel.spacing.x=unit(1, "lines"),panel.spacing.y=unit(1, "lines"))



#p9<-ggplot(data=Binary_EVEs_matrix_count_table, aes(x=Events_species, y=Freq, fill=consensus_genomic_structure)) +
#  geom_bar(stat="identity", width=.5,color="black",size=0.1)  +coord_flip() +
#  xlab("") +
#  theme_bw()+
#  theme(axis.title.y=element_blank(),
#        axis.text.y=element_blank(),
#        axis.ticks.y=element_blank())+
# theme(panel.grid.minor.y=element_blank(),
#       panel.grid.major.y=element_blank())+scale_fill_manual(values =c("dsDNA"="#0072DB","dsRNA"="#DE9100","ssDNA"="#00BBDE99","ssRNA"="#DBC900","Unknown"="grey")) +
# theme(legend.position="none")


#p10<-ggplot(data=Binary_dEVEs_matrix_count_table, aes(x=Events_species, y=Freq, fill=consensus_genomic_structure)) +
#  geom_bar(stat="identity", width=.5,color="black",size=0.1)  +coord_flip() +
#  xlab("") +
#theme_bw()+
#  theme(axis.title.y=element_blank(),
#        axis.text.y=element_blank(),
#        axis.ticks.y=element_blank())+
#  theme(panel.grid.minor.y=element_blank(),
#        panel.grid.major.y=element_blank())+scale_fill_manual(values =c("dsDNA"="#0072DB","dsRNA"="#DE9100","ssDNA"="#00BBDE99","ssRNA"="#DBC900","Unknown"="grey")) +
#  theme(legend.position="none")


#Change names
p$data$label[p$data$label=="Ganaspis_ganaspis"]<-"Ganaspis_sp"
p$data$label[p$data$label=="Copidosoma_sp"]<-"Copidosoma_aretas"
library(ggplot2)

cowplot::plot_grid(p, p9,p10,p11,ncol=4,align = "h",axis = "bt", rel_widths=c(0.5,0.3,0.15,0.1))

#35/50
##########

ggsave("/Users/bguinet/Desktop/Papier_scientifique/Main_phylogeny_figure2.pdf", width = 35, height = 50, units = "in",limitsize = FALSE)




#############Now we will create the dataframe with events in a 0/1 matrix
detachAllPackages()
library(dplyr)
library(phylogram)
library(tidyr)

anoleData<-read.csv("/Users/bguinet/Documents/All_informations_table_fam_busco_life_style.csv",sep=";",header=T)
anoleData$lifecycle1[anoleData$Species_name=="Megastigmus_dorsalis"]<-"ectoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Megastigmus_stigmatizans"]<-"ectoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Cecidostiba_semifascia"]<-"ectoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Cecidostiba_fungosa"]<-"ectoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Eretmocerus_eremicus"]<-"endoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Ormyrus_pomaceus"]<-"ectoparasitoide"
anoleData$lifecycle1[anoleData$Species_name=="Ormyrus_nitidulus"]<-"ectoparasitoide"

tree2=read.dendrogram(file = "/Users/bguinet/Desktop/Papier_scientifique/Revbayes_anestral_states_estimates/data/TmeTree_mcmc_10_test.nwk")

clade_order <- order.dendrogram(as.dendrogram(tree2))
clade_name <- labels(as.dendrogram(tree2))
clade_position <- data.frame(clade_name, clade_order)
clade_position <- clade_position[order(clade_position$clade_order),]

####get color depending of lifstyle
row.names(anoleData) <- anoleData$Species_name
new_order <- match(clade_position$clade_name,
                   row.names(anoleData))
anoleData <- anoleData[new_order,]
anoleData<- anoleData[!is.na(anoleData$Order),]
anoleData$New_query_bis<-row.names(anoleData)


anoleData<-anoleData[!anoleData$Species_name %in% c("Andricus_quercusramuli","Neuroterus_quercusbaccarum","Andricus_inflator",
                                                    "Andricus_curvator","Pseudoneuroterus_saliens","Trichopria_drosophilae","Andricus_grossulariae","Tribolium_castaneum","Annoplophora_glabripennis"),]

#EVEs matrix Events
Binary_table_alone_EVEs_event<-EVEs_alone_event_df_EVEs_count[c("species","Clusters","Event","consensus_genomic_structure")]
colnames(Binary_table_alone_EVEs_event)<-c("Events_species","Clustername","Events","consensus_genomic_structure")

if (Remove_controls=="yes"){
  write.table(Binary_table_alone_EVEs_event,"/Users/bguinet/Desktop/Papier_scientifique/Table_EVEs_events_all_alone_without_control.tab",sep=";")
}else{
  write.table(Binary_table_alone_EVEs_event,"/Users/bguinet/Desktop/Papier_scientifique/Table_EVEs_events_all_alone_without_cont.tab",sep=";")
}


Binary_table_shared_EVEs_event <-  EVEs_shared_event_df_EVEs_count[c("Events_species2","Clustername","Event","consensus_genomic_structure")]
colnames(Binary_table_shared_EVEs_event)<- c("Events_species","Clustername","Events","consensus_genomic_structure")


EVEs_matrix <-bind_rows(Binary_table_shared_EVEs_event[c("Events_species","Clustername","Events","consensus_genomic_structure")],Binary_table_alone_EVEs_event[c("Events_species","Clustername","Events","consensus_genomic_structure")])
EVEs_matrix$Event_numbers <- paste0("Event_",rownames(EVEs_matrix))


#dEVEs matrix Events
Binary_table_alone_dEVEs_event<-EVEs_alone_event_df_EVEs_count[EVEs_alone_event_df_EVEs_count$Nb_dEVEs>=1,][c("species","Clusters","Event",'consensus_genomic_structure')]
colnames(Binary_table_alone_dEVEs_event)<-c("Events_species","Clustername","Events",'consensus_genomic_structure')

if (Remove_controls=="yes"){
  write.table(Binary_table_alone_dEVEs_event,"/Users/bguinet/Desktop/Papier_scientifique/Table_dEVEs_events_all_alone_without_control.tab",sep=";")
}else{
  write.table(Binary_table_alone_dEVEs_event,"/Users/bguinet/Desktop/Papier_scientifique/Table_dEVEs_events_all_alone_without_cont.tab",sep=";")
}



Binary_table_shared_dEVEs_event <- EVEs_shared_event_df_EVEs_count[EVEs_shared_event_df_EVEs_count$Nb_dEVEs>=1,][c("Events_species2","Clustername","Event",'consensus_genomic_structure')]
colnames(Binary_table_shared_dEVEs_event)<- c("Events_species","Clustername","Events",'consensus_genomic_structure')

dEVEs_matrix <-bind_rows(Binary_table_shared_dEVEs_event[c("Events_species","Clustername","Events","consensus_genomic_structure")],Binary_table_alone_dEVEs_event[c("Events_species","Clustername","Events")])
dEVEs_matrix$Event_numbers <- paste0("Event_",rownames(dEVEs_matrix))


#EVEs and dEVEs matrix gene count ,
Binary_table_alone_EVEs_count<-EVEs_alone_event_df_EVEs_count[c("species","Clusters","Nb_EVEs","consensus_genomic_structure")]
colnames(Binary_table_alone_EVEs_count)<-c("Events_species","Clustername","Nb_EVEs","consensus_genomic_structure")

if (Remove_controls=="yes"){
  write.table(Binary_table_alone_EVEs_count,"/Users/bguinet/Desktop/Papier_scientifique/Table_EVEs_counts_all_alone_without_control.tab",sep=";")
}else{
  write.table(Binary_table_alone_EVEs_count,"/Users/bguinet/Desktop/Papier_scientifique/Table_EVEs_counts_all_alone_without_cont.tab",sep=";")
}


Binary_table_shared_EVEs_count <-  EVEs_shared_event_df_EVEs_count[c("Events_species2","Clustername","Nb_EVEs","consensus_genomic_structure")]
colnames(Binary_table_shared_EVEs_count)<- c("Events_species","Clustername","Nb_EVEs","consensus_genomic_structure")

EVEs_matrix_count <-bind_rows(Binary_table_shared_EVEs_count[c("Events_species","Clustername","Nb_EVEs","consensus_genomic_structure")],Binary_table_alone_EVEs_count[c("Events_species","Clustername","Nb_EVEs","consensus_genomic_structure")])

EVEs_matrix_count<-data.frame(EVEs_matrix_count[rep(seq_len(dim(EVEs_matrix_count)[1]), EVEs_matrix_count$Nb_EVEs), c(1,2,3,4), drop = FALSE], row.names=NULL)
EVEs_matrix_count$EVE_numbers <- paste0("EVE_",rownames(EVEs_matrix_count))


#dEVEs matrix count
Binary_table_alone_dEVEs_count<-EVEs_alone_event_df_EVEs_count[EVEs_alone_event_df_EVEs_count$Nb_dEVEs>=1,][c("species","Clusters","Nb_dEVEs","consensus_genomic_structure")]
colnames(Binary_table_alone_dEVEs_count)<-c("Events_species","Clustername","Nb_dEVEs","consensus_genomic_structure")

if (Remove_controls=="yes"){
  write.table(Binary_table_alone_dEVEs_count,"/Users/bguinet/Desktop/Papier_scientifique/Table_dEVEs_counts_all_alone_without_control.tab",sep=";")
}else{
  write.table(Binary_table_alone_dEVEs_count,"/Users/bguinet/Desktop/Papier_scientifique/Table_dEVEs_counts_all_alone_without_cont.tab",sep=";")
}


Binary_table_shared_dEVEs_count <-  EVEs_shared_event_df_EVEs_count[EVEs_shared_event_df_EVEs_count$Nb_dEVEs>=1,][c("Events_species2","Clustername","Nb_dEVEs","consensus_genomic_structure")]
colnames(Binary_table_shared_dEVEs_count)<- c("Events_species","Clustername","Nb_dEVEs","consensus_genomic_structure")

dEVEs_matrix_count <-bind_rows(Binary_table_shared_dEVEs_count[c("Events_species","Clustername","Nb_dEVEs","consensus_genomic_structure")],Binary_table_alone_dEVEs_count[c("Events_species","Clustername","Nb_dEVEs","consensus_genomic_structure")])
dEVEs_matrix_count<-data.frame(dEVEs_matrix_count[rep(seq_len(dim(dEVEs_matrix_count)[1]), dEVEs_matrix_count$Nb_dEVEs), c(1,2,3,4), drop = FALSE], row.names=NULL)
dEVEs_matrix_count$dEVE_numbers <- paste0("dEVE_",rownames(dEVEs_matrix_count))


##################################################
#Save informations into one unique matrix table ###
#################################################
EVEs_matrix_long <- EVEs_matrix %>%
  tidyr::separate_rows( Events_species, sep = "," ) %>%
  tidyr::separate_rows( Events, sep = "," )

colnames(EVEs_matrix_long)<-c("Species_name","Clustername","Event","Event_numbers","consensus_genomic_structure" )
Restricted_candidat_tosave<-Restricted_candidat %>% dplyr::left_join(EVEs_matrix_long,by=c("Species_name","Event"))

write.table(Restricted_candidat_tosave[c('Clustername.x','Species_name','Scaff_name','qstart','qend','target','evalue','consensus_genomic_structure','family','Domain_description','Gene.ontology..GO.','Gene.names','Protein.names','Consensus_function','Median_cov_depth_candidat','Median_cov_depth_BUSCO','FDR_pvalue_cov_median','Nb_repeat','count_eucaryote','x','Event','Event_numbers','Boot','Nsp','FDR_pvalue_dNdS','Pvalue_dNdS','Mean_dNdS','SE_dNdS','TPM_all','pseudogenized','ORF_len','ORF_perc','ORF_start','ORF_end','sequence')],"/Users/bguinet/Desktop/these/Table_summary_all_good_scaffols.txt",sep=";",row.names=FALSE)
t<-Restricted_candidat_tosave[!duplicated(Restricted_candidat_tosave$Event_numbers)]
as.data.frame(table(t$consensus_family))

#######

######################
#EVES matrix Events
######################
detachAllPackages()
library(dplyr)
library(tidyr)
Species_name_df<-anoleData[c('Species_name')]


colnames(Species_name_df)<-'Events_species'
Binary_EVEs_matrix <-as.data.frame(EVEs_matrix)[c("Events_species","Event_numbers","consensus_genomic_structure")] %>% separate_rows(Events_species,sep=",") %>%
  mutate(Val = 1) %>% type.convert(as.is = TRUE)

write.table(Binary_EVEs_matrix,"/Users/bguinet/Desktop/Papier_scientifique/Table_EVEs_events_all.tab",sep=";")

#Get number of event summary

max(table(Binary_EVEs_matrix[c("Events_species","Val")]))
median(table(Binary_EVEs_matrix[c("Events_species","Val")]))

Binary_EVEs_matrix_save   <-  Binary_EVEs_matrix

Binary_EVEs_matrix<-tidyr::pivot_wider(Binary_EVEs_matrix, names_from = Event_numbers, values_from = Val,
                                       values_fn = list(Val = ~any(. > 0)), values_fill = FALSE)%>% right_join(Species_name_df)

#Replace TRUE an FALSE
Binary_EVEs_matrix<-as.data.frame(Binary_EVEs_matrix)
Binary_EVEs_matrix[Binary_EVEs_matrix=="FALSE"]<-0
Binary_EVEs_matrix[Binary_EVEs_matrix=="TRUE"]<-1
Binary_EVEs_matrix[Binary_EVEs_matrix=="na"]<-0
Binary_EVEs_matrix[is.na(Binary_EVEs_matrix)]<-0


library(MASS)

if (Only_mulptiple_EVE=="yes"){
  write.table(Binary_EVEs_matrix,"/Users/bguinet/Desktop/these/Mapping_states_analysis/Matrix_Events/Matrix_Events_dsDNA_EVEs_multiple.txt",sep=";",row.names = F)
}else{
  write.table(Binary_EVEs_matrix,"/Users/bguinet/Desktop/these/Mapping_states_analysis/Matrix_Events/Matrix_Events_dsDNA_EVEs.txt",sep=";",row.names = F)
}

######################
#dEVES matrix Events
######################
detachAllPackages()
library(dplyr)
library(tidyr)

Species_name_df<-select(anoleData,'Species_name')
colnames(Species_name_df)<-'Events_species'
Binary_dEVEs_matrix <-as.data.frame(dEVEs_matrix)[c("Events_species","Event_numbers","consensus_genomic_structure")] %>% separate_rows(Events_species,sep=",") %>%
  mutate(Val = 1) %>% type.convert(as.is = TRUE)

write.table(Binary_dEVEs_matrix,"/Users/bguinet/Desktop/Papier_scientifique/Table_dEVEs_events_all.tab",sep=";")

Binary_dEVEs_matrix_save   <-  Binary_dEVEs_matrix

Binary_dEVEs_matrix<-tidyr::pivot_wider(Binary_dEVEs_matrix, names_from = Event_numbers, values_from = Val,
                                        values_fn = list(Val = ~any(. > 0)), values_fill = FALSE)%>%
  right_join(Species_name_df)

#Replace TRUE an FALSE
Binary_dEVEs_matrix<-as.data.frame(Binary_dEVEs_matrix)
Binary_dEVEs_matrix[Binary_dEVEs_matrix=="FALSE"]<-0
Binary_dEVEs_matrix[Binary_dEVEs_matrix=="TRUE"]<-1
Binary_dEVEs_matrix[Binary_dEVEs_matrix=="na"]<-0
Binary_dEVEs_matrix[is.na(Binary_dEVEs_matrix)]<-0

if (Only_mulptiple_EVE=="yes"){
  write.table(Binary_dEVEs_matrix,"/Users/bguinet/Desktop/these/Mapping_states_analysis/Matrix_Events/Matrix_Events_dsDNA_dEVEs_multiple.txt",sep=";",row.names = F)
}else{
  write.table(Binary_dEVEs_matrix,"/Users/bguinet/Desktop/these/Mapping_states_analysis/Matrix_Events/Matrix_Events_dsDNA_dEVEs.txt",sep=";",row.names = F)
}
######################
#EVES matrix count
######################
detachAllPackages()
library(dplyr)
library(tidyr)

Species_name_df<-anoleData[c('Species_name')]
colnames(Species_name_df)<-'Events_species'
Binary_EVEs_matrix_count <-as.data.frame(EVEs_matrix_count)[c("Events_species","EVE_numbers","consensus_genomic_structure")] %>% separate_rows(Events_species,sep=",") %>%
  mutate(Val = 1) %>% type.convert(as.is = TRUE)

write.table(Binary_EVEs_matrix_count,"/Users/bguinet/Desktop/Papier_scientifique/Table_EVEs_counts_all.tab",sep=";")

Binary_EVEs_matrix_count<-tidyr::pivot_wider(Binary_EVEs_matrix_count, names_from = EVE_numbers, values_from = Val,
                                             values_fn = list(Val = ~any(. > 0)), values_fill = FALSE)%>%
  right_join(Species_name_df)

#Replace TRUE an FALSE
Binary_EVEs_matrix_count<-as.data.frame(Binary_EVEs_matrix_count)
Binary_EVEs_matrix_count[Binary_EVEs_matrix_count=="FALSE"]<-0
Binary_EVEs_matrix_count[Binary_EVEs_matrix_count=="TRUE"]<-1
Binary_EVEs_matrix_count[Binary_EVEs_matrix_count=="na"]<-0
Binary_EVEs_matrix_count[is.na(Binary_EVEs_matrix_count)]<-0

if (Remove_controls=="yes"){
  write.table(Binary_EVEs_matrix_count,"/Users/bguinet/Desktop/these/Mapping_states_analysis/Matrix_Events/Matrix_Events_dsDNA_EVEs_count_without_controls.txt",sep=";",row.names = F)
}else{
  write.table(Binary_EVEs_matrix_count,"/Users/bguinet/Desktop/these/Mapping_states_analysis/Matrix_Events/Matrix_Events_dsDNA_EVEs_count.txt",sep=";",row.names = F)
}

######################
#dEVES matrix count
######################
detachAllPackages()
library(dplyr)
library(tidyr)

Species_name_df<-anoleData[c('Species_name')]
colnames(Species_name_df)<-'Events_species'
Binary_dEVEs_matrix_count <-as.data.frame(dEVEs_matrix_count)[c("Events_species","dEVE_numbers","consensus_genomic_structure")] %>% separate_rows(Events_species,sep=",") %>%
  mutate(Val = 1) %>% type.convert(as.is = TRUE)

write.table(Binary_dEVEs_matrix_count,"/Users/bguinet/Desktop/Papier_scientifique/Table_dEVEs_counts_all.tab",sep=";")

Binary_dEVEs_matrix_count<-tidyr::pivot_wider(Binary_dEVEs_matrix_count, names_from = dEVE_numbers, values_from = Val,
                                              values_fn = list(Val = ~any(. > 0)), values_fill = FALSE)%>%
  right_join(Species_name_df)

#Replace TRUE an FALSE
Binary_dEVEs_matrix_count<-as.data.frame(Binary_dEVEs_matrix_count)
Binary_dEVEs_matrix_count[Binary_dEVEs_matrix_count=="FALSE"]<-0
Binary_dEVEs_matrix_count[Binary_dEVEs_matrix_count=="TRUE"]<-1
Binary_dEVEs_matrix_count[Binary_dEVEs_matrix_count=="na"]<-0
Binary_dEVEs_matrix_count[is.na(Binary_dEVEs_matrix_count)]<-0

if (Remove_controls=="yes"){
  write.table(Binary_dEVEs_matrix_count,"/Users/bguinet/Desktop/these/Mapping_states_analysis/Matrix_Events/Matrix_Events_dsDNA_dEVEs_count_without_controls.txt",sep=";",row.names = F)
}else{
  write.table(Binary_dEVEs_matrix_count,"/Users/bguinet/Desktop/these/Mapping_states_analysis/Matrix_Events/Matrix_Events_dsDNA_dEVEs_count.txt",sep=";",row.names = F)
}
