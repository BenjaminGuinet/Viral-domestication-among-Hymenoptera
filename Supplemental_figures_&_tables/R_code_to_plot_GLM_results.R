---
title: "Analysis of nb EVEs"
output:
  pdf_document: default
  html_notebook: default
---


# New theme paramters
myt <- ttheme_default(
         # Use hjust and x to left justify the text
         # Alternate the row fill colours
                 core = list(fg_params=list(hjust = 1, x=1),
                             bg_params=list(fill=c("grey", "white"))),

         # Change column header to white text and red background
                 colhead = list(fg_params=list(col="white"),
                                bg_params=list(fill="black"))
 )


```{r}
#Open the data
data <- read.table("/Users/bguinet/Desktop/Data_jags.txt", sep=";")
data$Lifestyle<-as.character(data$Lifestyle)
```


In the following, we will write a generalized linear model using a Bayesian approach in which several response variables (number of EVEs, number of dEVEs, number of event EVEs & number of event dEVEs) will be a function of the explanatory variables:

- 1: Lifestyle comprising three modalities (Free-living, Endo-parasitoid and Ecto-parasitoid).
- 2: Branch size (counted in millions of years).

We previously included an interaction element between branch sizes and lifestyles in the model, the credibility intervals were too large so we eliminated them from the model and added them as additive effect.

**Model : Nb_EVEs_Events ~ Lifestyle + branch_length**

Since the results showed an overdispersion due to a too large representation of zero, we opted for a zero-inflated negative binomial model. The main motivation for zero-inflated count models is that our data present overdispersion and excess zeros which correspond to leaves of the hymenopteran tree that are too old for the algorythm to detect integration events. Thus we use the functio brm from brms package to write a negative binomial Zero-inflated count models providing a way of modeling the excess zeros in addition to allowing for overdispersion.




```{r}

# corrected dEVEs rates dsDNA only A
library(bayestestR)
library(ggstatsplot)
library(bayesplot)
library(ggpubr)
library(gridExtra)
library(grid)
library(ggprism)
library(rstatix)

##
##################################
# ## dEVEs Events all posteriors ##
##################################


Table_dsDNA_A_rate_corrected<- read.table("/Users/bguinet/Desktop/output_dsDNA_A_rate_corrected_all/Posteriors_dEVEs_Events_ALL.txt",sep=";",h=T)


Table_dsDNA_A_rate_corrected$b_value2_coef_inf<- exp(as.numeric(Table_dsDNA_A_rate_corrected$b_value2))
Table_dsDNA_A_rate_corrected$b_value3_coef_inf<- exp(as.numeric(Table_dsDNA_A_rate_corrected$b_value3))
Table_dsDNA_A_rate_corrected<-Table_dsDNA_A_rate_corrected[!Table_dsDNA_A_rate_corrected$b_value2_coef_inf=="Inf",]
Table_dsDNA_A_rate_corrected<-Table_dsDNA_A_rate_corrected[!Table_dsDNA_A_rate_corrected$b_value3_coef_inf=="Inf",]


Table_dsDNA_A_rate_corrected_endo <- as.data.frame(Table_dsDNA_A_rate_corrected$b_value2)

Table_dsDNA_A_rate_corrected_endo$Parameter <- "Endo-parasitoid"
colnames(Table_dsDNA_A_rate_corrected_endo)<- c("coef","Parameter")
Table_dsDNA_A_rate_corrected_endo$coef2 <- as.numeric(Table_dsDNA_A_rate_corrected_endo$coef)
Table_dsDNA_A_rate_corrected_endo$coef <- as.numeric(exp(Table_dsDNA_A_rate_corrected_endo$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_A_rate_corrected_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_A_rate_corrected_endo$rope <- posterior_tab*100
Table_dsDNA_A_rate_corrected_endo$med <- median(Table_dsDNA_A_rate_corrected_endo$coef)
Table_dsDNA_A_rate_corrected_endo$pd <- pd(Table_dsDNA_A_rate_corrected_endo$coef2)

Table_dsDNA_A_rate_corrected_ecto <- as.data.frame(Table_dsDNA_A_rate_corrected$b_value3)
Table_dsDNA_A_rate_corrected_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_dsDNA_A_rate_corrected_ecto)<- c("coef","Parameter")
Table_dsDNA_A_rate_corrected_ecto$coef2 <- as.numeric(Table_dsDNA_A_rate_corrected_ecto$coef)
Table_dsDNA_A_rate_corrected_ecto$coef <- as.numeric(exp(Table_dsDNA_A_rate_corrected_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_A_rate_corrected_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_A_rate_corrected_ecto$rope <- posterior_tab*100
Table_dsDNA_A_rate_corrected_ecto$med <- median(Table_dsDNA_A_rate_corrected_ecto$coef)
Table_dsDNA_A_rate_corrected_ecto$pd <- pd(Table_dsDNA_A_rate_corrected_ecto$coef2)



Table_dsDNA_A_rate_corrected_ecto_endo<-rbind(Table_dsDNA_A_rate_corrected_endo,Table_dsDNA_A_rate_corrected_ecto)

#Create summary table
statistics_Table_dsDNA_A_dEVEs_EVents= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("dsDNA_EVEs_Events",median(Table_dsDNA_A_rate_corrected_endo$coef),ci(Table_dsDNA_A_rate_corrected_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_A_rate_corrected_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_A_rate_corrected_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_dsDNA_A_rate_corrected_endo$pd[1]),
                             Ecto_parasitoid=c("dsDNA_EVEs_Events",median(Table_dsDNA_A_rate_corrected_ecto$coef),ci(Table_dsDNA_A_rate_corrected_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_A_rate_corrected_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_A_rate_corrected_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_dsDNA_A_rate_corrected_ecto$pd[1])))


colnames(statistics_Table_dsDNA_A_dEVEs_EVents) <- statistics_Table_dsDNA_A_dEVEs_EVents[1,]
statistics_Table_dsDNA_A_dEVEs_EVents<-statistics_Table_dsDNA_A_dEVEs_EVents[-1,]


Table_dsDNA_A_rate_corrected_ecto_endo$Parameter2 <- paste0(Table_dsDNA_A_rate_corrected_ecto_endo$Parameter,"\n %pd = ",round(Table_dsDNA_A_rate_corrected_ecto_endo$pd*100,2) ,"% ")

stat.test_dsDNA_A_rate_corrected<-Table_dsDNA_A_rate_corrected_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_dsDNA_A_rate_corrected <- stat.test_dsDNA_A_rate_corrected %>% add_xy_position(x = "Parameter2")
stat.test_dsDNA_A_rate_corrected$Prop <- paste0(round((nrow(Table_dsDNA_A_rate_corrected[Table_dsDNA_A_rate_corrected$b_value2 > Table_dsDNA_A_rate_corrected$b_value3,]) / nrow(Table_dsDNA_A_rate_corrected))*100,digits = 2),"%")

boxplot_dsDNA_A_dEVEs_rate_corrected_Events <- ggplot(Table_dsDNA_A_rate_corrected_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="dEvents",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_dsDNA_A_rate_corrected, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 9, size=10
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,11))


##
##################################
# ## dEVEs Counts all posteriors ##
##################################

Table_dsDNA_A_rate_corrected<- read.table("/Users/bguinet/Desktop/output_dsDNA_A_rate_corrected_all/Posteriors_dEVEs_Counts_ALL.txt",sep=";",h=T)


Table_dsDNA_A_rate_corrected$b_value2_coef_inf<- exp(as.numeric(Table_dsDNA_A_rate_corrected$b_value2))
Table_dsDNA_A_rate_corrected$b_value3_coef_inf<- exp(as.numeric(Table_dsDNA_A_rate_corrected$b_value3))
Table_dsDNA_A_rate_corrected<-Table_dsDNA_A_rate_corrected[!Table_dsDNA_A_rate_corrected$b_value2_coef_inf=="Inf",]
Table_dsDNA_A_rate_corrected<-Table_dsDNA_A_rate_corrected[!Table_dsDNA_A_rate_corrected$b_value3_coef_inf=="Inf",]


Table_dsDNA_A_rate_corrected_endo <- as.data.frame(Table_dsDNA_A_rate_corrected$b_value2)

Table_dsDNA_A_rate_corrected_endo$Parameter <- "Endo-parasitoid"
colnames(Table_dsDNA_A_rate_corrected_endo)<- c("coef","Parameter")
Table_dsDNA_A_rate_corrected_endo$coef2 <- as.numeric(Table_dsDNA_A_rate_corrected_endo$coef)
Table_dsDNA_A_rate_corrected_endo$coef <- as.numeric(exp(Table_dsDNA_A_rate_corrected_endo$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_A_rate_corrected_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_A_rate_corrected_endo$rope <- posterior_tab*100
Table_dsDNA_A_rate_corrected_endo$med <- median(Table_dsDNA_A_rate_corrected_endo$coef)
Table_dsDNA_A_rate_corrected_endo$pd <- pd(Table_dsDNA_A_rate_corrected_endo$coef2)

Table_dsDNA_A_rate_corrected_ecto <- as.data.frame(Table_dsDNA_A_rate_corrected$b_value3)
Table_dsDNA_A_rate_corrected_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_dsDNA_A_rate_corrected_ecto)<- c("coef","Parameter")
Table_dsDNA_A_rate_corrected_ecto$coef2 <- as.numeric(Table_dsDNA_A_rate_corrected_ecto$coef)
Table_dsDNA_A_rate_corrected_ecto$coef <- as.numeric(exp(Table_dsDNA_A_rate_corrected_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_A_rate_corrected_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_A_rate_corrected_ecto$rope <- posterior_tab*100
Table_dsDNA_A_rate_corrected_ecto$med <- median(Table_dsDNA_A_rate_corrected_ecto$coef)
Table_dsDNA_A_rate_corrected_ecto$pd <- pd(Table_dsDNA_A_rate_corrected_ecto$coef2)



Table_dsDNA_A_rate_corrected_ecto_endo<-rbind(Table_dsDNA_A_rate_corrected_endo,Table_dsDNA_A_rate_corrected_ecto)

#Create summary table
statistics_Table_dsDNA_A_dEVEs_Counts= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("dsDNA_EVEs_Counts",median(Table_dsDNA_A_rate_corrected_endo$coef),ci(Table_dsDNA_A_rate_corrected_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_A_rate_corrected_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_A_rate_corrected_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_dsDNA_A_rate_corrected_endo$pd[1]),
                             Ecto_parasitoid=c("dsDNA_EVEs_Counts",median(Table_dsDNA_A_rate_corrected_ecto$coef),ci(Table_dsDNA_A_rate_corrected_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_A_rate_corrected_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_A_rate_corrected_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_dsDNA_A_rate_corrected_ecto$pd[1])))



colnames(statistics_Table_dsDNA_A_dEVEs_Counts) <- statistics_Table_dsDNA_A_dEVEs_Counts[1,]
statistics_Table_dsDNA_A_dEVEs_Counts<-statistics_Table_dsDNA_A_dEVEs_Counts[-1,]


Table_dsDNA_A_rate_corrected_ecto_endo$Parameter2 <- paste0(Table_dsDNA_A_rate_corrected_ecto_endo$Parameter,"\n %pd = ",round(Table_dsDNA_A_rate_corrected_ecto_endo$pd*100,2) ,"% ")

stat.test_dsDNA_A_rate_corrected<-Table_dsDNA_A_rate_corrected_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_dsDNA_A_rate_corrected <- stat.test_dsDNA_A_rate_corrected %>% add_xy_position(x = "Parameter2")
stat.test_dsDNA_A_rate_corrected$Prop <- paste0(round((nrow(Table_dsDNA_A_rate_corrected[Table_dsDNA_A_rate_corrected$b_value2 > Table_dsDNA_A_rate_corrected$b_value3,]) / nrow(Table_dsDNA_A_rate_corrected))*100,digits = 2),"%")

boxplot_dsDNA_A_dEVEs_rate_corrected_Counts <- ggplot(Table_dsDNA_A_rate_corrected_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="dEVEs",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_dsDNA_A_rate_corrected, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 9, size=10
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,11))

## Combine All 2 plots

violin_plot_dsDNA_A_rate_corrected<-grid.arrange(arrangeGrob(
                         boxplot_dsDNA_A_dEVEs_rate_corrected_Counts  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank()) ,
                         boxplot_dsDNA_A_dEVEs_rate_corrected_Events  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank()),
                         nrow = 1,
                         ncol=2,
                         left = textGrob("Relative rate", rot = 90, vjust = 1)),nrow=1,top=textGrob("(relative number) dsDNA | A",
                                                                                                                 gp=gpar(fontsize=20,font=8)))



##
##################################
# ## dEVEs Events all posteriors ##
##################################


Table_dsDNA_rate_corrected<- read.table("/Users/bguinet/Desktop/output_dsDNA_rate_corrected_all/Posteriors_dEVEs_Events_ALL.txt",sep=";",h=T)


Table_dsDNA_rate_corrected$b_value2_coef_inf<- exp(as.numeric(Table_dsDNA_rate_corrected$b_value2))
Table_dsDNA_rate_corrected$b_value3_coef_inf<- exp(as.numeric(Table_dsDNA_rate_corrected$b_value3))
Table_dsDNA_rate_corrected<-Table_dsDNA_rate_corrected[!Table_dsDNA_rate_corrected$b_value2_coef_inf=="Inf",]
Table_dsDNA_rate_corrected<-Table_dsDNA_rate_corrected[!Table_dsDNA_rate_corrected$b_value3_coef_inf=="Inf",]


Table_dsDNA_rate_corrected_endo <- as.data.frame(Table_dsDNA_rate_corrected$b_value2)

Table_dsDNA_rate_corrected_endo$Parameter <- "Endo-parasitoid"
colnames(Table_dsDNA_rate_corrected_endo)<- c("coef","Parameter")
Table_dsDNA_rate_corrected_endo$coef2 <- as.numeric(Table_dsDNA_rate_corrected_endo$coef)
Table_dsDNA_rate_corrected_endo$coef <- as.numeric(exp(Table_dsDNA_rate_corrected_endo$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_rate_corrected_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_rate_corrected_endo$rope <- posterior_tab*100
Table_dsDNA_rate_corrected_endo$med <- median(Table_dsDNA_rate_corrected_endo$coef)
Table_dsDNA_rate_corrected_endo$pd <- pd(Table_dsDNA_rate_corrected_endo$coef2)

Table_dsDNA_rate_corrected_ecto <- as.data.frame(Table_dsDNA_rate_corrected$b_value3)
Table_dsDNA_rate_corrected_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_dsDNA_rate_corrected_ecto)<- c("coef","Parameter")
Table_dsDNA_rate_corrected_ecto$coef2 <- as.numeric(Table_dsDNA_rate_corrected_ecto$coef)
Table_dsDNA_rate_corrected_ecto$coef <- as.numeric(exp(Table_dsDNA_rate_corrected_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_rate_corrected_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_rate_corrected_ecto$rope <- posterior_tab*100
Table_dsDNA_rate_corrected_ecto$med <- median(Table_dsDNA_rate_corrected_ecto$coef)
Table_dsDNA_rate_corrected_ecto$pd <- pd(Table_dsDNA_rate_corrected_ecto$coef2)



Table_dsDNA_rate_corrected_ecto_endo<-rbind(Table_dsDNA_rate_corrected_endo,Table_dsDNA_rate_corrected_ecto)

#Create summary table
statistics_Table_dsDNA_A_dEVEs_EVents= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("dsDNA_EVEs_Events",median(Table_dsDNA_rate_corrected_endo$coef),ci(Table_dsDNA_rate_corrected_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_rate_corrected_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_rate_corrected_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_dsDNA_rate_corrected_endo$pd[1]),
                             Ecto_parasitoid=c("dsDNA_EVEs_Events",median(Table_dsDNA_rate_corrected_ecto$coef),ci(Table_dsDNA_rate_corrected_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_rate_corrected_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_rate_corrected_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_dsDNA_rate_corrected_ecto$pd[1])))



colnames(statistics_Table_dsDNA_A_dEVEs_EVents) <- statistics_Table_dsDNA_A_dEVEs_EVents[1,]
statistics_Table_dsDNA_A_dEVEs_EVents<-statistics_Table_dsDNA_A_dEVEs_EVents[-1,]


Table_dsDNA_rate_corrected_ecto_endo$Parameter2 <- paste0(Table_dsDNA_rate_corrected_ecto_endo$Parameter,"\n %pd = ",round(Table_dsDNA_rate_corrected_ecto_endo$pd*100,2) ,"% ")

stat.test_dsDNA_rate_corrected<-Table_dsDNA_rate_corrected_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_dsDNA_rate_corrected <- stat.test_dsDNA_rate_corrected %>% add_xy_position(x = "Parameter2")
stat.test_dsDNA_rate_corrected$Prop <- paste0(round((nrow(Table_dsDNA_rate_corrected[Table_dsDNA_rate_corrected$b_value2 > Table_dsDNA_rate_corrected$b_value3,]) / nrow(Table_dsDNA_rate_corrected))*100,digits = 2),"%")

boxplot_dsDNA_A_dEVEs_rate_corrected_Events <- ggplot(Table_dsDNA_rate_corrected_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="dEvents",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_dsDNA_rate_corrected, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 9, size=10
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,11))


##
##################################
# ## dEVEs Counts all posteriors ##
##################################

Table_dsDNA_rate_corrected<- read.table("/Users/bguinet/Desktop/output_dsDNA_rate_corrected_all/Posteriors_dEVEs_Counts_ALL.txt",sep=";",h=T)


Table_dsDNA_rate_corrected$b_value2_coef_inf<- exp(as.numeric(Table_dsDNA_rate_corrected$b_value2))
Table_dsDNA_rate_corrected$b_value3_coef_inf<- exp(as.numeric(Table_dsDNA_rate_corrected$b_value3))
Table_dsDNA_rate_corrected<-Table_dsDNA_rate_corrected[!Table_dsDNA_rate_corrected$b_value2_coef_inf=="Inf",]
Table_dsDNA_rate_corrected<-Table_dsDNA_rate_corrected[!Table_dsDNA_rate_corrected$b_value3_coef_inf=="Inf",]


Table_dsDNA_rate_corrected_endo <- as.data.frame(Table_dsDNA_rate_corrected$b_value2)

Table_dsDNA_rate_corrected_endo$Parameter <- "Endo-parasitoid"
colnames(Table_dsDNA_rate_corrected_endo)<- c("coef","Parameter")
Table_dsDNA_rate_corrected_endo$coef2 <- as.numeric(Table_dsDNA_rate_corrected_endo$coef)
Table_dsDNA_rate_corrected_endo$coef <- as.numeric(exp(Table_dsDNA_rate_corrected_endo$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_rate_corrected_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_rate_corrected_endo$rope <- posterior_tab*100
Table_dsDNA_rate_corrected_endo$med <- median(Table_dsDNA_rate_corrected_endo$coef)
Table_dsDNA_rate_corrected_endo$pd <- pd(Table_dsDNA_rate_corrected_endo$coef2)

Table_dsDNA_rate_corrected_ecto <- as.data.frame(Table_dsDNA_rate_corrected$b_value3)
Table_dsDNA_rate_corrected_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_dsDNA_rate_corrected_ecto)<- c("coef","Parameter")
Table_dsDNA_rate_corrected_ecto$coef2 <- as.numeric(Table_dsDNA_rate_corrected_ecto$coef)
Table_dsDNA_rate_corrected_ecto$coef <- as.numeric(exp(Table_dsDNA_rate_corrected_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_rate_corrected_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_rate_corrected_ecto$rope <- posterior_tab*100
Table_dsDNA_rate_corrected_ecto$med <- median(Table_dsDNA_rate_corrected_ecto$coef)
Table_dsDNA_rate_corrected_ecto$pd <- pd(Table_dsDNA_rate_corrected_ecto$coef2)



Table_dsDNA_rate_corrected_ecto_endo<-rbind(Table_dsDNA_rate_corrected_endo,Table_dsDNA_rate_corrected_ecto)

#Create summary table
statistics_Table_dsDNA_A_dEVEs_Counts= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("dsDNA_EVEs_Counts",median(Table_dsDNA_rate_corrected_endo$coef),ci(Table_dsDNA_rate_corrected_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_rate_corrected_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_rate_corrected_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_dsDNA_rate_corrected_endo$pd[1]),
                             Ecto_parasitoid=c("dsDNA_EVEs_Counts",median(Table_dsDNA_rate_corrected_ecto$coef),ci(Table_dsDNA_rate_corrected_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_rate_corrected_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_rate_corrected_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_dsDNA_rate_corrected_ecto$pd[1])))



colnames(statistics_Table_dsDNA_A_dEVEs_Counts) <- statistics_Table_dsDNA_A_dEVEs_Counts[1,]
statistics_Table_dsDNA_A_dEVEs_Counts<-statistics_Table_dsDNA_A_dEVEs_Counts[-1,]


Table_dsDNA_rate_corrected_ecto_endo$Parameter2 <- paste0(Table_dsDNA_rate_corrected_ecto_endo$Parameter,"\n %pd = ",round(Table_dsDNA_rate_corrected_ecto_endo$pd*100,2) ,"% ")

stat.test_dsDNA_rate_corrected<-Table_dsDNA_rate_corrected_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_dsDNA_rate_corrected <- stat.test_dsDNA_rate_corrected %>% add_xy_position(x = "Parameter2")
stat.test_dsDNA_rate_corrected$Prop <- paste0(round((nrow(Table_dsDNA_rate_corrected[Table_dsDNA_rate_corrected$b_value2 > Table_dsDNA_rate_corrected$b_value3,]) / nrow(Table_dsDNA_rate_corrected))*100,digits = 2),"%")

boxplot_dsDNA_A_dEVEs_rate_corrected_Counts <- ggplot(Table_dsDNA_rate_corrected_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="dEVEs",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_dsDNA_rate_corrected, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 9, size=10
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,11))

## Combine All 3 plots

violin_plot_dsDNA_rate_corrected<-grid.arrange(arrangeGrob(
                         boxplot_dsDNA_A_dEVEs_rate_corrected_Counts  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank()) ,
                         boxplot_dsDNA_A_dEVEs_rate_corrected_Events  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank()),
                         nrow = 1,
                         ncol=2,
                         left = textGrob("Relative rate", rot = 90, vjust = 1)),nrow=1,top=textGrob("(relative number) dsDNA",
                                                                                                                 gp=gpar(fontsize=20,font=8)))





#######

# corrected dEVEs rates dsDNA without controls
library(bayestestR)
library(ggstatsplot)
library(bayesplot)
library(ggpubr)
library(gridExtra)
library(grid)
library(ggprism)
library(rstatix)

##
##################################
# ## dEVEs Events all posteriors ##
##################################


Table_dsDNA_without_controls_rate_corrected<- read.table("/Users/bguinet/Desktop/output_dsDNA_without_controls_rate_corrected_all/Posteriors_dEVEs_Events_ALL.txt",sep=";",h=T)


Table_dsDNA_without_controls_rate_corrected$b_value2_coef_inf<- exp(as.numeric(Table_dsDNA_without_controls_rate_corrected$b_value2))
Table_dsDNA_without_controls_rate_corrected$b_value3_coef_inf<- exp(as.numeric(Table_dsDNA_without_controls_rate_corrected$b_value3))
Table_dsDNA_without_controls_rate_corrected<-Table_dsDNA_without_controls_rate_corrected[!Table_dsDNA_without_controls_rate_corrected$b_value2_coef_inf=="Inf",]
Table_dsDNA_without_controls_rate_corrected<-Table_dsDNA_without_controls_rate_corrected[!Table_dsDNA_without_controls_rate_corrected$b_value3_coef_inf=="Inf",]


Table_dsDNA_without_controls_rate_corrected_endo <- as.data.frame(Table_dsDNA_without_controls_rate_corrected$b_value2)

Table_dsDNA_without_controls_rate_corrected_endo$Parameter <- "Endo-parasitoid"
colnames(Table_dsDNA_without_controls_rate_corrected_endo)<- c("coef","Parameter")
Table_dsDNA_without_controls_rate_corrected_endo$coef2 <- as.numeric(Table_dsDNA_without_controls_rate_corrected_endo$coef)
Table_dsDNA_without_controls_rate_corrected_endo$coef <- as.numeric(exp(Table_dsDNA_without_controls_rate_corrected_endo$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_without_controls_rate_corrected_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_without_controls_rate_corrected_endo$rope <- posterior_tab*100
Table_dsDNA_without_controls_rate_corrected_endo$med <- median(Table_dsDNA_without_controls_rate_corrected_endo$coef)
Table_dsDNA_without_controls_rate_corrected_endo$pd <- pd(Table_dsDNA_without_controls_rate_corrected_endo$coef2)

Table_dsDNA_without_controls_rate_corrected_ecto <- as.data.frame(Table_dsDNA_without_controls_rate_corrected$b_value3)
Table_dsDNA_without_controls_rate_corrected_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_dsDNA_without_controls_rate_corrected_ecto)<- c("coef","Parameter")
Table_dsDNA_without_controls_rate_corrected_ecto$coef2 <- as.numeric(Table_dsDNA_without_controls_rate_corrected_ecto$coef)
Table_dsDNA_without_controls_rate_corrected_ecto$coef <- as.numeric(exp(Table_dsDNA_without_controls_rate_corrected_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_without_controls_rate_corrected_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_without_controls_rate_corrected_ecto$rope <- posterior_tab*100
Table_dsDNA_without_controls_rate_corrected_ecto$med <- median(Table_dsDNA_without_controls_rate_corrected_ecto$coef)
Table_dsDNA_without_controls_rate_corrected_ecto$pd <- pd(Table_dsDNA_without_controls_rate_corrected_ecto$coef2)



Table_dsDNA_without_controls_rate_corrected_ecto_endo<-rbind(Table_dsDNA_without_controls_rate_corrected_endo,Table_dsDNA_without_controls_rate_corrected_ecto)

#Create summary table
statistics_Table_dsDNA_A_dEVEs_EVents= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("dsDNA_EVEs_Events",median(Table_dsDNA_without_controls_rate_corrected_endo$coef),ci(Table_dsDNA_without_controls_rate_corrected_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_without_controls_rate_corrected_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_without_controls_rate_corrected_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_dsDNA_without_controls_rate_corrected_endo$pd[1]),
                             Ecto_parasitoid=c("dsDNA_EVEs_Events",median(Table_dsDNA_without_controls_rate_corrected_ecto$coef),ci(Table_dsDNA_without_controls_rate_corrected_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_without_controls_rate_corrected_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_without_controls_rate_corrected_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_dsDNA_without_controls_rate_corrected_ecto$pd[1])))



colnames(statistics_Table_dsDNA_A_dEVEs_EVents) <- statistics_Table_dsDNA_A_dEVEs_EVents[1,]
statistics_Table_dsDNA_A_dEVEs_EVents<-statistics_Table_dsDNA_A_dEVEs_EVents[-1,]


Table_dsDNA_without_controls_rate_corrected_ecto_endo$Parameter2 <- paste0(Table_dsDNA_without_controls_rate_corrected_ecto_endo$Parameter,"\n %pd = ",round(Table_dsDNA_without_controls_rate_corrected_ecto_endo$pd*100,2) ,"% ")

stat.test_dsDNA_without_controls_rate_corrected<-Table_dsDNA_without_controls_rate_corrected_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_dsDNA_without_controls_rate_corrected <- stat.test_dsDNA_without_controls_rate_corrected %>% add_xy_position(x = "Parameter2")
stat.test_dsDNA_without_controls_rate_corrected$Prop <- paste0(round((nrow(Table_dsDNA_without_controls_rate_corrected[Table_dsDNA_without_controls_rate_corrected$b_value2 > Table_dsDNA_without_controls_rate_corrected$b_value3,]) / nrow(Table_dsDNA_without_controls_rate_corrected))*100,digits = 2),"%")

boxplot_dsDNA_A_dEVEs_rate_corrected_Events <- ggplot(Table_dsDNA_without_controls_rate_corrected_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="dEvents",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_dsDNA_without_controls_rate_corrected, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 9, size=10
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,11))


##
##################################
# ## dEVEs Counts all posteriors ##
##################################

Table_dsDNA_without_controls_rate_corrected<- read.table("/Users/bguinet/Desktop/output_dsDNA_without_controls_rate_corrected_all/Posteriors_dEVEs_Counts_ALL.txt",sep=";",h=T)


Table_dsDNA_without_controls_rate_corrected$b_value2_coef_inf<- exp(as.numeric(Table_dsDNA_without_controls_rate_corrected$b_value2))
Table_dsDNA_without_controls_rate_corrected$b_value3_coef_inf<- exp(as.numeric(Table_dsDNA_without_controls_rate_corrected$b_value3))
Table_dsDNA_without_controls_rate_corrected<-Table_dsDNA_without_controls_rate_corrected[!Table_dsDNA_without_controls_rate_corrected$b_value2_coef_inf=="Inf",]
Table_dsDNA_without_controls_rate_corrected<-Table_dsDNA_without_controls_rate_corrected[!Table_dsDNA_without_controls_rate_corrected$b_value3_coef_inf=="Inf",]


Table_dsDNA_without_controls_rate_corrected_endo <- as.data.frame(Table_dsDNA_without_controls_rate_corrected$b_value2)

Table_dsDNA_without_controls_rate_corrected_endo$Parameter <- "Endo-parasitoid"
colnames(Table_dsDNA_without_controls_rate_corrected_endo)<- c("coef","Parameter")
Table_dsDNA_without_controls_rate_corrected_endo$coef2 <- as.numeric(Table_dsDNA_without_controls_rate_corrected_endo$coef)
Table_dsDNA_without_controls_rate_corrected_endo$coef <- as.numeric(exp(Table_dsDNA_without_controls_rate_corrected_endo$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_without_controls_rate_corrected_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_without_controls_rate_corrected_endo$rope <- posterior_tab*100
Table_dsDNA_without_controls_rate_corrected_endo$med <- median(Table_dsDNA_without_controls_rate_corrected_endo$coef)
Table_dsDNA_without_controls_rate_corrected_endo$pd <- pd(Table_dsDNA_without_controls_rate_corrected_endo$coef2)

Table_dsDNA_without_controls_rate_corrected_ecto <- as.data.frame(Table_dsDNA_without_controls_rate_corrected$b_value3)
Table_dsDNA_without_controls_rate_corrected_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_dsDNA_without_controls_rate_corrected_ecto)<- c("coef","Parameter")
Table_dsDNA_without_controls_rate_corrected_ecto$coef2 <- as.numeric(Table_dsDNA_without_controls_rate_corrected_ecto$coef)
Table_dsDNA_without_controls_rate_corrected_ecto$coef <- as.numeric(exp(Table_dsDNA_without_controls_rate_corrected_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_without_controls_rate_corrected_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_without_controls_rate_corrected_ecto$rope <- posterior_tab*100
Table_dsDNA_without_controls_rate_corrected_ecto$med <- median(Table_dsDNA_without_controls_rate_corrected_ecto$coef)
Table_dsDNA_without_controls_rate_corrected_ecto$pd <- pd(Table_dsDNA_without_controls_rate_corrected_ecto$coef2)



Table_dsDNA_without_controls_rate_corrected_ecto_endo<-rbind(Table_dsDNA_without_controls_rate_corrected_endo,Table_dsDNA_without_controls_rate_corrected_ecto)

#Create summary table
statistics_Table_dsDNA_A_dEVEs_Counts= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("dsDNA_EVEs_Counts",median(Table_dsDNA_without_controls_rate_corrected_endo$coef),ci(Table_dsDNA_without_controls_rate_corrected_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_without_controls_rate_corrected_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_without_controls_rate_corrected_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_dsDNA_without_controls_rate_corrected_endo$pd[1]),
                             Ecto_parasitoid=c("dsDNA_EVEs_Counts",median(Table_dsDNA_without_controls_rate_corrected_ecto$coef),ci(Table_dsDNA_without_controls_rate_corrected_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_without_controls_rate_corrected_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_without_controls_rate_corrected_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_dsDNA_without_controls_rate_corrected_ecto$pd[1])))



colnames(statistics_Table_dsDNA_A_dEVEs_Counts) <- statistics_Table_dsDNA_A_dEVEs_Counts[1,]
statistics_Table_dsDNA_A_dEVEs_Counts<-statistics_Table_dsDNA_A_dEVEs_Counts[-1,]


Table_dsDNA_without_controls_rate_corrected_ecto_endo$Parameter2 <- paste0(Table_dsDNA_without_controls_rate_corrected_ecto_endo$Parameter,"\n %pd = ",round(Table_dsDNA_without_controls_rate_corrected_ecto_endo$pd*100,2) ,"% ")

stat.test_dsDNA_without_controls_rate_corrected<-Table_dsDNA_without_controls_rate_corrected_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_dsDNA_without_controls_rate_corrected <- stat.test_dsDNA_without_controls_rate_corrected %>% add_xy_position(x = "Parameter2")
stat.test_dsDNA_without_controls_rate_corrected$Prop <- paste0(round((nrow(Table_dsDNA_without_controls_rate_corrected[Table_dsDNA_without_controls_rate_corrected$b_value2 > Table_dsDNA_without_controls_rate_corrected$b_value3,]) / nrow(Table_dsDNA_without_controls_rate_corrected))*100,digits = 2),"%")

boxplot_dsDNA_A_dEVEs_rate_corrected_Counts <- ggplot(Table_dsDNA_without_controls_rate_corrected_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="dEVEs",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_dsDNA_without_controls_rate_corrected, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 9, size=10
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,11))

## Combine All 2 plots

violin_plot_dsDNA_without_controls_rate_corrected<-grid.arrange(arrangeGrob(
                         boxplot_dsDNA_A_dEVEs_rate_corrected_Counts  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank()) ,
                         boxplot_dsDNA_A_dEVEs_rate_corrected_Events  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank()),
                         nrow = 1,
                         ncol=2,
                         left = textGrob("Relative rate", rot = 90, vjust = 1)),nrow=1,top=textGrob("(relative number) dsDNA | without controls ",
                                                                                                                 gp=gpar(fontsize=20,font=8)))





#######

# corrected dEVEs rates dsDNA without controls A
library(bayestestR)
library(ggstatsplot)
library(bayesplot)
library(ggpubr)
library(gridExtra)
library(grid)
library(ggprism)
library(rstatix)

##
##################################
# ## dEVEs Events all posteriors ##
##################################


Table_dsDNA_A_without_controls_rate_corrected<- read.table("/Users/bguinet/Desktop/output_dsDNA_A_without_controls_rate_corrected_all/Posteriors_dEVEs_Events_ALL.txt",sep=";",h=T)


Table_dsDNA_A_without_controls_rate_corrected$b_value2_coef_inf<- exp(as.numeric(Table_dsDNA_A_without_controls_rate_corrected$b_value2))
Table_dsDNA_A_without_controls_rate_corrected$b_value3_coef_inf<- exp(as.numeric(Table_dsDNA_A_without_controls_rate_corrected$b_value3))
Table_dsDNA_A_without_controls_rate_corrected<-Table_dsDNA_A_without_controls_rate_corrected[!Table_dsDNA_A_without_controls_rate_corrected$b_value2_coef_inf=="Inf",]
Table_dsDNA_A_without_controls_rate_corrected<-Table_dsDNA_A_without_controls_rate_corrected[!Table_dsDNA_A_without_controls_rate_corrected$b_value3_coef_inf=="Inf",]


Table_dsDNA_A_without_controls_rate_corrected_endo <- as.data.frame(Table_dsDNA_A_without_controls_rate_corrected$b_value2)

Table_dsDNA_A_without_controls_rate_corrected_endo$Parameter <- "Endo-parasitoid"
colnames(Table_dsDNA_A_without_controls_rate_corrected_endo)<- c("coef","Parameter")
Table_dsDNA_A_without_controls_rate_corrected_endo$coef2 <- as.numeric(Table_dsDNA_A_without_controls_rate_corrected_endo$coef)
Table_dsDNA_A_without_controls_rate_corrected_endo$coef <- as.numeric(exp(Table_dsDNA_A_without_controls_rate_corrected_endo$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_A_without_controls_rate_corrected_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_A_without_controls_rate_corrected_endo$rope <- posterior_tab*100
Table_dsDNA_A_without_controls_rate_corrected_endo$med <- median(Table_dsDNA_A_without_controls_rate_corrected_endo$coef)
Table_dsDNA_A_without_controls_rate_corrected_endo$pd <- pd(Table_dsDNA_A_without_controls_rate_corrected_endo$coef2)

Table_dsDNA_A_without_controls_rate_corrected_ecto <- as.data.frame(Table_dsDNA_A_without_controls_rate_corrected$b_value3)
Table_dsDNA_A_without_controls_rate_corrected_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_dsDNA_A_without_controls_rate_corrected_ecto)<- c("coef","Parameter")
Table_dsDNA_A_without_controls_rate_corrected_ecto$coef2 <- as.numeric(Table_dsDNA_A_without_controls_rate_corrected_ecto$coef)
Table_dsDNA_A_without_controls_rate_corrected_ecto$coef <- as.numeric(exp(Table_dsDNA_A_without_controls_rate_corrected_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_A_without_controls_rate_corrected_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_A_without_controls_rate_corrected_ecto$rope <- posterior_tab*100
Table_dsDNA_A_without_controls_rate_corrected_ecto$med <- median(Table_dsDNA_A_without_controls_rate_corrected_ecto$coef)
Table_dsDNA_A_without_controls_rate_corrected_ecto$pd <- pd(Table_dsDNA_A_without_controls_rate_corrected_ecto$coef2)



Table_dsDNA_A_without_controls_rate_corrected_ecto_endo<-rbind(Table_dsDNA_A_without_controls_rate_corrected_endo,Table_dsDNA_A_without_controls_rate_corrected_ecto)

#Create summary table
statistics_Table_dsDNA_A_dEVEs_EVents= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("dsDNA_EVEs_Events",median(Table_dsDNA_A_without_controls_rate_corrected_endo$coef),ci(Table_dsDNA_A_without_controls_rate_corrected_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_A_without_controls_rate_corrected_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_A_without_controls_rate_corrected_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_dsDNA_A_without_controls_rate_corrected_endo$pd[1]),
                             Ecto_parasitoid=c("dsDNA_EVEs_Events",median(Table_dsDNA_A_without_controls_rate_corrected_ecto$coef),ci(Table_dsDNA_A_without_controls_rate_corrected_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_A_without_controls_rate_corrected_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_A_without_controls_rate_corrected_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_dsDNA_A_without_controls_rate_corrected_ecto$pd[1])))



colnames(statistics_Table_dsDNA_A_dEVEs_EVents) <- statistics_Table_dsDNA_A_dEVEs_EVents[1,]
statistics_Table_dsDNA_A_dEVEs_EVents<-statistics_Table_dsDNA_A_dEVEs_EVents[-1,]


Table_dsDNA_A_without_controls_rate_corrected_ecto_endo$Parameter2 <- paste0(Table_dsDNA_A_without_controls_rate_corrected_ecto_endo$Parameter,"\n %pd = ",round(Table_dsDNA_A_without_controls_rate_corrected_ecto_endo$pd*100,2) ,"% ")

stat.test_dsDNA_A_without_controls_rate_corrected<-Table_dsDNA_A_without_controls_rate_corrected_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_dsDNA_A_without_controls_rate_corrected <- stat.test_dsDNA_A_without_controls_rate_corrected %>% add_xy_position(x = "Parameter2")
stat.test_dsDNA_A_without_controls_rate_corrected$Prop <- paste0(round((nrow(Table_dsDNA_A_without_controls_rate_corrected[Table_dsDNA_A_without_controls_rate_corrected$b_value2 > Table_dsDNA_A_without_controls_rate_corrected$b_value3,]) / nrow(Table_dsDNA_A_without_controls_rate_corrected))*100,digits = 2),"%")

boxplot_dsDNA_A_dEVEs_rate_corrected_Events <- ggplot(Table_dsDNA_A_without_controls_rate_corrected_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="dEvents",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_dsDNA_A_without_controls_rate_corrected, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 9, size=10
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,11))


##
##################################
# ## dEVEs Counts all posteriors ##
##################################

Table_dsDNA_A_without_controls_rate_corrected<- read.table("/Users/bguinet/Desktop/output_dsDNA_A_without_controls_rate_corrected_all/Posteriors_dEVEs_Counts_ALL.txt",sep=";",h=T)


Table_dsDNA_A_without_controls_rate_corrected$b_value2_coef_inf<- exp(as.numeric(Table_dsDNA_A_without_controls_rate_corrected$b_value2))
Table_dsDNA_A_without_controls_rate_corrected$b_value3_coef_inf<- exp(as.numeric(Table_dsDNA_A_without_controls_rate_corrected$b_value3))
Table_dsDNA_A_without_controls_rate_corrected<-Table_dsDNA_A_without_controls_rate_corrected[!Table_dsDNA_A_without_controls_rate_corrected$b_value2_coef_inf=="Inf",]
Table_dsDNA_A_without_controls_rate_corrected<-Table_dsDNA_A_without_controls_rate_corrected[!Table_dsDNA_A_without_controls_rate_corrected$b_value3_coef_inf=="Inf",]


Table_dsDNA_A_without_controls_rate_corrected_endo <- as.data.frame(Table_dsDNA_A_without_controls_rate_corrected$b_value2)

Table_dsDNA_A_without_controls_rate_corrected_endo$Parameter <- "Endo-parasitoid"
colnames(Table_dsDNA_A_without_controls_rate_corrected_endo)<- c("coef","Parameter")
Table_dsDNA_A_without_controls_rate_corrected_endo$coef2 <- as.numeric(Table_dsDNA_A_without_controls_rate_corrected_endo$coef)
Table_dsDNA_A_without_controls_rate_corrected_endo$coef <- as.numeric(exp(Table_dsDNA_A_without_controls_rate_corrected_endo$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_A_without_controls_rate_corrected_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_A_without_controls_rate_corrected_endo$rope <- posterior_tab*100
Table_dsDNA_A_without_controls_rate_corrected_endo$med <- median(Table_dsDNA_A_without_controls_rate_corrected_endo$coef)
Table_dsDNA_A_without_controls_rate_corrected_endo$pd <- pd(Table_dsDNA_A_without_controls_rate_corrected_endo$coef2)

Table_dsDNA_A_without_controls_rate_corrected_ecto <- as.data.frame(Table_dsDNA_A_without_controls_rate_corrected$b_value3)
Table_dsDNA_A_without_controls_rate_corrected_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_dsDNA_A_without_controls_rate_corrected_ecto)<- c("coef","Parameter")
Table_dsDNA_A_without_controls_rate_corrected_ecto$coef2 <- as.numeric(Table_dsDNA_A_without_controls_rate_corrected_ecto$coef)
Table_dsDNA_A_without_controls_rate_corrected_ecto$coef <- as.numeric(exp(Table_dsDNA_A_without_controls_rate_corrected_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_A_without_controls_rate_corrected_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_A_without_controls_rate_corrected_ecto$rope <- posterior_tab*100
Table_dsDNA_A_without_controls_rate_corrected_ecto$med <- median(Table_dsDNA_A_without_controls_rate_corrected_ecto$coef)
Table_dsDNA_A_without_controls_rate_corrected_ecto$pd <- pd(Table_dsDNA_A_without_controls_rate_corrected_ecto$coef2)



Table_dsDNA_A_without_controls_rate_corrected_ecto_endo<-rbind(Table_dsDNA_A_without_controls_rate_corrected_endo,Table_dsDNA_A_without_controls_rate_corrected_ecto)

#Create summary table
statistics_Table_dsDNA_A_dEVEs_Counts= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("dsDNA_EVEs_Counts",median(Table_dsDNA_A_without_controls_rate_corrected_endo$coef),ci(Table_dsDNA_A_without_controls_rate_corrected_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_A_without_controls_rate_corrected_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_A_without_controls_rate_corrected_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_dsDNA_A_without_controls_rate_corrected_endo$pd[1]),
                             Ecto_parasitoid=c("dsDNA_EVEs_Counts",median(Table_dsDNA_A_without_controls_rate_corrected_ecto$coef),ci(Table_dsDNA_A_without_controls_rate_corrected_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_A_without_controls_rate_corrected_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_A_without_controls_rate_corrected_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_dsDNA_A_without_controls_rate_corrected_ecto$pd[1])))



colnames(statistics_Table_dsDNA_A_dEVEs_Counts) <- statistics_Table_dsDNA_A_dEVEs_Counts[1,]
statistics_Table_dsDNA_A_dEVEs_Counts<-statistics_Table_dsDNA_A_dEVEs_Counts[-1,]


Table_dsDNA_A_without_controls_rate_corrected_ecto_endo$Parameter2 <- paste0(Table_dsDNA_A_without_controls_rate_corrected_ecto_endo$Parameter,"\n %pd = ",round(Table_dsDNA_A_without_controls_rate_corrected_ecto_endo$pd*100,2) ,"% ")

stat.test_dsDNA_A_without_controls_rate_corrected<-Table_dsDNA_A_without_controls_rate_corrected_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_dsDNA_A_without_controls_rate_corrected <- stat.test_dsDNA_A_without_controls_rate_corrected %>% add_xy_position(x = "Parameter2")
stat.test_dsDNA_A_without_controls_rate_corrected$Prop <- paste0(round((nrow(Table_dsDNA_A_without_controls_rate_corrected[Table_dsDNA_A_without_controls_rate_corrected$b_value2 > Table_dsDNA_A_without_controls_rate_corrected$b_value3,]) / nrow(Table_dsDNA_A_without_controls_rate_corrected))*100,digits = 2),"%")

boxplot_dsDNA_A_dEVEs_rate_corrected_Counts <- ggplot(Table_dsDNA_A_without_controls_rate_corrected_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="dEVEs",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_dsDNA_A_without_controls_rate_corrected, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 9, size=10
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,11))

## Combine All 2 plots

violin_plot_dsDNA_A_without_controls_rate_corrected<-grid.arrange(arrangeGrob(
                         boxplot_dsDNA_A_dEVEs_rate_corrected_Counts  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank()) ,
                         boxplot_dsDNA_A_dEVEs_rate_corrected_Events  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank()),
                         nrow = 1,
                         ncol=2,
                         left = textGrob("Relative rate", rot = 90, vjust = 1)),nrow=1,top=textGrob("(relative number) dsDNA | without controls | A ",
                                                                                                                 gp=gpar(fontsize=20,font=8)))



#### Gather all

library(cowplot)
allviolin_relative<-plot_grid(violin_plot_dsDNA_rate_corrected,
                     violin_plot_dsDNA_without_controls_rate_corrected,
                     violin_plot_dsDNA_A_rate_corrected,
                     violin_plot_dsDNA_A_without_controls_rate_corrected,
                     ncol=2,nrow=2,labels=LETTERS[1:4],label_size = 16)






#20*17


###############################################################################


#ALL A-B-C-D with controls

library(bayestestR)
library(ggstatsplot)
library(bayesplot)
library(ggpubr)
library(gridExtra)
library(grid)
library(ggprism)
library(rstatix)


##
##################################
# ## EVEs Events all posteriors ##
##################################
Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_dsDNA_ssDNA_dsRNA_ssRNA_All/Posteriors_EVEs_Events_ALL.txt",sep=";",h=T)

Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL$b_value2_coef_inf<- exp(as.numeric(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL$b_value2))
Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL$b_value3_coef_inf<- exp(as.numeric(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL$b_value3))
Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL<-Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL[!Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL$b_value2_coef_inf=="Inf",]
Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL<-Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL[!Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL$b_value3_coef_inf=="Inf",]


Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_endo <- as.data.frame(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL$b_value2)
Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_endo)<- c("coef","Parameter")
Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_endo$coef2 <- as.numeric(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_endo$coef)
Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_endo$coef <- as.numeric(exp(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_endo$rope <- posterior_tab*100
Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_endo$med <- median(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_endo$coef)
Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_endo$pd <- pd(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_endo$coef2)

Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_ecto <- as.data.frame(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL$b_value3)
Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_ecto)<- c("coef","Parameter")
Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_ecto$coef2 <- as.numeric(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_ecto$coef)
Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_ecto$coef <- as.numeric(exp(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_ecto$rope <- posterior_tab*100
Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_ecto$med <- median(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_ecto$coef)
Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_ecto$pd <- pd(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_ecto$coef2)



Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_ecto_endo<-rbind(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_endo,Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_ecto)

#Create summary table
statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("dsDNA_EVEs_Events",median(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_endo$coef),ci(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("dsDNA_EVEs_Events",median(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_ecto$coef),ci(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_ecto$pd[1])))



colnames(statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events) <- statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events[1,]
statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events<-statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events[-1,]


Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_ecto_endo$Parameter2 <- paste0(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL<-Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL <- stat.test_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL %>% add_xy_position(x = "Parameter2")
stat.test_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL$Prop <- paste0(round((nrow(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL[Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL$b_value2 > Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL$b_value3,]) / nrow(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL))*100,digits = 2),"%")

boxplot_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events <- ggplot(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="Events",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,15))


##################################
## dEVEs Events all posteriors
##################################

Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_dsDNA_ssDNA_dsRNA_ssRNA_All/Posteriors_dEVEs_Events_ALL.txt",sep=";",h=T)

Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL$b_value2_coef_inf<- exp(as.numeric(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL$b_value2))
Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL$b_value3_coef_inf<- exp(as.numeric(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL$b_value3))
Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL<-Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL[!Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL$b_value2_coef_inf=="Inf",]
Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL<-Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL[!Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL$b_value3_coef_inf=="Inf",]


Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_endo <- as.data.frame(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL$b_value2)
Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_endo)<- c("coef","Parameter")
Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_endo$coef2 <- as.numeric(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_endo$coef)
Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_endo$coef <- as.numeric(exp(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_endo$rope <- posterior_tab*100
Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_endo$med <- median(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_endo$coef)
Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_endo$pd <- pd(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_endo$coef2)

Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_ecto <- as.data.frame(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL$b_value3)
Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_ecto)<- c("coef","Parameter")
Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_ecto$coef2 <- as.numeric(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_ecto$coef)
Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_ecto$coef <- as.numeric(exp(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_ecto$rope <- posterior_tab*100
Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_ecto$med <- median(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_ecto$coef)
Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_ecto$pd <- pd(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_ecto$coef2)



Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_ecto_endo<-rbind(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_endo,Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_ecto)

#Create summary table
statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("dsDNA_dEVEs_Events",median(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_endo$coef),ci(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("dsDNA_dEVEs_Events",median(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_ecto$coef),ci(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_ecto$pd[1])))



colnames(statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events) <- statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events[1,]
statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events<-statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events[-1,]


Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_ecto_endo$Parameter2 <- paste0(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL<-Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL <- stat.test_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL %>% add_xy_position(x = "Parameter2")
stat.test_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL$Prop <- paste0(round((nrow(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL[Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL$b_value2 > Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL$b_value3,]) / nrow(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL))*100,digits = 2),"%")

boxplot_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events <- ggplot(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter),trim = T,adjust = 2)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="dEvents",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,15))

##################################
## EVEs Counts all posteriors####
##################################


Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_dsDNA_ssDNA_dsRNA_ssRNA_All/Posteriors_EVEs_Counts_ALL.txt",sep=";",h=T)

Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL$b_value2_coef_inf<- exp(as.numeric(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL$b_value2))
Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL$b_value3_coef_inf<- exp(as.numeric(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL$b_value3))
Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL<-Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL[!Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL$b_value2_coef_inf=="Inf",]
Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL<-Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL[!Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL$b_value3_coef_inf=="Inf",]


Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_endo <- as.data.frame(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL$b_value2)
Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_endo)<- c("coef","Parameter")
Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_endo$coef2 <- as.numeric(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_endo$coef)
Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_endo$coef <- as.numeric(exp(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_endo$rope <- posterior_tab*100
Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_endo$med <- median(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_endo$coef)
Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_endo$pd <- pd(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_endo$coef2)

Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_ecto <- as.data.frame(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL$b_value3)
Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_ecto)<- c("coef","Parameter")
Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_ecto$coef2 <- as.numeric(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_ecto$coef)
Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_ecto$coef <- as.numeric(exp(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_ecto$rope <- posterior_tab*100
Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_ecto$med <- median(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_ecto$coef)
Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_ecto$pd <- pd(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_ecto$coef2)



Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_ecto_endo<-rbind(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_endo,Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_ecto)

#Create summary table
statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("dsDNA_EVEs_Counts",median(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_endo$coef),ci(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("dsDNA_EVEs_Counts",median(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_ecto$coef),ci(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_ecto$pd[1])))



colnames(statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts) <- statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts[1,]
statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts<-statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts[-1,]


Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_ecto_endo$Parameter2 <- paste0(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL<-Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL <- stat.test_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL %>% add_xy_position(x = "Parameter2")
stat.test_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL$Prop <- paste0(round((nrow(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL[Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL$b_value2 > Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL$b_value3,]) / nrow(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL))*100,digits = 2),"%")

boxplot_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts <- ggplot(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="EVEs",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,15))


################################
## dEVEs Counts all posteriors #
################################
Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_dsDNA_ssDNA_dsRNA_ssRNA_All/Posteriors_dEVEs_Counts_ALL.txt",sep=";",h=T)

Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL$b_value2_coef_inf<- exp(as.numeric(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL$b_value2))
Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL$b_value3_coef_inf<- exp(as.numeric(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL$b_value3))
Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL<-Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL[!Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL$b_value2_coef_inf=="Inf",]
Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL<-Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL[!Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL$b_value3_coef_inf=="Inf",]


Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_endo <- as.data.frame(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL$b_value2)
Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_endo)<- c("coef","Parameter")
Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_endo$coef2 <- as.numeric(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_endo$coef)
Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_endo$coef <- as.numeric(exp(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_endo$rope <- posterior_tab*100
Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_endo$med <- median(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_endo$coef)
Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_endo$pd <- pd(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_endo$coef2)

Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_ecto <- as.data.frame(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL$b_value3)
Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_ecto)<- c("coef","Parameter")
Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_ecto$coef2 <- as.numeric(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_ecto$coef)
Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_ecto$coef <- as.numeric(exp(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_ecto$rope <- posterior_tab*100
Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_ecto$med <- median(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_ecto$coef)
Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_ecto$pd <- pd(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_ecto$coef2)



Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_ecto_endo<-rbind(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_endo,Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_ecto)

#Create summary table
statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("dsDNA_dEVEs_Counts",median(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_endo$coef),ci(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("dsDNA_dEVEs_Counts",median(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_ecto$coef),ci(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_ecto$pd[1])))



colnames(statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts) <- statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts[1,]
statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts<-statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts[-1,]


Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_ecto_endo$Parameter2 <- paste0(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL<-Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL <- stat.test_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL %>% add_xy_position(x = "Parameter2")
stat.test_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL$Prop <- paste0(round((nrow(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL[Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL$b_value2 > Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL$b_value3,]) / nrow(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL))*100,digits = 2),"%")

boxplot_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts <- ggplot(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="dEVEs",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,15))



## Combine All 4 plots

violin_plot_dsDNA_ssDNA_dsRNA_ssRNA<-grid.arrange(arrangeGrob(
                         boxplot_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank()) ,
                         boxplot_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank()),
                         boxplot_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank()) ,
                         boxplot_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank()),
                         nrow = 1,
                         ncol=4,
                         left = textGrob("Relative rate", rot = 90, vjust = 1)),nrow=1,top=textGrob("All viruses | A-B-C-D",
                                                                                                                 gp=gpar(fontsize=20,font=8)))



#Add posterior coefficient comparison between endo and ecto
statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events<-as.data.frame(statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events)
statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events$Endo_sup_Ecto<-round((nrow(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL[Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL$b_value2 > Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL$b_value3,]) / nrow(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events_ALL))*100,digits = 2)

statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events<-as.data.frame(statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events)
statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events$Endo_sup_Ecto<-round((nrow(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL[Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL$b_value2 > Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL$b_value3,]) / nrow(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events_ALL))*100,digits = 2)

statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts<-as.data.frame(statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts)
statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts$Endo_sup_Ecto<-round((nrow(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL[Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL$b_value2 > Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL$b_value3,]) / nrow(Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts_ALL))*100,digits = 2)


statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts<-as.data.frame(statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts)
statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts$Endo_sup_Ecto<-round((nrow(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL[Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL$b_value2 > Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL$b_value3,]) / nrow(Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts_ALL))*100,digits = 2)


## Combine all statistics tables

statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA<-as.data.frame(rbind(statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events,statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events,statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts,statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts))

statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA$Lifestyle<- gsub("\\..*","",rownames(statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA))

library(dplyr)
statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA$Median<-as.numeric(statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA$Median)
statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA$CI_low<-as.numeric(statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA$CI_low)
statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA$CI_high<-as.numeric(statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA$CI_high)
statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA$ROPE_Percentage<-as.numeric(statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA$ROPE_Percentage)*100
statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA<-statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA %>% mutate_if(is.numeric, round, digits=3)

statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA$Type <-c("EVEs_Events","EVEs_Events","dEVEs_Events","dEVEs_Events","EVEs_Numbers","EVEs_Numbers","dEVEs_Numbers","dEVEs_Numbers")

table_dsDNA_ssDNA_dsRNA_ssRNA<-statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA
table_dsDNA_ssDNA_dsRNA_ssRNA[nrow(table_dsDNA_ssDNA_dsRNA_ssRNA)+1,] <- NA






#ALL A with controls

library(bayestestR)
library(ggstatsplot)
library(bayesplot)
library(ggpubr)
library(gridExtra)
library(grid)



##
##################################
# ## EVEs Events all posteriors ##
##################################
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results2/output_dsDNA_ssDNA_dsRNA_ssRNA_A_All/Posteriors_EVEs_Events_ALL.txt",sep=";",h=T)

Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL$b_value2_coef_inf<- exp(as.numeric(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL$b_value2))
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL$b_value3_coef_inf<- exp(as.numeric(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL$b_value3))
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL<-Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL[!Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL$b_value2_coef_inf=="Inf",]
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL<-Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL[!Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL$b_value3_coef_inf=="Inf",]


Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_endo <- as.data.frame(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL$b_value2)
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_endo)<- c("coef","Parameter")
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_endo$coef2 <- as.numeric(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_endo$coef)
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_endo$coef <- as.numeric(exp(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_endo$rope <- posterior_tab*100
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_endo$med <- median(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_endo$coef)
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_endo$pd <- pd(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_endo$coef2)

Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_ecto <- as.data.frame(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL$b_value3)
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_ecto)<- c("coef","Parameter")
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_ecto$coef2 <- as.numeric(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_ecto$coef)
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_ecto$coef <- as.numeric(exp(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_ecto$rope <- posterior_tab*100
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_ecto$med <- median(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_ecto$coef)
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_ecto$pd <- pd(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_ecto$coef2)



Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_ecto_endo<-rbind(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_endo,Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_ecto)

#Create summary table
statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("dsDNA_EVEs_Events",median(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_endo$coef),ci(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("dsDNA_EVEs_Events",median(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_ecto$coef),ci(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_ecto$pd[1])))



colnames(statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events) <- statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events[1,]
statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events<-statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events[-1,]


Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_ecto_endo$Parameter2 <- paste0(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL<-Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL <- stat.test_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL %>% add_xy_position(x = "Parameter2")
stat.test_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL$Prop <- paste0(round((nrow(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL[Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL$b_value2 > Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL$b_value3,]) / nrow(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL))*100,digits = 2),"%")

boxplot_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events <- ggplot(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="Events",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,15))


##################################
## dEVEs Events all posteriors
##################################

Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results2/output_dsDNA_ssDNA_dsRNA_ssRNA_A_All/Posteriors_dEVEs_Events_ALL.txt",sep=";",h=T)

Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL$b_value2_coef_inf<- exp(as.numeric(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL$b_value2))
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL$b_value3_coef_inf<- exp(as.numeric(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL$b_value3))
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL<-Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL[!Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL$b_value2_coef_inf=="Inf",]
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL<-Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL[!Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL$b_value3_coef_inf=="Inf",]


Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_endo <- as.data.frame(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL$b_value2)
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_endo)<- c("coef","Parameter")
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_endo$coef2 <- as.numeric(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_endo$coef)
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_endo$coef <- as.numeric(exp(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_endo$rope <- posterior_tab*100
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_endo$med <- median(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_endo$coef)
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_endo$pd <- pd(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_endo$coef2)

Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_ecto <- as.data.frame(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL$b_value3)
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_ecto)<- c("coef","Parameter")
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_ecto$coef2 <- as.numeric(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_ecto$coef)
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_ecto$coef <- as.numeric(exp(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_ecto$rope <- posterior_tab*100
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_ecto$med <- median(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_ecto$coef)
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_ecto$pd <- pd(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_ecto$coef2)



Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_ecto_endo<-rbind(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_endo,Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_ecto)

#Create summary table
statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("dsDNA_dEVEs_Events",median(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_endo$coef),ci(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("dsDNA_dEVEs_Events",median(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_ecto$coef),ci(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_ecto$pd[1])))



colnames(statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events) <- statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events[1,]
statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events<-statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events[-1,]


Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_ecto_endo$Parameter2 <- paste0(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL<-Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL <- stat.test_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL %>% add_xy_position(x = "Parameter2")
stat.test_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL$Prop <- paste0(round((nrow(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL[Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL$b_value2 > Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL$b_value3,]) / nrow(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL))*100,digits = 2),"%")

boxplot_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events <- ggplot(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter),trim = T,adjust = 2)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="dEvents",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,15))

##################################
## EVEs Counts all posteriors####
##################################


Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results2/output_dsDNA_ssDNA_dsRNA_ssRNA_A_All/Posteriors_EVEs_Counts_ALL.txt",sep=";",h=T)

Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL$b_value2_coef_inf<- exp(as.numeric(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL$b_value2))
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL$b_value3_coef_inf<- exp(as.numeric(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL$b_value3))
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL<-Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL[!Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL$b_value2_coef_inf=="Inf",]
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL<-Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL[!Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL$b_value3_coef_inf=="Inf",]


Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_endo <- as.data.frame(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL$b_value2)
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_endo)<- c("coef","Parameter")
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_endo$coef2 <- as.numeric(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_endo$coef)
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_endo$coef <- as.numeric(exp(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_endo$rope <- posterior_tab*100
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_endo$med <- median(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_endo$coef)
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_endo$pd <- pd(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_endo$coef2)

Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_ecto <- as.data.frame(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL$b_value3)
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_ecto)<- c("coef","Parameter")
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_ecto$coef2 <- as.numeric(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_ecto$coef)
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_ecto$coef <- as.numeric(exp(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_ecto$rope <- posterior_tab*100
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_ecto$med <- median(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_ecto$coef)
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_ecto$pd <- pd(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_ecto$coef2)



Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_ecto_endo<-rbind(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_endo,Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_ecto)

#Create summary table
statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("dsDNA_EVEs_Counts",median(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_endo$coef),ci(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("dsDNA_EVEs_Counts",median(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_ecto$coef),ci(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_ecto$pd[1])))



colnames(statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts) <- statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts[1,]
statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts<-statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts[-1,]


Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_ecto_endo$Parameter2 <- paste0(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL<-Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL <- stat.test_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL %>% add_xy_position(x = "Parameter2")
stat.test_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL$Prop <- paste0(round((nrow(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL[Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL$b_value2 > Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL$b_value3,]) / nrow(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL))*100,digits = 2),"%")

boxplot_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts <- ggplot(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="EVEs",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,15))


################################
## dEVEs Counts all posteriors #
################################
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results2/output_dsDNA_ssDNA_dsRNA_ssRNA_A_All/Posteriors_dEVEs_Counts_ALL.txt",sep=";",h=T)

Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL$b_value2_coef_inf<- exp(as.numeric(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL$b_value2))
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL$b_value3_coef_inf<- exp(as.numeric(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL$b_value3))
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL<-Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL[!Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL$b_value2_coef_inf=="Inf",]
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL<-Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL[!Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL$b_value3_coef_inf=="Inf",]


Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_endo <- as.data.frame(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL$b_value2)
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_endo)<- c("coef","Parameter")
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_endo$coef2 <- as.numeric(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_endo$coef)
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_endo$coef <- as.numeric(exp(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_endo$rope <- posterior_tab*100
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_endo$med <- median(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_endo$coef)
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_endo$pd <- pd(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_endo$coef2)

Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_ecto <- as.data.frame(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL$b_value3)
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_ecto)<- c("coef","Parameter")
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_ecto$coef2 <- as.numeric(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_ecto$coef)
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_ecto$coef <- as.numeric(exp(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_ecto$rope <- posterior_tab*100
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_ecto$med <- median(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_ecto$coef)
Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_ecto$pd <- pd(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_ecto$coef2)



Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_ecto_endo<-rbind(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_endo,Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_ecto)

#Create summary table
statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("dsDNA_dEVEs_Counts",median(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_endo$coef),ci(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("dsDNA_dEVEs_Counts",median(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_ecto$coef),ci(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_ecto$pd[1])))



colnames(statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts) <- statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts[1,]
statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts<-statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts[-1,]


Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_ecto_endo$Parameter2 <- paste0(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL<-Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL <- stat.test_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL %>% add_xy_position(x = "Parameter2")
stat.test_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL$Prop <- paste0(round((nrow(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL[Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL$b_value2 > Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL$b_value3,]) / nrow(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL))*100,digits = 2),"%")

boxplot_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts <- ggplot(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="dEVEs",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,15))


## Combine All 4 plots

violin_plot_dsDNA_ssDNA_dsRNA_ssRNA_A<-grid.arrange(arrangeGrob(
                         boxplot_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank()) ,
                         boxplot_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank()),
                         boxplot_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank()) ,
                         boxplot_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank()),
                         nrow = 1,
                         ncol=4,
                         left = textGrob("Relative rate", rot = 90, vjust = 1)),nrow=1,top=textGrob("All viruses | A",
                                                                                                                 gp=gpar(fontsize=20,font=8)))


#Add posterior coefficient comparison between endo and ecto
statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events<-as.data.frame(statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events)
statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events$Endo_sup_Ecto<-round((nrow(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL[Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL$b_value2 > Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL$b_value3,]) / nrow(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events_ALL))*100,digits = 2)

statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events<-as.data.frame(statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events)
statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events$Endo_sup_Ecto<-round((nrow(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL[Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL$b_value2 > Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL$b_value3,]) / nrow(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events_ALL))*100,digits = 2)

statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts<-as.data.frame(statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts)
statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts$Endo_sup_Ecto<-round((nrow(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL[Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL$b_value2 > Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL$b_value3,]) / nrow(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts_ALL))*100,digits = 2)


statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts<-as.data.frame(statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts)
statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts$Endo_sup_Ecto<-round((nrow(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL[Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL$b_value2 > Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL$b_value3,]) / nrow(Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts_ALL))*100,digits = 2)



## Combine all statistics tables

statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A<-as.data.frame(rbind(statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events,statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events,statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts,statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts))

statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A$Lifestyle<- gsub("\\..*","",rownames(statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A))

library(dplyr)
statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A$Median<-as.numeric(statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A$Median)
statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A$CI_low<-as.numeric(statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A$CI_low)
statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A$CI_high<-as.numeric(statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A$CI_high)
statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A$ROPE_Percentage<-as.numeric(statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A$ROPE_Percentage)*100
statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A<-statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A %>% mutate_if(is.numeric, round, digits=3)

statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A$Type <-c("EVEs_Events","EVEs_Events","dEVEs_Events","dEVEs_Events","EVEs_Numbers","EVEs_Numbers","dEVEs_Numbers","dEVEs_Numbers")

table_dsDNA_ssDNA_dsRNA_ssRNA_A<-statistics_Table_dsDNA_ssDNA_dsRNA_ssRNA_A
table_dsDNA_ssDNA_dsRNA_ssRNA_A[nrow(table_dsDNA_ssDNA_dsRNA_ssRNA_A)+1,] <- NA





#dsDNA A-B-C-D with controls



library(bayestestR)
library(ggstatsplot)
library(bayesplot)
library(ggpubr)
library(gridExtra)
library(grid)



##################################
# ## EVEs Events all posteriors ##
##################################
Table_dsDNA_EVEs_Events_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_dsDNA_All/Posteriors_EVEs_Events_ALL.txt",sep=";",h=T)

Table_dsDNA_EVEs_Events_ALL$b_value2_coef_inf<- exp(as.numeric(Table_dsDNA_EVEs_Events_ALL$b_value2))
Table_dsDNA_EVEs_Events_ALL$b_value3_coef_inf<- exp(as.numeric(Table_dsDNA_EVEs_Events_ALL$b_value3))
Table_dsDNA_EVEs_Events_ALL<-Table_dsDNA_EVEs_Events_ALL[!Table_dsDNA_EVEs_Events_ALL$b_value2_coef_inf=="Inf",]
Table_dsDNA_EVEs_Events_ALL<-Table_dsDNA_EVEs_Events_ALL[!Table_dsDNA_EVEs_Events_ALL$b_value3_coef_inf=="Inf",]


Table_dsDNA_EVEs_Events_ALL_endo <- as.data.frame(Table_dsDNA_EVEs_Events_ALL$b_value2)
Table_dsDNA_EVEs_Events_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_dsDNA_EVEs_Events_ALL_endo)<- c("coef","Parameter")
Table_dsDNA_EVEs_Events_ALL_endo$coef2 <- as.numeric(Table_dsDNA_EVEs_Events_ALL_endo$coef)
Table_dsDNA_EVEs_Events_ALL_endo$coef <- as.numeric(exp(Table_dsDNA_EVEs_Events_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_EVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_EVEs_Events_ALL_endo$rope <- posterior_tab*100
Table_dsDNA_EVEs_Events_ALL_endo$med <- median(Table_dsDNA_EVEs_Events_ALL_endo$coef)
Table_dsDNA_EVEs_Events_ALL_endo$pd <- pd(Table_dsDNA_EVEs_Events_ALL_endo$coef2)

Table_dsDNA_EVEs_Events_ALL_ecto <- as.data.frame(Table_dsDNA_EVEs_Events_ALL$b_value3)
Table_dsDNA_EVEs_Events_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_dsDNA_EVEs_Events_ALL_ecto)<- c("coef","Parameter")
Table_dsDNA_EVEs_Events_ALL_ecto$coef2 <- as.numeric(Table_dsDNA_EVEs_Events_ALL_ecto$coef)
Table_dsDNA_EVEs_Events_ALL_ecto$coef <- as.numeric(exp(Table_dsDNA_EVEs_Events_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_EVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_EVEs_Events_ALL_ecto$rope <- posterior_tab*100
Table_dsDNA_EVEs_Events_ALL_ecto$med <- median(Table_dsDNA_EVEs_Events_ALL_ecto$coef)
Table_dsDNA_EVEs_Events_ALL_ecto$pd <- pd(Table_dsDNA_EVEs_Events_ALL_ecto$coef2)



Table_dsDNA_EVEs_Events_ALL_ecto_endo<-rbind(Table_dsDNA_EVEs_Events_ALL_endo,Table_dsDNA_EVEs_Events_ALL_ecto)

#Create summary table
statistics_Table_dsDNA_EVEs_Events= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("dsDNA_EVEs_Events",median(Table_dsDNA_EVEs_Events_ALL_endo$coef),ci(Table_dsDNA_EVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_EVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_EVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_dsDNA_EVEs_Events_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("dsDNA_EVEs_Events",median(Table_dsDNA_EVEs_Events_ALL_ecto$coef),ci(Table_dsDNA_EVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_EVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_EVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_dsDNA_EVEs_Events_ALL_ecto$pd[1])))



colnames(statistics_Table_dsDNA_EVEs_Events) <- statistics_Table_dsDNA_EVEs_Events[1,]
statistics_Table_dsDNA_EVEs_Events<-statistics_Table_dsDNA_EVEs_Events[-1,]


Table_dsDNA_EVEs_Events_ALL_ecto_endo$Parameter2 <- paste0(Table_dsDNA_EVEs_Events_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_dsDNA_EVEs_Events_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_dsDNA_EVEs_Events_ALL<-Table_dsDNA_EVEs_Events_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_dsDNA_EVEs_Events_ALL <- stat.test_dsDNA_EVEs_Events_ALL %>% add_xy_position(x = "Parameter2")
stat.test_dsDNA_EVEs_Events_ALL$Prop <- paste0(round((nrow(Table_dsDNA_EVEs_Events_ALL[Table_dsDNA_EVEs_Events_ALL$b_value2 > Table_dsDNA_EVEs_Events_ALL$b_value3,]) / nrow(Table_dsDNA_EVEs_Events_ALL))*100,digits = 2),"%")



boxplot_dsDNA_EVEs_Events <- ggplot(Table_dsDNA_EVEs_Events_ALL_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="Events",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_dsDNA_EVEs_Events_ALL, label ="PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,15))


##################################
## dEVEs Events all posteriors
##################################

Table_dsDNA_dEVEs_Events_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_dsDNA_All/Posteriors_dEVEs_Events_ALL.txt",sep=";",h=T)

Table_dsDNA_dEVEs_Events_ALL$b_value2_coef_inf<- exp(as.numeric(Table_dsDNA_dEVEs_Events_ALL$b_value2))
Table_dsDNA_dEVEs_Events_ALL$b_value3_coef_inf<- exp(as.numeric(Table_dsDNA_dEVEs_Events_ALL$b_value3))
Table_dsDNA_dEVEs_Events_ALL<-Table_dsDNA_dEVEs_Events_ALL[!Table_dsDNA_dEVEs_Events_ALL$b_value2_coef_inf=="Inf",]
Table_dsDNA_dEVEs_Events_ALL<-Table_dsDNA_dEVEs_Events_ALL[!Table_dsDNA_dEVEs_Events_ALL$b_value3_coef_inf=="Inf",]


Table_dsDNA_dEVEs_Events_ALL_endo <- as.data.frame(Table_dsDNA_dEVEs_Events_ALL$b_value2)
Table_dsDNA_dEVEs_Events_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_dsDNA_dEVEs_Events_ALL_endo)<- c("coef","Parameter")
Table_dsDNA_dEVEs_Events_ALL_endo$coef2 <- as.numeric(Table_dsDNA_dEVEs_Events_ALL_endo$coef)
Table_dsDNA_dEVEs_Events_ALL_endo$coef <- as.numeric(exp(Table_dsDNA_dEVEs_Events_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_dEVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_dEVEs_Events_ALL_endo$rope <- posterior_tab*100
Table_dsDNA_dEVEs_Events_ALL_endo$med <- median(Table_dsDNA_dEVEs_Events_ALL_endo$coef)
Table_dsDNA_dEVEs_Events_ALL_endo$pd <- pd(Table_dsDNA_dEVEs_Events_ALL_endo$coef2)

Table_dsDNA_dEVEs_Events_ALL_ecto <- as.data.frame(Table_dsDNA_dEVEs_Events_ALL$b_value3)
Table_dsDNA_dEVEs_Events_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_dsDNA_dEVEs_Events_ALL_ecto)<- c("coef","Parameter")
Table_dsDNA_dEVEs_Events_ALL_ecto$coef2 <- as.numeric(Table_dsDNA_dEVEs_Events_ALL_ecto$coef)
Table_dsDNA_dEVEs_Events_ALL_ecto$coef <- as.numeric(exp(Table_dsDNA_dEVEs_Events_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_dEVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_dEVEs_Events_ALL_ecto$rope <- posterior_tab*100
Table_dsDNA_dEVEs_Events_ALL_ecto$med <- median(Table_dsDNA_dEVEs_Events_ALL_ecto$coef)
Table_dsDNA_dEVEs_Events_ALL_ecto$pd <- pd(Table_dsDNA_dEVEs_Events_ALL_ecto$coef2)



Table_dsDNA_dEVEs_Events_ALL_ecto_endo<-rbind(Table_dsDNA_dEVEs_Events_ALL_endo,Table_dsDNA_dEVEs_Events_ALL_ecto)

#Create summary table
statistics_Table_dsDNA_dEVEs_Events= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("dsDNA_dEVEs_Events",median(Table_dsDNA_dEVEs_Events_ALL_endo$coef),ci(Table_dsDNA_dEVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_dEVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_dEVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_dsDNA_dEVEs_Events_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("dsDNA_dEVEs_Events",median(Table_dsDNA_dEVEs_Events_ALL_ecto$coef),ci(Table_dsDNA_dEVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_dEVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_dEVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_dsDNA_dEVEs_Events_ALL_ecto$pd[1])))



colnames(statistics_Table_dsDNA_dEVEs_Events) <- statistics_Table_dsDNA_dEVEs_Events[1,]
statistics_Table_dsDNA_dEVEs_Events<-statistics_Table_dsDNA_dEVEs_Events[-1,]


Table_dsDNA_dEVEs_Events_ALL_ecto_endo$Parameter2 <- paste0(Table_dsDNA_dEVEs_Events_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_dsDNA_dEVEs_Events_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_dsDNA_dEVEs_Events_ALL<-Table_dsDNA_dEVEs_Events_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_dsDNA_dEVEs_Events_ALL <- stat.test_dsDNA_dEVEs_Events_ALL %>% add_xy_position(x = "Parameter2")
stat.test_dsDNA_dEVEs_Events_ALL$Prop <- paste0(round((nrow(Table_dsDNA_dEVEs_Events_ALL[Table_dsDNA_dEVEs_Events_ALL$b_value2 > Table_dsDNA_dEVEs_Events_ALL$b_value3,]) / nrow(Table_dsDNA_dEVEs_Events_ALL))*100,digits = 2),"%")

boxplot_dsDNA_dEVEs_Events <- ggplot(Table_dsDNA_dEVEs_Events_ALL_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter),trim = T,adjust = 2)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="dEvents",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_dsDNA_dEVEs_Events_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,15))

##################################
## EVEs Counts all posteriors####
##################################


Table_dsDNA_EVEs_Counts_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_dsDNA_All/Posteriors_EVEs_Counts_ALL.txt",sep=";",h=T)

Table_dsDNA_EVEs_Counts_ALL$b_value2_coef_inf<- exp(as.numeric(Table_dsDNA_EVEs_Counts_ALL$b_value2))
Table_dsDNA_EVEs_Counts_ALL$b_value3_coef_inf<- exp(as.numeric(Table_dsDNA_EVEs_Counts_ALL$b_value3))
Table_dsDNA_EVEs_Counts_ALL<-Table_dsDNA_EVEs_Counts_ALL[!Table_dsDNA_EVEs_Counts_ALL$b_value2_coef_inf=="Inf",]
Table_dsDNA_EVEs_Counts_ALL<-Table_dsDNA_EVEs_Counts_ALL[!Table_dsDNA_EVEs_Counts_ALL$b_value3_coef_inf=="Inf",]


Table_dsDNA_EVEs_Counts_ALL_endo <- as.data.frame(Table_dsDNA_EVEs_Counts_ALL$b_value2)
Table_dsDNA_EVEs_Counts_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_dsDNA_EVEs_Counts_ALL_endo)<- c("coef","Parameter")
Table_dsDNA_EVEs_Counts_ALL_endo$coef2 <- as.numeric(Table_dsDNA_EVEs_Counts_ALL_endo$coef)
Table_dsDNA_EVEs_Counts_ALL_endo$coef <- as.numeric(exp(Table_dsDNA_EVEs_Counts_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_EVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_EVEs_Counts_ALL_endo$rope <- posterior_tab*100
Table_dsDNA_EVEs_Counts_ALL_endo$med <- median(Table_dsDNA_EVEs_Counts_ALL_endo$coef)
Table_dsDNA_EVEs_Counts_ALL_endo$pd <- pd(Table_dsDNA_EVEs_Counts_ALL_endo$coef2)

Table_dsDNA_EVEs_Counts_ALL_ecto <- as.data.frame(Table_dsDNA_EVEs_Counts_ALL$b_value3)
Table_dsDNA_EVEs_Counts_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_dsDNA_EVEs_Counts_ALL_ecto)<- c("coef","Parameter")
Table_dsDNA_EVEs_Counts_ALL_ecto$coef2 <- as.numeric(Table_dsDNA_EVEs_Counts_ALL_ecto$coef)
Table_dsDNA_EVEs_Counts_ALL_ecto$coef <- as.numeric(exp(Table_dsDNA_EVEs_Counts_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_EVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_EVEs_Counts_ALL_ecto$rope <- posterior_tab*100
Table_dsDNA_EVEs_Counts_ALL_ecto$med <- median(Table_dsDNA_EVEs_Counts_ALL_ecto$coef)
Table_dsDNA_EVEs_Counts_ALL_ecto$pd <- pd(Table_dsDNA_EVEs_Counts_ALL_ecto$coef2)



Table_dsDNA_EVEs_Counts_ALL_ecto_endo<-rbind(Table_dsDNA_EVEs_Counts_ALL_endo,Table_dsDNA_EVEs_Counts_ALL_ecto)

#Create summary table
statistics_Table_dsDNA_EVEs_Counts= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("dsDNA_EVEs_Counts",median(Table_dsDNA_EVEs_Counts_ALL_endo$coef),ci(Table_dsDNA_EVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_EVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_EVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_dsDNA_EVEs_Counts_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("dsDNA_EVEs_Counts",median(Table_dsDNA_EVEs_Counts_ALL_ecto$coef),ci(Table_dsDNA_EVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_EVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_EVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_dsDNA_EVEs_Counts_ALL_ecto$pd[1])))



colnames(statistics_Table_dsDNA_EVEs_Counts) <- statistics_Table_dsDNA_EVEs_Counts[1,]
statistics_Table_dsDNA_EVEs_Counts<-statistics_Table_dsDNA_EVEs_Counts[-1,]


Table_dsDNA_EVEs_Counts_ALL_ecto_endo$Parameter2 <- paste0(Table_dsDNA_EVEs_Counts_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_dsDNA_EVEs_Counts_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_dsDNA_EVEs_Counts_ALL<-Table_dsDNA_EVEs_Counts_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_dsDNA_EVEs_Counts_ALL <- stat.test_dsDNA_EVEs_Counts_ALL %>% add_xy_position(x = "Parameter2")
stat.test_dsDNA_EVEs_Counts_ALL$Prop <- paste0(round((nrow(Table_dsDNA_EVEs_Counts_ALL[Table_dsDNA_EVEs_Counts_ALL$b_value2 > Table_dsDNA_EVEs_Counts_ALL$b_value3,]) / nrow(Table_dsDNA_EVEs_Counts_ALL))*100,digits = 2),"%")

boxplot_dsDNA_EVEs_Counts <- ggplot(Table_dsDNA_EVEs_Counts_ALL_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="EVEs",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_dsDNA_EVEs_Counts_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,15))


################################
## dEVEs Counts all posteriors #
################################
Table_dsDNA_dEVEs_Counts_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_dsDNA_All/Posteriors_dEVEs_Counts_ALL.txt",sep=";",h=T)

Table_dsDNA_dEVEs_Counts_ALL$b_value2_coef_inf<- exp(as.numeric(Table_dsDNA_dEVEs_Counts_ALL$b_value2))
Table_dsDNA_dEVEs_Counts_ALL$b_value3_coef_inf<- exp(as.numeric(Table_dsDNA_dEVEs_Counts_ALL$b_value3))
Table_dsDNA_dEVEs_Counts_ALL<-Table_dsDNA_dEVEs_Counts_ALL[!Table_dsDNA_dEVEs_Counts_ALL$b_value2_coef_inf=="Inf",]
Table_dsDNA_dEVEs_Counts_ALL<-Table_dsDNA_dEVEs_Counts_ALL[!Table_dsDNA_dEVEs_Counts_ALL$b_value3_coef_inf=="Inf",]


Table_dsDNA_dEVEs_Counts_ALL_endo <- as.data.frame(Table_dsDNA_dEVEs_Counts_ALL$b_value2)
Table_dsDNA_dEVEs_Counts_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_dsDNA_dEVEs_Counts_ALL_endo)<- c("coef","Parameter")
Table_dsDNA_dEVEs_Counts_ALL_endo$coef2 <- as.numeric(Table_dsDNA_dEVEs_Counts_ALL_endo$coef)
Table_dsDNA_dEVEs_Counts_ALL_endo$coef <- as.numeric(exp(Table_dsDNA_dEVEs_Counts_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_dEVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_dEVEs_Counts_ALL_endo$rope <- posterior_tab*100
Table_dsDNA_dEVEs_Counts_ALL_endo$med <- median(Table_dsDNA_dEVEs_Counts_ALL_endo$coef)
Table_dsDNA_dEVEs_Counts_ALL_endo$pd <- pd(Table_dsDNA_dEVEs_Counts_ALL_endo$coef2)

Table_dsDNA_dEVEs_Counts_ALL_ecto <- as.data.frame(Table_dsDNA_dEVEs_Counts_ALL$b_value3)
Table_dsDNA_dEVEs_Counts_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_dsDNA_dEVEs_Counts_ALL_ecto)<- c("coef","Parameter")
Table_dsDNA_dEVEs_Counts_ALL_ecto$coef2 <- as.numeric(Table_dsDNA_dEVEs_Counts_ALL_ecto$coef)
Table_dsDNA_dEVEs_Counts_ALL_ecto$coef <- as.numeric(exp(Table_dsDNA_dEVEs_Counts_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_dEVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_dEVEs_Counts_ALL_ecto$rope <- posterior_tab*100
Table_dsDNA_dEVEs_Counts_ALL_ecto$med <- median(Table_dsDNA_dEVEs_Counts_ALL_ecto$coef)
Table_dsDNA_dEVEs_Counts_ALL_ecto$pd <- pd(Table_dsDNA_dEVEs_Counts_ALL_ecto$coef2)



Table_dsDNA_dEVEs_Counts_ALL_ecto_endo<-rbind(Table_dsDNA_dEVEs_Counts_ALL_endo,Table_dsDNA_dEVEs_Counts_ALL_ecto)

#Create summary table
statistics_Table_dsDNA_dEVEs_Counts= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("dsDNA_dEVEs_Counts",median(Table_dsDNA_dEVEs_Counts_ALL_endo$coef),ci(Table_dsDNA_dEVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_dEVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_dEVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_dsDNA_dEVEs_Counts_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("dsDNA_dEVEs_Counts",median(Table_dsDNA_dEVEs_Counts_ALL_ecto$coef),ci(Table_dsDNA_dEVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_dEVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_dEVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_dsDNA_dEVEs_Counts_ALL_ecto$pd[1])))



colnames(statistics_Table_dsDNA_dEVEs_Counts) <- statistics_Table_dsDNA_dEVEs_Counts[1,]
statistics_Table_dsDNA_dEVEs_Counts<-statistics_Table_dsDNA_dEVEs_Counts[-1,]


Table_dsDNA_dEVEs_Counts_ALL_ecto_endo$Parameter2 <- paste0(Table_dsDNA_dEVEs_Counts_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_dsDNA_dEVEs_Counts_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_dsDNA_dEVEs_Counts_ALL<-Table_dsDNA_dEVEs_Counts_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_dsDNA_dEVEs_Counts_ALL <- stat.test_dsDNA_dEVEs_Counts_ALL %>% add_xy_position(x = "Parameter2")
stat.test_dsDNA_dEVEs_Counts_ALL$Prop <- paste0(round((nrow(Table_dsDNA_dEVEs_Counts_ALL[Table_dsDNA_dEVEs_Counts_ALL$b_value2 > Table_dsDNA_dEVEs_Counts_ALL$b_value3,]) / nrow(Table_dsDNA_dEVEs_Counts_ALL))*100,digits = 2),"%")

boxplot_dsDNA_dEVEs_Counts <- ggplot(Table_dsDNA_dEVEs_Counts_ALL_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="dEVEs",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_dsDNA_dEVEs_Counts_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,15))


## Combine All 4 plots

violin_plot_dsDNA<-grid.arrange(arrangeGrob(
                             boxplot_dsDNA_EVEs_Counts  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank()) ,
                             boxplot_dsDNA_EVEs_Events  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank()),
                             boxplot_dsDNA_dEVEs_Counts  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank()) ,
                             boxplot_dsDNA_dEVEs_Events  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank()),
                         nrow = 1,
                         ncol=4,
                         left = textGrob("Relative rate", rot = 90, vjust = 1)),nrow=1,top=textGrob("dsDNA | A-B-C-D",
                                                                                                                 gp=gpar(fontsize=20,font=8)))

violin_plot_dsDNA_publication<-grid.arrange(arrangeGrob(

                                  boxplot_dsDNA_EVEs_Counts  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank()) ,
                             boxplot_dsDNA_EVEs_Events  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank()) ,
                             boxplot_dsDNA_dEVEs_Counts  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank()),
                             boxplot_dsDNA_dEVEs_Events  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank()),
                         nrow = 2,
                         ncol=2,
                         left = textGrob("Relative rate", rot = 90, vjust = 1)),nrow=1)#,top=textGrob("dsDNA | A-B-C-D",
                                                                                                                # gp=gpar(fontsize=20,font=8)))


ggsave(filename="/Users/bguinet/Desktop/Papier_scientifique/dsDNA_violin_plot.pdf",violin_plot_dsDNA_publication, width = 250, height = 350, units = "mm",device="pdf")


#Add posterior coefficient comparison between endo and ecto
statistics_Table_dsDNA_EVEs_Events<-as.data.frame(statistics_Table_dsDNA_EVEs_Events)
statistics_Table_dsDNA_EVEs_Events$Endo_sup_Ecto<-round((nrow(Table_dsDNA_EVEs_Events_ALL[Table_dsDNA_EVEs_Events_ALL$b_value2 > Table_dsDNA_EVEs_Events_ALL$b_value3,]) / nrow(Table_dsDNA_EVEs_Events_ALL))*100,digits = 2)

statistics_Table_dsDNA_dEVEs_Events<-as.data.frame(statistics_Table_dsDNA_dEVEs_Events)
statistics_Table_dsDNA_dEVEs_Events$Endo_sup_Ecto<-round((nrow(Table_dsDNA_dEVEs_Events_ALL[Table_dsDNA_dEVEs_Events_ALL$b_value2 > Table_dsDNA_dEVEs_Events_ALL$b_value3,]) / nrow(Table_dsDNA_dEVEs_Events_ALL))*100,digits = 2)

statistics_Table_dsDNA_EVEs_Counts<-as.data.frame(statistics_Table_dsDNA_EVEs_Counts)
statistics_Table_dsDNA_EVEs_Counts$Endo_sup_Ecto<-round((nrow(Table_dsDNA_EVEs_Counts_ALL[Table_dsDNA_EVEs_Counts_ALL$b_value2 > Table_dsDNA_EVEs_Counts_ALL$b_value3,]) / nrow(Table_dsDNA_EVEs_Counts_ALL))*100,digits = 2)


statistics_Table_dsDNA_dEVEs_Counts<-as.data.frame(statistics_Table_dsDNA_dEVEs_Counts)
statistics_Table_dsDNA_dEVEs_Counts$Endo_sup_Ecto<-round((nrow(Table_dsDNA_dEVEs_Counts_ALL[Table_dsDNA_dEVEs_Counts_ALL$b_value2 > Table_dsDNA_dEVEs_Counts_ALL$b_value3,]) / nrow(Table_dsDNA_dEVEs_Counts_ALL))*100,digits = 2)


## Combine all statistics tables

statistics_Table_dsDNA<-as.data.frame(rbind(statistics_Table_dsDNA_EVEs_Events,statistics_Table_dsDNA_dEVEs_Events,statistics_Table_dsDNA_EVEs_Counts,statistics_Table_dsDNA_dEVEs_Counts))

statistics_Table_dsDNA$Lifestyle<- gsub("\\..*","",rownames(statistics_Table_dsDNA))

library(dplyr)
statistics_Table_dsDNA$Median<-as.numeric(statistics_Table_dsDNA$Median)
statistics_Table_dsDNA$CI_low<-as.numeric(statistics_Table_dsDNA$CI_low)
statistics_Table_dsDNA$CI_high<-as.numeric(statistics_Table_dsDNA$CI_high)
statistics_Table_dsDNA$ROPE_Percentage<-as.numeric(statistics_Table_dsDNA$ROPE_Percentage)*100
statistics_Table_dsDNA<-statistics_Table_dsDNA %>% mutate_if(is.numeric, round, digits=3)

statistics_Table_dsDNA$Type <-c("EVEs_Events","EVEs_Events","dEVEs_Events","dEVEs_Events","EVEs_Numbers","EVEs_Numbers","dEVEs_Numbers","dEVEs_Numbers")

table_dsDNA<-statistics_Table_dsDNA
table_dsDNA[nrow(table_dsDNA)+1,] <- NA






#dsDNA A with controls


library(bayestestR)
library(ggstatsplot)
library(bayesplot)
library(ggpubr)
library(gridExtra)
library(grid)



##################################
# ## EVEs Events all posteriors ##
##################################
Table_dsDNA_A_EVEs_Events_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_dsDNA_A_All/Posteriors_EVEs_Events_ALL.txt",sep=";",h=T)

Table_dsDNA_A_EVEs_Events_ALL$b_value2_coef_inf<- exp(as.numeric(Table_dsDNA_A_EVEs_Events_ALL$b_value2))
Table_dsDNA_A_EVEs_Events_ALL$b_value3_coef_inf<- exp(as.numeric(Table_dsDNA_A_EVEs_Events_ALL$b_value3))
Table_dsDNA_A_EVEs_Events_ALL<-Table_dsDNA_A_EVEs_Events_ALL[!Table_dsDNA_A_EVEs_Events_ALL$b_value2_coef_inf=="Inf",]
Table_dsDNA_A_EVEs_Events_ALL<-Table_dsDNA_A_EVEs_Events_ALL[!Table_dsDNA_A_EVEs_Events_ALL$b_value3_coef_inf=="Inf",]


Table_dsDNA_A_EVEs_Events_ALL_endo <- as.data.frame(Table_dsDNA_A_EVEs_Events_ALL$b_value2)
Table_dsDNA_A_EVEs_Events_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_dsDNA_A_EVEs_Events_ALL_endo)<- c("coef","Parameter")
Table_dsDNA_A_EVEs_Events_ALL_endo$coef2 <- as.numeric(Table_dsDNA_A_EVEs_Events_ALL_endo$coef)
Table_dsDNA_A_EVEs_Events_ALL_endo$coef <- as.numeric(exp(Table_dsDNA_A_EVEs_Events_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_A_EVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_A_EVEs_Events_ALL_endo$rope <- posterior_tab*100
Table_dsDNA_A_EVEs_Events_ALL_endo$med <- median(Table_dsDNA_A_EVEs_Events_ALL_endo$coef)
Table_dsDNA_A_EVEs_Events_ALL_endo$pd <- pd(Table_dsDNA_A_EVEs_Events_ALL_endo$coef2)

Table_dsDNA_A_EVEs_Events_ALL_ecto <- as.data.frame(Table_dsDNA_A_EVEs_Events_ALL$b_value3)
Table_dsDNA_A_EVEs_Events_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_dsDNA_A_EVEs_Events_ALL_ecto)<- c("coef","Parameter")
Table_dsDNA_A_EVEs_Events_ALL_ecto$coef2 <- as.numeric(Table_dsDNA_A_EVEs_Events_ALL_ecto$coef)
Table_dsDNA_A_EVEs_Events_ALL_ecto$coef <- as.numeric(exp(Table_dsDNA_A_EVEs_Events_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_A_EVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_A_EVEs_Events_ALL_ecto$rope <- posterior_tab*100
Table_dsDNA_A_EVEs_Events_ALL_ecto$med <- median(Table_dsDNA_A_EVEs_Events_ALL_ecto$coef)
Table_dsDNA_A_EVEs_Events_ALL_ecto$pd <- pd(Table_dsDNA_A_EVEs_Events_ALL_ecto$coef2)



Table_dsDNA_A_EVEs_Events_ALL_ecto_endo<-rbind(Table_dsDNA_A_EVEs_Events_ALL_endo,Table_dsDNA_A_EVEs_Events_ALL_ecto)

#Create summary table
statistics_Table_dsDNA_A_EVEs_Events= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("dsDNA_A_EVEs_Events",median(Table_dsDNA_A_EVEs_Events_ALL_endo$coef),ci(Table_dsDNA_A_EVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_A_EVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_A_EVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_dsDNA_A_EVEs_Events_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("dsDNA_A_EVEs_Events",median(Table_dsDNA_A_EVEs_Events_ALL_ecto$coef),ci(Table_dsDNA_A_EVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_A_EVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_A_EVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_dsDNA_A_EVEs_Events_ALL_ecto$pd[1])))



colnames(statistics_Table_dsDNA_A_EVEs_Events) <- statistics_Table_dsDNA_A_EVEs_Events[1,]
statistics_Table_dsDNA_A_EVEs_Events<-statistics_Table_dsDNA_A_EVEs_Events[-1,]


Table_dsDNA_A_EVEs_Events_ALL_ecto_endo$Parameter2 <- paste0(Table_dsDNA_A_EVEs_Events_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_dsDNA_A_EVEs_Events_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_dsDNA_A_EVEs_Events_ALL<-Table_dsDNA_A_EVEs_Events_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_dsDNA_A_EVEs_Events_ALL <- stat.test_dsDNA_A_EVEs_Events_ALL %>% add_xy_position(x = "Parameter2")
stat.test_dsDNA_A_EVEs_Events_ALL$Prop <- paste0(round((nrow(Table_dsDNA_A_EVEs_Events_ALL[Table_dsDNA_A_EVEs_Events_ALL$b_value2 > Table_dsDNA_A_EVEs_Events_ALL$b_value3,]) / nrow(Table_dsDNA_A_EVEs_Events_ALL))*100,digits = 2),"%")

boxplot_dsDNA_A_EVEs_Events <- ggplot(Table_dsDNA_A_EVEs_Events_ALL_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="Events",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_dsDNA_A_EVEs_Events_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,15))


##################################
## dEVEs Events all posteriors
##################################

Table_dsDNA_A_dEVEs_Events_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_dsDNA_A_All/Posteriors_dEVEs_Events_ALL.txt",sep=";",h=T)

Table_dsDNA_A_dEVEs_Events_ALL$b_value2_coef_inf<- exp(as.numeric(Table_dsDNA_A_dEVEs_Events_ALL$b_value2))
Table_dsDNA_A_dEVEs_Events_ALL$b_value3_coef_inf<- exp(as.numeric(Table_dsDNA_A_dEVEs_Events_ALL$b_value3))
Table_dsDNA_A_dEVEs_Events_ALL<-Table_dsDNA_A_dEVEs_Events_ALL[!Table_dsDNA_A_dEVEs_Events_ALL$b_value2_coef_inf=="Inf",]
Table_dsDNA_A_dEVEs_Events_ALL<-Table_dsDNA_A_dEVEs_Events_ALL[!Table_dsDNA_A_dEVEs_Events_ALL$b_value3_coef_inf=="Inf",]


Table_dsDNA_A_dEVEs_Events_ALL_endo <- as.data.frame(Table_dsDNA_A_dEVEs_Events_ALL$b_value2)
Table_dsDNA_A_dEVEs_Events_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_dsDNA_A_dEVEs_Events_ALL_endo)<- c("coef","Parameter")
Table_dsDNA_A_dEVEs_Events_ALL_endo$coef2 <- as.numeric(Table_dsDNA_A_dEVEs_Events_ALL_endo$coef)
Table_dsDNA_A_dEVEs_Events_ALL_endo$coef <- as.numeric(exp(Table_dsDNA_A_dEVEs_Events_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_A_dEVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_A_dEVEs_Events_ALL_endo$rope <- posterior_tab*100
Table_dsDNA_A_dEVEs_Events_ALL_endo$med <- median(Table_dsDNA_A_dEVEs_Events_ALL_endo$coef)
Table_dsDNA_A_dEVEs_Events_ALL_endo$pd <- pd(Table_dsDNA_A_dEVEs_Events_ALL_endo$coef2)

Table_dsDNA_A_dEVEs_Events_ALL_ecto <- as.data.frame(Table_dsDNA_A_dEVEs_Events_ALL$b_value3)
Table_dsDNA_A_dEVEs_Events_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_dsDNA_A_dEVEs_Events_ALL_ecto)<- c("coef","Parameter")
Table_dsDNA_A_dEVEs_Events_ALL_ecto$coef2 <- as.numeric(Table_dsDNA_A_dEVEs_Events_ALL_ecto$coef)
Table_dsDNA_A_dEVEs_Events_ALL_ecto$coef <- as.numeric(exp(Table_dsDNA_A_dEVEs_Events_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_A_dEVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_A_dEVEs_Events_ALL_ecto$rope <- posterior_tab*100
Table_dsDNA_A_dEVEs_Events_ALL_ecto$med <- median(Table_dsDNA_A_dEVEs_Events_ALL_ecto$coef)
Table_dsDNA_A_dEVEs_Events_ALL_ecto$pd <- pd(Table_dsDNA_A_dEVEs_Events_ALL_ecto$coef2)



Table_dsDNA_A_dEVEs_Events_ALL_ecto_endo<-rbind(Table_dsDNA_A_dEVEs_Events_ALL_endo,Table_dsDNA_A_dEVEs_Events_ALL_ecto)

#Create summary table
statistics_Table_dsDNA_A_dEVEs_Events= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("dsDNA_A_dEVEs_Events",median(Table_dsDNA_A_dEVEs_Events_ALL_endo$coef),ci(Table_dsDNA_A_dEVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_A_dEVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_A_dEVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_dsDNA_A_dEVEs_Events_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("dsDNA_A_dEVEs_Events",median(Table_dsDNA_A_dEVEs_Events_ALL_ecto$coef),ci(Table_dsDNA_A_dEVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_A_dEVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_A_dEVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_dsDNA_A_dEVEs_Events_ALL_ecto$pd[1])))



colnames(statistics_Table_dsDNA_A_dEVEs_Events) <- statistics_Table_dsDNA_A_dEVEs_Events[1,]
statistics_Table_dsDNA_A_dEVEs_Events<-statistics_Table_dsDNA_A_dEVEs_Events[-1,]


Table_dsDNA_A_dEVEs_Events_ALL_ecto_endo$Parameter2 <- paste0(Table_dsDNA_A_dEVEs_Events_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_dsDNA_A_dEVEs_Events_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_dsDNA_A_dEVEs_Events_ALL<-Table_dsDNA_A_dEVEs_Events_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_dsDNA_A_dEVEs_Events_ALL <- stat.test_dsDNA_A_dEVEs_Events_ALL %>% add_xy_position(x = "Parameter2")
stat.test_dsDNA_A_dEVEs_Events_ALL$Prop <- paste0(round((nrow(Table_dsDNA_A_dEVEs_Events_ALL[Table_dsDNA_A_dEVEs_Events_ALL$b_value2 > Table_dsDNA_A_dEVEs_Events_ALL$b_value3,]) / nrow(Table_dsDNA_A_dEVEs_Events_ALL))*100,digits = 2),"%")

boxplot_dsDNA_A_dEVEs_Events <- ggplot(Table_dsDNA_A_dEVEs_Events_ALL_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter),trim = T,adjust = 2)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="dEvents",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_dsDNA_A_dEVEs_Events_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,15))

##################################
## EVEs Counts all posteriors####
##################################


Table_dsDNA_A_EVEs_Counts_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_dsDNA_A_All/Posteriors_EVEs_Counts_ALL.txt",sep=";",h=T)

Table_dsDNA_A_EVEs_Counts_ALL$b_value2_coef_inf<- exp(as.numeric(Table_dsDNA_A_EVEs_Counts_ALL$b_value2))
Table_dsDNA_A_EVEs_Counts_ALL$b_value3_coef_inf<- exp(as.numeric(Table_dsDNA_A_EVEs_Counts_ALL$b_value3))
Table_dsDNA_A_EVEs_Counts_ALL<-Table_dsDNA_A_EVEs_Counts_ALL[!Table_dsDNA_A_EVEs_Counts_ALL$b_value2_coef_inf=="Inf",]
Table_dsDNA_A_EVEs_Counts_ALL<-Table_dsDNA_A_EVEs_Counts_ALL[!Table_dsDNA_A_EVEs_Counts_ALL$b_value3_coef_inf=="Inf",]


Table_dsDNA_A_EVEs_Counts_ALL_endo <- as.data.frame(Table_dsDNA_A_EVEs_Counts_ALL$b_value2)
Table_dsDNA_A_EVEs_Counts_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_dsDNA_A_EVEs_Counts_ALL_endo)<- c("coef","Parameter")
Table_dsDNA_A_EVEs_Counts_ALL_endo$coef2 <- as.numeric(Table_dsDNA_A_EVEs_Counts_ALL_endo$coef)
Table_dsDNA_A_EVEs_Counts_ALL_endo$coef <- as.numeric(exp(Table_dsDNA_A_EVEs_Counts_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_A_EVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_A_EVEs_Counts_ALL_endo$rope <- posterior_tab*100
Table_dsDNA_A_EVEs_Counts_ALL_endo$med <- median(Table_dsDNA_A_EVEs_Counts_ALL_endo$coef)
Table_dsDNA_A_EVEs_Counts_ALL_endo$pd <- pd(Table_dsDNA_A_EVEs_Counts_ALL_endo$coef2)

Table_dsDNA_A_EVEs_Counts_ALL_ecto <- as.data.frame(Table_dsDNA_A_EVEs_Counts_ALL$b_value3)
Table_dsDNA_A_EVEs_Counts_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_dsDNA_A_EVEs_Counts_ALL_ecto)<- c("coef","Parameter")
Table_dsDNA_A_EVEs_Counts_ALL_ecto$coef2 <- as.numeric(Table_dsDNA_A_EVEs_Counts_ALL_ecto$coef)
Table_dsDNA_A_EVEs_Counts_ALL_ecto$coef <- as.numeric(exp(Table_dsDNA_A_EVEs_Counts_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_A_EVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_A_EVEs_Counts_ALL_ecto$rope <- posterior_tab*100
Table_dsDNA_A_EVEs_Counts_ALL_ecto$med <- median(Table_dsDNA_A_EVEs_Counts_ALL_ecto$coef)
Table_dsDNA_A_EVEs_Counts_ALL_ecto$pd <- pd(Table_dsDNA_A_EVEs_Counts_ALL_ecto$coef2)



Table_dsDNA_A_EVEs_Counts_ALL_ecto_endo<-rbind(Table_dsDNA_A_EVEs_Counts_ALL_endo,Table_dsDNA_A_EVEs_Counts_ALL_ecto)

#Create summary table
statistics_Table_dsDNA_A_EVEs_Counts= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("dsDNA_A_EVEs_Counts",median(Table_dsDNA_A_EVEs_Counts_ALL_endo$coef),ci(Table_dsDNA_A_EVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_A_EVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_A_EVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_dsDNA_A_EVEs_Counts_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("dsDNA_A_EVEs_Counts",median(Table_dsDNA_A_EVEs_Counts_ALL_ecto$coef),ci(Table_dsDNA_A_EVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_A_EVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_A_EVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_dsDNA_A_EVEs_Counts_ALL_ecto$pd[1])))



colnames(statistics_Table_dsDNA_A_EVEs_Counts) <- statistics_Table_dsDNA_A_EVEs_Counts[1,]
statistics_Table_dsDNA_A_EVEs_Counts<-statistics_Table_dsDNA_A_EVEs_Counts[-1,]


Table_dsDNA_A_EVEs_Counts_ALL_ecto_endo$Parameter2 <- paste0(Table_dsDNA_A_EVEs_Counts_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_dsDNA_A_EVEs_Counts_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_dsDNA_A_EVEs_Counts_ALL<-Table_dsDNA_A_EVEs_Counts_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_dsDNA_A_EVEs_Counts_ALL <- stat.test_dsDNA_A_EVEs_Counts_ALL %>% add_xy_position(x = "Parameter2")
stat.test_dsDNA_A_EVEs_Counts_ALL$Prop <- paste0(round((nrow(Table_dsDNA_A_EVEs_Counts_ALL[Table_dsDNA_A_EVEs_Counts_ALL$b_value2 > Table_dsDNA_A_EVEs_Counts_ALL$b_value3,]) / nrow(Table_dsDNA_A_EVEs_Counts_ALL))*100,digits = 2),"%")

boxplot_dsDNA_A_EVEs_Counts <- ggplot(Table_dsDNA_A_EVEs_Counts_ALL_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="EVEs",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_dsDNA_A_EVEs_Counts_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,15))


################################
## dEVEs Counts all posteriors #
################################
Table_dsDNA_A_dEVEs_Counts_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_dsDNA_A_All/Posteriors_dEVEs_Counts_ALL.txt",sep=";",h=T)

Table_dsDNA_A_dEVEs_Counts_ALL$b_value2_coef_inf<- exp(as.numeric(Table_dsDNA_A_dEVEs_Counts_ALL$b_value2))
Table_dsDNA_A_dEVEs_Counts_ALL$b_value3_coef_inf<- exp(as.numeric(Table_dsDNA_A_dEVEs_Counts_ALL$b_value3))
Table_dsDNA_A_dEVEs_Counts_ALL<-Table_dsDNA_A_dEVEs_Counts_ALL[!Table_dsDNA_A_dEVEs_Counts_ALL$b_value2_coef_inf=="Inf",]
Table_dsDNA_A_dEVEs_Counts_ALL<-Table_dsDNA_A_dEVEs_Counts_ALL[!Table_dsDNA_A_dEVEs_Counts_ALL$b_value3_coef_inf=="Inf",]


Table_dsDNA_A_dEVEs_Counts_ALL_endo <- as.data.frame(Table_dsDNA_A_dEVEs_Counts_ALL$b_value2)
Table_dsDNA_A_dEVEs_Counts_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_dsDNA_A_dEVEs_Counts_ALL_endo)<- c("coef","Parameter")
Table_dsDNA_A_dEVEs_Counts_ALL_endo$coef2 <- as.numeric(Table_dsDNA_A_dEVEs_Counts_ALL_endo$coef)
Table_dsDNA_A_dEVEs_Counts_ALL_endo$coef <- as.numeric(exp(Table_dsDNA_A_dEVEs_Counts_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_A_dEVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_A_dEVEs_Counts_ALL_endo$rope <- posterior_tab*100
Table_dsDNA_A_dEVEs_Counts_ALL_endo$med <- median(Table_dsDNA_A_dEVEs_Counts_ALL_endo$coef)
Table_dsDNA_A_dEVEs_Counts_ALL_endo$pd <- pd(Table_dsDNA_A_dEVEs_Counts_ALL_endo$coef2)

Table_dsDNA_A_dEVEs_Counts_ALL_ecto <- as.data.frame(Table_dsDNA_A_dEVEs_Counts_ALL$b_value3)
Table_dsDNA_A_dEVEs_Counts_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_dsDNA_A_dEVEs_Counts_ALL_ecto)<- c("coef","Parameter")
Table_dsDNA_A_dEVEs_Counts_ALL_ecto$coef2 <- as.numeric(Table_dsDNA_A_dEVEs_Counts_ALL_ecto$coef)
Table_dsDNA_A_dEVEs_Counts_ALL_ecto$coef <- as.numeric(exp(Table_dsDNA_A_dEVEs_Counts_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_A_dEVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_A_dEVEs_Counts_ALL_ecto$rope <- posterior_tab*100
Table_dsDNA_A_dEVEs_Counts_ALL_ecto$med <- median(Table_dsDNA_A_dEVEs_Counts_ALL_ecto$coef)
Table_dsDNA_A_dEVEs_Counts_ALL_ecto$pd <- pd(Table_dsDNA_A_dEVEs_Counts_ALL_ecto$coef2)



Table_dsDNA_A_dEVEs_Counts_ALL_ecto_endo<-rbind(Table_dsDNA_A_dEVEs_Counts_ALL_endo,Table_dsDNA_A_dEVEs_Counts_ALL_ecto)

#Create summary table
statistics_Table_dsDNA_A_dEVEs_Counts= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("dsDNA_A_dEVEs_Counts",median(Table_dsDNA_A_dEVEs_Counts_ALL_endo$coef),ci(Table_dsDNA_A_dEVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_A_dEVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_A_dEVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_dsDNA_A_dEVEs_Counts_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("dsDNA_A_dEVEs_Counts",median(Table_dsDNA_A_dEVEs_Counts_ALL_ecto$coef),ci(Table_dsDNA_A_dEVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_A_dEVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_A_dEVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_dsDNA_A_dEVEs_Counts_ALL_ecto$pd[1])))



colnames(statistics_Table_dsDNA_A_dEVEs_Counts) <- statistics_Table_dsDNA_A_dEVEs_Counts[1,]
statistics_Table_dsDNA_A_dEVEs_Counts<-statistics_Table_dsDNA_A_dEVEs_Counts[-1,]


Table_dsDNA_A_dEVEs_Counts_ALL_ecto_endo$Parameter2 <- paste0(Table_dsDNA_A_dEVEs_Counts_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_dsDNA_A_dEVEs_Counts_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_dsDNA_A_dEVEs_Counts_ALL<-Table_dsDNA_A_dEVEs_Counts_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_dsDNA_A_dEVEs_Counts_ALL <- stat.test_dsDNA_A_dEVEs_Counts_ALL %>% add_xy_position(x = "Parameter2")
stat.test_dsDNA_A_dEVEs_Counts_ALL$Prop <- paste0(round((nrow(Table_dsDNA_A_dEVEs_Counts_ALL[Table_dsDNA_A_dEVEs_Counts_ALL$b_value2 > Table_dsDNA_A_dEVEs_Counts_ALL$b_value3,]) / nrow(Table_dsDNA_A_dEVEs_Counts_ALL))*100,digits = 2),"%")

boxplot_dsDNA_A_dEVEs_Counts <- ggplot(Table_dsDNA_A_dEVEs_Counts_ALL_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="dEVEs",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_dsDNA_A_dEVEs_Counts_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,15))

## Combine All 4 plots

violin_plot_dsDNA_A<-grid.arrange(arrangeGrob(
                             boxplot_dsDNA_A_EVEs_Counts  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank()) ,
                             boxplot_dsDNA_A_EVEs_Events  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank()),
                            boxplot_dsDNA_A_dEVEs_Counts  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank()) ,
                             boxplot_dsDNA_A_dEVEs_Events  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank()),
                         nrow = 1,
                         ncol=4,
                         left = textGrob("Relative rate", rot = 90, vjust = 1)),nrow=1,top=textGrob("dsDNA | A",
                                                                                                                 gp=gpar(fontsize=20,font=8)))


#Add posterior coefficient comparison between endo and ecto
statistics_Table_dsDNA_A_EVEs_Events<-as.data.frame(statistics_Table_dsDNA_A_EVEs_Events)
statistics_Table_dsDNA_A_EVEs_Events$Endo_sup_Ecto<-round((nrow(Table_dsDNA_A_EVEs_Events_ALL[Table_dsDNA_A_EVEs_Events_ALL$b_value2 > Table_dsDNA_A_EVEs_Events_ALL$b_value3,]) / nrow(Table_dsDNA_A_EVEs_Events_ALL))*100,digits = 2)

statistics_Table_dsDNA_A_dEVEs_Events<-as.data.frame(statistics_Table_dsDNA_A_dEVEs_Events)
statistics_Table_dsDNA_A_dEVEs_Events$Endo_sup_Ecto<-round((nrow(Table_dsDNA_A_dEVEs_Events_ALL[Table_dsDNA_A_dEVEs_Events_ALL$b_value2 > Table_dsDNA_A_dEVEs_Events_ALL$b_value3,]) / nrow(Table_dsDNA_A_dEVEs_Events_ALL))*100,digits = 2)

statistics_Table_dsDNA_A_EVEs_Counts<-as.data.frame(statistics_Table_dsDNA_A_EVEs_Counts)
statistics_Table_dsDNA_A_EVEs_Counts$Endo_sup_Ecto<-round((nrow(Table_dsDNA_A_EVEs_Counts_ALL[Table_dsDNA_A_EVEs_Counts_ALL$b_value2 > Table_dsDNA_A_EVEs_Counts_ALL$b_value3,]) / nrow(Table_dsDNA_A_EVEs_Counts_ALL))*100,digits = 2)


statistics_Table_dsDNA_A_dEVEs_Counts<-as.data.frame(statistics_Table_dsDNA_A_dEVEs_Counts)
statistics_Table_dsDNA_A_dEVEs_Counts$Endo_sup_Ecto<-round((nrow(Table_dsDNA_A_dEVEs_Counts_ALL[Table_dsDNA_A_dEVEs_Counts_ALL$b_value2 > Table_dsDNA_A_dEVEs_Counts_ALL$b_value3,]) / nrow(Table_dsDNA_A_dEVEs_Counts_ALL))*100,digits = 2)


## Combine all statistics tables

statistics_Table_dsDNA_A<-as.data.frame(rbind(statistics_Table_dsDNA_A_EVEs_Events,statistics_Table_dsDNA_A_dEVEs_Events,statistics_Table_dsDNA_A_EVEs_Counts,statistics_Table_dsDNA_A_dEVEs_Counts))

statistics_Table_dsDNA_A$Lifestyle<- gsub("\\..*","",rownames(statistics_Table_dsDNA_A))

library(dplyr)
statistics_Table_dsDNA_A$Median<-as.numeric(statistics_Table_dsDNA_A$Median)
statistics_Table_dsDNA_A$CI_low<-as.numeric(statistics_Table_dsDNA_A$CI_low)
statistics_Table_dsDNA_A$CI_high<-as.numeric(statistics_Table_dsDNA_A$CI_high)
statistics_Table_dsDNA_A$ROPE_Percentage<-as.numeric(statistics_Table_dsDNA_A$ROPE_Percentage)*100
statistics_Table_dsDNA_A<-statistics_Table_dsDNA_A %>% mutate_if(is.numeric, round, digits=3)

statistics_Table_dsDNA_A$Type <-c("EVEs_Events","EVEs_Events","dEVEs_Events","dEVEs_Events","EVEs_Numbers","EVEs_Numbers","dEVEs_Numbers","dEVEs_Numbers")

table_dsDNA_A<-statistics_Table_dsDNA_A
table_dsDNA_A[nrow(table_dsDNA_A)+1,] <- NA




#dsDNA A-B-C-D without controls



library(bayestestR)
library(ggstatsplot)
library(bayesplot)
library(ggpubr)
library(gridExtra)
library(grid)




##################################
# ## EVEs Events all posteriors ##
##################################
Table_dsDNA_without_controls_EVEs_Events_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_dsDNA_without_controls_All/Posteriors_EVEs_Events_ALL.txt",sep=";",h=T)

Table_dsDNA_without_controls_EVEs_Events_ALL$b_value2_coef_inf<- exp(as.numeric(Table_dsDNA_without_controls_EVEs_Events_ALL$b_value2))
Table_dsDNA_without_controls_EVEs_Events_ALL$b_value3_coef_inf<- exp(as.numeric(Table_dsDNA_without_controls_EVEs_Events_ALL$b_value3))
Table_dsDNA_without_controls_EVEs_Events_ALL<-Table_dsDNA_without_controls_EVEs_Events_ALL[!Table_dsDNA_without_controls_EVEs_Events_ALL$b_value2_coef_inf=="Inf",]
Table_dsDNA_without_controls_EVEs_Events_ALL<-Table_dsDNA_without_controls_EVEs_Events_ALL[!Table_dsDNA_without_controls_EVEs_Events_ALL$b_value3_coef_inf=="Inf",]


Table_dsDNA_without_controls_EVEs_Events_ALL_endo <- as.data.frame(Table_dsDNA_without_controls_EVEs_Events_ALL$b_value2)
Table_dsDNA_without_controls_EVEs_Events_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_dsDNA_without_controls_EVEs_Events_ALL_endo)<- c("coef","Parameter")
Table_dsDNA_without_controls_EVEs_Events_ALL_endo$coef2 <- as.numeric(Table_dsDNA_without_controls_EVEs_Events_ALL_endo$coef)
Table_dsDNA_without_controls_EVEs_Events_ALL_endo$coef <- as.numeric(exp(Table_dsDNA_without_controls_EVEs_Events_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_without_controls_EVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_without_controls_EVEs_Events_ALL_endo$rope <- posterior_tab*100
Table_dsDNA_without_controls_EVEs_Events_ALL_endo$med <- median(Table_dsDNA_without_controls_EVEs_Events_ALL_endo$coef)
Table_dsDNA_without_controls_EVEs_Events_ALL_endo$pd <- pd(Table_dsDNA_without_controls_EVEs_Events_ALL_endo$coef2)

Table_dsDNA_without_controls_EVEs_Events_ALL_ecto <- as.data.frame(Table_dsDNA_without_controls_EVEs_Events_ALL$b_value3)
Table_dsDNA_without_controls_EVEs_Events_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_dsDNA_without_controls_EVEs_Events_ALL_ecto)<- c("coef","Parameter")
Table_dsDNA_without_controls_EVEs_Events_ALL_ecto$coef2 <- as.numeric(Table_dsDNA_without_controls_EVEs_Events_ALL_ecto$coef)
Table_dsDNA_without_controls_EVEs_Events_ALL_ecto$coef <- as.numeric(exp(Table_dsDNA_without_controls_EVEs_Events_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_without_controls_EVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_without_controls_EVEs_Events_ALL_ecto$rope <- posterior_tab*100
Table_dsDNA_without_controls_EVEs_Events_ALL_ecto$med <- median(Table_dsDNA_without_controls_EVEs_Events_ALL_ecto$coef)
Table_dsDNA_without_controls_EVEs_Events_ALL_ecto$pd <- pd(Table_dsDNA_without_controls_EVEs_Events_ALL_ecto$coef2)



Table_dsDNA_without_controls_EVEs_Events_ALL_ecto_endo<-rbind(Table_dsDNA_without_controls_EVEs_Events_ALL_endo,Table_dsDNA_without_controls_EVEs_Events_ALL_ecto)

#Create summary table
statistics_Table_dsDNA_without_controls_EVEs_Events= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("dsDNA_without_controls_EVEs_Events",median(Table_dsDNA_without_controls_EVEs_Events_ALL_endo$coef),ci(Table_dsDNA_without_controls_EVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_without_controls_EVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_without_controls_EVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_dsDNA_without_controls_EVEs_Events_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("dsDNA_without_controls_EVEs_Events",median(Table_dsDNA_without_controls_EVEs_Events_ALL_ecto$coef),ci(Table_dsDNA_without_controls_EVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_without_controls_EVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_without_controls_EVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_dsDNA_without_controls_EVEs_Events_ALL_ecto$pd[1])))



colnames(statistics_Table_dsDNA_without_controls_EVEs_Events) <- statistics_Table_dsDNA_without_controls_EVEs_Events[1,]
statistics_Table_dsDNA_without_controls_EVEs_Events<-statistics_Table_dsDNA_without_controls_EVEs_Events[-1,]


Table_dsDNA_without_controls_EVEs_Events_ALL_ecto_endo$Parameter2 <- paste0(Table_dsDNA_without_controls_EVEs_Events_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_dsDNA_without_controls_EVEs_Events_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_dsDNA_without_controls_EVEs_Events_ALL<-Table_dsDNA_without_controls_EVEs_Events_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_dsDNA_without_controls_EVEs_Events_ALL <- stat.test_dsDNA_without_controls_EVEs_Events_ALL %>% add_xy_position(x = "Parameter2")
stat.test_dsDNA_without_controls_EVEs_Events_ALL$Prop <- paste0(round((nrow(Table_dsDNA_without_controls_EVEs_Events_ALL[Table_dsDNA_without_controls_EVEs_Events_ALL$b_value2 > Table_dsDNA_without_controls_EVEs_Events_ALL$b_value3,]) / nrow(Table_dsDNA_without_controls_EVEs_Events_ALL))*100,digits = 2),"%")

boxplot_dsDNA_without_controls_EVEs_Events <- ggplot(Table_dsDNA_without_controls_EVEs_Events_ALL_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="Events",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_dsDNA_without_controls_EVEs_Events_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,15))


##################################
## dEVEs Events all posteriors
##################################

Table_dsDNA_without_controls_dEVEs_Events_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_dsDNA_without_controls_All/Posteriors_dEVEs_Events_ALL.txt",sep=";",h=T)

Table_dsDNA_without_controls_dEVEs_Events_ALL$b_value2_coef_inf<- exp(as.numeric(Table_dsDNA_without_controls_dEVEs_Events_ALL$b_value2))
Table_dsDNA_without_controls_dEVEs_Events_ALL$b_value3_coef_inf<- exp(as.numeric(Table_dsDNA_without_controls_dEVEs_Events_ALL$b_value3))
Table_dsDNA_without_controls_dEVEs_Events_ALL<-Table_dsDNA_without_controls_dEVEs_Events_ALL[!Table_dsDNA_without_controls_dEVEs_Events_ALL$b_value2_coef_inf=="Inf",]
Table_dsDNA_without_controls_dEVEs_Events_ALL<-Table_dsDNA_without_controls_dEVEs_Events_ALL[!Table_dsDNA_without_controls_dEVEs_Events_ALL$b_value3_coef_inf=="Inf",]


Table_dsDNA_without_controls_dEVEs_Events_ALL_endo <- as.data.frame(Table_dsDNA_without_controls_dEVEs_Events_ALL$b_value2)
Table_dsDNA_without_controls_dEVEs_Events_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_dsDNA_without_controls_dEVEs_Events_ALL_endo)<- c("coef","Parameter")
Table_dsDNA_without_controls_dEVEs_Events_ALL_endo$coef2 <- as.numeric(Table_dsDNA_without_controls_dEVEs_Events_ALL_endo$coef)
Table_dsDNA_without_controls_dEVEs_Events_ALL_endo$coef <- as.numeric(exp(Table_dsDNA_without_controls_dEVEs_Events_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_without_controls_dEVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_without_controls_dEVEs_Events_ALL_endo$rope <- posterior_tab*100
Table_dsDNA_without_controls_dEVEs_Events_ALL_endo$med <- median(Table_dsDNA_without_controls_dEVEs_Events_ALL_endo$coef)
Table_dsDNA_without_controls_dEVEs_Events_ALL_endo$pd <- pd(Table_dsDNA_without_controls_dEVEs_Events_ALL_endo$coef2)

Table_dsDNA_without_controls_dEVEs_Events_ALL_ecto <- as.data.frame(Table_dsDNA_without_controls_dEVEs_Events_ALL$b_value3)
Table_dsDNA_without_controls_dEVEs_Events_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_dsDNA_without_controls_dEVEs_Events_ALL_ecto)<- c("coef","Parameter")
Table_dsDNA_without_controls_dEVEs_Events_ALL_ecto$coef2 <- as.numeric(Table_dsDNA_without_controls_dEVEs_Events_ALL_ecto$coef)
Table_dsDNA_without_controls_dEVEs_Events_ALL_ecto$coef <- as.numeric(exp(Table_dsDNA_without_controls_dEVEs_Events_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_without_controls_dEVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_without_controls_dEVEs_Events_ALL_ecto$rope <- posterior_tab*100
Table_dsDNA_without_controls_dEVEs_Events_ALL_ecto$med <- median(Table_dsDNA_without_controls_dEVEs_Events_ALL_ecto$coef)
Table_dsDNA_without_controls_dEVEs_Events_ALL_ecto$pd <- pd(Table_dsDNA_without_controls_dEVEs_Events_ALL_ecto$coef2)



Table_dsDNA_without_controls_dEVEs_Events_ALL_ecto_endo<-rbind(Table_dsDNA_without_controls_dEVEs_Events_ALL_endo,Table_dsDNA_without_controls_dEVEs_Events_ALL_ecto)

#Create summary table
statistics_Table_dsDNA_without_controls_dEVEs_Events= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("dsDNA_without_controls_dEVEs_Events",median(Table_dsDNA_without_controls_dEVEs_Events_ALL_endo$coef),ci(Table_dsDNA_without_controls_dEVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_without_controls_dEVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_without_controls_dEVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_dsDNA_without_controls_dEVEs_Events_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("dsDNA_without_controls_dEVEs_Events",median(Table_dsDNA_without_controls_dEVEs_Events_ALL_ecto$coef),ci(Table_dsDNA_without_controls_dEVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_without_controls_dEVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_without_controls_dEVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_dsDNA_without_controls_dEVEs_Events_ALL_ecto$pd[1])))



colnames(statistics_Table_dsDNA_without_controls_dEVEs_Events) <- statistics_Table_dsDNA_without_controls_dEVEs_Events[1,]
statistics_Table_dsDNA_without_controls_dEVEs_Events<-statistics_Table_dsDNA_without_controls_dEVEs_Events[-1,]


Table_dsDNA_without_controls_dEVEs_Events_ALL_ecto_endo$Parameter2 <- paste0(Table_dsDNA_without_controls_dEVEs_Events_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_dsDNA_without_controls_dEVEs_Events_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_dsDNA_without_controls_dEVEs_Events_ALL<-Table_dsDNA_without_controls_dEVEs_Events_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_dsDNA_without_controls_dEVEs_Events_ALL <- stat.test_dsDNA_without_controls_dEVEs_Events_ALL %>% add_xy_position(x = "Parameter2")
stat.test_dsDNA_without_controls_dEVEs_Events_ALL$Prop <- paste0(round((nrow(Table_dsDNA_without_controls_dEVEs_Events_ALL[Table_dsDNA_without_controls_dEVEs_Events_ALL$b_value2 > Table_dsDNA_without_controls_dEVEs_Events_ALL$b_value3,]) / nrow(Table_dsDNA_without_controls_dEVEs_Events_ALL))*100,digits = 2),"%")

boxplot_dsDNA_without_controls_dEVEs_Events <- ggplot(Table_dsDNA_without_controls_dEVEs_Events_ALL_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter),trim = T,adjust = 2)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="dEvents",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_dsDNA_without_controls_dEVEs_Events_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,15))

##################################
## EVEs Counts all posteriors####
##################################


Table_dsDNA_without_controls_EVEs_Counts_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_dsDNA_without_controls_All/Posteriors_EVEs_Counts_ALL.txt",sep=";",h=T)

Table_dsDNA_without_controls_EVEs_Counts_ALL$b_value2_coef_inf<- exp(as.numeric(Table_dsDNA_without_controls_EVEs_Counts_ALL$b_value2))
Table_dsDNA_without_controls_EVEs_Counts_ALL$b_value3_coef_inf<- exp(as.numeric(Table_dsDNA_without_controls_EVEs_Counts_ALL$b_value3))
Table_dsDNA_without_controls_EVEs_Counts_ALL<-Table_dsDNA_without_controls_EVEs_Counts_ALL[!Table_dsDNA_without_controls_EVEs_Counts_ALL$b_value2_coef_inf=="Inf",]
Table_dsDNA_without_controls_EVEs_Counts_ALL<-Table_dsDNA_without_controls_EVEs_Counts_ALL[!Table_dsDNA_without_controls_EVEs_Counts_ALL$b_value3_coef_inf=="Inf",]


Table_dsDNA_without_controls_EVEs_Counts_ALL_endo <- as.data.frame(Table_dsDNA_without_controls_EVEs_Counts_ALL$b_value2)
Table_dsDNA_without_controls_EVEs_Counts_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_dsDNA_without_controls_EVEs_Counts_ALL_endo)<- c("coef","Parameter")
Table_dsDNA_without_controls_EVEs_Counts_ALL_endo$coef2 <- as.numeric(Table_dsDNA_without_controls_EVEs_Counts_ALL_endo$coef)
Table_dsDNA_without_controls_EVEs_Counts_ALL_endo$coef <- as.numeric(exp(Table_dsDNA_without_controls_EVEs_Counts_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_without_controls_EVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_without_controls_EVEs_Counts_ALL_endo$rope <- posterior_tab*100
Table_dsDNA_without_controls_EVEs_Counts_ALL_endo$med <- median(Table_dsDNA_without_controls_EVEs_Counts_ALL_endo$coef)
Table_dsDNA_without_controls_EVEs_Counts_ALL_endo$pd <- pd(Table_dsDNA_without_controls_EVEs_Counts_ALL_endo$coef2)

Table_dsDNA_without_controls_EVEs_Counts_ALL_ecto <- as.data.frame(Table_dsDNA_without_controls_EVEs_Counts_ALL$b_value3)
Table_dsDNA_without_controls_EVEs_Counts_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_dsDNA_without_controls_EVEs_Counts_ALL_ecto)<- c("coef","Parameter")
Table_dsDNA_without_controls_EVEs_Counts_ALL_ecto$coef2 <- as.numeric(Table_dsDNA_without_controls_EVEs_Counts_ALL_ecto$coef)
Table_dsDNA_without_controls_EVEs_Counts_ALL_ecto$coef <- as.numeric(exp(Table_dsDNA_without_controls_EVEs_Counts_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_without_controls_EVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_without_controls_EVEs_Counts_ALL_ecto$rope <- posterior_tab*100
Table_dsDNA_without_controls_EVEs_Counts_ALL_ecto$med <- median(Table_dsDNA_without_controls_EVEs_Counts_ALL_ecto$coef)
Table_dsDNA_without_controls_EVEs_Counts_ALL_ecto$pd <- pd(Table_dsDNA_without_controls_EVEs_Counts_ALL_ecto$coef2)



Table_dsDNA_without_controls_EVEs_Counts_ALL_ecto_endo<-rbind(Table_dsDNA_without_controls_EVEs_Counts_ALL_endo,Table_dsDNA_without_controls_EVEs_Counts_ALL_ecto)

#Create summary table
statistics_Table_dsDNA_without_controls_EVEs_Counts= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("dsDNA_without_controls_EVEs_Counts",median(Table_dsDNA_without_controls_EVEs_Counts_ALL_endo$coef),ci(Table_dsDNA_without_controls_EVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_without_controls_EVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_without_controls_EVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_dsDNA_without_controls_EVEs_Counts_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("dsDNA_without_controls_EVEs_Counts",median(Table_dsDNA_without_controls_EVEs_Counts_ALL_ecto$coef),ci(Table_dsDNA_without_controls_EVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_without_controls_EVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_without_controls_EVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_dsDNA_without_controls_EVEs_Counts_ALL_ecto$pd[1])))



colnames(statistics_Table_dsDNA_without_controls_EVEs_Counts) <- statistics_Table_dsDNA_without_controls_EVEs_Counts[1,]
statistics_Table_dsDNA_without_controls_EVEs_Counts<-statistics_Table_dsDNA_without_controls_EVEs_Counts[-1,]


Table_dsDNA_without_controls_EVEs_Counts_ALL_ecto_endo$Parameter2 <- paste0(Table_dsDNA_without_controls_EVEs_Counts_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_dsDNA_without_controls_EVEs_Counts_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_dsDNA_without_controls_EVEs_Counts_ALL<-Table_dsDNA_without_controls_EVEs_Counts_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_dsDNA_without_controls_EVEs_Counts_ALL <- stat.test_dsDNA_without_controls_EVEs_Counts_ALL %>% add_xy_position(x = "Parameter2")
stat.test_dsDNA_without_controls_EVEs_Counts_ALL$Prop <- paste0(round((nrow(Table_dsDNA_without_controls_EVEs_Counts_ALL[Table_dsDNA_without_controls_EVEs_Counts_ALL$b_value2 > Table_dsDNA_without_controls_EVEs_Counts_ALL$b_value3,]) / nrow(Table_dsDNA_without_controls_EVEs_Counts_ALL))*100,digits = 2),"%")

boxplot_dsDNA_without_controls_EVEs_Counts <- ggplot(Table_dsDNA_without_controls_EVEs_Counts_ALL_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="EVEs",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_dsDNA_without_controls_EVEs_Counts_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,15))


################################
## dEVEs Counts all posteriors #
################################
Table_dsDNA_without_controls_dEVEs_Counts_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_dsDNA_without_controls_All/Posteriors_dEVEs_Counts_ALL.txt",sep=";",h=T)

Table_dsDNA_without_controls_dEVEs_Counts_ALL$b_value2_coef_inf<- exp(as.numeric(Table_dsDNA_without_controls_dEVEs_Counts_ALL$b_value2))
Table_dsDNA_without_controls_dEVEs_Counts_ALL$b_value3_coef_inf<- exp(as.numeric(Table_dsDNA_without_controls_dEVEs_Counts_ALL$b_value3))
Table_dsDNA_without_controls_dEVEs_Counts_ALL<-Table_dsDNA_without_controls_dEVEs_Counts_ALL[!Table_dsDNA_without_controls_dEVEs_Counts_ALL$b_value2_coef_inf=="Inf",]
Table_dsDNA_without_controls_dEVEs_Counts_ALL<-Table_dsDNA_without_controls_dEVEs_Counts_ALL[!Table_dsDNA_without_controls_dEVEs_Counts_ALL$b_value3_coef_inf=="Inf",]


Table_dsDNA_without_controls_dEVEs_Counts_ALL_endo <- as.data.frame(Table_dsDNA_without_controls_dEVEs_Counts_ALL$b_value2)
Table_dsDNA_without_controls_dEVEs_Counts_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_dsDNA_without_controls_dEVEs_Counts_ALL_endo)<- c("coef","Parameter")
Table_dsDNA_without_controls_dEVEs_Counts_ALL_endo$coef2 <- as.numeric(Table_dsDNA_without_controls_dEVEs_Counts_ALL_endo$coef)
Table_dsDNA_without_controls_dEVEs_Counts_ALL_endo$coef <- as.numeric(exp(Table_dsDNA_without_controls_dEVEs_Counts_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_without_controls_dEVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_without_controls_dEVEs_Counts_ALL_endo$rope <- posterior_tab*100
Table_dsDNA_without_controls_dEVEs_Counts_ALL_endo$med <- median(Table_dsDNA_without_controls_dEVEs_Counts_ALL_endo$coef)
Table_dsDNA_without_controls_dEVEs_Counts_ALL_endo$pd <- pd(Table_dsDNA_without_controls_dEVEs_Counts_ALL_endo$coef2)

Table_dsDNA_without_controls_dEVEs_Counts_ALL_ecto <- as.data.frame(Table_dsDNA_without_controls_dEVEs_Counts_ALL$b_value3)
Table_dsDNA_without_controls_dEVEs_Counts_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_dsDNA_without_controls_dEVEs_Counts_ALL_ecto)<- c("coef","Parameter")
Table_dsDNA_without_controls_dEVEs_Counts_ALL_ecto$coef2 <- as.numeric(Table_dsDNA_without_controls_dEVEs_Counts_ALL_ecto$coef)
Table_dsDNA_without_controls_dEVEs_Counts_ALL_ecto$coef <- as.numeric(exp(Table_dsDNA_without_controls_dEVEs_Counts_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_without_controls_dEVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_without_controls_dEVEs_Counts_ALL_ecto$rope <- posterior_tab*100
Table_dsDNA_without_controls_dEVEs_Counts_ALL_ecto$med <- median(Table_dsDNA_without_controls_dEVEs_Counts_ALL_ecto$coef)
Table_dsDNA_without_controls_dEVEs_Counts_ALL_ecto$pd <- pd(Table_dsDNA_without_controls_dEVEs_Counts_ALL_ecto$coef2)



Table_dsDNA_without_controls_dEVEs_Counts_ALL_ecto_endo<-rbind(Table_dsDNA_without_controls_dEVEs_Counts_ALL_endo,Table_dsDNA_without_controls_dEVEs_Counts_ALL_ecto)

#Create summary table
statistics_Table_dsDNA_without_controls_dEVEs_Counts= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("dsDNA_without_controls_dEVEs_Counts",median(Table_dsDNA_without_controls_dEVEs_Counts_ALL_endo$coef),ci(Table_dsDNA_without_controls_dEVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_without_controls_dEVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_without_controls_dEVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_dsDNA_without_controls_dEVEs_Counts_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("dsDNA_without_controls_dEVEs_Counts",median(Table_dsDNA_without_controls_dEVEs_Counts_ALL_ecto$coef),ci(Table_dsDNA_without_controls_dEVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_without_controls_dEVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_without_controls_dEVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_dsDNA_without_controls_dEVEs_Counts_ALL_ecto$pd[1])))



colnames(statistics_Table_dsDNA_without_controls_dEVEs_Counts) <- statistics_Table_dsDNA_without_controls_dEVEs_Counts[1,]
statistics_Table_dsDNA_without_controls_dEVEs_Counts<-statistics_Table_dsDNA_without_controls_dEVEs_Counts[-1,]


Table_dsDNA_without_controls_dEVEs_Counts_ALL_ecto_endo$Parameter2 <- paste0(Table_dsDNA_without_controls_dEVEs_Counts_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_dsDNA_without_controls_dEVEs_Counts_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_dsDNA_without_controls_dEVEs_Counts_ALL<-Table_dsDNA_without_controls_dEVEs_Counts_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_dsDNA_without_controls_dEVEs_Counts_ALL <- stat.test_dsDNA_without_controls_dEVEs_Counts_ALL %>% add_xy_position(x = "Parameter2")
stat.test_dsDNA_without_controls_dEVEs_Counts_ALL$Prop <- paste0(round((nrow(Table_dsDNA_without_controls_dEVEs_Counts_ALL[Table_dsDNA_without_controls_dEVEs_Counts_ALL$b_value2 > Table_dsDNA_without_controls_dEVEs_Counts_ALL$b_value3,]) / nrow(Table_dsDNA_without_controls_dEVEs_Counts_ALL))*100,digits = 2),"%")

boxplot_dsDNA_without_controls_dEVEs_Counts <- ggplot(Table_dsDNA_without_controls_dEVEs_Counts_ALL_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="dEVEs",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_dsDNA_without_controls_dEVEs_Counts_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,15))


## Combine All 4 plots

violin_plot_dsDNA_without_controls<-grid.arrange(arrangeGrob(
                             boxplot_dsDNA_without_controls_EVEs_Counts  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank()) ,
                             boxplot_dsDNA_without_controls_EVEs_Events  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank()),
                             boxplot_dsDNA_without_controls_dEVEs_Counts  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank()) ,
                             boxplot_dsDNA_without_controls_dEVEs_Events  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank()),
                         nrow = 1,
                         ncol=4,
                         left = textGrob("Relative rate", rot = 90, vjust = 1)),nrow=1,top=textGrob("dsDNA_no_controls | A-B-C-D",
                                                                                                                 gp=gpar(fontsize=20,font=8)))



#Add posterior coefficient comparison between endo and ecto
statistics_Table_dsDNA_without_controls_EVEs_Events<-as.data.frame(statistics_Table_dsDNA_without_controls_EVEs_Events)
statistics_Table_dsDNA_without_controls_EVEs_Events$Endo_sup_Ecto<-round((nrow(Table_dsDNA_without_controls_EVEs_Events_ALL[Table_dsDNA_without_controls_EVEs_Events_ALL$b_value2 > Table_dsDNA_without_controls_EVEs_Events_ALL$b_value3,]) / nrow(Table_dsDNA_without_controls_EVEs_Events_ALL))*100,digits = 2)

statistics_Table_dsDNA_without_controls_dEVEs_Events<-as.data.frame(statistics_Table_dsDNA_without_controls_dEVEs_Events)
statistics_Table_dsDNA_without_controls_dEVEs_Events$Endo_sup_Ecto<-round((nrow(Table_dsDNA_without_controls_dEVEs_Events_ALL[Table_dsDNA_without_controls_dEVEs_Events_ALL$b_value2 > Table_dsDNA_without_controls_dEVEs_Events_ALL$b_value3,]) / nrow(Table_dsDNA_without_controls_dEVEs_Events_ALL))*100,digits = 2)

statistics_Table_dsDNA_without_controls_EVEs_Counts<-as.data.frame(statistics_Table_dsDNA_without_controls_EVEs_Counts)
statistics_Table_dsDNA_without_controls_EVEs_Counts$Endo_sup_Ecto<-round((nrow(Table_dsDNA_without_controls_EVEs_Counts_ALL[Table_dsDNA_without_controls_EVEs_Counts_ALL$b_value2 > Table_dsDNA_without_controls_EVEs_Counts_ALL$b_value3,]) / nrow(Table_dsDNA_without_controls_EVEs_Counts_ALL))*100,digits = 2)


statistics_Table_dsDNA_without_controls_dEVEs_Counts<-as.data.frame(statistics_Table_dsDNA_without_controls_dEVEs_Counts)
statistics_Table_dsDNA_without_controls_dEVEs_Counts$Endo_sup_Ecto<-round((nrow(Table_dsDNA_without_controls_dEVEs_Counts_ALL[Table_dsDNA_without_controls_dEVEs_Counts_ALL$b_value2 > Table_dsDNA_without_controls_dEVEs_Counts_ALL$b_value3,]) / nrow(Table_dsDNA_without_controls_dEVEs_Counts_ALL))*100,digits = 2)


## Combine all statistics tables

statistics_Table_dsDNA_without_controls<-as.data.frame(rbind(statistics_Table_dsDNA_without_controls_EVEs_Events,statistics_Table_dsDNA_without_controls_dEVEs_Events,statistics_Table_dsDNA_without_controls_EVEs_Counts,statistics_Table_dsDNA_without_controls_dEVEs_Counts))

statistics_Table_dsDNA_without_controls$Lifestyle<- gsub("\\..*","",rownames(statistics_Table_dsDNA_without_controls))

library(dplyr)
statistics_Table_dsDNA_without_controls$Median<-as.numeric(statistics_Table_dsDNA_without_controls$Median)
statistics_Table_dsDNA_without_controls$CI_low<-as.numeric(statistics_Table_dsDNA_without_controls$CI_low)
statistics_Table_dsDNA_without_controls$CI_high<-as.numeric(statistics_Table_dsDNA_without_controls$CI_high)
statistics_Table_dsDNA_without_controls$ROPE_Percentage<-as.numeric(statistics_Table_dsDNA_without_controls$ROPE_Percentage)*100
statistics_Table_dsDNA_without_controls<-statistics_Table_dsDNA_without_controls %>% mutate_if(is.numeric, round, digits=3)

statistics_Table_dsDNA_without_controls$Type <-c("EVEs_Events","EVEs_Events","dEVEs_Events","dEVEs_Events","EVEs_Numbers","EVEs_Numbers","dEVEs_Numbers","dEVEs_Numbers")

table_dsDNA_without_controls<-statistics_Table_dsDNA_without_controls
table_dsDNA_without_controls[nrow(table_dsDNA_without_controls)+1,] <- NA



#dsDNA A without controls


library(bayestestR)
library(ggstatsplot)
library(bayesplot)
library(ggpubr)
library(gridExtra)
library(grid)



##################################
# ## EVEs Events all posteriors ##
##################################
Table_dsDNA_without_controls_A_EVEs_Events_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_dsDNA_without_controls_A_All/Posteriors_EVEs_Events_ALL.txt",sep=";",h=T)

Table_dsDNA_without_controls_A_EVEs_Events_ALL$b_value2_coef_inf<- exp(as.numeric(Table_dsDNA_without_controls_A_EVEs_Events_ALL$b_value2))
Table_dsDNA_without_controls_A_EVEs_Events_ALL$b_value3_coef_inf<- exp(as.numeric(Table_dsDNA_without_controls_A_EVEs_Events_ALL$b_value3))
Table_dsDNA_without_controls_A_EVEs_Events_ALL<-Table_dsDNA_without_controls_A_EVEs_Events_ALL[!Table_dsDNA_without_controls_A_EVEs_Events_ALL$b_value2_coef_inf=="Inf",]
Table_dsDNA_without_controls_A_EVEs_Events_ALL<-Table_dsDNA_without_controls_A_EVEs_Events_ALL[!Table_dsDNA_without_controls_A_EVEs_Events_ALL$b_value3_coef_inf=="Inf",]


Table_dsDNA_without_controls_A_EVEs_Events_ALL_endo <- as.data.frame(Table_dsDNA_without_controls_A_EVEs_Events_ALL$b_value2)
Table_dsDNA_without_controls_A_EVEs_Events_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_dsDNA_without_controls_A_EVEs_Events_ALL_endo)<- c("coef","Parameter")
Table_dsDNA_without_controls_A_EVEs_Events_ALL_endo$coef2 <- as.numeric(Table_dsDNA_without_controls_A_EVEs_Events_ALL_endo$coef)
Table_dsDNA_without_controls_A_EVEs_Events_ALL_endo$coef <- as.numeric(exp(Table_dsDNA_without_controls_A_EVEs_Events_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_without_controls_A_EVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_without_controls_A_EVEs_Events_ALL_endo$rope <- posterior_tab*100
Table_dsDNA_without_controls_A_EVEs_Events_ALL_endo$med <- median(Table_dsDNA_without_controls_A_EVEs_Events_ALL_endo$coef)
Table_dsDNA_without_controls_A_EVEs_Events_ALL_endo$pd <- pd(Table_dsDNA_without_controls_A_EVEs_Events_ALL_endo$coef2)

Table_dsDNA_without_controls_A_EVEs_Events_ALL_ecto <- as.data.frame(Table_dsDNA_without_controls_A_EVEs_Events_ALL$b_value3)
Table_dsDNA_without_controls_A_EVEs_Events_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_dsDNA_without_controls_A_EVEs_Events_ALL_ecto)<- c("coef","Parameter")
Table_dsDNA_without_controls_A_EVEs_Events_ALL_ecto$coef2 <- as.numeric(Table_dsDNA_without_controls_A_EVEs_Events_ALL_ecto$coef)
Table_dsDNA_without_controls_A_EVEs_Events_ALL_ecto$coef <- as.numeric(exp(Table_dsDNA_without_controls_A_EVEs_Events_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_without_controls_A_EVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_without_controls_A_EVEs_Events_ALL_ecto$rope <- posterior_tab*100
Table_dsDNA_without_controls_A_EVEs_Events_ALL_ecto$med <- median(Table_dsDNA_without_controls_A_EVEs_Events_ALL_ecto$coef)
Table_dsDNA_without_controls_A_EVEs_Events_ALL_ecto$pd <- pd(Table_dsDNA_without_controls_A_EVEs_Events_ALL_ecto$coef2)



Table_dsDNA_without_controls_A_EVEs_Events_ALL_ecto_endo<-rbind(Table_dsDNA_without_controls_A_EVEs_Events_ALL_endo,Table_dsDNA_without_controls_A_EVEs_Events_ALL_ecto)

#Create summary table
statistics_Table_dsDNA_without_controls_A_EVEs_Events= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("dsDNA_without_controls_A_EVEs_Events",median(Table_dsDNA_without_controls_A_EVEs_Events_ALL_endo$coef),ci(Table_dsDNA_without_controls_A_EVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_without_controls_A_EVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_without_controls_A_EVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_dsDNA_without_controls_A_EVEs_Events_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("dsDNA_without_controls_A_EVEs_Events",median(Table_dsDNA_without_controls_A_EVEs_Events_ALL_ecto$coef),ci(Table_dsDNA_without_controls_A_EVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_without_controls_A_EVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_without_controls_A_EVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_dsDNA_without_controls_A_EVEs_Events_ALL_ecto$pd[1])))



colnames(statistics_Table_dsDNA_without_controls_A_EVEs_Events) <- statistics_Table_dsDNA_without_controls_A_EVEs_Events[1,]
statistics_Table_dsDNA_without_controls_A_EVEs_Events<-statistics_Table_dsDNA_without_controls_A_EVEs_Events[-1,]


Table_dsDNA_without_controls_A_EVEs_Events_ALL_ecto_endo$Parameter2 <- paste0(Table_dsDNA_without_controls_A_EVEs_Events_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_dsDNA_without_controls_A_EVEs_Events_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_dsDNA_without_controls_A_EVEs_Events_ALL<-Table_dsDNA_without_controls_A_EVEs_Events_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_dsDNA_without_controls_A_EVEs_Events_ALL <- stat.test_dsDNA_without_controls_A_EVEs_Events_ALL %>% add_xy_position(x = "Parameter2")
stat.test_dsDNA_without_controls_A_EVEs_Events_ALL$Prop <- paste0(round((nrow(Table_dsDNA_without_controls_A_EVEs_Events_ALL[Table_dsDNA_without_controls_A_EVEs_Events_ALL$b_value2 > Table_dsDNA_without_controls_A_EVEs_Events_ALL$b_value3,]) / nrow(Table_dsDNA_without_controls_A_EVEs_Events_ALL))*100,digits = 2),"%")

boxplot_dsDNA_without_controls_A_EVEs_Events <- ggplot(Table_dsDNA_without_controls_A_EVEs_Events_ALL_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="Events",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_dsDNA_without_controls_A_EVEs_Events_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,15))


##################################
## dEVEs Events all posteriors
##################################

Table_dsDNA_without_controls_A_dEVEs_Events_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_dsDNA_without_controls_A_All/Posteriors_dEVEs_Events_ALL.txt",sep=";",h=T)

Table_dsDNA_without_controls_A_dEVEs_Events_ALL$b_value2_coef_inf<- exp(as.numeric(Table_dsDNA_without_controls_A_dEVEs_Events_ALL$b_value2))
Table_dsDNA_without_controls_A_dEVEs_Events_ALL$b_value3_coef_inf<- exp(as.numeric(Table_dsDNA_without_controls_A_dEVEs_Events_ALL$b_value3))
Table_dsDNA_without_controls_A_dEVEs_Events_ALL<-Table_dsDNA_without_controls_A_dEVEs_Events_ALL[!Table_dsDNA_without_controls_A_dEVEs_Events_ALL$b_value2_coef_inf=="Inf",]
Table_dsDNA_without_controls_A_dEVEs_Events_ALL<-Table_dsDNA_without_controls_A_dEVEs_Events_ALL[!Table_dsDNA_without_controls_A_dEVEs_Events_ALL$b_value3_coef_inf=="Inf",]


Table_dsDNA_without_controls_A_dEVEs_Events_ALL_endo <- as.data.frame(Table_dsDNA_without_controls_A_dEVEs_Events_ALL$b_value2)
Table_dsDNA_without_controls_A_dEVEs_Events_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_dsDNA_without_controls_A_dEVEs_Events_ALL_endo)<- c("coef","Parameter")
Table_dsDNA_without_controls_A_dEVEs_Events_ALL_endo$coef2 <- as.numeric(Table_dsDNA_without_controls_A_dEVEs_Events_ALL_endo$coef)
Table_dsDNA_without_controls_A_dEVEs_Events_ALL_endo$coef <- as.numeric(exp(Table_dsDNA_without_controls_A_dEVEs_Events_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_without_controls_A_dEVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_without_controls_A_dEVEs_Events_ALL_endo$rope <- posterior_tab*100
Table_dsDNA_without_controls_A_dEVEs_Events_ALL_endo$med <- median(Table_dsDNA_without_controls_A_dEVEs_Events_ALL_endo$coef)
Table_dsDNA_without_controls_A_dEVEs_Events_ALL_endo$pd <- pd(Table_dsDNA_without_controls_A_dEVEs_Events_ALL_endo$coef2)

Table_dsDNA_without_controls_A_dEVEs_Events_ALL_ecto <- as.data.frame(Table_dsDNA_without_controls_A_dEVEs_Events_ALL$b_value3)
Table_dsDNA_without_controls_A_dEVEs_Events_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_dsDNA_without_controls_A_dEVEs_Events_ALL_ecto)<- c("coef","Parameter")
Table_dsDNA_without_controls_A_dEVEs_Events_ALL_ecto$coef2 <- as.numeric(Table_dsDNA_without_controls_A_dEVEs_Events_ALL_ecto$coef)
Table_dsDNA_without_controls_A_dEVEs_Events_ALL_ecto$coef <- as.numeric(exp(Table_dsDNA_without_controls_A_dEVEs_Events_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_without_controls_A_dEVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_without_controls_A_dEVEs_Events_ALL_ecto$rope <- posterior_tab*100
Table_dsDNA_without_controls_A_dEVEs_Events_ALL_ecto$med <- median(Table_dsDNA_without_controls_A_dEVEs_Events_ALL_ecto$coef)
Table_dsDNA_without_controls_A_dEVEs_Events_ALL_ecto$pd <- pd(Table_dsDNA_without_controls_A_dEVEs_Events_ALL_ecto$coef2)



Table_dsDNA_without_controls_A_dEVEs_Events_ALL_ecto_endo<-rbind(Table_dsDNA_without_controls_A_dEVEs_Events_ALL_endo,Table_dsDNA_without_controls_A_dEVEs_Events_ALL_ecto)

#Create summary table
statistics_Table_dsDNA_without_controls_A_dEVEs_Events= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("dsDNA_without_controls_A_dEVEs_Events",median(Table_dsDNA_without_controls_A_dEVEs_Events_ALL_endo$coef),ci(Table_dsDNA_without_controls_A_dEVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_without_controls_A_dEVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_without_controls_A_dEVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_dsDNA_without_controls_A_dEVEs_Events_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("dsDNA_without_controls_A_dEVEs_Events",median(Table_dsDNA_without_controls_A_dEVEs_Events_ALL_ecto$coef),ci(Table_dsDNA_without_controls_A_dEVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_without_controls_A_dEVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_without_controls_A_dEVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_dsDNA_without_controls_A_dEVEs_Events_ALL_ecto$pd[1])))



colnames(statistics_Table_dsDNA_without_controls_A_dEVEs_Events) <- statistics_Table_dsDNA_without_controls_A_dEVEs_Events[1,]
statistics_Table_dsDNA_without_controls_A_dEVEs_Events<-statistics_Table_dsDNA_without_controls_A_dEVEs_Events[-1,]


Table_dsDNA_without_controls_A_dEVEs_Events_ALL_ecto_endo$Parameter2 <- paste0(Table_dsDNA_without_controls_A_dEVEs_Events_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_dsDNA_without_controls_A_dEVEs_Events_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_dsDNA_without_controls_A_dEVEs_Events_ALL<-Table_dsDNA_without_controls_A_dEVEs_Events_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_dsDNA_without_controls_A_dEVEs_Events_ALL <- stat.test_dsDNA_without_controls_A_dEVEs_Events_ALL %>% add_xy_position(x = "Parameter2")
stat.test_dsDNA_without_controls_A_dEVEs_Events_ALL$Prop <- paste0(round((nrow(Table_dsDNA_without_controls_A_dEVEs_Events_ALL[Table_dsDNA_without_controls_A_dEVEs_Events_ALL$b_value2 > Table_dsDNA_without_controls_A_dEVEs_Events_ALL$b_value3,]) / nrow(Table_dsDNA_without_controls_A_dEVEs_Events_ALL))*100,digits = 2),"%")

boxplot_dsDNA_without_controls_A_dEVEs_Events <- ggplot(Table_dsDNA_without_controls_A_dEVEs_Events_ALL_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter),trim = T,adjust = 2)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="dEvents",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_dsDNA_without_controls_A_dEVEs_Events_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,15))

##################################
## EVEs Counts all posteriors####
##################################


Table_dsDNA_without_controls_A_EVEs_Counts_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_dsDNA_without_controls_A_All/Posteriors_EVEs_Counts_ALL.txt",sep=";",h=T)

Table_dsDNA_without_controls_A_EVEs_Counts_ALL$b_value2_coef_inf<- exp(as.numeric(Table_dsDNA_without_controls_A_EVEs_Counts_ALL$b_value2))
Table_dsDNA_without_controls_A_EVEs_Counts_ALL$b_value3_coef_inf<- exp(as.numeric(Table_dsDNA_without_controls_A_EVEs_Counts_ALL$b_value3))
Table_dsDNA_without_controls_A_EVEs_Counts_ALL<-Table_dsDNA_without_controls_A_EVEs_Counts_ALL[!Table_dsDNA_without_controls_A_EVEs_Counts_ALL$b_value2_coef_inf=="Inf",]
Table_dsDNA_without_controls_A_EVEs_Counts_ALL<-Table_dsDNA_without_controls_A_EVEs_Counts_ALL[!Table_dsDNA_without_controls_A_EVEs_Counts_ALL$b_value3_coef_inf=="Inf",]


Table_dsDNA_without_controls_A_EVEs_Counts_ALL_endo <- as.data.frame(Table_dsDNA_without_controls_A_EVEs_Counts_ALL$b_value2)
Table_dsDNA_without_controls_A_EVEs_Counts_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_dsDNA_without_controls_A_EVEs_Counts_ALL_endo)<- c("coef","Parameter")
Table_dsDNA_without_controls_A_EVEs_Counts_ALL_endo$coef2 <- as.numeric(Table_dsDNA_without_controls_A_EVEs_Counts_ALL_endo$coef)
Table_dsDNA_without_controls_A_EVEs_Counts_ALL_endo$coef <- as.numeric(exp(Table_dsDNA_without_controls_A_EVEs_Counts_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_without_controls_A_EVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_without_controls_A_EVEs_Counts_ALL_endo$rope <- posterior_tab*100
Table_dsDNA_without_controls_A_EVEs_Counts_ALL_endo$med <- median(Table_dsDNA_without_controls_A_EVEs_Counts_ALL_endo$coef)
Table_dsDNA_without_controls_A_EVEs_Counts_ALL_endo$pd <- pd(Table_dsDNA_without_controls_A_EVEs_Counts_ALL_endo$coef2)

Table_dsDNA_without_controls_A_EVEs_Counts_ALL_ecto <- as.data.frame(Table_dsDNA_without_controls_A_EVEs_Counts_ALL$b_value3)
Table_dsDNA_without_controls_A_EVEs_Counts_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_dsDNA_without_controls_A_EVEs_Counts_ALL_ecto)<- c("coef","Parameter")
Table_dsDNA_without_controls_A_EVEs_Counts_ALL_ecto$coef2 <- as.numeric(Table_dsDNA_without_controls_A_EVEs_Counts_ALL_ecto$coef)
Table_dsDNA_without_controls_A_EVEs_Counts_ALL_ecto$coef <- as.numeric(exp(Table_dsDNA_without_controls_A_EVEs_Counts_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_without_controls_A_EVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_without_controls_A_EVEs_Counts_ALL_ecto$rope <- posterior_tab*100
Table_dsDNA_without_controls_A_EVEs_Counts_ALL_ecto$med <- median(Table_dsDNA_without_controls_A_EVEs_Counts_ALL_ecto$coef)
Table_dsDNA_without_controls_A_EVEs_Counts_ALL_ecto$pd <- pd(Table_dsDNA_without_controls_A_EVEs_Counts_ALL_ecto$coef2)



Table_dsDNA_without_controls_A_EVEs_Counts_ALL_ecto_endo<-rbind(Table_dsDNA_without_controls_A_EVEs_Counts_ALL_endo,Table_dsDNA_without_controls_A_EVEs_Counts_ALL_ecto)

#Create summary table
statistics_Table_dsDNA_without_controls_A_EVEs_Counts= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("dsDNA_without_controls_A_EVEs_Counts",median(Table_dsDNA_without_controls_A_EVEs_Counts_ALL_endo$coef),ci(Table_dsDNA_without_controls_A_EVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_without_controls_A_EVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_without_controls_A_EVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_dsDNA_without_controls_A_EVEs_Counts_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("dsDNA_without_controls_A_EVEs_Counts",median(Table_dsDNA_without_controls_A_EVEs_Counts_ALL_ecto$coef),ci(Table_dsDNA_without_controls_A_EVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_without_controls_A_EVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_without_controls_A_EVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_dsDNA_without_controls_A_EVEs_Counts_ALL_ecto$pd[1])))



colnames(statistics_Table_dsDNA_without_controls_A_EVEs_Counts) <- statistics_Table_dsDNA_without_controls_A_EVEs_Counts[1,]
statistics_Table_dsDNA_without_controls_A_EVEs_Counts<-statistics_Table_dsDNA_without_controls_A_EVEs_Counts[-1,]


Table_dsDNA_without_controls_A_EVEs_Counts_ALL_ecto_endo$Parameter2 <- paste0(Table_dsDNA_without_controls_A_EVEs_Counts_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_dsDNA_without_controls_A_EVEs_Counts_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_dsDNA_without_controls_A_EVEs_Counts_ALL<-Table_dsDNA_without_controls_A_EVEs_Counts_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_dsDNA_without_controls_A_EVEs_Counts_ALL <- stat.test_dsDNA_without_controls_A_EVEs_Counts_ALL %>% add_xy_position(x = "Parameter2")
stat.test_dsDNA_without_controls_A_EVEs_Counts_ALL$Prop <- paste0(round((nrow(Table_dsDNA_without_controls_A_EVEs_Counts_ALL[Table_dsDNA_without_controls_A_EVEs_Counts_ALL$b_value2 > Table_dsDNA_without_controls_A_EVEs_Counts_ALL$b_value3,]) / nrow(Table_dsDNA_without_controls_A_EVEs_Counts_ALL))*100,digits = 2),"%")

boxplot_dsDNA_without_controls_A_EVEs_Counts <- ggplot(Table_dsDNA_without_controls_A_EVEs_Counts_ALL_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="EVEs",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_dsDNA_without_controls_A_EVEs_Counts_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,15))


################################
## dEVEs Counts all posteriors #
################################
Table_dsDNA_without_controls_A_dEVEs_Counts_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_dsDNA_without_controls_A_All/Posteriors_dEVEs_Counts_ALL.txt",sep=";",h=T)

Table_dsDNA_without_controls_A_dEVEs_Counts_ALL$b_value2_coef_inf<- exp(as.numeric(Table_dsDNA_without_controls_A_dEVEs_Counts_ALL$b_value2))
Table_dsDNA_without_controls_A_dEVEs_Counts_ALL$b_value3_coef_inf<- exp(as.numeric(Table_dsDNA_without_controls_A_dEVEs_Counts_ALL$b_value3))
Table_dsDNA_without_controls_A_dEVEs_Counts_ALL<-Table_dsDNA_without_controls_A_dEVEs_Counts_ALL[!Table_dsDNA_without_controls_A_dEVEs_Counts_ALL$b_value2_coef_inf=="Inf",]
Table_dsDNA_without_controls_A_dEVEs_Counts_ALL<-Table_dsDNA_without_controls_A_dEVEs_Counts_ALL[!Table_dsDNA_without_controls_A_dEVEs_Counts_ALL$b_value3_coef_inf=="Inf",]


Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_endo <- as.data.frame(Table_dsDNA_without_controls_A_dEVEs_Counts_ALL$b_value2)
Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_endo)<- c("coef","Parameter")
Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_endo$coef2 <- as.numeric(Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_endo$coef)
Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_endo$coef <- as.numeric(exp(Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_endo$rope <- posterior_tab*100
Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_endo$med <- median(Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_endo$coef)
Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_endo$pd <- pd(Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_endo$coef2)

Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_ecto <- as.data.frame(Table_dsDNA_without_controls_A_dEVEs_Counts_ALL$b_value3)
Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_ecto)<- c("coef","Parameter")
Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_ecto$coef2 <- as.numeric(Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_ecto$coef)
Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_ecto$coef <- as.numeric(exp(Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_ecto$rope <- posterior_tab*100
Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_ecto$med <- median(Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_ecto$coef)
Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_ecto$pd <- pd(Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_ecto$coef2)



Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_ecto_endo<-rbind(Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_endo,Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_ecto)

#Create summary table
statistics_Table_dsDNA_without_controls_A_dEVEs_Counts= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("dsDNA_without_controls_A_dEVEs_Counts",median(Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_endo$coef),ci(Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("dsDNA_without_controls_A_dEVEs_Counts",median(Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_ecto$coef),ci(Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_ecto$pd[1])))



colnames(statistics_Table_dsDNA_without_controls_A_dEVEs_Counts) <- statistics_Table_dsDNA_without_controls_A_dEVEs_Counts[1,]
statistics_Table_dsDNA_without_controls_A_dEVEs_Counts<-statistics_Table_dsDNA_without_controls_A_dEVEs_Counts[-1,]


Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_ecto_endo$Parameter2 <- paste0(Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_dsDNA_without_controls_A_dEVEs_Counts_ALL<-Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_dsDNA_without_controls_A_dEVEs_Counts_ALL <- stat.test_dsDNA_without_controls_A_dEVEs_Counts_ALL %>% add_xy_position(x = "Parameter2")
stat.test_dsDNA_without_controls_A_dEVEs_Counts_ALL$Prop <- paste0(round((nrow(Table_dsDNA_without_controls_A_dEVEs_Counts_ALL[Table_dsDNA_without_controls_A_dEVEs_Counts_ALL$b_value2 > Table_dsDNA_without_controls_A_dEVEs_Counts_ALL$b_value3,]) / nrow(Table_dsDNA_without_controls_A_dEVEs_Counts_ALL))*100,digits = 2),"%")

boxplot_dsDNA_without_controls_A_dEVEs_Counts <- ggplot(Table_dsDNA_without_controls_A_dEVEs_Counts_ALL_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="dEVEs",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_dsDNA_without_controls_A_dEVEs_Counts_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,15))




## Combine All 4 plots

violin_plot_dsDNA_without_controls_A<-grid.arrange(arrangeGrob(
                             boxplot_dsDNA_without_controls_A_EVEs_Counts  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank()) ,
                             boxplot_dsDNA_without_controls_A_EVEs_Events  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank()),
                             boxplot_dsDNA_without_controls_A_dEVEs_Counts  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank()) ,
                             boxplot_dsDNA_without_controls_A_dEVEs_Events  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank()),
                         nrow = 1,
                         ncol=4,
                         left = textGrob("Relative rate", rot = 90, vjust = 1)),nrow=1,top=textGrob("dsDNA_no_controls | A",
                                                                                                                 gp=gpar(fontsize=20,font=8)))






#Add posterior coefficient comparison between endo and ecto
statistics_Table_dsDNA_without_controls_A_EVEs_Events<-as.data.frame(statistics_Table_dsDNA_without_controls_A_EVEs_Events)
statistics_Table_dsDNA_without_controls_A_EVEs_Events$Endo_sup_Ecto<-round((nrow(Table_dsDNA_without_controls_A_EVEs_Events_ALL[Table_dsDNA_without_controls_A_EVEs_Events_ALL$b_value2 > Table_dsDNA_without_controls_A_EVEs_Events_ALL$b_value3,]) / nrow(Table_dsDNA_without_controls_A_EVEs_Events_ALL))*100,digits = 2)

statistics_Table_dsDNA_without_controls_A_dEVEs_Events<-as.data.frame(statistics_Table_dsDNA_without_controls_A_dEVEs_Events)
statistics_Table_dsDNA_without_controls_A_dEVEs_Events$Endo_sup_Ecto<-round((nrow(Table_dsDNA_without_controls_A_dEVEs_Events_ALL[Table_dsDNA_without_controls_A_dEVEs_Events_ALL$b_value2 > Table_dsDNA_without_controls_A_dEVEs_Events_ALL$b_value3,]) / nrow(Table_dsDNA_without_controls_A_dEVEs_Events_ALL))*100,digits = 2)

statistics_Table_dsDNA_without_controls_A_EVEs_Counts<-as.data.frame(statistics_Table_dsDNA_without_controls_A_EVEs_Counts)
statistics_Table_dsDNA_without_controls_A_EVEs_Counts$Endo_sup_Ecto<-round((nrow(Table_dsDNA_without_controls_A_EVEs_Counts_ALL[Table_dsDNA_without_controls_A_EVEs_Counts_ALL$b_value2 > Table_dsDNA_without_controls_A_EVEs_Counts_ALL$b_value3,]) / nrow(Table_dsDNA_without_controls_A_EVEs_Counts_ALL))*100,digits = 2)


statistics_Table_dsDNA_without_controls_A_dEVEs_Counts<-as.data.frame(statistics_Table_dsDNA_without_controls_A_dEVEs_Counts)
statistics_Table_dsDNA_without_controls_A_dEVEs_Counts$Endo_sup_Ecto<-round((nrow(Table_dsDNA_without_controls_A_dEVEs_Counts_ALL[Table_dsDNA_without_controls_A_dEVEs_Counts_ALL$b_value2 > Table_dsDNA_without_controls_A_dEVEs_Counts_ALL$b_value3,]) / nrow(Table_dsDNA_without_controls_A_dEVEs_Counts_ALL))*100,digits = 2)


## Combine all statistics tables

statistics_Table_dsDNA_without_controls_A<-as.data.frame(rbind(statistics_Table_dsDNA_without_controls_A_EVEs_Events,statistics_Table_dsDNA_without_controls_A_dEVEs_Events,statistics_Table_dsDNA_without_controls_A_EVEs_Counts,statistics_Table_dsDNA_without_controls_A_dEVEs_Counts))

statistics_Table_dsDNA_without_controls_A$Lifestyle<- gsub("\\..*","",rownames(statistics_Table_dsDNA_without_controls_A))

library(dplyr)
statistics_Table_dsDNA_without_controls_A$Median<-as.numeric(statistics_Table_dsDNA_without_controls_A$Median)
statistics_Table_dsDNA_without_controls_A$CI_low<-as.numeric(statistics_Table_dsDNA_without_controls_A$CI_low)
statistics_Table_dsDNA_without_controls_A$CI_high<-as.numeric(statistics_Table_dsDNA_without_controls_A$CI_high)
statistics_Table_dsDNA_without_controls_A$ROPE_Percentage<-as.numeric(statistics_Table_dsDNA_without_controls_A$ROPE_Percentage)*100
statistics_Table_dsDNA_without_controls_A<-statistics_Table_dsDNA_without_controls_A %>% mutate_if(is.numeric, round, digits=3)

statistics_Table_dsDNA_without_controls_A$Type <-c("EVEs_Events","EVEs_Events","dEVEs_Events","dEVEs_Events","EVEs_Numbers","EVEs_Numbers","dEVEs_Numbers","dEVEs_Numbers")

table_dsDNA_without_controls_A<-statistics_Table_dsDNA_without_controls_A
table_dsDNA_without_controls_A[nrow(table_dsDNA_without_controls_A)+1,] <- NA




#dsRNA A-B-C-D



library(bayestestR)
library(ggstatsplot)
library(bayesplot)
library(ggpubr)
library(gridExtra)
library(grid)


##################################
# ## EVEs Events all posteriors ##
##################################
Table_dsRNA_EVEs_Events_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_dsRNA_All/Posteriors_EVEs_Events_ALL.txt",sep=";",h=T)

Table_dsRNA_EVEs_Events_ALL$b_value2_coef_inf<- exp(as.numeric(Table_dsRNA_EVEs_Events_ALL$b_value2))
Table_dsRNA_EVEs_Events_ALL$b_value3_coef_inf<- exp(as.numeric(Table_dsRNA_EVEs_Events_ALL$b_value3))
Table_dsRNA_EVEs_Events_ALL<-Table_dsRNA_EVEs_Events_ALL[!Table_dsRNA_EVEs_Events_ALL$b_value2_coef_inf=="Inf",]
Table_dsRNA_EVEs_Events_ALL<-Table_dsRNA_EVEs_Events_ALL[!Table_dsRNA_EVEs_Events_ALL$b_value3_coef_inf=="Inf",]


Table_dsRNA_EVEs_Events_ALL_endo <- as.data.frame(Table_dsRNA_EVEs_Events_ALL$b_value2)
Table_dsRNA_EVEs_Events_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_dsRNA_EVEs_Events_ALL_endo)<- c("coef","Parameter")
Table_dsRNA_EVEs_Events_ALL_endo$coef2 <- as.numeric(Table_dsRNA_EVEs_Events_ALL_endo$coef)
Table_dsRNA_EVEs_Events_ALL_endo$coef <- as.numeric(exp(Table_dsRNA_EVEs_Events_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_dsRNA_EVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsRNA_EVEs_Events_ALL_endo$rope <- posterior_tab*100
Table_dsRNA_EVEs_Events_ALL_endo$med <- median(Table_dsRNA_EVEs_Events_ALL_endo$coef)
Table_dsRNA_EVEs_Events_ALL_endo$pd <- pd(Table_dsRNA_EVEs_Events_ALL_endo$coef2)

Table_dsRNA_EVEs_Events_ALL_ecto <- as.data.frame(Table_dsRNA_EVEs_Events_ALL$b_value3)
Table_dsRNA_EVEs_Events_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_dsRNA_EVEs_Events_ALL_ecto)<- c("coef","Parameter")
Table_dsRNA_EVEs_Events_ALL_ecto$coef2 <- as.numeric(Table_dsRNA_EVEs_Events_ALL_ecto$coef)
Table_dsRNA_EVEs_Events_ALL_ecto$coef <- as.numeric(exp(Table_dsRNA_EVEs_Events_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_dsRNA_EVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsRNA_EVEs_Events_ALL_ecto$rope <- posterior_tab*100
Table_dsRNA_EVEs_Events_ALL_ecto$med <- median(Table_dsRNA_EVEs_Events_ALL_ecto$coef)
Table_dsRNA_EVEs_Events_ALL_ecto$pd <- pd(Table_dsRNA_EVEs_Events_ALL_ecto$coef2)



Table_dsRNA_EVEs_Events_ALL_ecto_endo<-rbind(Table_dsRNA_EVEs_Events_ALL_endo,Table_dsRNA_EVEs_Events_ALL_ecto)

#Create summary table
statistics_Table_dsRNA_EVEs_Events= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("dsRNA_A_EVEs_Events",median(Table_dsRNA_EVEs_Events_ALL_endo$coef),ci(Table_dsRNA_EVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsRNA_EVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsRNA_EVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_dsRNA_EVEs_Events_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("dsRNA_A_EVEs_Events",median(Table_dsRNA_EVEs_Events_ALL_ecto$coef),ci(Table_dsRNA_EVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsRNA_EVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsRNA_EVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_dsRNA_EVEs_Events_ALL_ecto$pd[1])))



colnames(statistics_Table_dsRNA_EVEs_Events) <- statistics_Table_dsRNA_EVEs_Events[1,]
statistics_Table_dsRNA_EVEs_Events<-statistics_Table_dsRNA_EVEs_Events[-1,]


Table_dsRNA_EVEs_Events_ALL_ecto_endo$Parameter2 <- paste0(Table_dsRNA_EVEs_Events_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_dsRNA_EVEs_Events_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_dsRNA_EVEs_Events_ALL<-Table_dsRNA_EVEs_Events_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_dsRNA_EVEs_Events_ALL <- stat.test_dsRNA_EVEs_Events_ALL %>% add_xy_position(x = "Parameter2")
stat.test_dsRNA_EVEs_Events_ALL$Prop <- paste0(round((nrow(Table_dsRNA_EVEs_Events_ALL[Table_dsRNA_EVEs_Events_ALL$b_value2 > Table_dsRNA_EVEs_Events_ALL$b_value3,]) / nrow(Table_dsRNA_EVEs_Events_ALL))*100,digits = 2),"%")

boxplot_dsRNA_EVEs_Events <- ggplot(Table_dsRNA_EVEs_Events_ALL_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="Events",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_dsRNA_EVEs_Events_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,50))


##################################
## dEVEs Events all posteriors
##################################

Table_dsRNA_dEVEs_Events_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_dsRNA_All/Posteriors_dEVEs_Events_ALL.txt",sep=";",h=T)

Table_dsRNA_dEVEs_Events_ALL$b_value2_coef_inf<- exp(as.numeric(Table_dsRNA_dEVEs_Events_ALL$b_value2))
Table_dsRNA_dEVEs_Events_ALL$b_value3_coef_inf<- exp(as.numeric(Table_dsRNA_dEVEs_Events_ALL$b_value3))
Table_dsRNA_dEVEs_Events_ALL<-Table_dsRNA_dEVEs_Events_ALL[!Table_dsRNA_dEVEs_Events_ALL$b_value2_coef_inf=="Inf",]
Table_dsRNA_dEVEs_Events_ALL<-Table_dsRNA_dEVEs_Events_ALL[!Table_dsRNA_dEVEs_Events_ALL$b_value3_coef_inf=="Inf",]


Table_dsRNA_dEVEs_Events_ALL_endo <- as.data.frame(Table_dsRNA_dEVEs_Events_ALL$b_value2)
Table_dsRNA_dEVEs_Events_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_dsRNA_dEVEs_Events_ALL_endo)<- c("coef","Parameter")
Table_dsRNA_dEVEs_Events_ALL_endo$coef2 <- as.numeric(Table_dsRNA_dEVEs_Events_ALL_endo$coef)
Table_dsRNA_dEVEs_Events_ALL_endo$coef <- as.numeric(exp(Table_dsRNA_dEVEs_Events_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_dsRNA_dEVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsRNA_dEVEs_Events_ALL_endo$rope <- posterior_tab*100
Table_dsRNA_dEVEs_Events_ALL_endo$med <- median(Table_dsRNA_dEVEs_Events_ALL_endo$coef)
Table_dsRNA_dEVEs_Events_ALL_endo$pd <- pd(Table_dsRNA_dEVEs_Events_ALL_endo$coef2)

Table_dsRNA_dEVEs_Events_ALL_ecto <- as.data.frame(Table_dsRNA_dEVEs_Events_ALL$b_value3)
Table_dsRNA_dEVEs_Events_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_dsRNA_dEVEs_Events_ALL_ecto)<- c("coef","Parameter")
Table_dsRNA_dEVEs_Events_ALL_ecto$coef2 <- as.numeric(Table_dsRNA_dEVEs_Events_ALL_ecto$coef)
Table_dsRNA_dEVEs_Events_ALL_ecto$coef <- as.numeric(exp(Table_dsRNA_dEVEs_Events_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_dsRNA_dEVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsRNA_dEVEs_Events_ALL_ecto$rope <- posterior_tab*100
Table_dsRNA_dEVEs_Events_ALL_ecto$med <- median(Table_dsRNA_dEVEs_Events_ALL_ecto$coef)
Table_dsRNA_dEVEs_Events_ALL_ecto$pd <- pd(Table_dsRNA_dEVEs_Events_ALL_ecto$coef2)



Table_dsRNA_dEVEs_Events_ALL_ecto_endo<-rbind(Table_dsRNA_dEVEs_Events_ALL_endo,Table_dsRNA_dEVEs_Events_ALL_ecto)

#Create summary table
statistics_Table_dsRNA_dEVEs_Events= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("dsRNA_A_dEVEs_Events",median(Table_dsRNA_dEVEs_Events_ALL_endo$coef),ci(Table_dsRNA_dEVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsRNA_dEVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsRNA_dEVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_dsRNA_dEVEs_Events_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("dsRNA_A_dEVEs_Events",median(Table_dsRNA_dEVEs_Events_ALL_ecto$coef),ci(Table_dsRNA_dEVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsRNA_dEVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsRNA_dEVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_dsRNA_dEVEs_Events_ALL_ecto$pd[1])))



colnames(statistics_Table_dsRNA_dEVEs_Events) <- statistics_Table_dsRNA_dEVEs_Events[1,]
statistics_Table_dsRNA_dEVEs_Events<-statistics_Table_dsRNA_dEVEs_Events[-1,]


Table_dsRNA_dEVEs_Events_ALL_ecto_endo$Parameter2 <- paste0(Table_dsRNA_dEVEs_Events_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_dsRNA_dEVEs_Events_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_dsRNA_dEVEs_Events_ALL<-Table_dsRNA_dEVEs_Events_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_dsRNA_dEVEs_Events_ALL <- stat.test_dsRNA_dEVEs_Events_ALL %>% add_xy_position(x = "Parameter2")
stat.test_dsRNA_dEVEs_Events_ALL$Prop <- paste0(round((nrow(Table_dsRNA_dEVEs_Events_ALL[Table_dsRNA_dEVEs_Events_ALL$b_value2 > Table_dsRNA_dEVEs_Events_ALL$b_value3,]) / nrow(Table_dsRNA_dEVEs_Events_ALL))*100,digits = 2),"%")

boxplot_dsRNA_dEVEs_Events <- ggplot(Table_dsRNA_dEVEs_Events_ALL_ecto_endo[Table_dsRNA_dEVEs_Events_ALL_ecto_endo$coef < 1000,], aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter),trim = T,adjust = 2)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="dEvents",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_dsRNA_dEVEs_Events_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,50))

##################################
## EVEs Counts all posteriors####
##################################


Table_dsRNA_EVEs_Counts_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_dsRNA_All/Posteriors_EVEs_Counts_ALL.txt",sep=";",h=T)

Table_dsRNA_EVEs_Counts_ALL$b_value2_coef_inf<- exp(as.numeric(Table_dsRNA_EVEs_Counts_ALL$b_value2))
Table_dsRNA_EVEs_Counts_ALL$b_value3_coef_inf<- exp(as.numeric(Table_dsRNA_EVEs_Counts_ALL$b_value3))
Table_dsRNA_EVEs_Counts_ALL<-Table_dsRNA_EVEs_Counts_ALL[!Table_dsRNA_EVEs_Counts_ALL$b_value2_coef_inf=="Inf",]
Table_dsRNA_EVEs_Counts_ALL<-Table_dsRNA_EVEs_Counts_ALL[!Table_dsRNA_EVEs_Counts_ALL$b_value3_coef_inf=="Inf",]


Table_dsRNA_EVEs_Counts_ALL_endo <- as.data.frame(Table_dsRNA_EVEs_Counts_ALL$b_value2)
Table_dsRNA_EVEs_Counts_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_dsRNA_EVEs_Counts_ALL_endo)<- c("coef","Parameter")
Table_dsRNA_EVEs_Counts_ALL_endo$coef2 <- as.numeric(Table_dsRNA_EVEs_Counts_ALL_endo$coef)
Table_dsRNA_EVEs_Counts_ALL_endo$coef <- as.numeric(exp(Table_dsRNA_EVEs_Counts_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_dsRNA_EVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsRNA_EVEs_Counts_ALL_endo$rope <- posterior_tab*100
Table_dsRNA_EVEs_Counts_ALL_endo$med <- median(Table_dsRNA_EVEs_Counts_ALL_endo$coef)
Table_dsRNA_EVEs_Counts_ALL_endo$pd <- pd(Table_dsRNA_EVEs_Counts_ALL_endo$coef2)

Table_dsRNA_EVEs_Counts_ALL_ecto <- as.data.frame(Table_dsRNA_EVEs_Counts_ALL$b_value3)
Table_dsRNA_EVEs_Counts_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_dsRNA_EVEs_Counts_ALL_ecto)<- c("coef","Parameter")
Table_dsRNA_EVEs_Counts_ALL_ecto$coef2 <- as.numeric(Table_dsRNA_EVEs_Counts_ALL_ecto$coef)
Table_dsRNA_EVEs_Counts_ALL_ecto$coef <- as.numeric(exp(Table_dsRNA_EVEs_Counts_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_dsRNA_EVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsRNA_EVEs_Counts_ALL_ecto$rope <- posterior_tab*100
Table_dsRNA_EVEs_Counts_ALL_ecto$med <- median(Table_dsRNA_EVEs_Counts_ALL_ecto$coef)
Table_dsRNA_EVEs_Counts_ALL_ecto$pd <- pd(Table_dsRNA_EVEs_Counts_ALL_ecto$coef2)



Table_dsRNA_EVEs_Counts_ALL_ecto_endo<-rbind(Table_dsRNA_EVEs_Counts_ALL_endo,Table_dsRNA_EVEs_Counts_ALL_ecto)

#Create summary table
statistics_Table_dsRNA_EVEs_Counts= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("dsRNA_A_EVEs_Counts",median(Table_dsRNA_EVEs_Counts_ALL_endo$coef),ci(Table_dsRNA_EVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsRNA_EVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsRNA_EVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_dsRNA_EVEs_Counts_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("dsRNA_A_EVEs_Counts",median(Table_dsRNA_EVEs_Counts_ALL_ecto$coef),ci(Table_dsRNA_EVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsRNA_EVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsRNA_EVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_dsRNA_EVEs_Counts_ALL_ecto$pd[1])))



colnames(statistics_Table_dsRNA_EVEs_Counts) <- statistics_Table_dsRNA_EVEs_Counts[1,]
statistics_Table_dsRNA_EVEs_Counts<-statistics_Table_dsRNA_EVEs_Counts[-1,]


Table_dsRNA_EVEs_Counts_ALL_ecto_endo$Parameter2 <- paste0(Table_dsRNA_EVEs_Counts_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_dsRNA_EVEs_Counts_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_dsRNA_EVEs_Counts_ALL<-Table_dsRNA_EVEs_Counts_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_dsRNA_EVEs_Counts_ALL <- stat.test_dsRNA_EVEs_Counts_ALL %>% add_xy_position(x = "Parameter2")
stat.test_dsRNA_EVEs_Counts_ALL$Prop <- paste0(round((nrow(Table_dsRNA_EVEs_Counts_ALL[Table_dsRNA_EVEs_Counts_ALL$b_value2 > Table_dsRNA_EVEs_Counts_ALL$b_value3,]) / nrow(Table_dsRNA_EVEs_Counts_ALL))*100,digits = 2),"%")

boxplot_dsRNA_EVEs_Counts <- ggplot(Table_dsRNA_EVEs_Counts_ALL_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="EVEs",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_dsRNA_EVEs_Counts_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,50))


################################
## dEVEs Counts all posteriors #
################################
Table_dsRNA_dEVEs_Counts_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_dsRNA_All/Posteriors_dEVEs_Counts_ALL.txt",sep=";",h=T)

Table_dsRNA_dEVEs_Counts_ALL$b_value2_coef_inf<- exp(as.numeric(Table_dsRNA_dEVEs_Counts_ALL$b_value2))
Table_dsRNA_dEVEs_Counts_ALL$b_value3_coef_inf<- exp(as.numeric(Table_dsRNA_dEVEs_Counts_ALL$b_value3))
Table_dsRNA_dEVEs_Counts_ALL<-Table_dsRNA_dEVEs_Counts_ALL[!Table_dsRNA_dEVEs_Counts_ALL$b_value2_coef_inf=="Inf",]
Table_dsRNA_dEVEs_Counts_ALL<-Table_dsRNA_dEVEs_Counts_ALL[!Table_dsRNA_dEVEs_Counts_ALL$b_value3_coef_inf=="Inf",]


Table_dsRNA_dEVEs_Counts_ALL_endo <- as.data.frame(Table_dsRNA_dEVEs_Counts_ALL$b_value2)
Table_dsRNA_dEVEs_Counts_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_dsRNA_dEVEs_Counts_ALL_endo)<- c("coef","Parameter")
Table_dsRNA_dEVEs_Counts_ALL_endo$coef2 <- as.numeric(Table_dsRNA_dEVEs_Counts_ALL_endo$coef)
Table_dsRNA_dEVEs_Counts_ALL_endo$coef <- as.numeric(exp(Table_dsRNA_dEVEs_Counts_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_dsRNA_dEVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsRNA_dEVEs_Counts_ALL_endo$rope <- posterior_tab*100
Table_dsRNA_dEVEs_Counts_ALL_endo$med <- median(Table_dsRNA_dEVEs_Counts_ALL_endo$coef)
Table_dsRNA_dEVEs_Counts_ALL_endo$pd <- pd(Table_dsRNA_dEVEs_Counts_ALL_endo$coef2)

Table_dsRNA_dEVEs_Counts_ALL_ecto <- as.data.frame(Table_dsRNA_dEVEs_Counts_ALL$b_value3)
Table_dsRNA_dEVEs_Counts_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_dsRNA_dEVEs_Counts_ALL_ecto)<- c("coef","Parameter")
Table_dsRNA_dEVEs_Counts_ALL_ecto$coef2 <- as.numeric(Table_dsRNA_dEVEs_Counts_ALL_ecto$coef)
Table_dsRNA_dEVEs_Counts_ALL_ecto$coef <- as.numeric(exp(Table_dsRNA_dEVEs_Counts_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_dsRNA_dEVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsRNA_dEVEs_Counts_ALL_ecto$rope <- posterior_tab*100
Table_dsRNA_dEVEs_Counts_ALL_ecto$med <- median(Table_dsRNA_dEVEs_Counts_ALL_ecto$coef)
Table_dsRNA_dEVEs_Counts_ALL_ecto$pd <- pd(Table_dsRNA_dEVEs_Counts_ALL_ecto$coef2)



Table_dsRNA_dEVEs_Counts_ALL_ecto_endo<-rbind(Table_dsRNA_dEVEs_Counts_ALL_endo,Table_dsRNA_dEVEs_Counts_ALL_ecto)

#Create summary table
statistics_Table_dsRNA_dEVEs_Counts= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("dsRNA_A_dEVEs_Counts",median(Table_dsRNA_dEVEs_Counts_ALL_endo$coef),ci(Table_dsRNA_dEVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsRNA_dEVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsRNA_dEVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_dsRNA_dEVEs_Counts_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("dsRNA_A_dEVEs_Counts",median(Table_dsRNA_dEVEs_Counts_ALL_ecto$coef),ci(Table_dsRNA_dEVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsRNA_dEVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsRNA_dEVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_dsRNA_dEVEs_Counts_ALL_ecto$pd[1])))



colnames(statistics_Table_dsRNA_dEVEs_Counts) <- statistics_Table_dsRNA_dEVEs_Counts[1,]
statistics_Table_dsRNA_dEVEs_Counts<-statistics_Table_dsRNA_dEVEs_Counts[-1,]


Table_dsRNA_dEVEs_Counts_ALL_ecto_endo$Parameter2 <- paste0(Table_dsRNA_dEVEs_Counts_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_dsRNA_dEVEs_Counts_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_dsRNA_dEVEs_Counts_ALL<-Table_dsRNA_dEVEs_Counts_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_dsRNA_dEVEs_Counts_ALL <- stat.test_dsRNA_dEVEs_Counts_ALL %>% add_xy_position(x = "Parameter2")
stat.test_dsRNA_dEVEs_Counts_ALL$Prop <- paste0(round((nrow(Table_dsRNA_dEVEs_Counts_ALL[Table_dsRNA_dEVEs_Counts_ALL$b_value2 > Table_dsRNA_dEVEs_Counts_ALL$b_value3,]) / nrow(Table_dsRNA_dEVEs_Counts_ALL))*100,digits = 2),"%")

boxplot_dsRNA_dEVEs_Counts <- ggplot(Table_dsRNA_dEVEs_Counts_ALL_ecto_endo[Table_dsRNA_dEVEs_Counts_ALL_ecto_endo$coef < 1000,], aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="dEVEs",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_dsRNA_dEVEs_Counts_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,50))



## Combine All 4 plots

violin_plot_dsRNA<-grid.arrange(arrangeGrob(
                             boxplot_dsRNA_EVEs_Counts  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank()) ,
                             boxplot_dsRNA_EVEs_Events  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank()),
                             boxplot_dsRNA_dEVEs_Counts  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank()) ,
                             boxplot_dsRNA_dEVEs_Events  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank()),
                         nrow = 1,
                         ncol=4,
                         left = textGrob("Relative rate", rot = 90, vjust = 1)),nrow=1,top=textGrob("dsRNA | A-B-C-D",
                                                                                                                 gp=gpar(fontsize=20,font=8)))




#Add posterior coefficient comparison between endo and ecto
statistics_Table_dsRNA_EVEs_Events<-as.data.frame(statistics_Table_dsRNA_EVEs_Events)
statistics_Table_dsRNA_EVEs_Events$Endo_sup_Ecto<-round((nrow(Table_dsRNA_EVEs_Events_ALL[Table_dsRNA_EVEs_Events_ALL$b_value2 > Table_dsRNA_EVEs_Events_ALL$b_value3,]) / nrow(Table_dsRNA_EVEs_Events_ALL))*100,digits = 2)

statistics_Table_dsRNA_dEVEs_Events<-as.data.frame(statistics_Table_dsRNA_dEVEs_Events)
statistics_Table_dsRNA_dEVEs_Events$Endo_sup_Ecto<-round((nrow(Table_dsRNA_dEVEs_Events_ALL[Table_dsRNA_dEVEs_Events_ALL$b_value2 > Table_dsRNA_dEVEs_Events_ALL$b_value3,]) / nrow(Table_dsRNA_dEVEs_Events_ALL))*100,digits = 2)

statistics_Table_dsRNA_EVEs_Counts<-as.data.frame(statistics_Table_dsRNA_EVEs_Counts)
statistics_Table_dsRNA_EVEs_Counts$Endo_sup_Ecto<-round((nrow(Table_dsRNA_EVEs_Counts_ALL[Table_dsRNA_EVEs_Counts_ALL$b_value2 > Table_dsRNA_EVEs_Counts_ALL$b_value3,]) / nrow(Table_dsRNA_EVEs_Counts_ALL))*100,digits = 2)


statistics_Table_dsRNA_dEVEs_Counts<-as.data.frame(statistics_Table_dsRNA_dEVEs_Counts)
statistics_Table_dsRNA_dEVEs_Counts$Endo_sup_Ecto<-round((nrow(Table_dsRNA_dEVEs_Counts_ALL[Table_dsRNA_dEVEs_Counts_ALL$b_value2 > Table_dsRNA_dEVEs_Counts_ALL$b_value3,]) / nrow(Table_dsRNA_dEVEs_Counts_ALL))*100,digits = 2)



## Combine all statistics tables

statistics_Table_dsRNA<-as.data.frame(rbind(statistics_Table_dsRNA_EVEs_Events,statistics_Table_dsRNA_dEVEs_Events,statistics_Table_dsRNA_EVEs_Counts,statistics_Table_dsRNA_dEVEs_Counts))

statistics_Table_dsRNA$Lifestyle<- gsub("\\..*","",rownames(statistics_Table_dsRNA))

library(dplyr)
statistics_Table_dsRNA$Median<-as.numeric(statistics_Table_dsRNA$Median)
statistics_Table_dsRNA$CI_low<-as.numeric(statistics_Table_dsRNA$CI_low)
statistics_Table_dsRNA$CI_high<-as.numeric(statistics_Table_dsRNA$CI_high)
statistics_Table_dsRNA$ROPE_Percentage<-as.numeric(statistics_Table_dsRNA$ROPE_Percentage)*100
statistics_Table_dsRNA<-statistics_Table_dsRNA %>% mutate_if(is.numeric, round, digits=3)

statistics_Table_dsRNA$Type <-c("EVEs_Events","EVEs_Events","dEVEs_Events","dEVEs_Events","EVEs_Numbers","EVEs_Numbers","dEVEs_Numbers","dEVEs_Numbers")

table_dsRNA<-statistics_Table_dsRNA
table_dsRNA[nrow(table_dsRNA)+1,] <- NA




#dsRNA_A A

library(bayestestR)
library(ggstatsplot)
library(bayesplot)
library(ggpubr)
library(gridExtra)
library(grid)


##################################
# ## EVEs Events all posteriors ##
##################################
Table_dsRNA_A_EVEs_Events_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_dsRNA_A_All/Posteriors_EVEs_Events_ALL.txt",sep=";",h=T)

Table_dsRNA_A_EVEs_Events_ALL$b_value2_coef_inf<- exp(as.numeric(Table_dsRNA_A_EVEs_Events_ALL$b_value2))
Table_dsRNA_A_EVEs_Events_ALL$b_value3_coef_inf<- exp(as.numeric(Table_dsRNA_A_EVEs_Events_ALL$b_value3))
Table_dsRNA_A_EVEs_Events_ALL<-Table_dsRNA_A_EVEs_Events_ALL[!Table_dsRNA_A_EVEs_Events_ALL$b_value2_coef_inf=="Inf",]
Table_dsRNA_A_EVEs_Events_ALL<-Table_dsRNA_A_EVEs_Events_ALL[!Table_dsRNA_A_EVEs_Events_ALL$b_value3_coef_inf=="Inf",]


Table_dsRNA_A_EVEs_Events_ALL_endo <- as.data.frame(Table_dsRNA_A_EVEs_Events_ALL$b_value2)
Table_dsRNA_A_EVEs_Events_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_dsRNA_A_EVEs_Events_ALL_endo)<- c("coef","Parameter")
Table_dsRNA_A_EVEs_Events_ALL_endo$coef2 <- as.numeric(Table_dsRNA_A_EVEs_Events_ALL_endo$coef)
Table_dsRNA_A_EVEs_Events_ALL_endo$coef <- as.numeric(exp(Table_dsRNA_A_EVEs_Events_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_dsRNA_A_EVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsRNA_A_EVEs_Events_ALL_endo$rope <- posterior_tab*100
Table_dsRNA_A_EVEs_Events_ALL_endo$med <- median(Table_dsRNA_A_EVEs_Events_ALL_endo$coef)
Table_dsRNA_A_EVEs_Events_ALL_endo$pd <- pd(Table_dsRNA_A_EVEs_Events_ALL_endo$coef2)

Table_dsRNA_A_EVEs_Events_ALL_ecto <- as.data.frame(Table_dsRNA_A_EVEs_Events_ALL$b_value3)
Table_dsRNA_A_EVEs_Events_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_dsRNA_A_EVEs_Events_ALL_ecto)<- c("coef","Parameter")
Table_dsRNA_A_EVEs_Events_ALL_ecto$coef2 <- as.numeric(Table_dsRNA_A_EVEs_Events_ALL_ecto$coef)
Table_dsRNA_A_EVEs_Events_ALL_ecto$coef <- as.numeric(exp(Table_dsRNA_A_EVEs_Events_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_dsRNA_A_EVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsRNA_A_EVEs_Events_ALL_ecto$rope <- posterior_tab*100
Table_dsRNA_A_EVEs_Events_ALL_ecto$med <- median(Table_dsRNA_A_EVEs_Events_ALL_ecto$coef)
Table_dsRNA_A_EVEs_Events_ALL_ecto$pd <- pd(Table_dsRNA_A_EVEs_Events_ALL_ecto$coef2)



Table_dsRNA_A_EVEs_Events_ALL_ecto_endo<-rbind(Table_dsRNA_A_EVEs_Events_ALL_endo,Table_dsRNA_A_EVEs_Events_ALL_ecto)

#Create summary table
statistics_Table_dsRNA_A_EVEs_Events= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("dsRNA_A_EVEs_Events",median(Table_dsRNA_A_EVEs_Events_ALL_endo$coef),ci(Table_dsRNA_A_EVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsRNA_A_EVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsRNA_A_EVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_dsRNA_A_EVEs_Events_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("dsRNA_A_EVEs_Events",median(Table_dsRNA_A_EVEs_Events_ALL_ecto$coef),ci(Table_dsRNA_A_EVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsRNA_A_EVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsRNA_A_EVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_dsRNA_A_EVEs_Events_ALL_ecto$pd[1])))



colnames(statistics_Table_dsRNA_A_EVEs_Events) <- statistics_Table_dsRNA_A_EVEs_Events[1,]
statistics_Table_dsRNA_A_EVEs_Events<-statistics_Table_dsRNA_A_EVEs_Events[-1,]


Table_dsRNA_A_EVEs_Events_ALL_ecto_endo$Parameter2 <- paste0(Table_dsRNA_A_EVEs_Events_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_dsRNA_A_EVEs_Events_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_dsRNA_A_EVEs_Events_ALL<-Table_dsRNA_A_EVEs_Events_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_dsRNA_A_EVEs_Events_ALL <- stat.test_dsRNA_A_EVEs_Events_ALL %>% add_xy_position(x = "Parameter2")
stat.test_dsRNA_A_EVEs_Events_ALL$Prop <- paste0(round((nrow(Table_dsRNA_A_EVEs_Events_ALL[Table_dsRNA_A_EVEs_Events_ALL$b_value2 > Table_dsRNA_A_EVEs_Events_ALL$b_value3,]) / nrow(Table_dsRNA_A_EVEs_Events_ALL))*100,digits = 2),"%")

boxplot_dsRNA_A_EVEs_Events <- ggplot(Table_dsRNA_A_EVEs_Events_ALL_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="Events",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_dsRNA_A_EVEs_Events_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,50))


##################################
## dEVEs Events all posteriors
##################################

Table_dsRNA_A_dEVEs_Events_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_dsRNA_A_All/Posteriors_dEVEs_Events_ALL.txt",sep=";",h=T)

Table_dsRNA_A_dEVEs_Events_ALL$b_value2_coef_inf<- exp(as.numeric(Table_dsRNA_A_dEVEs_Events_ALL$b_value2))
Table_dsRNA_A_dEVEs_Events_ALL$b_value3_coef_inf<- exp(as.numeric(Table_dsRNA_A_dEVEs_Events_ALL$b_value3))
Table_dsRNA_A_dEVEs_Events_ALL<-Table_dsRNA_A_dEVEs_Events_ALL[!Table_dsRNA_A_dEVEs_Events_ALL$b_value2_coef_inf=="Inf",]
Table_dsRNA_A_dEVEs_Events_ALL<-Table_dsRNA_A_dEVEs_Events_ALL[!Table_dsRNA_A_dEVEs_Events_ALL$b_value3_coef_inf=="Inf",]


Table_dsRNA_A_dEVEs_Events_ALL_endo <- as.data.frame(Table_dsRNA_A_dEVEs_Events_ALL$b_value2)
Table_dsRNA_A_dEVEs_Events_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_dsRNA_A_dEVEs_Events_ALL_endo)<- c("coef","Parameter")
Table_dsRNA_A_dEVEs_Events_ALL_endo$coef2 <- as.numeric(Table_dsRNA_A_dEVEs_Events_ALL_endo$coef)
Table_dsRNA_A_dEVEs_Events_ALL_endo$coef <- as.numeric(exp(Table_dsRNA_A_dEVEs_Events_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_dsRNA_A_dEVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsRNA_A_dEVEs_Events_ALL_endo$rope <- posterior_tab*100
Table_dsRNA_A_dEVEs_Events_ALL_endo$med <- median(Table_dsRNA_A_dEVEs_Events_ALL_endo$coef)
Table_dsRNA_A_dEVEs_Events_ALL_endo$pd <- pd(Table_dsRNA_A_dEVEs_Events_ALL_endo$coef2)

Table_dsRNA_A_dEVEs_Events_ALL_ecto <- as.data.frame(Table_dsRNA_A_dEVEs_Events_ALL$b_value3)
Table_dsRNA_A_dEVEs_Events_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_dsRNA_A_dEVEs_Events_ALL_ecto)<- c("coef","Parameter")
Table_dsRNA_A_dEVEs_Events_ALL_ecto$coef2 <- as.numeric(Table_dsRNA_A_dEVEs_Events_ALL_ecto$coef)
Table_dsRNA_A_dEVEs_Events_ALL_ecto$coef <- as.numeric(exp(Table_dsRNA_A_dEVEs_Events_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_dsRNA_A_dEVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsRNA_A_dEVEs_Events_ALL_ecto$rope <- posterior_tab*100
Table_dsRNA_A_dEVEs_Events_ALL_ecto$med <- median(Table_dsRNA_A_dEVEs_Events_ALL_ecto$coef)
Table_dsRNA_A_dEVEs_Events_ALL_ecto$pd <- pd(Table_dsRNA_A_dEVEs_Events_ALL_ecto$coef2)



Table_dsRNA_A_dEVEs_Events_ALL_ecto_endo<-rbind(Table_dsRNA_A_dEVEs_Events_ALL_endo,Table_dsRNA_A_dEVEs_Events_ALL_ecto)

#Create summary table
statistics_Table_dsRNA_A_dEVEs_Events= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("dsRNA_A_dEVEs_Events",median(Table_dsRNA_A_dEVEs_Events_ALL_endo$coef),ci(Table_dsRNA_A_dEVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsRNA_A_dEVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsRNA_A_dEVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_dsRNA_A_dEVEs_Events_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("dsRNA_A_dEVEs_Events",median(Table_dsRNA_A_dEVEs_Events_ALL_ecto$coef),ci(Table_dsRNA_A_dEVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsRNA_A_dEVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsRNA_A_dEVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_dsRNA_A_dEVEs_Events_ALL_ecto$pd[1])))



colnames(statistics_Table_dsRNA_A_dEVEs_Events) <- statistics_Table_dsRNA_A_dEVEs_Events[1,]
statistics_Table_dsRNA_A_dEVEs_Events<-statistics_Table_dsRNA_A_dEVEs_Events[-1,]


Table_dsRNA_A_dEVEs_Events_ALL_ecto_endo$Parameter2 <- paste0(Table_dsRNA_A_dEVEs_Events_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_dsRNA_A_dEVEs_Events_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_dsRNA_A_dEVEs_Events_ALL<-Table_dsRNA_A_dEVEs_Events_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_dsRNA_A_dEVEs_Events_ALL <- stat.test_dsRNA_A_dEVEs_Events_ALL %>% add_xy_position(x = "Parameter2")
stat.test_dsRNA_A_dEVEs_Events_ALL$Prop <- paste0(round((nrow(Table_dsRNA_A_dEVEs_Events_ALL[Table_dsRNA_A_dEVEs_Events_ALL$b_value2 > Table_dsRNA_A_dEVEs_Events_ALL$b_value3,]) / nrow(Table_dsRNA_A_dEVEs_Events_ALL))*100,digits = 2),"%")

boxplot_dsRNA_A_dEVEs_Events <- ggplot(Table_dsRNA_A_dEVEs_Events_ALL_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter),trim = T,adjust = 2)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="dEvents",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_dsRNA_A_dEVEs_Events_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,50))


##################################
## EVEs Counts all posteriors####
##################################


Table_dsRNA_A_EVEs_Counts_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_dsRNA_A_All/Posteriors_EVEs_Counts_ALL.txt",sep=";",h=T)

Table_dsRNA_A_EVEs_Counts_ALL$b_value2_coef_inf<- exp(as.numeric(Table_dsRNA_A_EVEs_Counts_ALL$b_value2))
Table_dsRNA_A_EVEs_Counts_ALL$b_value3_coef_inf<- exp(as.numeric(Table_dsRNA_A_EVEs_Counts_ALL$b_value3))
Table_dsRNA_A_EVEs_Counts_ALL<-Table_dsRNA_A_EVEs_Counts_ALL[!Table_dsRNA_A_EVEs_Counts_ALL$b_value2_coef_inf=="Inf",]
Table_dsRNA_A_EVEs_Counts_ALL<-Table_dsRNA_A_EVEs_Counts_ALL[!Table_dsRNA_A_EVEs_Counts_ALL$b_value3_coef_inf=="Inf",]


Table_dsRNA_A_EVEs_Counts_ALL_endo <- as.data.frame(Table_dsRNA_A_EVEs_Counts_ALL$b_value2)
Table_dsRNA_A_EVEs_Counts_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_dsRNA_A_EVEs_Counts_ALL_endo)<- c("coef","Parameter")
Table_dsRNA_A_EVEs_Counts_ALL_endo$coef2 <- as.numeric(Table_dsRNA_A_EVEs_Counts_ALL_endo$coef)
Table_dsRNA_A_EVEs_Counts_ALL_endo$coef <- as.numeric(exp(Table_dsRNA_A_EVEs_Counts_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_dsRNA_A_EVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsRNA_A_EVEs_Counts_ALL_endo$rope <- posterior_tab*100
Table_dsRNA_A_EVEs_Counts_ALL_endo$med <- median(Table_dsRNA_A_EVEs_Counts_ALL_endo$coef)
Table_dsRNA_A_EVEs_Counts_ALL_endo$pd <- pd(Table_dsRNA_A_EVEs_Counts_ALL_endo$coef2)

Table_dsRNA_A_EVEs_Counts_ALL_ecto <- as.data.frame(Table_dsRNA_A_EVEs_Counts_ALL$b_value3)
Table_dsRNA_A_EVEs_Counts_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_dsRNA_A_EVEs_Counts_ALL_ecto)<- c("coef","Parameter")
Table_dsRNA_A_EVEs_Counts_ALL_ecto$coef2 <- as.numeric(Table_dsRNA_A_EVEs_Counts_ALL_ecto$coef)
Table_dsRNA_A_EVEs_Counts_ALL_ecto$coef <- as.numeric(exp(Table_dsRNA_A_EVEs_Counts_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_dsRNA_A_EVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsRNA_A_EVEs_Counts_ALL_ecto$rope <- posterior_tab*100
Table_dsRNA_A_EVEs_Counts_ALL_ecto$med <- median(Table_dsRNA_A_EVEs_Counts_ALL_ecto$coef)
Table_dsRNA_A_EVEs_Counts_ALL_ecto$pd <- pd(Table_dsRNA_A_EVEs_Counts_ALL_ecto$coef2)



Table_dsRNA_A_EVEs_Counts_ALL_ecto_endo<-rbind(Table_dsRNA_A_EVEs_Counts_ALL_endo,Table_dsRNA_A_EVEs_Counts_ALL_ecto)

#Create summary table
statistics_Table_dsRNA_A_EVEs_Counts= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("dsRNA_A_EVEs_Counts",median(Table_dsRNA_A_EVEs_Counts_ALL_endo$coef),ci(Table_dsRNA_A_EVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsRNA_A_EVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsRNA_A_EVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_dsRNA_A_EVEs_Counts_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("dsRNA_A_EVEs_Counts",median(Table_dsRNA_A_EVEs_Counts_ALL_ecto$coef),ci(Table_dsRNA_A_EVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsRNA_A_EVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsRNA_A_EVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_dsRNA_A_EVEs_Counts_ALL_ecto$pd[1])))



colnames(statistics_Table_dsRNA_A_EVEs_Counts) <- statistics_Table_dsRNA_A_EVEs_Counts[1,]
statistics_Table_dsRNA_A_EVEs_Counts<-statistics_Table_dsRNA_A_EVEs_Counts[-1,]


Table_dsRNA_A_EVEs_Counts_ALL_ecto_endo$Parameter2 <- paste0(Table_dsRNA_A_EVEs_Counts_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_dsRNA_A_EVEs_Counts_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_dsRNA_A_EVEs_Counts_ALL<-Table_dsRNA_A_EVEs_Counts_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_dsRNA_A_EVEs_Counts_ALL <- stat.test_dsRNA_A_EVEs_Counts_ALL %>% add_xy_position(x = "Parameter2")
stat.test_dsRNA_A_EVEs_Counts_ALL$Prop <- paste0(round((nrow(Table_dsRNA_A_EVEs_Counts_ALL[Table_dsRNA_A_EVEs_Counts_ALL$b_value2 > Table_dsRNA_A_EVEs_Counts_ALL$b_value3,]) / nrow(Table_dsRNA_A_EVEs_Counts_ALL))*100,digits = 2),"%")

boxplot_dsRNA_A_EVEs_Counts <- ggplot(Table_dsRNA_A_EVEs_Counts_ALL_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="EVEs",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_dsRNA_A_EVEs_Counts_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,50))


################################
## dEVEs Counts all posteriors #
################################
Table_dsRNA_A_dEVEs_Counts_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_dsRNA_A_All/Posteriors_dEVEs_Counts_ALL.txt",sep=";",h=T)

Table_dsRNA_A_dEVEs_Counts_ALL$b_value2_coef_inf<- exp(as.numeric(Table_dsRNA_A_dEVEs_Counts_ALL$b_value2))
Table_dsRNA_A_dEVEs_Counts_ALL$b_value3_coef_inf<- exp(as.numeric(Table_dsRNA_A_dEVEs_Counts_ALL$b_value3))
Table_dsRNA_A_dEVEs_Counts_ALL<-Table_dsRNA_A_dEVEs_Counts_ALL[!Table_dsRNA_A_dEVEs_Counts_ALL$b_value2_coef_inf=="Inf",]
Table_dsRNA_A_dEVEs_Counts_ALL<-Table_dsRNA_A_dEVEs_Counts_ALL[!Table_dsRNA_A_dEVEs_Counts_ALL$b_value3_coef_inf=="Inf",]


Table_dsRNA_A_dEVEs_Counts_ALL_endo <- as.data.frame(Table_dsRNA_A_dEVEs_Counts_ALL$b_value2)
Table_dsRNA_A_dEVEs_Counts_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_dsRNA_A_dEVEs_Counts_ALL_endo)<- c("coef","Parameter")
Table_dsRNA_A_dEVEs_Counts_ALL_endo$coef2 <- as.numeric(Table_dsRNA_A_dEVEs_Counts_ALL_endo$coef)
Table_dsRNA_A_dEVEs_Counts_ALL_endo$coef <- as.numeric(exp(Table_dsRNA_A_dEVEs_Counts_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_dsRNA_A_dEVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsRNA_A_dEVEs_Counts_ALL_endo$rope <- posterior_tab*100
Table_dsRNA_A_dEVEs_Counts_ALL_endo$med <- median(Table_dsRNA_A_dEVEs_Counts_ALL_endo$coef)
Table_dsRNA_A_dEVEs_Counts_ALL_endo$pd <- pd(Table_dsRNA_A_dEVEs_Counts_ALL_endo$coef2)

Table_dsRNA_A_dEVEs_Counts_ALL_ecto <- as.data.frame(Table_dsRNA_A_dEVEs_Counts_ALL$b_value3)
Table_dsRNA_A_dEVEs_Counts_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_dsRNA_A_dEVEs_Counts_ALL_ecto)<- c("coef","Parameter")
Table_dsRNA_A_dEVEs_Counts_ALL_ecto$coef2 <- as.numeric(Table_dsRNA_A_dEVEs_Counts_ALL_ecto$coef)
Table_dsRNA_A_dEVEs_Counts_ALL_ecto$coef <- as.numeric(exp(Table_dsRNA_A_dEVEs_Counts_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_dsRNA_A_dEVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_dsRNA_A_dEVEs_Counts_ALL_ecto$rope <- posterior_tab*100
Table_dsRNA_A_dEVEs_Counts_ALL_ecto$med <- median(Table_dsRNA_A_dEVEs_Counts_ALL_ecto$coef)
Table_dsRNA_A_dEVEs_Counts_ALL_ecto$pd <- pd(Table_dsRNA_A_dEVEs_Counts_ALL_ecto$coef2)



Table_dsRNA_A_dEVEs_Counts_ALL_ecto_endo<-rbind(Table_dsRNA_A_dEVEs_Counts_ALL_endo,Table_dsRNA_A_dEVEs_Counts_ALL_ecto)

#Create summary table
statistics_Table_dsRNA_A_dEVEs_Counts= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("dsRNA_A_dEVEs_Counts",median(Table_dsRNA_A_dEVEs_Counts_ALL_endo$coef),ci(Table_dsRNA_A_dEVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsRNA_A_dEVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsRNA_A_dEVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_dsRNA_A_dEVEs_Counts_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("dsRNA_A_dEVEs_Counts",median(Table_dsRNA_A_dEVEs_Counts_ALL_ecto$coef),ci(Table_dsRNA_A_dEVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_dsRNA_A_dEVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_dsRNA_A_dEVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_dsRNA_A_dEVEs_Counts_ALL_ecto$pd[1])))



colnames(statistics_Table_dsRNA_A_dEVEs_Counts) <- statistics_Table_dsRNA_A_dEVEs_Counts[1,]
statistics_Table_dsRNA_A_dEVEs_Counts<-statistics_Table_dsRNA_A_dEVEs_Counts[-1,]


Table_dsRNA_A_dEVEs_Counts_ALL_ecto_endo$Parameter2 <- paste0(Table_dsRNA_A_dEVEs_Counts_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_dsRNA_A_dEVEs_Counts_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_dsRNA_A_dEVEs_Counts_ALL<-Table_dsRNA_A_dEVEs_Counts_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_dsRNA_A_dEVEs_Counts_ALL <- stat.test_dsRNA_A_dEVEs_Counts_ALL %>% add_xy_position(x = "Parameter2")
stat.test_dsRNA_A_dEVEs_Counts_ALL$Prop <- paste0(round((nrow(Table_dsRNA_A_dEVEs_Counts_ALL[Table_dsRNA_A_dEVEs_Counts_ALL$b_value2 > Table_dsRNA_A_dEVEs_Counts_ALL$b_value3,]) / nrow(Table_dsRNA_A_dEVEs_Counts_ALL))*100,digits = 2),"%")

boxplot_dsRNA_A_dEVEs_Counts <- ggplot(Table_dsRNA_A_dEVEs_Counts_ALL_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="dEVEs",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_dsRNA_A_dEVEs_Counts_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,200))



## Combine All 4 plots

violin_plot_dsRNA_A<-grid.arrange(arrangeGrob(
                             boxplot_dsRNA_A_EVEs_Counts  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank()) ,
                             boxplot_dsRNA_A_EVEs_Events  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank()),
                             boxplot_dsRNA_A_dEVEs_Counts  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank()) ,
                             boxplot_dsRNA_A_dEVEs_Events  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank()),
                         nrow = 1,
                         ncol=4,
                         left = textGrob("Relative rate", rot = 90, vjust = 1)),nrow=1,top=textGrob("dsRNA | A",
                                                                                                                 gp=gpar(fontsize=20,font=8)))



#Add posterior coefficient comparison between endo and ecto
statistics_Table_dsRNA_A_EVEs_Events<-as.data.frame(statistics_Table_dsRNA_A_EVEs_Events)
statistics_Table_dsRNA_A_EVEs_Events$Endo_sup_Ecto<-round((nrow(Table_dsRNA_A_EVEs_Events_ALL[Table_dsRNA_A_EVEs_Events_ALL$b_value2 > Table_dsRNA_A_EVEs_Events_ALL$b_value3,]) / nrow(Table_dsRNA_A_EVEs_Events_ALL))*100,digits = 2)

statistics_Table_dsRNA_A_dEVEs_Events<-as.data.frame(statistics_Table_dsRNA_A_dEVEs_Events)
statistics_Table_dsRNA_A_dEVEs_Events$Endo_sup_Ecto<-round((nrow(Table_dsRNA_A_dEVEs_Events_ALL[Table_dsRNA_A_dEVEs_Events_ALL$b_value2 > Table_dsRNA_A_dEVEs_Events_ALL$b_value3,]) / nrow(Table_dsRNA_A_dEVEs_Events_ALL))*100,digits = 2)

statistics_Table_dsRNA_A_EVEs_Counts<-as.data.frame(statistics_Table_dsRNA_A_EVEs_Counts)
statistics_Table_dsRNA_A_EVEs_Counts$Endo_sup_Ecto<-round((nrow(Table_dsRNA_A_EVEs_Counts_ALL[Table_dsRNA_A_EVEs_Counts_ALL$b_value2 > Table_dsRNA_A_EVEs_Counts_ALL$b_value3,]) / nrow(Table_dsRNA_A_EVEs_Counts_ALL))*100,digits = 2)


statistics_Table_dsRNA_A_dEVEs_Counts<-as.data.frame(statistics_Table_dsRNA_A_dEVEs_Counts)
statistics_Table_dsRNA_A_dEVEs_Counts$Endo_sup_Ecto<-round((nrow(Table_dsRNA_A_dEVEs_Counts_ALL[Table_dsRNA_A_dEVEs_Counts_ALL$b_value2 > Table_dsRNA_A_dEVEs_Counts_ALL$b_value3,]) / nrow(Table_dsRNA_A_dEVEs_Counts_ALL))*100,digits = 2)


## Combine all statistics tables

statistics_Table_dsRNA_A<-as.data.frame(rbind(statistics_Table_dsRNA_A_EVEs_Events,statistics_Table_dsRNA_A_dEVEs_Events,statistics_Table_dsRNA_A_EVEs_Counts,statistics_Table_dsRNA_A_dEVEs_Counts))

statistics_Table_dsRNA_A$Lifestyle<- gsub("\\..*","",rownames(statistics_Table_dsRNA_A))

library(dplyr)
statistics_Table_dsRNA_A$Median<-as.numeric(statistics_Table_dsRNA_A$Median)
statistics_Table_dsRNA_A$CI_low<-as.numeric(statistics_Table_dsRNA_A$CI_low)
statistics_Table_dsRNA_A$CI_high<-as.numeric(statistics_Table_dsRNA_A$CI_high)
statistics_Table_dsRNA_A$ROPE_Percentage<-as.numeric(statistics_Table_dsRNA_A$ROPE_Percentage)*100
statistics_Table_dsRNA_A<-statistics_Table_dsRNA_A %>% mutate_if(is.numeric, round, digits=3)

statistics_Table_dsRNA_A$Type <-c("EVEs_Events","EVEs_Events","dEVEs_Events","dEVEs_Events","EVEs_Numbers","EVEs_Numbers","dEVEs_Numbers","dEVEs_Numbers")

table_dsRNA_A<-statistics_Table_dsRNA_A
table_dsRNA_A[nrow(table_dsRNA_A)+1,] <- NA





#ssRNA A-B-C-D


library(bayestestR)
library(ggstatsplot)
library(bayesplot)
library(ggpubr)
library(gridExtra)
library(grid)



##################################
# ## EVEs Events all posteriors ##
##################################
Table_ssRNA_EVEs_Events_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_ssRNA_All/Posteriors_EVEs_Events_ALL.txt",sep=";",h=T)

Table_ssRNA_EVEs_Events_ALL$b_value2_coef_inf<- exp(as.numeric(Table_ssRNA_EVEs_Events_ALL$b_value2))
Table_ssRNA_EVEs_Events_ALL$b_value3_coef_inf<- exp(as.numeric(Table_ssRNA_EVEs_Events_ALL$b_value3))
Table_ssRNA_EVEs_Events_ALL<-Table_ssRNA_EVEs_Events_ALL[!Table_ssRNA_EVEs_Events_ALL$b_value2_coef_inf=="Inf",]
Table_ssRNA_EVEs_Events_ALL<-Table_ssRNA_EVEs_Events_ALL[!Table_ssRNA_EVEs_Events_ALL$b_value3_coef_inf=="Inf",]


Table_ssRNA_EVEs_Events_ALL_endo <- as.data.frame(Table_ssRNA_EVEs_Events_ALL$b_value2)
Table_ssRNA_EVEs_Events_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_ssRNA_EVEs_Events_ALL_endo)<- c("coef","Parameter")
Table_ssRNA_EVEs_Events_ALL_endo$coef2 <- as.numeric(Table_ssRNA_EVEs_Events_ALL_endo$coef)
Table_ssRNA_EVEs_Events_ALL_endo$coef <- as.numeric(exp(Table_ssRNA_EVEs_Events_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_ssRNA_EVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_ssRNA_EVEs_Events_ALL_endo$rope <- posterior_tab*100
Table_ssRNA_EVEs_Events_ALL_endo$med <- median(Table_ssRNA_EVEs_Events_ALL_endo$coef)
Table_ssRNA_EVEs_Events_ALL_endo$pd <- pd(Table_ssRNA_EVEs_Events_ALL_endo$coef2)

Table_ssRNA_EVEs_Events_ALL_ecto <- as.data.frame(Table_ssRNA_EVEs_Events_ALL$b_value3)
Table_ssRNA_EVEs_Events_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_ssRNA_EVEs_Events_ALL_ecto)<- c("coef","Parameter")
Table_ssRNA_EVEs_Events_ALL_ecto$coef2 <- as.numeric(Table_ssRNA_EVEs_Events_ALL_ecto$coef)
Table_ssRNA_EVEs_Events_ALL_ecto$coef <- as.numeric(exp(Table_ssRNA_EVEs_Events_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_ssRNA_EVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_ssRNA_EVEs_Events_ALL_ecto$rope <- posterior_tab*100
Table_ssRNA_EVEs_Events_ALL_ecto$med <- median(Table_ssRNA_EVEs_Events_ALL_ecto$coef)
Table_ssRNA_EVEs_Events_ALL_ecto$pd <- pd(Table_ssRNA_EVEs_Events_ALL_ecto$coef2)



Table_ssRNA_EVEs_Events_ALL_ecto_endo<-rbind(Table_ssRNA_EVEs_Events_ALL_endo,Table_ssRNA_EVEs_Events_ALL_ecto)

#Create summary table
statistics_Table_ssRNA_EVEs_Events= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("ssRNA_EVEs_Events",median(Table_ssRNA_EVEs_Events_ALL_endo$coef),ci(Table_ssRNA_EVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_ssRNA_EVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_ssRNA_EVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_ssRNA_EVEs_Events_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("ssRNA_EVEs_Events",median(Table_ssRNA_EVEs_Events_ALL_ecto$coef),ci(Table_ssRNA_EVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_ssRNA_EVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_ssRNA_EVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_ssRNA_EVEs_Events_ALL_ecto$pd[1])))



colnames(statistics_Table_ssRNA_EVEs_Events) <- statistics_Table_ssRNA_EVEs_Events[1,]
statistics_Table_ssRNA_EVEs_Events<-statistics_Table_ssRNA_EVEs_Events[-1,]


Table_ssRNA_EVEs_Events_ALL_ecto_endo$Parameter2 <- paste0(Table_ssRNA_EVEs_Events_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_ssRNA_EVEs_Events_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_ssRNA_EVEs_Events_ALL<-Table_ssRNA_EVEs_Events_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_ssRNA_EVEs_Events_ALL <- stat.test_ssRNA_EVEs_Events_ALL %>% add_xy_position(x = "Parameter2")
stat.test_ssRNA_EVEs_Events_ALL$Prop <- paste0(round((nrow(Table_ssRNA_EVEs_Events_ALL[Table_ssRNA_EVEs_Events_ALL$b_value2 > Table_ssRNA_EVEs_Events_ALL$b_value3,]) / nrow(Table_ssRNA_EVEs_Events_ALL))*100,digits = 2),"%")

boxplot_ssRNA_EVEs_Events <- ggplot(Table_ssRNA_EVEs_Events_ALL_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="Events",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_ssRNA_EVEs_Events_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,15))


##################################
## dEVEs Events all posteriors
##################################

Table_ssRNA_dEVEs_Events_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_ssRNA_All/Posteriors_dEVEs_Events_ALL.txt",sep=";",h=T)

Table_ssRNA_dEVEs_Events_ALL$b_value2_coef_inf<- exp(as.numeric(Table_ssRNA_dEVEs_Events_ALL$b_value2))
Table_ssRNA_dEVEs_Events_ALL$b_value3_coef_inf<- exp(as.numeric(Table_ssRNA_dEVEs_Events_ALL$b_value3))
Table_ssRNA_dEVEs_Events_ALL<-Table_ssRNA_dEVEs_Events_ALL[!Table_ssRNA_dEVEs_Events_ALL$b_value2_coef_inf=="Inf",]
Table_ssRNA_dEVEs_Events_ALL<-Table_ssRNA_dEVEs_Events_ALL[!Table_ssRNA_dEVEs_Events_ALL$b_value3_coef_inf=="Inf",]


Table_ssRNA_dEVEs_Events_ALL_endo <- as.data.frame(Table_ssRNA_dEVEs_Events_ALL$b_value2)
Table_ssRNA_dEVEs_Events_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_ssRNA_dEVEs_Events_ALL_endo)<- c("coef","Parameter")
Table_ssRNA_dEVEs_Events_ALL_endo$coef2 <- as.numeric(Table_ssRNA_dEVEs_Events_ALL_endo$coef)
Table_ssRNA_dEVEs_Events_ALL_endo$coef <- as.numeric(exp(Table_ssRNA_dEVEs_Events_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_ssRNA_dEVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_ssRNA_dEVEs_Events_ALL_endo$rope <- posterior_tab*100
Table_ssRNA_dEVEs_Events_ALL_endo$med <- median(Table_ssRNA_dEVEs_Events_ALL_endo$coef)
Table_ssRNA_dEVEs_Events_ALL_endo$pd <- pd(Table_ssRNA_dEVEs_Events_ALL_endo$coef2)

Table_ssRNA_dEVEs_Events_ALL_ecto <- as.data.frame(Table_ssRNA_dEVEs_Events_ALL$b_value3)
Table_ssRNA_dEVEs_Events_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_ssRNA_dEVEs_Events_ALL_ecto)<- c("coef","Parameter")
Table_ssRNA_dEVEs_Events_ALL_ecto$coef2 <- as.numeric(Table_ssRNA_dEVEs_Events_ALL_ecto$coef)
Table_ssRNA_dEVEs_Events_ALL_ecto$coef <- as.numeric(exp(Table_ssRNA_dEVEs_Events_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_ssRNA_dEVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_ssRNA_dEVEs_Events_ALL_ecto$rope <- posterior_tab*100
Table_ssRNA_dEVEs_Events_ALL_ecto$med <- median(Table_ssRNA_dEVEs_Events_ALL_ecto$coef)
Table_ssRNA_dEVEs_Events_ALL_ecto$pd <- pd(Table_ssRNA_dEVEs_Events_ALL_ecto$coef2)



Table_ssRNA_dEVEs_Events_ALL_ecto_endo<-rbind(Table_ssRNA_dEVEs_Events_ALL_endo,Table_ssRNA_dEVEs_Events_ALL_ecto)

#Create summary table
statistics_Table_ssRNA_dEVEs_Events= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("ssRNA_dEVEs_Events",median(Table_ssRNA_dEVEs_Events_ALL_endo$coef),ci(Table_ssRNA_dEVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_ssRNA_dEVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_ssRNA_dEVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_ssRNA_dEVEs_Events_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("ssRNA_dEVEs_Events",median(Table_ssRNA_dEVEs_Events_ALL_ecto$coef),ci(Table_ssRNA_dEVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_ssRNA_dEVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_ssRNA_dEVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_ssRNA_dEVEs_Events_ALL_ecto$pd[1])))



colnames(statistics_Table_ssRNA_dEVEs_Events) <- statistics_Table_ssRNA_dEVEs_Events[1,]
statistics_Table_ssRNA_dEVEs_Events<-statistics_Table_ssRNA_dEVEs_Events[-1,]


Table_ssRNA_dEVEs_Events_ALL_ecto_endo$Parameter2 <- paste0(Table_ssRNA_dEVEs_Events_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_ssRNA_dEVEs_Events_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_ssRNA_dEVEs_Events_ALL<-Table_ssRNA_dEVEs_Events_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_ssRNA_dEVEs_Events_ALL <- stat.test_ssRNA_dEVEs_Events_ALL %>% add_xy_position(x = "Parameter2")
stat.test_ssRNA_dEVEs_Events_ALL$Prop <- paste0(round((nrow(Table_ssRNA_dEVEs_Events_ALL[Table_ssRNA_dEVEs_Events_ALL$b_value2 > Table_ssRNA_dEVEs_Events_ALL$b_value3,]) / nrow(Table_ssRNA_dEVEs_Events_ALL))*100,digits = 2),"%")

boxplot_ssRNA_dEVEs_Events <- ggplot(Table_ssRNA_dEVEs_Events_ALL_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter),trim = T,adjust = 2)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="dEvents",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_ssRNA_dEVEs_Events_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,15))

##################################
## EVEs Counts all posteriors####
##################################


Table_ssRNA_EVEs_Counts_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_ssRNA_All/Posteriors_EVEs_Counts_ALL.txt",sep=";",h=T)

Table_ssRNA_EVEs_Counts_ALL$b_value2_coef_inf<- exp(as.numeric(Table_ssRNA_EVEs_Counts_ALL$b_value2))
Table_ssRNA_EVEs_Counts_ALL$b_value3_coef_inf<- exp(as.numeric(Table_ssRNA_EVEs_Counts_ALL$b_value3))
Table_ssRNA_EVEs_Counts_ALL<-Table_ssRNA_EVEs_Counts_ALL[!Table_ssRNA_EVEs_Counts_ALL$b_value2_coef_inf=="Inf",]
Table_ssRNA_EVEs_Counts_ALL<-Table_ssRNA_EVEs_Counts_ALL[!Table_ssRNA_EVEs_Counts_ALL$b_value3_coef_inf=="Inf",]


Table_ssRNA_EVEs_Counts_ALL_endo <- as.data.frame(Table_ssRNA_EVEs_Counts_ALL$b_value2)
Table_ssRNA_EVEs_Counts_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_ssRNA_EVEs_Counts_ALL_endo)<- c("coef","Parameter")
Table_ssRNA_EVEs_Counts_ALL_endo$coef2 <- as.numeric(Table_ssRNA_EVEs_Counts_ALL_endo$coef)
Table_ssRNA_EVEs_Counts_ALL_endo$coef <- as.numeric(exp(Table_ssRNA_EVEs_Counts_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_ssRNA_EVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_ssRNA_EVEs_Counts_ALL_endo$rope <- posterior_tab*100
Table_ssRNA_EVEs_Counts_ALL_endo$med <- median(Table_ssRNA_EVEs_Counts_ALL_endo$coef)
Table_ssRNA_EVEs_Counts_ALL_endo$pd <- pd(Table_ssRNA_EVEs_Counts_ALL_endo$coef2)

Table_ssRNA_EVEs_Counts_ALL_ecto <- as.data.frame(Table_ssRNA_EVEs_Counts_ALL$b_value3)
Table_ssRNA_EVEs_Counts_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_ssRNA_EVEs_Counts_ALL_ecto)<- c("coef","Parameter")
Table_ssRNA_EVEs_Counts_ALL_ecto$coef2 <- as.numeric(Table_ssRNA_EVEs_Counts_ALL_ecto$coef)
Table_ssRNA_EVEs_Counts_ALL_ecto$coef <- as.numeric(exp(Table_ssRNA_EVEs_Counts_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_ssRNA_EVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_ssRNA_EVEs_Counts_ALL_ecto$rope <- posterior_tab*100
Table_ssRNA_EVEs_Counts_ALL_ecto$med <- median(Table_ssRNA_EVEs_Counts_ALL_ecto$coef)
Table_ssRNA_EVEs_Counts_ALL_ecto$pd <- pd(Table_ssRNA_EVEs_Counts_ALL_ecto$coef2)



Table_ssRNA_EVEs_Counts_ALL_ecto_endo<-rbind(Table_ssRNA_EVEs_Counts_ALL_endo,Table_ssRNA_EVEs_Counts_ALL_ecto)

#Create summary table
statistics_Table_ssRNA_EVEs_Counts= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("ssRNA_EVEs_Counts",median(Table_ssRNA_EVEs_Counts_ALL_endo$coef),ci(Table_ssRNA_EVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_ssRNA_EVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_ssRNA_EVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_ssRNA_EVEs_Counts_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("ssRNA_EVEs_Counts",median(Table_ssRNA_EVEs_Counts_ALL_ecto$coef),ci(Table_ssRNA_EVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_ssRNA_EVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_ssRNA_EVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_ssRNA_EVEs_Counts_ALL_ecto$pd[1])))



colnames(statistics_Table_ssRNA_EVEs_Counts) <- statistics_Table_ssRNA_EVEs_Counts[1,]
statistics_Table_ssRNA_EVEs_Counts<-statistics_Table_ssRNA_EVEs_Counts[-1,]


Table_ssRNA_EVEs_Counts_ALL_ecto_endo$Parameter2 <- paste0(Table_ssRNA_EVEs_Counts_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_ssRNA_EVEs_Counts_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_ssRNA_EVEs_Counts_ALL<-Table_ssRNA_EVEs_Counts_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_ssRNA_EVEs_Counts_ALL <- stat.test_ssRNA_EVEs_Counts_ALL %>% add_xy_position(x = "Parameter2")
stat.test_ssRNA_EVEs_Counts_ALL$Prop <- paste0(round((nrow(Table_ssRNA_EVEs_Counts_ALL[Table_ssRNA_EVEs_Counts_ALL$b_value2 > Table_ssRNA_EVEs_Counts_ALL$b_value3,]) / nrow(Table_ssRNA_EVEs_Counts_ALL))*100,digits = 2),"%")

boxplot_ssRNA_EVEs_Counts <- ggplot(Table_ssRNA_EVEs_Counts_ALL_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="EVEs",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_ssRNA_EVEs_Counts_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,15))


################################
## dEVEs Counts all posteriors #
################################
Table_ssRNA_dEVEs_Counts_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_ssRNA_All/Posteriors_dEVEs_Counts_ALL.txt",sep=";",h=T)

Table_ssRNA_dEVEs_Counts_ALL$b_value2_coef_inf<- exp(as.numeric(Table_ssRNA_dEVEs_Counts_ALL$b_value2))
Table_ssRNA_dEVEs_Counts_ALL$b_value3_coef_inf<- exp(as.numeric(Table_ssRNA_dEVEs_Counts_ALL$b_value3))
Table_ssRNA_dEVEs_Counts_ALL<-Table_ssRNA_dEVEs_Counts_ALL[!Table_ssRNA_dEVEs_Counts_ALL$b_value2_coef_inf=="Inf",]
Table_ssRNA_dEVEs_Counts_ALL<-Table_ssRNA_dEVEs_Counts_ALL[!Table_ssRNA_dEVEs_Counts_ALL$b_value3_coef_inf=="Inf",]


Table_ssRNA_dEVEs_Counts_ALL_endo <- as.data.frame(Table_ssRNA_dEVEs_Counts_ALL$b_value2)
Table_ssRNA_dEVEs_Counts_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_ssRNA_dEVEs_Counts_ALL_endo)<- c("coef","Parameter")
Table_ssRNA_dEVEs_Counts_ALL_endo$coef2 <- as.numeric(Table_ssRNA_dEVEs_Counts_ALL_endo$coef)
Table_ssRNA_dEVEs_Counts_ALL_endo$coef <- as.numeric(exp(Table_ssRNA_dEVEs_Counts_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_ssRNA_dEVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_ssRNA_dEVEs_Counts_ALL_endo$rope <- posterior_tab*100
Table_ssRNA_dEVEs_Counts_ALL_endo$med <- median(Table_ssRNA_dEVEs_Counts_ALL_endo$coef)
Table_ssRNA_dEVEs_Counts_ALL_endo$pd <- pd(Table_ssRNA_dEVEs_Counts_ALL_endo$coef2)

Table_ssRNA_dEVEs_Counts_ALL_ecto <- as.data.frame(Table_ssRNA_dEVEs_Counts_ALL$b_value3)
Table_ssRNA_dEVEs_Counts_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_ssRNA_dEVEs_Counts_ALL_ecto)<- c("coef","Parameter")
Table_ssRNA_dEVEs_Counts_ALL_ecto$coef2 <- as.numeric(Table_ssRNA_dEVEs_Counts_ALL_ecto$coef)
Table_ssRNA_dEVEs_Counts_ALL_ecto$coef <- as.numeric(exp(Table_ssRNA_dEVEs_Counts_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_ssRNA_dEVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_ssRNA_dEVEs_Counts_ALL_ecto$rope <- posterior_tab*100
Table_ssRNA_dEVEs_Counts_ALL_ecto$med <- median(Table_ssRNA_dEVEs_Counts_ALL_ecto$coef)
Table_ssRNA_dEVEs_Counts_ALL_ecto$pd <- pd(Table_ssRNA_dEVEs_Counts_ALL_ecto$coef2)



Table_ssRNA_dEVEs_Counts_ALL_ecto_endo<-rbind(Table_ssRNA_dEVEs_Counts_ALL_endo,Table_ssRNA_dEVEs_Counts_ALL_ecto)

#Create summary table
statistics_Table_ssRNA_dEVEs_Counts= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("ssRNA_dEVEs_Counts",median(Table_ssRNA_dEVEs_Counts_ALL_endo$coef),ci(Table_ssRNA_dEVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_ssRNA_dEVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_ssRNA_dEVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_ssRNA_dEVEs_Counts_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("ssRNA_dEVEs_Counts",median(Table_ssRNA_dEVEs_Counts_ALL_ecto$coef),ci(Table_ssRNA_dEVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_ssRNA_dEVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_ssRNA_dEVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_ssRNA_dEVEs_Counts_ALL_ecto$pd[1])))



colnames(statistics_Table_ssRNA_dEVEs_Counts) <- statistics_Table_ssRNA_dEVEs_Counts[1,]
statistics_Table_ssRNA_dEVEs_Counts<-statistics_Table_ssRNA_dEVEs_Counts[-1,]


Table_ssRNA_dEVEs_Counts_ALL_ecto_endo$Parameter2 <- paste0(Table_ssRNA_dEVEs_Counts_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_ssRNA_dEVEs_Counts_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_ssRNA_dEVEs_Counts_ALL<-Table_ssRNA_dEVEs_Counts_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_ssRNA_dEVEs_Counts_ALL <- stat.test_ssRNA_dEVEs_Counts_ALL %>% add_xy_position(x = "Parameter2")
stat.test_ssRNA_dEVEs_Counts_ALL$Prop <- paste0(round((nrow(Table_ssRNA_dEVEs_Counts_ALL[Table_ssRNA_dEVEs_Counts_ALL$b_value2 > Table_ssRNA_dEVEs_Counts_ALL$b_value3,]) / nrow(Table_ssRNA_dEVEs_Counts_ALL))*100,digits = 2),"%")

boxplot_ssRNA_dEVEs_Counts <- ggplot(Table_ssRNA_dEVEs_Counts_ALL_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="dEVEs",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_ssRNA_dEVEs_Counts_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,15))


## Combine All 4 plots

violin_plot_ssRNA<-grid.arrange(arrangeGrob(
                             boxplot_ssRNA_EVEs_Counts  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank()) ,
                             boxplot_ssRNA_EVEs_Events  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank()),
                             boxplot_ssRNA_dEVEs_Counts  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank()) ,
                             boxplot_ssRNA_dEVEs_Events  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank()),
                         nrow = 1,
                         ncol=4,
                         left = textGrob("Relative rate", rot = 90, vjust = 1)),nrow=1,top=textGrob("ssRNA | A-B-C-D",
                                                                                                                 gp=gpar(fontsize=20,font=8)))






#Add posterior coefficient comparison between endo and ecto
statistics_Table_ssRNA_EVEs_Events<-as.data.frame(statistics_Table_ssRNA_EVEs_Events)
statistics_Table_ssRNA_EVEs_Events$Endo_sup_Ecto<-round((nrow(Table_ssRNA_EVEs_Events_ALL[Table_ssRNA_EVEs_Events_ALL$b_value2 > Table_ssRNA_EVEs_Events_ALL$b_value3,]) / nrow(Table_ssRNA_EVEs_Events_ALL))*100,digits = 2)

statistics_Table_ssRNA_dEVEs_Events<-as.data.frame(statistics_Table_ssRNA_dEVEs_Events)
statistics_Table_ssRNA_dEVEs_Events$Endo_sup_Ecto<-round((nrow(Table_ssRNA_dEVEs_Events_ALL[Table_ssRNA_dEVEs_Events_ALL$b_value2 > Table_ssRNA_dEVEs_Events_ALL$b_value3,]) / nrow(Table_ssRNA_dEVEs_Events_ALL))*100,digits = 2)

statistics_Table_ssRNA_EVEs_Counts<-as.data.frame(statistics_Table_ssRNA_EVEs_Counts)
statistics_Table_ssRNA_EVEs_Counts$Endo_sup_Ecto<-round((nrow(Table_ssRNA_EVEs_Counts_ALL[Table_ssRNA_EVEs_Counts_ALL$b_value2 > Table_ssRNA_EVEs_Counts_ALL$b_value3,]) / nrow(Table_ssRNA_EVEs_Counts_ALL))*100,digits = 2)


statistics_Table_ssRNA_dEVEs_Counts<-as.data.frame(statistics_Table_ssRNA_dEVEs_Counts)
statistics_Table_ssRNA_dEVEs_Counts$Endo_sup_Ecto<-round((nrow(Table_ssRNA_dEVEs_Counts_ALL[Table_ssRNA_dEVEs_Counts_ALL$b_value2 > Table_ssRNA_dEVEs_Counts_ALL$b_value3,]) / nrow(Table_ssRNA_dEVEs_Counts_ALL))*100,digits = 2)


## Combine all statistics tables

statistics_Table_ssRNA<-as.data.frame(rbind(statistics_Table_ssRNA_EVEs_Events,statistics_Table_ssRNA_dEVEs_Events,statistics_Table_ssRNA_EVEs_Counts,statistics_Table_ssRNA_dEVEs_Counts))

statistics_Table_ssRNA$Lifestyle<- gsub("\\..*","",rownames(statistics_Table_ssRNA))

library(dplyr)
statistics_Table_ssRNA$Median<-as.numeric(statistics_Table_ssRNA$Median)
statistics_Table_ssRNA$CI_low<-as.numeric(statistics_Table_ssRNA$CI_low)
statistics_Table_ssRNA$CI_high<-as.numeric(statistics_Table_ssRNA$CI_high)
statistics_Table_ssRNA$ROPE_Percentage<-as.numeric(statistics_Table_ssRNA$ROPE_Percentage)*100
statistics_Table_ssRNA<-statistics_Table_ssRNA %>% mutate_if(is.numeric, round, digits=3)

statistics_Table_ssRNA$Type <-c("EVEs_Events","EVEs_Events","dEVEs_Events","dEVEs_Events","EVEs_Numbers","EVEs_Numbers","dEVEs_Numbers","dEVEs_Numbers")

table_ssRNA<-statistics_Table_ssRNA
table_ssRNA[nrow(table_ssRNA)+1,] <- NA



#ssRNA_A A


library(bayestestR)
library(ggstatsplot)
library(bayesplot)
library(ggpubr)
library(gridExtra)
library(grid)



##################################
# ## EVEs Events all posteriors ##
##################################
Table_ssRNA_A_EVEs_Events_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_ssRNA_A_All/Posteriors_EVEs_Events_ALL.txt",sep=";",h=T)

Table_ssRNA_A_EVEs_Events_ALL$b_value2_coef_inf<- exp(as.numeric(Table_ssRNA_A_EVEs_Events_ALL$b_value2))
Table_ssRNA_A_EVEs_Events_ALL$b_value3_coef_inf<- exp(as.numeric(Table_ssRNA_A_EVEs_Events_ALL$b_value3))
Table_ssRNA_A_EVEs_Events_ALL<-Table_ssRNA_A_EVEs_Events_ALL[!Table_ssRNA_A_EVEs_Events_ALL$b_value2_coef_inf=="Inf",]
Table_ssRNA_A_EVEs_Events_ALL<-Table_ssRNA_A_EVEs_Events_ALL[!Table_ssRNA_A_EVEs_Events_ALL$b_value3_coef_inf=="Inf",]


Table_ssRNA_A_EVEs_Events_ALL_endo <- as.data.frame(Table_ssRNA_A_EVEs_Events_ALL$b_value2)
Table_ssRNA_A_EVEs_Events_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_ssRNA_A_EVEs_Events_ALL_endo)<- c("coef","Parameter")
Table_ssRNA_A_EVEs_Events_ALL_endo$coef2 <- as.numeric(Table_ssRNA_A_EVEs_Events_ALL_endo$coef)
Table_ssRNA_A_EVEs_Events_ALL_endo$coef <- as.numeric(exp(Table_ssRNA_A_EVEs_Events_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_ssRNA_A_EVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_ssRNA_A_EVEs_Events_ALL_endo$rope <- posterior_tab*100
Table_ssRNA_A_EVEs_Events_ALL_endo$med <- median(Table_ssRNA_A_EVEs_Events_ALL_endo$coef)
Table_ssRNA_A_EVEs_Events_ALL_endo$pd <- pd(Table_ssRNA_A_EVEs_Events_ALL_endo$coef2)

Table_ssRNA_A_EVEs_Events_ALL_ecto <- as.data.frame(Table_ssRNA_A_EVEs_Events_ALL$b_value3)
Table_ssRNA_A_EVEs_Events_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_ssRNA_A_EVEs_Events_ALL_ecto)<- c("coef","Parameter")
Table_ssRNA_A_EVEs_Events_ALL_ecto$coef2 <- as.numeric(Table_ssRNA_A_EVEs_Events_ALL_ecto$coef)
Table_ssRNA_A_EVEs_Events_ALL_ecto$coef <- as.numeric(exp(Table_ssRNA_A_EVEs_Events_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_ssRNA_A_EVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_ssRNA_A_EVEs_Events_ALL_ecto$rope <- posterior_tab*100
Table_ssRNA_A_EVEs_Events_ALL_ecto$med <- median(Table_ssRNA_A_EVEs_Events_ALL_ecto$coef)
Table_ssRNA_A_EVEs_Events_ALL_ecto$pd <- pd(Table_ssRNA_A_EVEs_Events_ALL_ecto$coef2)



Table_ssRNA_A_EVEs_Events_ALL_ecto_endo<-rbind(Table_ssRNA_A_EVEs_Events_ALL_endo,Table_ssRNA_A_EVEs_Events_ALL_ecto)

#Create summary table
statistics_Table_ssRNA_A_EVEs_Events= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("ssRNA_EVEs_Events",median(Table_ssRNA_A_EVEs_Events_ALL_endo$coef),ci(Table_ssRNA_A_EVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_ssRNA_A_EVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_ssRNA_A_EVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_ssRNA_A_EVEs_Events_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("ssRNA_EVEs_Events",median(Table_ssRNA_A_EVEs_Events_ALL_ecto$coef),ci(Table_ssRNA_A_EVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_ssRNA_A_EVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_ssRNA_A_EVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_ssRNA_A_EVEs_Events_ALL_ecto$pd[1])))



colnames(statistics_Table_ssRNA_A_EVEs_Events) <- statistics_Table_ssRNA_A_EVEs_Events[1,]
statistics_Table_ssRNA_A_EVEs_Events<-statistics_Table_ssRNA_A_EVEs_Events[-1,]


Table_ssRNA_A_EVEs_Events_ALL_ecto_endo$Parameter2 <- paste0(Table_ssRNA_A_EVEs_Events_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_ssRNA_A_EVEs_Events_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_ssRNA_A_EVEs_Events_ALL<-Table_ssRNA_A_EVEs_Events_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_ssRNA_A_EVEs_Events_ALL <- stat.test_ssRNA_A_EVEs_Events_ALL %>% add_xy_position(x = "Parameter2")
stat.test_ssRNA_A_EVEs_Events_ALL$Prop <- paste0(round((nrow(Table_ssRNA_A_EVEs_Events_ALL[Table_ssRNA_A_EVEs_Events_ALL$b_value2 > Table_ssRNA_A_EVEs_Events_ALL$b_value3,]) / nrow(Table_ssRNA_A_EVEs_Events_ALL))*100,digits = 2),"%")

boxplot_ssRNA_A_EVEs_Events <- ggplot(Table_ssRNA_A_EVEs_Events_ALL_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="Events",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_ssRNA_A_EVEs_Events_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,15))


##################################
## dEVEs Events all posteriors
##################################

Table_ssRNA_A_dEVEs_Events_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_ssRNA_A_All/Posteriors_dEVEs_Events_ALL.txt",sep=";",h=T)

Table_ssRNA_A_dEVEs_Events_ALL$b_value2_coef_inf<- exp(as.numeric(Table_ssRNA_A_dEVEs_Events_ALL$b_value2))
Table_ssRNA_A_dEVEs_Events_ALL$b_value3_coef_inf<- exp(as.numeric(Table_ssRNA_A_dEVEs_Events_ALL$b_value3))
Table_ssRNA_A_dEVEs_Events_ALL<-Table_ssRNA_A_dEVEs_Events_ALL[!Table_ssRNA_A_dEVEs_Events_ALL$b_value2_coef_inf=="Inf",]
Table_ssRNA_A_dEVEs_Events_ALL<-Table_ssRNA_A_dEVEs_Events_ALL[!Table_ssRNA_A_dEVEs_Events_ALL$b_value3_coef_inf=="Inf",]


Table_ssRNA_A_dEVEs_Events_ALL_endo <- as.data.frame(Table_ssRNA_A_dEVEs_Events_ALL$b_value2)
Table_ssRNA_A_dEVEs_Events_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_ssRNA_A_dEVEs_Events_ALL_endo)<- c("coef","Parameter")
Table_ssRNA_A_dEVEs_Events_ALL_endo$coef2 <- as.numeric(Table_ssRNA_A_dEVEs_Events_ALL_endo$coef)
Table_ssRNA_A_dEVEs_Events_ALL_endo$coef <- as.numeric(exp(Table_ssRNA_A_dEVEs_Events_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_ssRNA_A_dEVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_ssRNA_A_dEVEs_Events_ALL_endo$rope <- posterior_tab*100
Table_ssRNA_A_dEVEs_Events_ALL_endo$med <- median(Table_ssRNA_A_dEVEs_Events_ALL_endo$coef)
Table_ssRNA_A_dEVEs_Events_ALL_endo$pd <- pd(Table_ssRNA_A_dEVEs_Events_ALL_endo$coef2)

Table_ssRNA_A_dEVEs_Events_ALL_ecto <- as.data.frame(Table_ssRNA_A_dEVEs_Events_ALL$b_value3)
Table_ssRNA_A_dEVEs_Events_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_ssRNA_A_dEVEs_Events_ALL_ecto)<- c("coef","Parameter")
Table_ssRNA_A_dEVEs_Events_ALL_ecto$coef2 <- as.numeric(Table_ssRNA_A_dEVEs_Events_ALL_ecto$coef)
Table_ssRNA_A_dEVEs_Events_ALL_ecto$coef <- as.numeric(exp(Table_ssRNA_A_dEVEs_Events_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_ssRNA_A_dEVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_ssRNA_A_dEVEs_Events_ALL_ecto$rope <- posterior_tab*100
Table_ssRNA_A_dEVEs_Events_ALL_ecto$med <- median(Table_ssRNA_A_dEVEs_Events_ALL_ecto$coef)
Table_ssRNA_A_dEVEs_Events_ALL_ecto$pd <- pd(Table_ssRNA_A_dEVEs_Events_ALL_ecto$coef2)



Table_ssRNA_A_dEVEs_Events_ALL_ecto_endo<-rbind(Table_ssRNA_A_dEVEs_Events_ALL_endo,Table_ssRNA_A_dEVEs_Events_ALL_ecto)

#Create summary table
statistics_Table_ssRNA_A_dEVEs_Events= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("ssRNA_dEVEs_Events",median(Table_ssRNA_A_dEVEs_Events_ALL_endo$coef),ci(Table_ssRNA_A_dEVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_ssRNA_A_dEVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_ssRNA_A_dEVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_ssRNA_A_dEVEs_Events_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("ssRNA_dEVEs_Events",median(Table_ssRNA_A_dEVEs_Events_ALL_ecto$coef),ci(Table_ssRNA_A_dEVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_ssRNA_A_dEVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_ssRNA_A_dEVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_ssRNA_A_dEVEs_Events_ALL_ecto$pd[1])))



colnames(statistics_Table_ssRNA_A_dEVEs_Events) <- statistics_Table_ssRNA_A_dEVEs_Events[1,]
statistics_Table_ssRNA_A_dEVEs_Events<-statistics_Table_ssRNA_A_dEVEs_Events[-1,]


Table_ssRNA_A_dEVEs_Events_ALL_ecto_endo$Parameter2 <- paste0(Table_ssRNA_A_dEVEs_Events_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_ssRNA_A_dEVEs_Events_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_ssRNA_A_dEVEs_Events_ALL<-Table_ssRNA_A_dEVEs_Events_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_ssRNA_A_dEVEs_Events_ALL <- stat.test_ssRNA_A_dEVEs_Events_ALL %>% add_xy_position(x = "Parameter2")
stat.test_ssRNA_A_dEVEs_Events_ALL$Prop <- paste0(round((nrow(Table_ssRNA_A_dEVEs_Events_ALL[Table_ssRNA_A_dEVEs_Events_ALL$b_value2 > Table_ssRNA_A_dEVEs_Events_ALL$b_value3,]) / nrow(Table_ssRNA_A_dEVEs_Events_ALL))*100,digits = 2),"%")

boxplot_ssRNA_A_dEVEs_Events <- ggplot(Table_ssRNA_A_dEVEs_Events_ALL_ecto_endo[Table_ssRNA_A_dEVEs_Events_ALL_ecto_endo$coef<1000,], aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter),trim = T,adjust = 2)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="dEvents",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_ssRNA_A_dEVEs_Events_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,15))

##################################
## EVEs Counts all posteriors####
##################################


Table_ssRNA_A_EVEs_Counts_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_ssRNA_A_All/Posteriors_EVEs_Counts_ALL.txt",sep=";",h=T)

Table_ssRNA_A_EVEs_Counts_ALL$b_value2_coef_inf<- exp(as.numeric(Table_ssRNA_A_EVEs_Counts_ALL$b_value2))
Table_ssRNA_A_EVEs_Counts_ALL$b_value3_coef_inf<- exp(as.numeric(Table_ssRNA_A_EVEs_Counts_ALL$b_value3))
Table_ssRNA_A_EVEs_Counts_ALL<-Table_ssRNA_A_EVEs_Counts_ALL[!Table_ssRNA_A_EVEs_Counts_ALL$b_value2_coef_inf=="Inf",]
Table_ssRNA_A_EVEs_Counts_ALL<-Table_ssRNA_A_EVEs_Counts_ALL[!Table_ssRNA_A_EVEs_Counts_ALL$b_value3_coef_inf=="Inf",]


Table_ssRNA_A_EVEs_Counts_ALL_endo <- as.data.frame(Table_ssRNA_A_EVEs_Counts_ALL$b_value2)
Table_ssRNA_A_EVEs_Counts_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_ssRNA_A_EVEs_Counts_ALL_endo)<- c("coef","Parameter")
Table_ssRNA_A_EVEs_Counts_ALL_endo$coef2 <- as.numeric(Table_ssRNA_A_EVEs_Counts_ALL_endo$coef)
Table_ssRNA_A_EVEs_Counts_ALL_endo$coef <- as.numeric(exp(Table_ssRNA_A_EVEs_Counts_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_ssRNA_A_EVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_ssRNA_A_EVEs_Counts_ALL_endo$rope <- posterior_tab*100
Table_ssRNA_A_EVEs_Counts_ALL_endo$med <- median(Table_ssRNA_A_EVEs_Counts_ALL_endo$coef)
Table_ssRNA_A_EVEs_Counts_ALL_endo$pd <- pd(Table_ssRNA_A_EVEs_Counts_ALL_endo$coef2)

Table_ssRNA_A_EVEs_Counts_ALL_ecto <- as.data.frame(Table_ssRNA_A_EVEs_Counts_ALL$b_value3)
Table_ssRNA_A_EVEs_Counts_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_ssRNA_A_EVEs_Counts_ALL_ecto)<- c("coef","Parameter")
Table_ssRNA_A_EVEs_Counts_ALL_ecto$coef2 <- as.numeric(Table_ssRNA_A_EVEs_Counts_ALL_ecto$coef)
Table_ssRNA_A_EVEs_Counts_ALL_ecto$coef <- as.numeric(exp(Table_ssRNA_A_EVEs_Counts_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_ssRNA_A_EVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_ssRNA_A_EVEs_Counts_ALL_ecto$rope <- posterior_tab*100
Table_ssRNA_A_EVEs_Counts_ALL_ecto$med <- median(Table_ssRNA_A_EVEs_Counts_ALL_ecto$coef)
Table_ssRNA_A_EVEs_Counts_ALL_ecto$pd <- pd(Table_ssRNA_A_EVEs_Counts_ALL_ecto$coef2)



Table_ssRNA_A_EVEs_Counts_ALL_ecto_endo<-rbind(Table_ssRNA_A_EVEs_Counts_ALL_endo,Table_ssRNA_A_EVEs_Counts_ALL_ecto)

#Create summary table
statistics_Table_ssRNA_A_EVEs_Counts= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("ssRNA_EVEs_Counts",median(Table_ssRNA_A_EVEs_Counts_ALL_endo$coef),ci(Table_ssRNA_A_EVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_ssRNA_A_EVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_ssRNA_A_EVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_ssRNA_A_EVEs_Counts_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("ssRNA_EVEs_Counts",median(Table_ssRNA_A_EVEs_Counts_ALL_ecto$coef),ci(Table_ssRNA_A_EVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_ssRNA_A_EVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_ssRNA_A_EVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_ssRNA_A_EVEs_Counts_ALL_ecto$pd[1])))



colnames(statistics_Table_ssRNA_A_EVEs_Counts) <- statistics_Table_ssRNA_A_EVEs_Counts[1,]
statistics_Table_ssRNA_A_EVEs_Counts<-statistics_Table_ssRNA_A_EVEs_Counts[-1,]


Table_ssRNA_A_EVEs_Counts_ALL_ecto_endo$Parameter2 <- paste0(Table_ssRNA_A_EVEs_Counts_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_ssRNA_A_EVEs_Counts_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_ssRNA_A_EVEs_Counts_ALL<-Table_ssRNA_A_EVEs_Counts_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_ssRNA_A_EVEs_Counts_ALL <- stat.test_ssRNA_A_EVEs_Counts_ALL %>% add_xy_position(x = "Parameter2")
stat.test_ssRNA_A_EVEs_Counts_ALL$Prop <- paste0(round((nrow(Table_ssRNA_A_EVEs_Counts_ALL[Table_ssRNA_A_EVEs_Counts_ALL$b_value2 > Table_ssRNA_A_EVEs_Counts_ALL$b_value3,]) / nrow(Table_ssRNA_A_EVEs_Counts_ALL))*100,digits = 2),"%")

boxplot_ssRNA_A_EVEs_Counts <- ggplot(Table_ssRNA_A_EVEs_Counts_ALL_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="EVEs",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_ssRNA_A_EVEs_Counts_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,15))


################################
## dEVEs Counts all posteriors #
################################
Table_ssRNA_A_dEVEs_Counts_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_ssRNA_A_All/Posteriors_dEVEs_Counts_ALL.txt",sep=";",h=T)

Table_ssRNA_A_dEVEs_Counts_ALL$b_value2_coef_inf<- exp(as.numeric(Table_ssRNA_A_dEVEs_Counts_ALL$b_value2))
Table_ssRNA_A_dEVEs_Counts_ALL$b_value3_coef_inf<- exp(as.numeric(Table_ssRNA_A_dEVEs_Counts_ALL$b_value3))
Table_ssRNA_A_dEVEs_Counts_ALL<-Table_ssRNA_A_dEVEs_Counts_ALL[!Table_ssRNA_A_dEVEs_Counts_ALL$b_value2_coef_inf=="Inf",]
Table_ssRNA_A_dEVEs_Counts_ALL<-Table_ssRNA_A_dEVEs_Counts_ALL[!Table_ssRNA_A_dEVEs_Counts_ALL$b_value3_coef_inf=="Inf",]


Table_ssRNA_A_dEVEs_Counts_ALL_endo <- as.data.frame(Table_ssRNA_A_dEVEs_Counts_ALL$b_value2)
Table_ssRNA_A_dEVEs_Counts_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_ssRNA_A_dEVEs_Counts_ALL_endo)<- c("coef","Parameter")
Table_ssRNA_A_dEVEs_Counts_ALL_endo$coef2 <- as.numeric(Table_ssRNA_A_dEVEs_Counts_ALL_endo$coef)
Table_ssRNA_A_dEVEs_Counts_ALL_endo$coef <- as.numeric(exp(Table_ssRNA_A_dEVEs_Counts_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_ssRNA_A_dEVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_ssRNA_A_dEVEs_Counts_ALL_endo$rope <- posterior_tab*100
Table_ssRNA_A_dEVEs_Counts_ALL_endo$med <- median(Table_ssRNA_A_dEVEs_Counts_ALL_endo$coef)
Table_ssRNA_A_dEVEs_Counts_ALL_endo$pd <- pd(Table_ssRNA_A_dEVEs_Counts_ALL_endo$coef2)

Table_ssRNA_A_dEVEs_Counts_ALL_ecto <- as.data.frame(Table_ssRNA_A_dEVEs_Counts_ALL$b_value3)
Table_ssRNA_A_dEVEs_Counts_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_ssRNA_A_dEVEs_Counts_ALL_ecto)<- c("coef","Parameter")
Table_ssRNA_A_dEVEs_Counts_ALL_ecto$coef2 <- as.numeric(Table_ssRNA_A_dEVEs_Counts_ALL_ecto$coef)
Table_ssRNA_A_dEVEs_Counts_ALL_ecto$coef <- as.numeric(exp(Table_ssRNA_A_dEVEs_Counts_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_ssRNA_A_dEVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_ssRNA_A_dEVEs_Counts_ALL_ecto$rope <- posterior_tab*100
Table_ssRNA_A_dEVEs_Counts_ALL_ecto$med <- median(Table_ssRNA_A_dEVEs_Counts_ALL_ecto$coef)
Table_ssRNA_A_dEVEs_Counts_ALL_ecto$pd <- pd(Table_ssRNA_A_dEVEs_Counts_ALL_ecto$coef2)



Table_ssRNA_A_dEVEs_Counts_ALL_ecto_endo<-rbind(Table_ssRNA_A_dEVEs_Counts_ALL_endo,Table_ssRNA_A_dEVEs_Counts_ALL_ecto)

#Create summary table
statistics_Table_ssRNA_A_dEVEs_Counts= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("ssRNA_dEVEs_Counts",median(Table_ssRNA_A_dEVEs_Counts_ALL_endo$coef),ci(Table_ssRNA_A_dEVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_ssRNA_A_dEVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_ssRNA_A_dEVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_ssRNA_A_dEVEs_Counts_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("ssRNA_dEVEs_Counts",median(Table_ssRNA_A_dEVEs_Counts_ALL_ecto$coef),ci(Table_ssRNA_A_dEVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_ssRNA_A_dEVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_ssRNA_A_dEVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_ssRNA_A_dEVEs_Counts_ALL_ecto$pd[1])))



colnames(statistics_Table_ssRNA_A_dEVEs_Counts) <- statistics_Table_ssRNA_A_dEVEs_Counts[1,]
statistics_Table_ssRNA_A_dEVEs_Counts<-statistics_Table_ssRNA_A_dEVEs_Counts[-1,]


Table_ssRNA_A_dEVEs_Counts_ALL_ecto_endo$Parameter2 <- paste0(Table_ssRNA_A_dEVEs_Counts_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_ssRNA_A_dEVEs_Counts_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_ssRNA_A_dEVEs_Counts_ALL<-Table_ssRNA_A_dEVEs_Counts_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_ssRNA_A_dEVEs_Counts_ALL <- stat.test_ssRNA_A_dEVEs_Counts_ALL %>% add_xy_position(x = "Parameter2")
stat.test_ssRNA_A_dEVEs_Counts_ALL$Prop <- paste0(round((nrow(Table_ssRNA_A_dEVEs_Counts_ALL[Table_ssRNA_A_dEVEs_Counts_ALL$b_value2 > Table_ssRNA_A_dEVEs_Counts_ALL$b_value3,]) / nrow(Table_ssRNA_A_dEVEs_Counts_ALL))*100,digits = 2),"%")

boxplot_ssRNA_A_dEVEs_Counts <- ggplot(Table_ssRNA_A_dEVEs_Counts_ALL_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="dEVEs",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_ssRNA_A_dEVEs_Counts_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,15))




## Combine All 4 plots

violin_plot_ssRNA_A<-grid.arrange(arrangeGrob(
                             boxplot_ssRNA_A_EVEs_Counts  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank()) ,
                             boxplot_ssRNA_A_EVEs_Events  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank()),
                             boxplot_ssRNA_A_dEVEs_Counts  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank()) ,
                             boxplot_ssRNA_A_dEVEs_Events  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank()),
                         nrow = 1,
                         ncol=4,
                         left = textGrob("Relative rate", rot = 90, vjust = 1)),nrow=1,top=textGrob("ssRNA | A",
                                                                                                                 gp=gpar(fontsize=20,font=8)))




#Add posterior coefficient comparison between endo and ecto
statistics_Table_ssRNA_A_EVEs_Events<-as.data.frame(statistics_Table_ssRNA_A_EVEs_Events)
statistics_Table_ssRNA_A_EVEs_Events$Endo_sup_Ecto<-round((nrow(Table_ssRNA_A_EVEs_Events_ALL[Table_ssRNA_A_EVEs_Events_ALL$b_value2 > Table_ssRNA_A_EVEs_Events_ALL$b_value3,]) / nrow(Table_ssRNA_A_EVEs_Events_ALL))*100,digits = 2)

statistics_Table_ssRNA_A_dEVEs_Events<-as.data.frame(statistics_Table_ssRNA_A_dEVEs_Events)
statistics_Table_ssRNA_A_dEVEs_Events$Endo_sup_Ecto<-round((nrow(Table_ssRNA_A_dEVEs_Events_ALL[Table_ssRNA_A_dEVEs_Events_ALL$b_value2 > Table_ssRNA_A_dEVEs_Events_ALL$b_value3,]) / nrow(Table_ssRNA_A_dEVEs_Events_ALL))*100,digits = 2)

statistics_Table_ssRNA_A_EVEs_Counts<-as.data.frame(statistics_Table_ssRNA_A_EVEs_Counts)
statistics_Table_ssRNA_A_EVEs_Counts$Endo_sup_Ecto<-round((nrow(Table_ssRNA_A_EVEs_Counts_ALL[Table_ssRNA_A_EVEs_Counts_ALL$b_value2 > Table_ssRNA_A_EVEs_Counts_ALL$b_value3,]) / nrow(Table_ssRNA_A_EVEs_Counts_ALL))*100,digits = 2)


statistics_Table_ssRNA_A_dEVEs_Counts<-as.data.frame(statistics_Table_ssRNA_A_dEVEs_Counts)
statistics_Table_ssRNA_A_dEVEs_Counts$Endo_sup_Ecto<-round((nrow(Table_ssRNA_A_dEVEs_Counts_ALL[Table_ssRNA_A_dEVEs_Counts_ALL$b_value2 > Table_ssRNA_A_dEVEs_Counts_ALL$b_value3,]) / nrow(Table_ssRNA_A_dEVEs_Counts_ALL))*100,digits = 2)


## Combine all statistics tables

statistics_Table_ssRNA_A<-as.data.frame(rbind(statistics_Table_ssRNA_A_EVEs_Events,statistics_Table_ssRNA_A_dEVEs_Events,statistics_Table_ssRNA_A_EVEs_Counts,statistics_Table_ssRNA_A_dEVEs_Counts))

statistics_Table_ssRNA_A$Lifestyle<- gsub("\\..*","",rownames(statistics_Table_ssRNA_A))

library(dplyr)
statistics_Table_ssRNA_A$Median<-as.numeric(statistics_Table_ssRNA_A$Median)
statistics_Table_ssRNA_A$CI_low<-as.numeric(statistics_Table_ssRNA_A$CI_low)
statistics_Table_ssRNA_A$CI_high<-as.numeric(statistics_Table_ssRNA_A$CI_high)
statistics_Table_ssRNA_A$ROPE_Percentage<-as.numeric(statistics_Table_ssRNA_A$ROPE_Percentage)*100
statistics_Table_ssRNA_A<-statistics_Table_ssRNA_A %>% mutate_if(is.numeric, round, digits=3)

statistics_Table_ssRNA_A$Type <-c("EVEs_Events","EVEs_Events","dEVEs_Events","dEVEs_Events","EVEs_Numbers","EVEs_Numbers","dEVEs_Numbers","dEVEs_Numbers")

table_ssRNA_A<-statistics_Table_ssRNA_A
table_ssRNA_A[nrow(table_ssRNA_A)+1,] <- NA




#ssDNA A-B-C-D


library(bayestestR)
library(ggstatsplot)
library(bayesplot)
library(ggpubr)
library(gridExtra)
library(grid)
library(rstatix)



##################################
# ## EVEs Events all posteriors ##
##################################
Table_ssDNA_EVEs_Events_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_ssDNA_All/Posteriors_EVEs_Events_ALL.txt",sep=";",h=T)

Table_ssDNA_EVEs_Events_ALL$b_value2_coef_inf<- exp(as.numeric(Table_ssDNA_EVEs_Events_ALL$b_value2))
Table_ssDNA_EVEs_Events_ALL$b_value3_coef_inf<- exp(as.numeric(Table_ssDNA_EVEs_Events_ALL$b_value3))
Table_ssDNA_EVEs_Events_ALL<-Table_ssDNA_EVEs_Events_ALL[!Table_ssDNA_EVEs_Events_ALL$b_value2_coef_inf=="Inf",]
Table_ssDNA_EVEs_Events_ALL<-Table_ssDNA_EVEs_Events_ALL[!Table_ssDNA_EVEs_Events_ALL$b_value3_coef_inf=="Inf",]


Table_ssDNA_EVEs_Events_ALL_endo <- as.data.frame(Table_ssDNA_EVEs_Events_ALL$b_value2)
Table_ssDNA_EVEs_Events_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_ssDNA_EVEs_Events_ALL_endo)<- c("coef","Parameter")
Table_ssDNA_EVEs_Events_ALL_endo$coef2 <- as.numeric(Table_ssDNA_EVEs_Events_ALL_endo$coef)
Table_ssDNA_EVEs_Events_ALL_endo$coef <- as.numeric(exp(Table_ssDNA_EVEs_Events_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_ssDNA_EVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_ssDNA_EVEs_Events_ALL_endo$rope <- posterior_tab*100
Table_ssDNA_EVEs_Events_ALL_endo$med <- median(Table_ssDNA_EVEs_Events_ALL_endo$coef)
Table_ssDNA_EVEs_Events_ALL_endo$pd <- pd(Table_ssDNA_EVEs_Events_ALL_endo$coef2)

Table_ssDNA_EVEs_Events_ALL_ecto <- as.data.frame(Table_ssDNA_EVEs_Events_ALL$b_value3)
Table_ssDNA_EVEs_Events_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_ssDNA_EVEs_Events_ALL_ecto)<- c("coef","Parameter")
Table_ssDNA_EVEs_Events_ALL_ecto$coef2 <- as.numeric(Table_ssDNA_EVEs_Events_ALL_ecto$coef)
Table_ssDNA_EVEs_Events_ALL_ecto$coef <- as.numeric(exp(Table_ssDNA_EVEs_Events_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_ssDNA_EVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_ssDNA_EVEs_Events_ALL_ecto$rope <- posterior_tab*100
Table_ssDNA_EVEs_Events_ALL_ecto$med <- median(Table_ssDNA_EVEs_Events_ALL_ecto$coef)
Table_ssDNA_EVEs_Events_ALL_ecto$pd <- pd(Table_ssDNA_EVEs_Events_ALL_ecto$coef2)



Table_ssDNA_EVEs_Events_ALL_ecto_endo<-rbind(Table_ssDNA_EVEs_Events_ALL_endo,Table_ssDNA_EVEs_Events_ALL_ecto)

#Create summary table
statistics_Table_ssDNA_EVEs_Events= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("ssRNA_EVEs_Events",median(Table_ssDNA_EVEs_Events_ALL_endo$coef),ci(Table_ssDNA_EVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_ssDNA_EVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_ssDNA_EVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_ssDNA_EVEs_Events_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("ssRNA_EVEs_Events",median(Table_ssDNA_EVEs_Events_ALL_ecto$coef),ci(Table_ssDNA_EVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_ssDNA_EVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_ssDNA_EVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_ssDNA_EVEs_Events_ALL_ecto$pd[1])))



colnames(statistics_Table_ssDNA_EVEs_Events) <- statistics_Table_ssDNA_EVEs_Events[1,]
statistics_Table_ssDNA_EVEs_Events<-statistics_Table_ssDNA_EVEs_Events[-1,]


Table_ssDNA_EVEs_Events_ALL_ecto_endo$Parameter2 <- paste0(Table_ssDNA_EVEs_Events_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_ssDNA_EVEs_Events_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_ssDNA_EVEs_Events_ALL<-Table_ssDNA_EVEs_Events_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_ssDNA_EVEs_Events_ALL <- stat.test_ssDNA_EVEs_Events_ALL %>% add_xy_position(x = "Parameter2")
stat.test_ssDNA_EVEs_Events_ALL$Prop <- paste0(round((nrow(Table_ssDNA_EVEs_Events_ALL[Table_ssDNA_EVEs_Events_ALL$b_value2 > Table_ssDNA_EVEs_Events_ALL$b_value3,]) / nrow(Table_ssDNA_EVEs_Events_ALL))*100,digits = 2),"%")

boxplot_ssDNA_EVEs_Events <- ggplot(Table_ssDNA_EVEs_Events_ALL_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="Events",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_ssDNA_EVEs_Events_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,15))


##################################
## dEVEs Events all posteriors
##################################

Table_ssDNA_dEVEs_Events_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_ssDNA_All/Posteriors_dEVEs_Events_ALL.txt",sep=";",h=T)

Table_ssDNA_dEVEs_Events_ALL$b_value2_coef_inf<- exp(as.numeric(Table_ssDNA_dEVEs_Events_ALL$b_value2))
Table_ssDNA_dEVEs_Events_ALL$b_value3_coef_inf<- exp(as.numeric(Table_ssDNA_dEVEs_Events_ALL$b_value3))
Table_ssDNA_dEVEs_Events_ALL<-Table_ssDNA_dEVEs_Events_ALL[!Table_ssDNA_dEVEs_Events_ALL$b_value2_coef_inf=="Inf",]
Table_ssDNA_dEVEs_Events_ALL<-Table_ssDNA_dEVEs_Events_ALL[!Table_ssDNA_dEVEs_Events_ALL$b_value3_coef_inf=="Inf",]


Table_ssDNA_dEVEs_Events_ALL_endo <- as.data.frame(Table_ssDNA_dEVEs_Events_ALL$b_value2)
Table_ssDNA_dEVEs_Events_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_ssDNA_dEVEs_Events_ALL_endo)<- c("coef","Parameter")
Table_ssDNA_dEVEs_Events_ALL_endo$coef2 <- as.numeric(Table_ssDNA_dEVEs_Events_ALL_endo$coef)
Table_ssDNA_dEVEs_Events_ALL_endo$coef <- as.numeric(exp(Table_ssDNA_dEVEs_Events_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_ssDNA_dEVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_ssDNA_dEVEs_Events_ALL_endo$rope <- posterior_tab*100
Table_ssDNA_dEVEs_Events_ALL_endo$med <- median(Table_ssDNA_dEVEs_Events_ALL_endo$coef)
Table_ssDNA_dEVEs_Events_ALL_endo$pd <- pd(Table_ssDNA_dEVEs_Events_ALL_endo$coef2)

Table_ssDNA_dEVEs_Events_ALL_ecto <- as.data.frame(Table_ssDNA_dEVEs_Events_ALL$b_value3)
Table_ssDNA_dEVEs_Events_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_ssDNA_dEVEs_Events_ALL_ecto)<- c("coef","Parameter")
Table_ssDNA_dEVEs_Events_ALL_ecto$coef2 <- as.numeric(Table_ssDNA_dEVEs_Events_ALL_ecto$coef)
Table_ssDNA_dEVEs_Events_ALL_ecto$coef <- as.numeric(exp(Table_ssDNA_dEVEs_Events_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_ssDNA_dEVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_ssDNA_dEVEs_Events_ALL_ecto$rope <- posterior_tab*100
Table_ssDNA_dEVEs_Events_ALL_ecto$med <- median(Table_ssDNA_dEVEs_Events_ALL_ecto$coef)
Table_ssDNA_dEVEs_Events_ALL_ecto$pd <- pd(Table_ssDNA_dEVEs_Events_ALL_ecto$coef2)



Table_ssDNA_dEVEs_Events_ALL_ecto_endo<-rbind(Table_ssDNA_dEVEs_Events_ALL_endo,Table_ssDNA_dEVEs_Events_ALL_ecto)

#Create summary table
statistics_Table_ssDNA_dEVEs_Events= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("ssRNA_dEVEs_Events",median(Table_ssDNA_dEVEs_Events_ALL_endo$coef),ci(Table_ssDNA_dEVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_ssDNA_dEVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_ssDNA_dEVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_ssDNA_dEVEs_Events_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("ssRNA_dEVEs_Events",median(Table_ssDNA_dEVEs_Events_ALL_ecto$coef),ci(Table_ssDNA_dEVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_ssDNA_dEVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_ssDNA_dEVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_ssDNA_dEVEs_Events_ALL_ecto$pd[1])))



colnames(statistics_Table_ssDNA_dEVEs_Events) <- statistics_Table_ssDNA_dEVEs_Events[1,]
statistics_Table_ssDNA_dEVEs_Events<-statistics_Table_ssDNA_dEVEs_Events[-1,]


Table_ssDNA_dEVEs_Events_ALL_ecto_endo$Parameter2 <- paste0(Table_ssDNA_dEVEs_Events_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_ssDNA_dEVEs_Events_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_ssDNA_dEVEs_Events_ALL<-Table_ssDNA_dEVEs_Events_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_ssDNA_dEVEs_Events_ALL <- stat.test_ssDNA_dEVEs_Events_ALL %>% add_xy_position(x = "Parameter2")
stat.test_ssDNA_dEVEs_Events_ALL$Prop <- paste0(round((nrow(Table_ssDNA_dEVEs_Events_ALL[Table_ssDNA_dEVEs_Events_ALL$b_value2 > Table_ssDNA_dEVEs_Events_ALL$b_value3,]) / nrow(Table_ssDNA_dEVEs_Events_ALL))*100,digits = 2),"%")

boxplot_ssDNA_dEVEs_Events <- ggplot(Table_ssDNA_dEVEs_Events_ALL_ecto_endo[Table_ssDNA_dEVEs_Events_ALL_ecto_endo$coef<1000,], aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter),trim = T,adjust = 2)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="dEvents",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_ssDNA_dEVEs_Events_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,15))

##################################
## EVEs Counts all posteriors####
##################################


Table_ssDNA_EVEs_Counts_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_ssDNA_All/Posteriors_EVEs_Counts_ALL.txt",sep=";",h=T)

Table_ssDNA_EVEs_Counts_ALL$b_value2_coef_inf<- exp(as.numeric(Table_ssDNA_EVEs_Counts_ALL$b_value2))
Table_ssDNA_EVEs_Counts_ALL$b_value3_coef_inf<- exp(as.numeric(Table_ssDNA_EVEs_Counts_ALL$b_value3))
Table_ssDNA_EVEs_Counts_ALL<-Table_ssDNA_EVEs_Counts_ALL[!Table_ssDNA_EVEs_Counts_ALL$b_value2_coef_inf=="Inf",]
Table_ssDNA_EVEs_Counts_ALL<-Table_ssDNA_EVEs_Counts_ALL[!Table_ssDNA_EVEs_Counts_ALL$b_value3_coef_inf=="Inf",]


Table_ssDNA_EVEs_Counts_ALL_endo <- as.data.frame(Table_ssDNA_EVEs_Counts_ALL$b_value2)
Table_ssDNA_EVEs_Counts_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_ssDNA_EVEs_Counts_ALL_endo)<- c("coef","Parameter")
Table_ssDNA_EVEs_Counts_ALL_endo$coef2 <- as.numeric(Table_ssDNA_EVEs_Counts_ALL_endo$coef)
Table_ssDNA_EVEs_Counts_ALL_endo$coef <- as.numeric(exp(Table_ssDNA_EVEs_Counts_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_ssDNA_EVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_ssDNA_EVEs_Counts_ALL_endo$rope <- posterior_tab*100
Table_ssDNA_EVEs_Counts_ALL_endo$med <- median(Table_ssDNA_EVEs_Counts_ALL_endo$coef)
Table_ssDNA_EVEs_Counts_ALL_endo$pd <- pd(Table_ssDNA_EVEs_Counts_ALL_endo$coef2)

Table_ssDNA_EVEs_Counts_ALL_ecto <- as.data.frame(Table_ssDNA_EVEs_Counts_ALL$b_value3)
Table_ssDNA_EVEs_Counts_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_ssDNA_EVEs_Counts_ALL_ecto)<- c("coef","Parameter")
Table_ssDNA_EVEs_Counts_ALL_ecto$coef2 <- as.numeric(Table_ssDNA_EVEs_Counts_ALL_ecto$coef)
Table_ssDNA_EVEs_Counts_ALL_ecto$coef <- as.numeric(exp(Table_ssDNA_EVEs_Counts_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_ssDNA_EVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_ssDNA_EVEs_Counts_ALL_ecto$rope <- posterior_tab*100
Table_ssDNA_EVEs_Counts_ALL_ecto$med <- median(Table_ssDNA_EVEs_Counts_ALL_ecto$coef)
Table_ssDNA_EVEs_Counts_ALL_ecto$pd <- pd(Table_ssDNA_EVEs_Counts_ALL_ecto$coef2)



Table_ssDNA_EVEs_Counts_ALL_ecto_endo<-rbind(Table_ssDNA_EVEs_Counts_ALL_endo,Table_ssDNA_EVEs_Counts_ALL_ecto)

#Create summary table
statistics_Table_ssDNA_EVEs_Counts= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("ssRNA_EVEs_Counts",median(Table_ssDNA_EVEs_Counts_ALL_endo$coef),ci(Table_ssDNA_EVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_ssDNA_EVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_ssDNA_EVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_ssDNA_EVEs_Counts_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("ssRNA_EVEs_Counts",median(Table_ssDNA_EVEs_Counts_ALL_ecto$coef),ci(Table_ssDNA_EVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_ssDNA_EVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_ssDNA_EVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_ssDNA_EVEs_Counts_ALL_ecto$pd[1])))



colnames(statistics_Table_ssDNA_EVEs_Counts) <- statistics_Table_ssDNA_EVEs_Counts[1,]
statistics_Table_ssDNA_EVEs_Counts<-statistics_Table_ssDNA_EVEs_Counts[-1,]


Table_ssDNA_EVEs_Counts_ALL_ecto_endo$Parameter2 <- paste0(Table_ssDNA_EVEs_Counts_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_ssDNA_EVEs_Counts_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_ssDNA_EVEs_Counts_ALL<-Table_ssDNA_EVEs_Counts_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_ssDNA_EVEs_Counts_ALL <- stat.test_ssDNA_EVEs_Counts_ALL %>% add_xy_position(x = "Parameter2")
stat.test_ssDNA_EVEs_Counts_ALL$Prop <- paste0(round((nrow(Table_ssDNA_EVEs_Counts_ALL[Table_ssDNA_EVEs_Counts_ALL$b_value2 > Table_ssDNA_EVEs_Counts_ALL$b_value3,]) / nrow(Table_ssDNA_EVEs_Counts_ALL))*100,digits = 2),"%")

boxplot_ssDNA_EVEs_Counts <- ggplot(Table_ssDNA_EVEs_Counts_ALL_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="EVEs",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_ssDNA_EVEs_Counts_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,15))


################################
## dEVEs Counts all posteriors #
################################
Table_ssDNA_dEVEs_Counts_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_ssDNA_All/Posteriors_dEVEs_Counts_ALL.txt",sep=";",h=T)

Table_ssDNA_dEVEs_Counts_ALL$b_value2_coef_inf<- exp(as.numeric(Table_ssDNA_dEVEs_Counts_ALL$b_value2))
Table_ssDNA_dEVEs_Counts_ALL$b_value3_coef_inf<- exp(as.numeric(Table_ssDNA_dEVEs_Counts_ALL$b_value3))
Table_ssDNA_dEVEs_Counts_ALL<-Table_ssDNA_dEVEs_Counts_ALL[!Table_ssDNA_dEVEs_Counts_ALL$b_value2_coef_inf=="Inf",]
Table_ssDNA_dEVEs_Counts_ALL<-Table_ssDNA_dEVEs_Counts_ALL[!Table_ssDNA_dEVEs_Counts_ALL$b_value3_coef_inf=="Inf",]


Table_ssDNA_dEVEs_Counts_ALL_endo <- as.data.frame(Table_ssDNA_dEVEs_Counts_ALL$b_value2)
Table_ssDNA_dEVEs_Counts_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_ssDNA_dEVEs_Counts_ALL_endo)<- c("coef","Parameter")
Table_ssDNA_dEVEs_Counts_ALL_endo$coef2 <- as.numeric(Table_ssDNA_dEVEs_Counts_ALL_endo$coef)
Table_ssDNA_dEVEs_Counts_ALL_endo$coef <- as.numeric(exp(Table_ssDNA_dEVEs_Counts_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_ssDNA_dEVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_ssDNA_dEVEs_Counts_ALL_endo$rope <- posterior_tab*100
Table_ssDNA_dEVEs_Counts_ALL_endo$med <- median(Table_ssDNA_dEVEs_Counts_ALL_endo$coef)
Table_ssDNA_dEVEs_Counts_ALL_endo$pd <- pd(Table_ssDNA_dEVEs_Counts_ALL_endo$coef2)

Table_ssDNA_dEVEs_Counts_ALL_ecto <- as.data.frame(Table_ssDNA_dEVEs_Counts_ALL$b_value3)
Table_ssDNA_dEVEs_Counts_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_ssDNA_dEVEs_Counts_ALL_ecto)<- c("coef","Parameter")
Table_ssDNA_dEVEs_Counts_ALL_ecto$coef2 <- as.numeric(Table_ssDNA_dEVEs_Counts_ALL_ecto$coef)
Table_ssDNA_dEVEs_Counts_ALL_ecto$coef <- as.numeric(exp(Table_ssDNA_dEVEs_Counts_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_ssDNA_dEVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_ssDNA_dEVEs_Counts_ALL_ecto$rope <- posterior_tab*100
Table_ssDNA_dEVEs_Counts_ALL_ecto$med <- median(Table_ssDNA_dEVEs_Counts_ALL_ecto$coef)
Table_ssDNA_dEVEs_Counts_ALL_ecto$pd <- pd(Table_ssDNA_dEVEs_Counts_ALL_ecto$coef2)



Table_ssDNA_dEVEs_Counts_ALL_ecto_endo<-rbind(Table_ssDNA_dEVEs_Counts_ALL_endo,Table_ssDNA_dEVEs_Counts_ALL_ecto)

#Create summary table
statistics_Table_ssDNA_dEVEs_Counts= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("ssRNA_dEVEs_Counts",median(Table_ssDNA_dEVEs_Counts_ALL_endo$coef),ci(Table_ssDNA_dEVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_ssDNA_dEVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_ssDNA_dEVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_ssDNA_dEVEs_Counts_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("ssRNA_dEVEs_Counts",median(Table_ssDNA_dEVEs_Counts_ALL_ecto$coef),ci(Table_ssDNA_dEVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_ssDNA_dEVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_ssDNA_dEVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_ssDNA_dEVEs_Counts_ALL_ecto$pd[1])))



colnames(statistics_Table_ssDNA_dEVEs_Counts) <- statistics_Table_ssDNA_dEVEs_Counts[1,]
statistics_Table_ssDNA_dEVEs_Counts<-statistics_Table_ssDNA_dEVEs_Counts[-1,]


Table_ssDNA_dEVEs_Counts_ALL_ecto_endo$Parameter2 <- paste0(Table_ssDNA_dEVEs_Counts_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_ssDNA_dEVEs_Counts_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_ssDNA_dEVEs_Counts_ALL<-Table_ssDNA_dEVEs_Counts_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_ssDNA_dEVEs_Counts_ALL <- stat.test_ssDNA_dEVEs_Counts_ALL %>% add_xy_position(x = "Parameter2")
stat.test_ssDNA_dEVEs_Counts_ALL$Prop <- paste0(round((nrow(Table_ssDNA_dEVEs_Counts_ALL[Table_ssDNA_dEVEs_Counts_ALL$b_value2 > Table_ssDNA_dEVEs_Counts_ALL$b_value3,]) / nrow(Table_ssDNA_dEVEs_Counts_ALL))*100,digits = 2),"%")

boxplot_ssDNA_dEVEs_Counts <- ggplot(Table_ssDNA_dEVEs_Counts_ALL_ecto_endo[Table_ssDNA_dEVEs_Counts_ALL_ecto_endo$coef<1000,], aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="dEVEs",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_ssDNA_dEVEs_Counts_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,15))



## Combine All 4 plots

violin_plot_ssDNA<-grid.arrange(arrangeGrob(
                             boxplot_ssDNA_EVEs_Counts  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank()) ,
                             boxplot_ssDNA_EVEs_Events  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank()),
                             boxplot_ssDNA_dEVEs_Counts  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank()) ,
                             boxplot_ssDNA_dEVEs_Events  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank()),
                         nrow = 1,
                         ncol=4,
                         left = textGrob("Relative rate", rot = 90, vjust = 1)),nrow=1,top=textGrob("ssDNA | A-B-C-D",
                                                                                                                 gp=gpar(fontsize=20,font=8)))




#Add posterior coefficient comparison between endo and ecto
statistics_Table_ssDNA_EVEs_Events<-as.data.frame(statistics_Table_ssDNA_EVEs_Events)
statistics_Table_ssDNA_EVEs_Events$Endo_sup_Ecto<-round((nrow(Table_ssDNA_EVEs_Events_ALL[Table_ssDNA_EVEs_Events_ALL$b_value2 > Table_ssDNA_EVEs_Events_ALL$b_value3,]) / nrow(Table_ssDNA_EVEs_Events_ALL))*100,digits = 2)

statistics_Table_ssDNA_dEVEs_Events<-as.data.frame(statistics_Table_ssDNA_dEVEs_Events)
statistics_Table_ssDNA_dEVEs_Events$Endo_sup_Ecto<-round((nrow(Table_ssDNA_dEVEs_Events_ALL[Table_ssDNA_dEVEs_Events_ALL$b_value2 > Table_ssDNA_dEVEs_Events_ALL$b_value3,]) / nrow(Table_ssDNA_dEVEs_Events_ALL))*100,digits = 2)

statistics_Table_ssDNA_EVEs_Counts<-as.data.frame(statistics_Table_ssDNA_EVEs_Counts)
statistics_Table_ssDNA_EVEs_Counts$Endo_sup_Ecto<-round((nrow(Table_ssDNA_EVEs_Counts_ALL[Table_ssDNA_EVEs_Counts_ALL$b_value2 > Table_ssDNA_EVEs_Counts_ALL$b_value3,]) / nrow(Table_ssDNA_EVEs_Counts_ALL))*100,digits = 2)


statistics_Table_ssDNA_dEVEs_Counts<-as.data.frame(statistics_Table_ssDNA_dEVEs_Counts)
statistics_Table_ssDNA_dEVEs_Counts$Endo_sup_Ecto<-round((nrow(Table_ssDNA_dEVEs_Counts_ALL[Table_ssDNA_dEVEs_Counts_ALL$b_value2 > Table_ssDNA_dEVEs_Counts_ALL$b_value3,]) / nrow(Table_ssDNA_dEVEs_Counts_ALL))*100,digits = 2)


## Combine all statistics tables

statistics_Table_ssDNA<-as.data.frame(rbind(statistics_Table_ssDNA_EVEs_Events,statistics_Table_ssDNA_dEVEs_Events,statistics_Table_ssDNA_EVEs_Counts,statistics_Table_ssDNA_dEVEs_Counts))

statistics_Table_ssDNA$Lifestyle<- gsub("\\..*","",rownames(statistics_Table_ssDNA))

library(dplyr)
statistics_Table_ssDNA$Median<-as.numeric(statistics_Table_ssDNA$Median)
statistics_Table_ssDNA$CI_low<-as.numeric(statistics_Table_ssDNA$CI_low)
statistics_Table_ssDNA$CI_high<-as.numeric(statistics_Table_ssDNA$CI_high)
statistics_Table_ssDNA$ROPE_Percentage<-as.numeric(statistics_Table_ssDNA$ROPE_Percentage)*100
statistics_Table_ssDNA<-statistics_Table_ssDNA %>% mutate_if(is.numeric, round, digits=3)

statistics_Table_ssDNA$Type <-c("EVEs_Events","EVEs_Events","dEVEs_Events","dEVEs_Events","EVEs_Numbers","EVEs_Numbers","dEVEs_Numbers","dEVEs_Numbers")

table_ssDNA<-statistics_Table_ssDNA
table_ssDNA[nrow(table_ssDNA)+1,] <- NA






#ssDNA_A A


library(bayestestR)
library(ggstatsplot)
library(bayesplot)
library(ggpubr)
library(gridExtra)
library(grid)



##################################
# ## EVEs Events all posteriors ##
##################################
Table_ssDNA_A_EVEs_Events_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_ssDNA_A_All/Posteriors_EVEs_Events_ALL.txt",sep=";",h=T)

Table_ssDNA_A_EVEs_Events_ALL$b_value2_coef_inf<- exp(as.numeric(Table_ssDNA_A_EVEs_Events_ALL$b_value2))
Table_ssDNA_A_EVEs_Events_ALL$b_value3_coef_inf<- exp(as.numeric(Table_ssDNA_A_EVEs_Events_ALL$b_value3))
Table_ssDNA_A_EVEs_Events_ALL<-Table_ssDNA_A_EVEs_Events_ALL[!Table_ssDNA_A_EVEs_Events_ALL$b_value2_coef_inf=="Inf",]
Table_ssDNA_A_EVEs_Events_ALL<-Table_ssDNA_A_EVEs_Events_ALL[!Table_ssDNA_A_EVEs_Events_ALL$b_value3_coef_inf=="Inf",]


Table_ssDNA_A_EVEs_Events_ALL_endo <- as.data.frame(Table_ssDNA_A_EVEs_Events_ALL$b_value2)
Table_ssDNA_A_EVEs_Events_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_ssDNA_A_EVEs_Events_ALL_endo)<- c("coef","Parameter")
Table_ssDNA_A_EVEs_Events_ALL_endo$coef2 <- as.numeric(Table_ssDNA_A_EVEs_Events_ALL_endo$coef)
Table_ssDNA_A_EVEs_Events_ALL_endo$coef <- as.numeric(exp(Table_ssDNA_A_EVEs_Events_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_ssDNA_A_EVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_ssDNA_A_EVEs_Events_ALL_endo$rope <- posterior_tab*100
Table_ssDNA_A_EVEs_Events_ALL_endo$med <- median(Table_ssDNA_A_EVEs_Events_ALL_endo$coef)
Table_ssDNA_A_EVEs_Events_ALL_endo$pd <- pd(Table_ssDNA_A_EVEs_Events_ALL_endo$coef2)

Table_ssDNA_A_EVEs_Events_ALL_ecto <- as.data.frame(Table_ssDNA_A_EVEs_Events_ALL$b_value3)
Table_ssDNA_A_EVEs_Events_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_ssDNA_A_EVEs_Events_ALL_ecto)<- c("coef","Parameter")
Table_ssDNA_A_EVEs_Events_ALL_ecto$coef2 <- as.numeric(Table_ssDNA_A_EVEs_Events_ALL_ecto$coef)
Table_ssDNA_A_EVEs_Events_ALL_ecto$coef <- as.numeric(exp(Table_ssDNA_A_EVEs_Events_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_ssDNA_A_EVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_ssDNA_A_EVEs_Events_ALL_ecto$rope <- posterior_tab*100
Table_ssDNA_A_EVEs_Events_ALL_ecto$med <- median(Table_ssDNA_A_EVEs_Events_ALL_ecto$coef)
Table_ssDNA_A_EVEs_Events_ALL_ecto$pd <- pd(Table_ssDNA_A_EVEs_Events_ALL_ecto$coef2)



Table_ssDNA_A_EVEs_Events_ALL_ecto_endo<-rbind(Table_ssDNA_A_EVEs_Events_ALL_endo,Table_ssDNA_A_EVEs_Events_ALL_ecto)

#Create summary table
statistics_Table_ssDNA_A_EVEs_Events= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("ssRNA_EVEs_Events",median(Table_ssDNA_A_EVEs_Events_ALL_endo$coef),ci(Table_ssDNA_A_EVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_ssDNA_A_EVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_ssDNA_A_EVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_ssDNA_A_EVEs_Events_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("ssRNA_EVEs_Events",median(Table_ssDNA_A_EVEs_Events_ALL_ecto$coef),ci(Table_ssDNA_A_EVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_ssDNA_A_EVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_ssDNA_A_EVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_ssDNA_A_EVEs_Events_ALL_ecto$pd[1])))



colnames(statistics_Table_ssDNA_A_EVEs_Events) <- statistics_Table_ssDNA_A_EVEs_Events[1,]
statistics_Table_ssDNA_A_EVEs_Events<-statistics_Table_ssDNA_A_EVEs_Events[-1,]


Table_ssDNA_A_EVEs_Events_ALL_ecto_endo$Parameter2 <- paste0(Table_ssDNA_A_EVEs_Events_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_ssDNA_A_EVEs_Events_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_ssDNA_A_EVEs_Events_ALL<-Table_ssDNA_A_EVEs_Events_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_ssDNA_A_EVEs_Events_ALL <- stat.test_ssDNA_A_EVEs_Events_ALL %>% add_xy_position(x = "Parameter2")
stat.test_ssDNA_A_EVEs_Events_ALL$Prop <- paste0(round((nrow(Table_ssDNA_A_EVEs_Events_ALL[Table_ssDNA_A_EVEs_Events_ALL$b_value2 > Table_ssDNA_A_EVEs_Events_ALL$b_value3,]) / nrow(Table_ssDNA_A_EVEs_Events_ALL))*100,digits = 2),"%")

boxplot_ssDNA_A_EVEs_Events <- ggplot(Table_ssDNA_A_EVEs_Events_ALL_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="Events",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_ssDNA_A_EVEs_Events_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,15))


##################################
## dEVEs Events all posteriors
##################################

Table_ssDNA_A_dEVEs_Events_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_ssDNA_A_All/Posteriors_dEVEs_Events_ALL.txt",sep=";",h=T)

Table_ssDNA_A_dEVEs_Events_ALL$b_value2_coef_inf<- exp(as.numeric(Table_ssDNA_A_dEVEs_Events_ALL$b_value2))
Table_ssDNA_A_dEVEs_Events_ALL$b_value3_coef_inf<- exp(as.numeric(Table_ssDNA_A_dEVEs_Events_ALL$b_value3))
Table_ssDNA_A_dEVEs_Events_ALL<-Table_ssDNA_A_dEVEs_Events_ALL[!Table_ssDNA_A_dEVEs_Events_ALL$b_value2_coef_inf=="Inf",]
Table_ssDNA_A_dEVEs_Events_ALL<-Table_ssDNA_A_dEVEs_Events_ALL[!Table_ssDNA_A_dEVEs_Events_ALL$b_value3_coef_inf=="Inf",]


Table_ssDNA_A_dEVEs_Events_ALL_endo <- as.data.frame(Table_ssDNA_A_dEVEs_Events_ALL$b_value2)
Table_ssDNA_A_dEVEs_Events_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_ssDNA_A_dEVEs_Events_ALL_endo)<- c("coef","Parameter")
Table_ssDNA_A_dEVEs_Events_ALL_endo$coef2 <- as.numeric(Table_ssDNA_A_dEVEs_Events_ALL_endo$coef)
Table_ssDNA_A_dEVEs_Events_ALL_endo$coef <- as.numeric(exp(Table_ssDNA_A_dEVEs_Events_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_ssDNA_A_dEVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_ssDNA_A_dEVEs_Events_ALL_endo$rope <- posterior_tab*100
Table_ssDNA_A_dEVEs_Events_ALL_endo$med <- median(Table_ssDNA_A_dEVEs_Events_ALL_endo$coef)
Table_ssDNA_A_dEVEs_Events_ALL_endo$pd <- pd(Table_ssDNA_A_dEVEs_Events_ALL_endo$coef2)

Table_ssDNA_A_dEVEs_Events_ALL_ecto <- as.data.frame(Table_ssDNA_A_dEVEs_Events_ALL$b_value3)
Table_ssDNA_A_dEVEs_Events_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_ssDNA_A_dEVEs_Events_ALL_ecto)<- c("coef","Parameter")
Table_ssDNA_A_dEVEs_Events_ALL_ecto$coef2 <- as.numeric(Table_ssDNA_A_dEVEs_Events_ALL_ecto$coef)
Table_ssDNA_A_dEVEs_Events_ALL_ecto$coef <- as.numeric(exp(Table_ssDNA_A_dEVEs_Events_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_ssDNA_A_dEVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_ssDNA_A_dEVEs_Events_ALL_ecto$rope <- posterior_tab*100
Table_ssDNA_A_dEVEs_Events_ALL_ecto$med <- median(Table_ssDNA_A_dEVEs_Events_ALL_ecto$coef)
Table_ssDNA_A_dEVEs_Events_ALL_ecto$pd <- pd(Table_ssDNA_A_dEVEs_Events_ALL_ecto$coef2)



Table_ssDNA_A_dEVEs_Events_ALL_ecto_endo<-rbind(Table_ssDNA_A_dEVEs_Events_ALL_endo,Table_ssDNA_A_dEVEs_Events_ALL_ecto)

#Create summary table
statistics_Table_ssDNA_A_dEVEs_Events= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("ssRNA_dEVEs_Events",median(Table_ssDNA_A_dEVEs_Events_ALL_endo$coef),ci(Table_ssDNA_A_dEVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_ssDNA_A_dEVEs_Events_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_ssDNA_A_dEVEs_Events_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_ssDNA_A_dEVEs_Events_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("ssRNA_dEVEs_Events",median(Table_ssDNA_A_dEVEs_Events_ALL_ecto$coef),ci(Table_ssDNA_A_dEVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_ssDNA_A_dEVEs_Events_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_ssDNA_A_dEVEs_Events_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_ssDNA_A_dEVEs_Events_ALL_ecto$pd[1])))



colnames(statistics_Table_ssDNA_A_dEVEs_Events) <- statistics_Table_ssDNA_A_dEVEs_Events[1,]
statistics_Table_ssDNA_A_dEVEs_Events<-statistics_Table_ssDNA_A_dEVEs_Events[-1,]


Table_ssDNA_A_dEVEs_Events_ALL_ecto_endo$Parameter2 <- paste0(Table_ssDNA_A_dEVEs_Events_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_ssDNA_A_dEVEs_Events_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_ssDNA_A_dEVEs_Events_ALL<-Table_ssDNA_A_dEVEs_Events_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_ssDNA_A_dEVEs_Events_ALL <- stat.test_ssDNA_A_dEVEs_Events_ALL %>% add_xy_position(x = "Parameter2")
stat.test_ssDNA_A_dEVEs_Events_ALL$Prop <- paste0(round((nrow(Table_ssDNA_A_dEVEs_Events_ALL[Table_ssDNA_A_dEVEs_Events_ALL$b_value2 > Table_ssDNA_A_dEVEs_Events_ALL$b_value3,]) / nrow(Table_ssDNA_A_dEVEs_Events_ALL))*100,digits = 2),"%")

boxplot_ssDNA_A_dEVEs_Events <- ggplot(Table_ssDNA_A_dEVEs_Events_ALL_ecto_endo[Table_ssDNA_A_dEVEs_Events_ALL_ecto_endo$coef<1000,], aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter),trim = T,adjust = 2)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="dEvents",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_ssDNA_A_dEVEs_Events_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,15))

##################################
## EVEs Counts all posteriors####
##################################


Table_ssDNA_A_EVEs_Counts_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_ssDNA_A_All/Posteriors_EVEs_Counts_ALL.txt",sep=";",h=T)

Table_ssDNA_A_EVEs_Counts_ALL$b_value2_coef_inf<- exp(as.numeric(Table_ssDNA_A_EVEs_Counts_ALL$b_value2))
Table_ssDNA_A_EVEs_Counts_ALL$b_value3_coef_inf<- exp(as.numeric(Table_ssDNA_A_EVEs_Counts_ALL$b_value3))
Table_ssDNA_A_EVEs_Counts_ALL<-Table_ssDNA_A_EVEs_Counts_ALL[!Table_ssDNA_A_EVEs_Counts_ALL$b_value2_coef_inf=="Inf",]
Table_ssDNA_A_EVEs_Counts_ALL<-Table_ssDNA_A_EVEs_Counts_ALL[!Table_ssDNA_A_EVEs_Counts_ALL$b_value3_coef_inf=="Inf",]


Table_ssDNA_A_EVEs_Counts_ALL_endo <- as.data.frame(Table_ssDNA_A_EVEs_Counts_ALL$b_value2)
Table_ssDNA_A_EVEs_Counts_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_ssDNA_A_EVEs_Counts_ALL_endo)<- c("coef","Parameter")
Table_ssDNA_A_EVEs_Counts_ALL_endo$coef2 <- as.numeric(Table_ssDNA_A_EVEs_Counts_ALL_endo$coef)
Table_ssDNA_A_EVEs_Counts_ALL_endo$coef <- as.numeric(exp(Table_ssDNA_A_EVEs_Counts_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_ssDNA_A_EVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_ssDNA_A_EVEs_Counts_ALL_endo$rope <- posterior_tab*100
Table_ssDNA_A_EVEs_Counts_ALL_endo$med <- median(Table_ssDNA_A_EVEs_Counts_ALL_endo$coef)
Table_ssDNA_A_EVEs_Counts_ALL_endo$pd <- pd(Table_ssDNA_A_EVEs_Counts_ALL_endo$coef2)

Table_ssDNA_A_EVEs_Counts_ALL_ecto <- as.data.frame(Table_ssDNA_A_EVEs_Counts_ALL$b_value3)
Table_ssDNA_A_EVEs_Counts_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_ssDNA_A_EVEs_Counts_ALL_ecto)<- c("coef","Parameter")
Table_ssDNA_A_EVEs_Counts_ALL_ecto$coef2 <- as.numeric(Table_ssDNA_A_EVEs_Counts_ALL_ecto$coef)
Table_ssDNA_A_EVEs_Counts_ALL_ecto$coef <- as.numeric(exp(Table_ssDNA_A_EVEs_Counts_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_ssDNA_A_EVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_ssDNA_A_EVEs_Counts_ALL_ecto$rope <- posterior_tab*100
Table_ssDNA_A_EVEs_Counts_ALL_ecto$med <- median(Table_ssDNA_A_EVEs_Counts_ALL_ecto$coef)
Table_ssDNA_A_EVEs_Counts_ALL_ecto$pd <- pd(Table_ssDNA_A_EVEs_Counts_ALL_ecto$coef2)



Table_ssDNA_A_EVEs_Counts_ALL_ecto_endo<-rbind(Table_ssDNA_A_EVEs_Counts_ALL_endo,Table_ssDNA_A_EVEs_Counts_ALL_ecto)

#Create summary table
statistics_Table_ssDNA_A_EVEs_Counts= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("ssRNA_EVEs_Counts",median(Table_ssDNA_A_EVEs_Counts_ALL_endo$coef),ci(Table_ssDNA_A_EVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_ssDNA_A_EVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_ssDNA_A_EVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_ssDNA_A_EVEs_Counts_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("ssRNA_EVEs_Counts",median(Table_ssDNA_A_EVEs_Counts_ALL_ecto$coef),ci(Table_ssDNA_A_EVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_ssDNA_A_EVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_ssDNA_A_EVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_ssDNA_A_EVEs_Counts_ALL_ecto$pd[1])))



colnames(statistics_Table_ssDNA_A_EVEs_Counts) <- statistics_Table_ssDNA_A_EVEs_Counts[1,]
statistics_Table_ssDNA_A_EVEs_Counts<-statistics_Table_ssDNA_A_EVEs_Counts[-1,]


Table_ssDNA_A_EVEs_Counts_ALL_ecto_endo$Parameter2 <- paste0(Table_ssDNA_A_EVEs_Counts_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_ssDNA_A_EVEs_Counts_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_ssDNA_A_EVEs_Counts_ALL<-Table_ssDNA_A_EVEs_Counts_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_ssDNA_A_EVEs_Counts_ALL <- stat.test_ssDNA_A_EVEs_Counts_ALL %>% add_xy_position(x = "Parameter2")
stat.test_ssDNA_A_EVEs_Counts_ALL$Prop <- paste0(round((nrow(Table_ssDNA_A_EVEs_Counts_ALL[Table_ssDNA_A_EVEs_Counts_ALL$b_value2 > Table_ssDNA_A_EVEs_Counts_ALL$b_value3,]) / nrow(Table_ssDNA_A_EVEs_Counts_ALL))*100,digits = 2),"%")

boxplot_ssDNA_A_EVEs_Counts <- ggplot(Table_ssDNA_A_EVEs_Counts_ALL_ecto_endo, aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="EVEs",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_ssDNA_A_EVEs_Counts_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,15))


################################
## dEVEs Counts all posteriors #
################################
Table_ssDNA_A_dEVEs_Counts_ALL<- read.table("/Users/bguinet/Desktop/Papier_scientifique/Gathered_results5/output_ssDNA_A_All/Posteriors_dEVEs_Counts_ALL.txt",sep=";",h=T)

Table_ssDNA_A_dEVEs_Counts_ALL$b_value2_coef_inf<- exp(as.numeric(Table_ssDNA_A_dEVEs_Counts_ALL$b_value2))
Table_ssDNA_A_dEVEs_Counts_ALL$b_value3_coef_inf<- exp(as.numeric(Table_ssDNA_A_dEVEs_Counts_ALL$b_value3))
Table_ssDNA_A_dEVEs_Counts_ALL<-Table_ssDNA_A_dEVEs_Counts_ALL[!Table_ssDNA_A_dEVEs_Counts_ALL$b_value2_coef_inf=="Inf",]
Table_ssDNA_A_dEVEs_Counts_ALL<-Table_ssDNA_A_dEVEs_Counts_ALL[!Table_ssDNA_A_dEVEs_Counts_ALL$b_value3_coef_inf=="Inf",]


Table_ssDNA_A_dEVEs_Counts_ALL_endo <- as.data.frame(Table_ssDNA_A_dEVEs_Counts_ALL$b_value2)
Table_ssDNA_A_dEVEs_Counts_ALL_endo$Parameter <- "Endo-parasitoid"
colnames(Table_ssDNA_A_dEVEs_Counts_ALL_endo)<- c("coef","Parameter")
Table_ssDNA_A_dEVEs_Counts_ALL_endo$coef2 <- as.numeric(Table_ssDNA_A_dEVEs_Counts_ALL_endo$coef)
Table_ssDNA_A_dEVEs_Counts_ALL_endo$coef <- as.numeric(exp(Table_ssDNA_A_dEVEs_Counts_ALL_endo$coef))
posterior_tab<- as.data.frame(rope(Table_ssDNA_A_dEVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_ssDNA_A_dEVEs_Counts_ALL_endo$rope <- posterior_tab*100
Table_ssDNA_A_dEVEs_Counts_ALL_endo$med <- median(Table_ssDNA_A_dEVEs_Counts_ALL_endo$coef)
Table_ssDNA_A_dEVEs_Counts_ALL_endo$pd <- pd(Table_ssDNA_A_dEVEs_Counts_ALL_endo$coef2)

Table_ssDNA_A_dEVEs_Counts_ALL_ecto <- as.data.frame(Table_ssDNA_A_dEVEs_Counts_ALL$b_value3)
Table_ssDNA_A_dEVEs_Counts_ALL_ecto$Parameter <- "Ecto-parasitoid"
colnames(Table_ssDNA_A_dEVEs_Counts_ALL_ecto)<- c("coef","Parameter")
Table_ssDNA_A_dEVEs_Counts_ALL_ecto$coef2 <- as.numeric(Table_ssDNA_A_dEVEs_Counts_ALL_ecto$coef)
Table_ssDNA_A_dEVEs_Counts_ALL_ecto$coef <- as.numeric(exp(Table_ssDNA_A_dEVEs_Counts_ALL_ecto$coef))
posterior_tab<- as.data.frame(rope(Table_ssDNA_A_dEVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1))$ROPE_Percentage
Table_ssDNA_A_dEVEs_Counts_ALL_ecto$rope <- posterior_tab*100
Table_ssDNA_A_dEVEs_Counts_ALL_ecto$med <- median(Table_ssDNA_A_dEVEs_Counts_ALL_ecto$coef)
Table_ssDNA_A_dEVEs_Counts_ALL_ecto$pd <- pd(Table_ssDNA_A_dEVEs_Counts_ALL_ecto$coef2)



Table_ssDNA_A_dEVEs_Counts_ALL_ecto_endo<-rbind(Table_ssDNA_A_dEVEs_Counts_ALL_endo,Table_ssDNA_A_dEVEs_Counts_ALL_ecto)

#Create summary table
statistics_Table_ssDNA_A_dEVEs_Counts= t(data.frame(Type= c('Type','Median',"CI_low","CI_high","ROPE_Percentage","pd"),
                             Endo_parasitoid=c("ssRNA_dEVEs_Counts",median(Table_ssDNA_A_dEVEs_Counts_ALL_endo$coef),ci(Table_ssDNA_A_dEVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_ssDNA_A_dEVEs_Counts_ALL_endo$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_ssDNA_A_dEVEs_Counts_ALL_endo$coef,range=c(0,1), ci = 1)$ROPE_Percentage , Table_ssDNA_A_dEVEs_Counts_ALL_endo$pd[1]),
                             Ecto_parasitoid=c("ssRNA_dEVEs_Counts",median(Table_ssDNA_A_dEVEs_Counts_ALL_ecto$coef),ci(Table_ssDNA_A_dEVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_low,ci(Table_ssDNA_A_dEVEs_Counts_ALL_ecto$coef, method = "HDI",ci=0.89)$CI_high,rope(Table_ssDNA_A_dEVEs_Counts_ALL_ecto$coef,range=c(0,1), ci = 1)$ROPE_Percentage, Table_ssDNA_A_dEVEs_Counts_ALL_ecto$pd[1])))



colnames(statistics_Table_ssDNA_A_dEVEs_Counts) <- statistics_Table_ssDNA_A_dEVEs_Counts[1,]
statistics_Table_ssDNA_A_dEVEs_Counts<-statistics_Table_ssDNA_A_dEVEs_Counts[-1,]


Table_ssDNA_A_dEVEs_Counts_ALL_ecto_endo$Parameter2 <- paste0(Table_ssDNA_A_dEVEs_Counts_ALL_ecto_endo$Parameter,"\n %pd = ",round(Table_ssDNA_A_dEVEs_Counts_ALL_ecto_endo$pd*100,2) ,"% ")

stat.test_ssDNA_A_dEVEs_Counts_ALL<-Table_ssDNA_A_dEVEs_Counts_ALL_ecto_endo %>%
  t_test(coef ~ Parameter2) %>%
  add_significance()

stat.test_ssDNA_A_dEVEs_Counts_ALL <- stat.test_ssDNA_A_dEVEs_Counts_ALL %>% add_xy_position(x = "Parameter2")
stat.test_ssDNA_A_dEVEs_Counts_ALL$Prop <- paste0(round((nrow(Table_ssDNA_A_dEVEs_Counts_ALL[Table_ssDNA_A_dEVEs_Counts_ALL$b_value2 > Table_ssDNA_A_dEVEs_Counts_ALL$b_value3,]) / nrow(Table_ssDNA_A_dEVEs_Counts_ALL))*100,digits = 2),"%")

boxplot_ssDNA_A_dEVEs_Counts <- ggplot(Table_ssDNA_A_dEVEs_Counts_ALL_ecto_endo[Table_ssDNA_A_dEVEs_Counts_ALL_ecto_endo$coef<1000,], aes(x=Parameter2, y=coef)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0, ymax=1, alpha=0.4)+ geom_hline(yintercept=1, color = '#1FA1D0')+
  geom_violin(aes(fill=Parameter))+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="dEVEs",y="Relative rate", x ="")+
  scale_fill_manual(values = c("#FFC20A","#10A401"))+   stat_pvalue_manual(
  stat.test_ssDNA_A_dEVEs_Counts_ALL, label = "PMCMC = {Prop}",
  vjust = -1, bracket.nudge.y = 1,y.position = 10, size=5
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(name="Relative rate",n.breaks = 20)  + theme(plot.title = element_text(size=22)) + theme(legend.position="none")+coord_cartesian(ylim = c(0,15))



#Add posterior coefficient comparison between endo and ecto
statistics_Table_ssDNA_A_EVEs_Events<-as.data.frame(statistics_Table_ssDNA_A_EVEs_Events)
statistics_Table_ssDNA_A_EVEs_Events$Endo_sup_Ecto<-round((nrow(Table_ssDNA_A_EVEs_Events_ALL[Table_ssDNA_A_EVEs_Events_ALL$b_value2 > Table_ssDNA_A_EVEs_Events_ALL$b_value3,]) / nrow(Table_ssDNA_A_EVEs_Events_ALL))*100,digits = 2)

statistics_Table_ssDNA_A_dEVEs_Events<-as.data.frame(statistics_Table_ssDNA_A_dEVEs_Events)
statistics_Table_ssDNA_A_dEVEs_Events$Endo_sup_Ecto<-round((nrow(Table_ssDNA_A_dEVEs_Events_ALL[Table_ssDNA_A_dEVEs_Events_ALL$b_value2 > Table_ssDNA_A_dEVEs_Events_ALL$b_value3,]) / nrow(Table_ssDNA_A_dEVEs_Events_ALL))*100,digits = 2)

statistics_Table_ssDNA_A_EVEs_Counts<-as.data.frame(statistics_Table_ssDNA_A_EVEs_Counts)
statistics_Table_ssDNA_A_EVEs_Counts$Endo_sup_Ecto<-round((nrow(Table_ssDNA_A_EVEs_Counts_ALL[Table_ssDNA_A_EVEs_Counts_ALL$b_value2 > Table_ssDNA_A_EVEs_Counts_ALL$b_value3,]) / nrow(Table_ssDNA_A_EVEs_Counts_ALL))*100,digits = 2)


statistics_Table_ssDNA_A_dEVEs_Counts<-as.data.frame(statistics_Table_ssDNA_A_dEVEs_Counts)
statistics_Table_ssDNA_A_dEVEs_Counts$Endo_sup_Ecto<-round((nrow(Table_ssDNA_A_dEVEs_Counts_ALL[Table_ssDNA_A_dEVEs_Counts_ALL$b_value2 > Table_ssDNA_A_dEVEs_Counts_ALL$b_value3,]) / nrow(Table_ssDNA_A_dEVEs_Counts_ALL))*100,digits = 2)

## Combine All 4 plots

violin_plot_ssDNA_A<-grid.arrange(arrangeGrob(
                             boxplot_ssDNA_A_EVEs_Counts  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank()) ,
                             boxplot_ssDNA_A_EVEs_Events  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank()),
                             boxplot_ssDNA_A_dEVEs_Counts  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank()) ,
                             boxplot_ssDNA_A_dEVEs_Events  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank()),
                         nrow = 1,
                         ncol=4,
                         left = textGrob("Relative rate", rot = 90, vjust = 1)),nrow=1,top=textGrob("ssDNA | A",
                                                                                                                 gp=gpar(fontsize=20,font=8)))



## Combine all statistics tables

statistics_Table_ssDNA_A<-as.data.frame(rbind(statistics_Table_ssDNA_A_EVEs_Events,statistics_Table_ssDNA_A_dEVEs_Events,statistics_Table_ssDNA_A_EVEs_Counts,statistics_Table_ssDNA_A_dEVEs_Counts))

statistics_Table_ssDNA_A$Lifestyle<- gsub("\\..*","",rownames(statistics_Table_ssDNA_A))

library(dplyr)
statistics_Table_ssDNA_A$Median<-as.numeric(statistics_Table_ssDNA_A$Median)
statistics_Table_ssDNA_A$CI_low<-as.numeric(statistics_Table_ssDNA_A$CI_low)
statistics_Table_ssDNA_A$CI_high<-as.numeric(statistics_Table_ssDNA_A$CI_high)
statistics_Table_ssDNA_A$ROPE_Percentage<-as.numeric(statistics_Table_ssDNA_A$ROPE_Percentage)*100
statistics_Table_ssDNA_A<-statistics_Table_ssDNA_A %>% mutate_if(is.numeric, round, digits=3)

statistics_Table_ssDNA_A$Type <-c("EVEs_Events","EVEs_Events","dEVEs_Events","dEVEs_Events","EVEs_Numbers","EVEs_Numbers","dEVEs_Numbers","dEVEs_Numbers")

table_ssDNA_A<-statistics_Table_ssDNA_A
table_ssDNA_A[nrow(table_ssDNA_A)+1,] <- NA




# All EVEs violin plots

violin_plot_EVEs<-grid.arrange(arrangeGrob(
                             boxplot_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Counts+ theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank())+ labs(title="A"),
                             boxplot_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Counts + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank())+ labs(title="") ,
                             boxplot_dsDNA_EVEs_Counts  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank())+ labs(title="E") ,
                             boxplot_dsDNA_A_EVEs_Counts  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank())+ labs(title=""),
                             boxplot_dsDNA_without_controls_EVEs_Counts + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank())+ labs(title="I"),
                             boxplot_dsDNA_without_controls_A_EVEs_Counts + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank())+ labs(title=""),
                             boxplot_ssDNA_EVEs_Counts  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank())+ labs(title="M") ,
                             boxplot_ssDNA_A_EVEs_Counts  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank())+ labs(title=""),
                             boxplot_dsRNA_EVEs_Counts  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank())+ labs(title="Q"),
                             boxplot_dsRNA_A_EVEs_Counts  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank())+ labs(title=""),
                             boxplot_ssRNA_EVEs_Counts + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank())+ labs(title="U"),
                             boxplot_ssRNA_A_EVEs_Counts + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank())+ labs(title=""),
                         nrow = 6,
                         ncol=2,
                         left = textGrob("Relative rate", rot = 90, vjust = 1)),
                         nrow=1,top=textGrob("EVEs",gp=gpar(fontsize=20,font=8)))

violin_plot_Events<-grid.arrange(arrangeGrob(
                             boxplot_dsDNA_ssDNA_dsRNA_ssRNA_EVEs_Events+ theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank())+ labs(title="B") ,
                             boxplot_dsDNA_ssDNA_dsRNA_ssRNA_A_EVEs_Events + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank())+ labs(title="") ,
                             boxplot_dsDNA_EVEs_Events  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank())+ labs(title="F") ,
                             boxplot_dsDNA_A_EVEs_Events  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank())+ labs(title=""),
                             boxplot_dsDNA_without_controls_EVEs_Events + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank())+ labs(title="J"),
                             boxplot_dsDNA_without_controls_A_EVEs_Events + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank())+ labs(title=""),
                             boxplot_ssDNA_EVEs_Events  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank())+ labs(title="N") ,
                             boxplot_ssDNA_A_EVEs_Events  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank())+ labs(title=""),
                             boxplot_dsRNA_EVEs_Events  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank())+ labs(title="R"),
                             boxplot_dsRNA_A_EVEs_Events  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank())+ labs(title=""),
                             boxplot_ssRNA_EVEs_Events + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank())+ labs(title="V"),
                             boxplot_ssRNA_A_EVEs_Events + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank())+ labs(title=""),
                         nrow = 6,
                         ncol=2,
                         left = textGrob("Relative rate", rot = 90, vjust = 1)),nrow=1,top=textGrob("Events",gp=gpar(fontsize=20,font=8)))


violin_plot_dEVEs<-grid.arrange(arrangeGrob(
                             boxplot_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Counts+ theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank())+ labs(title="C") ,
                             boxplot_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Counts + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank())+ labs(title="") ,
                             boxplot_dsDNA_dEVEs_Counts  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank())+ labs(title="G") ,
                             boxplot_dsDNA_A_dEVEs_Counts  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank())+ labs(title=""),
                             boxplot_dsDNA_without_controls_dEVEs_Counts + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank())+ labs(title="K"),
                             boxplot_dsDNA_without_controls_A_dEVEs_Counts + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank())+ labs(title=""),
                             boxplot_ssDNA_dEVEs_Counts  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank())+ labs(title="O") ,
                             boxplot_ssDNA_A_dEVEs_Counts  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank())+ labs(title=""),
                             boxplot_dsRNA_dEVEs_Counts  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank())+ labs(title="S"),
                             boxplot_dsRNA_A_dEVEs_Counts  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank())+ labs(title=""),
                             boxplot_ssRNA_dEVEs_Counts + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank())+ labs(title="W"),
                             boxplot_ssRNA_A_dEVEs_Counts + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank())+ labs(title=""),
                         nrow = 6,
                         ncol=2,
                         left = textGrob("Relative rate", rot = 90, vjust = 1)),nrow=1,top=textGrob("dEVEs",gp=gpar(fontsize=20,font=8)))

violin_plot_dEvents<-grid.arrange(arrangeGrob(
                             boxplot_dsDNA_ssDNA_dsRNA_ssRNA_dEVEs_Events+ theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank())+ labs(title="D") ,
                             boxplot_dsDNA_ssDNA_dsRNA_ssRNA_A_dEVEs_Events + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank())+ labs(title="") ,
                             boxplot_dsDNA_dEVEs_Events  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank())+ labs(title="H") ,
                             boxplot_dsDNA_A_dEVEs_Events  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank())+ labs(title=""),
                             boxplot_dsDNA_without_controls_dEVEs_Events + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank())+ labs(title="L"),
                             boxplot_dsDNA_without_controls_A_dEVEs_Events + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank())+ labs(title=""),
                             boxplot_ssDNA_dEVEs_Events  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank())+ labs(title="P") ,
                             boxplot_ssDNA_A_dEVEs_Events  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank())+ labs(title=""),
                             boxplot_dsRNA_dEVEs_Events  + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank())+ labs(title="T"),
                             boxplot_dsRNA_A_dEVEs_Events  + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank())+ labs(title=""),
                             boxplot_ssRNA_dEVEs_Events + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank())+ labs(title="X"),
                             boxplot_ssRNA_A_dEVEs_Events + theme(axis.title.x=element_blank())+ theme(axis.title.y=element_blank())+ labs(title=""),
                         nrow = 6,
                         ncol=2,
                         left = textGrob("Relative rate", rot = 90, vjust = 1)),nrow=1,top=textGrob("dEvents",gp=gpar(fontsize=20,font=8)))



```
#ENDS

#Gather all violing plot together

library(cowplot)
#plot_grid(violin_plot_dsDNA,violin_plot_dsDNA_A,violin_plot_dsDNA_without_controls, violin_plot_dsDNA_without_controls_A,ncol=2,nrow=2, labels=LETTERS[1:4])


#allviolin<-plot_grid(violin_plot_dsDNA_ssDNA_dsRNA_ssRNA,violin_plot_dsDNA_ssDNA_dsRNA_ssRNA_A,violin_plot_dsDNA,violin_plot_dsDNA_A,violin_plot_dsDNA_without_controls, violin_plot_dsDNA_without_controls_A,violin_plot_ssDNA,violin_plot_ssDNA_A,violin_plot_dsRNA,violin_plot_dsRNA_A,violin_plot_ssRNA,violin_plot_ssRNA_A,ncol=2,nrow=6,labels=LETTERS[1:12],label_size = 16)

allviolin<-plot_grid(violin_plot_EVEs,violin_plot_Events,violin_plot_dEVEs,violin_plot_dEvents, ncol=4,nrow=1,label_size = 16,scale = 0.9)


ggsave(filename="/Users/bguinet/Desktop/Papier_scientifique/All_violin_plots_new_order.pdf",allviolin, width = 33.08, height = 46.76, units = "in",device="pdf")



# size 33.08 , 46.76

#save with

#Gather all statistics tables

table_ssDNA$type<-"ssDNA"
table_ssDNA_A$type<-"ssDNA_A"
table_ssRNA$type<-"ssRNA"
table_ssRNA_A$type<-"ssRNA_A"
table_dsRNA$type<-"dsRNA"
table_dsRNA_A$type<-"dsRNA_A"
table_dsDNA$type<-"dsDNA"
table_dsDNA_A$type<-"dsDNA_A"
table_dsDNA_without_controls_A$type <- "dsDNA_A_without_controls"
table_dsDNA_ssDNA_dsRNA_ssRNA$type <-"dsDNA_ssDNA_dsRNA_ssRNA"
table_dsDNA_ssDNA_dsRNA_ssRNA_A$type <-"dsDNA_ssDNA_dsRNA_ssRNA_A"
table_dsDNA_without_controls$type <- "dsDNA_without_controls"

All_statistic_table<-rbind(table_dsDNA_ssDNA_dsRNA_ssRNA,table_dsDNA_ssDNA_dsRNA_ssRNA_A,table_ssDNA,table_ssDNA_A,table_ssRNA,table_ssRNA_A,table_dsRNA,table_dsRNA_A,table_dsDNA,table_dsDNA_A,table_dsDNA_without_controls_A,table_dsDNA_without_controls)

All_statistic_table$pd<-as.numeric(All_statistic_table$pd)*100

All_statistic_table$pvalue<- 2 * (1 - All_statistic_table$pd /100)


write.table(All_statistic_table,"/Users/bguinet/Desktop/Papier_scientifique/All_coef_statistics_tables2.txt",sep=";")
