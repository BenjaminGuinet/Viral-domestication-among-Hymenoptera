#!/usr/bin/env Rscript
#library(dplyr)

if (length(args)==0) {
  stop("Please, add a file with one specie name / line and the path /path/Project/", call.=FALSE)
}

##### Libraries #############################
library(brms)
#library(rstanarm)
library(rstan)
library(tidybayes)
library(dplyr)
library(reshape2)
library(data.table)
#library(ggplot2)
library(bayestestR)
library(bayesplot)

args <- commandArgs(trailingOnly = TRUE)

#Iteration number
iteration_number <- args[1]

print("Collecting all informations to build the models...")

#Load the tree data with EVEs/dEVEs counts and Branch lengths in MA
tree5<- read.table("/beegfs/data/bguinet/these/Statistics_analysis/GLM_baysian_analysis/Table_with_EVEs_candidates_and_phylogeny_informations_dsDNA.txt",sep=";",h=T)
colnames(tree5)<- c("node","parent", "branch.length","label", "index","Nb_EVEs_Events", "Nb_dEVEs_Events", "Nb_EVEs2", "Nb_dEVEs2",    "Node_age")

#Load all the mk states scenarios
#Read the states table
states_scenarios_table<- read.table("/beegfs/data/bguinet/these/Statistics_analysis/GLM_baysian_analysis/State_analysis/mk.states.txt",h=T)

#Change the format of the table
states_scenarios_table_i<-melt(as.data.table(states_scenarios_table[iteration_number,], keep.rownames = "Iteration"), id.vars = "Iteration")
states_scenarios_table_i$variable<-gsub("end_","",states_scenarios_table_i$variable)
states_scenarios_table_i$variable<- as.character(states_scenarios_table_i$variable)
states_scenarios_table_i<-states_scenarios_table_i[!is.na(states_scenarios_table_i$variable),]
states_scenarios_table_i<-tail(states_scenarios_table_i, -2)

#Load the MA scenario tree to get posterior probabilities of node states
Bigtable<- read.table("/beegfs/data/bguinet/these/Statistics_analysis/GLM_baysian_analysis/Table_with_MA_states.txt",sep=";",h=T)
Bigtable$index <- as.character(Bigtable$index)

#Merge informations
tab<-merge(x=states_scenarios_table_i,y=Bigtable,by.x=c("variable"),by.y=c("index"),ALL=T)
tab<-merge(x=tab,y=tree5,by.x=c("parent","node"),by.y=c("parent", "node"),ALL=T)

tab$value<-as.character(tab$value)
tab<-tab[!abs(tab$Node_age)>160,]
tab<-tab[!tab$anc_state_1_pp<0.9,]


print("All informations loaded succesfully")
print("\n")

print(paste0("Running EVEs_Events_model2 on iteration : ", iteration_number))


#Create dataframe for ttest

columns= c("Name","Pvalue")
df = data.frame(matrix(nrow = 0,
                           ncol = length(columns)))
colnames(df) = columns


#####
#Posterior predictiv check
####
#
# The proportion of zeros
#yrepnzb <- posterior_predict(mymodel)
#(prop_zero_test4 <- ppc_stat(y=tab$Nb_EVEs_Events, yrepnzb, stat=function(y) mean(y==0)))
# The max count
#(max_test_znb <- ppc_stat(y=tab$Nb_EVEs_Events, yrepnzb, stat="max"))


#Run the bayesian GLM model
#EVEs Event model
EVEs_Events_model = brm(Nb_EVEs_Events ~  value + branch.length.x  , data = tab,
                     family = zero_inflated_negbinomial(),
                     iter  = 10000, chains = 4,
                     seed  = 1234,
                     cores = 1,thin=5,backend = "rstan")


#Get summary model statistics
posterior_results_EVEs_Events<- describe_posterior(EVEs_Events_model,ci = 0.89)
posterior_results_EVEs_Events$Parameter<-paste0("Iteration_",iteration_number,":",posterior_results_EVEs_Events$Parameter)
#posterior_results_EVEs_Events$Bayesian_pvalue <-  2*(1- as.integer(gsub("%","",posterior_results_EVEs_Events$pd)) / 100)

#Save the posterior parameters
posteriors_EVEs_Events <- insight::get_parameters(EVEs_Events_model)
posteriors_EVEs_Events$Iteration <- paste0("Iteration_",iteration_number)
write.table(as.data.frame(posteriors_EVEs_Events),paste0("/beegfs/data/bguinet/these/Statistics_analysis/GLM_baysian_analysis/output_dsDNA/Posteriors_EVEs_Events_iteration_",iteration_number,".txt"),sep=";")

#pairwise test comparaison
posterior_EVEs_Events_model<-EVEs_Events_model %>%
  spread_draws(b_value2,b_value3)

df[nrow(df) + 1,] = c("Pvalue_EVEs_event",t.test(posterior_EVEs_Events_model$b_value2,posterior_EVEs_Events_model$b_value3, paired=TRUE,alternative = "greater")$p.value)

#Write results
posterior_results_EVEs_Events$Bayesian_pvalue <- 2*(1-posterior_results_EVEs_Events$pd*100/100)
print("Results:")
print(posterior_results_EVEs_Events)
write.table(as.data.frame(posterior_results_EVEs_Events),paste0("/beegfs/data/bguinet/these/Statistics_analysis/GLM_baysian_analysis/output_dsDNA/Table_EVEs_Events_iteration_",iteration_number,".txt"),sep=";")

#------------------------------


print(paste0("Running dEVEs_Events_model on iteration : ", iteration_number))
#dEVEs Event model
dEVEs_Events_model = brm(Nb_dEVEs_Events ~  value + branch.length.x  , data = tab,
                     family = zero_inflated_negbinomial(),
                     iter  = 10000, chains = 4,
                     seed  = 1234,
                     cores = 1,thin=5,backend = "rstan")


#Get summary model statistics
posterior_results_dEVEs_Events<- describe_posterior(dEVEs_Events_model,ci = 0.89)
posterior_results_dEVEs_Events$Parameter<-paste0("Iteration_",iteration_number,":",posterior_results_dEVEs_Events$Parameter)

#Save the posterior parameters
posteriors_dEVEs_Events <- insight::get_parameters(dEVEs_Events_model)
posteriors_dEVEs_Events$Iteration <- paste0("Iteration_",iteration_number)
write.table(as.data.frame(posteriors_dEVEs_Events),paste0("/beegfs/data/bguinet/these/Statistics_analysis/GLM_baysian_analysis/output_dsDNA/Posteriors_dEVEs_Events_iteration_",iteration_number,".txt"),sep=";")

#pairwise test comparaison
posterior_dEVEs_Events_model<-dEVEs_Events_model %>%
  spread_draws(b_value2,b_value3)

df[nrow(df) + 1,] = c("Pvalue_dEVEs_event",t.test(posterior_dEVEs_Events_model$b_value2,posterior_dEVEs_Events_model$b_value3, paired=TRUE,alternative = "greater")$p.value)

#Write results
posterior_results_dEVEs_Events$Bayesian_pvalue <- 2*(1-posterior_results_dEVEs_Events$pd*100/100)
print("Results:")
print(posterior_results_dEVEs_Events)

write.table(as.data.frame(posterior_results_dEVEs_Events),paste0("/beegfs/data/bguinet/these/Statistics_analysis/GLM_baysian_analysis/output_dsDNA/Table_dEVEs_Events_iteration_",iteration_number,".txt"),sep=";")

#-------------------------------
print(paste0("Running EVEs_Counts_model on iteration : ", iteration_number ))

#EVEs count model
EVEs_Counts_model = brm(Nb_EVEs2 ~  value + branch.length.x  , data = tab,
                     family = zero_inflated_negbinomial(),
                     iter  = 10000, chains = 4,
                     seed  = 1234,
                     cores = 1,thin=5,backend = "rstan")


#Get summary model statistics
posterior_results_EVEs_Counts<- describe_posterior(EVEs_Counts_model,ci = 0.89)
posterior_results_EVEs_Counts$Parameter<-paste0("Iteration_",iteration_number,":",posterior_results_EVEs_Counts$Parameter)

#Save the posterior parameters
posteriors_EVEs_Counts <- insight::get_parameters(EVEs_Counts_model)
posteriors_EVEs_Counts$Iteration <- paste0("Iteration_",iteration_number)
write.table(as.data.frame(posteriors_EVEs_Counts),paste0("/beegfs/data/bguinet/these/Statistics_analysis/GLM_baysian_analysis/output_dsDNA/Posteriors_EVEs_Counts_iteration_",iteration_number,".txt"),sep=";")

#pairwise test comparaison
posterior_EVEs_Counts_model<-EVEs_Counts_model %>%
  spread_draws(b_value2,b_value3)

df[nrow(df) + 1,] = c("Pvalue_EVEs_Counts",t.test(posterior_EVEs_Counts_model$b_value2,posterior_EVEs_Counts_model$b_value3, paired=TRUE,alternative = "greater")$p.value)


#Write results
posterior_results_EVEs_Counts$Bayesian_pvalue <- 2*(1-posterior_results_EVEs_Counts$pd*100/100)

print("Results:")
print(posterior_results_EVEs_Counts)

write.table(as.data.frame(posterior_results_EVEs_Counts),paste0("/beegfs/data/bguinet/these/Statistics_analysis/GLM_baysian_analysis/output_dsDNA/Table_EVEs_Counts_iteration_",iteration_number,".txt"),sep=";")

#---------------------------------
print(paste0("Running dEVEs_Counts_model on iteration : ",iteration_number))
#dEVEs count model
dEVEs_Counts_model = brm(Nb_dEVEs2 ~  value + branch.length.x  , data = tab,
                     family = zero_inflated_negbinomial(),
                     iter  = 10000, chains = 4,
                     seed  = 1234,
                     cores = 1,thin=5,backend = "rstan")


#Get summary model statistics
posterior_results_dEVEs_Counts<- describe_posterior(dEVEs_Counts_model,ci = 0.89)
posterior_results_dEVEs_Counts$Parameter<-paste0("Iteration_",iteration_number,":",posterior_results_dEVEs_Counts$Parameter)

#Save the posterior parameters
posteriors_dEVEs_Counts <- insight::get_parameters(dEVEs_Counts_model)
posteriors_dEVEs_Counts$Iteration <- paste0("Iteration_",iteration_number)
write.table(as.data.frame(posteriors_dEVEs_Counts),paste0("/beegfs/data/bguinet/these/Statistics_analysis/GLM_baysian_analysis/output_dsDNA/Posteriors_dEVEs_Counts_iteration_",iteration_number,".txt"),sep=";")

#pairwise test comparaison
posterior_dEVEs_Counts_model<-dEVEs_Counts_model %>%
  spread_draws(b_value2,b_value3)

df[nrow(df) + 1,] = c("Pvalue_dEVEs_Counts",t.test(posterior_dEVEs_Counts_model$b_value2,posterior_dEVEs_Counts_model$b_value3, paired=TRUE,alternative = "greater")$p.value)


#Write results
posterior_results_dEVEs_Counts$Bayesian_pvalue <- 2*(1-posterior_results_dEVEs_Counts$pd*100/100)

print("Results:")
print(posterior_results_dEVEs_Counts)
write.table(as.data.frame(posterior_results_dEVEs_Counts),paste0("/beegfs/data/bguinet/these/Statistics_analysis/GLM_baysian_analysis/output_dsDNA/Table_dEVEs_Counts_iteration_",iteration_number,".txt"),sep=";")

#________________
write.table(as.data.frame(df),paste0("/beegfs/data/bguinet/these/Statistics_analysis/GLM_baysian_analysis/output_dsDNA/Table_student_test_iteration_",iteration_number,".txt"),sep=";")
