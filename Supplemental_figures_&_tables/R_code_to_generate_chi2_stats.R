# Statistics to look for difference between number of Events among genomic structures
M <- as.table(rbind(c(155, 75, 401,1844), c(131,26,59,152))) # création d'une table 2 lignes/3 colonnes
dimnames(M) <- list(gender=c("Nb species","EVEs event"), Genome_structure=c("dsDNA","ssDNA", "dsRNA","ssRNA"))#
print(M)
(test <- chisq.test(M)) # affichage des résultats du test
#The null hypothesis that the proportions are the same in the two databases is very strongly rejected with P-value near 0.
chisq.test(M)$resi

# Statistics to look for difference between number of Events among genomic structures without controls
M_without_controls<- as.table(rbind(c(155, 75, 401,1844), c(127,26,59,152))) # création d'une table 2 lignes/3 colonnes
dimnames(M_without_controls) <- list(gender=c("Nb species","EVEs event"), Genome_structure=c("dsDNA","ssDNA", "dsRNA","ssRNA"))#
print(M_without_controls)
(test <- chisq.test(M_without_controls)) # affichage des résultats du test
#The null hypothesis that the proportions are the same in the two databases is very strongly rejected with P-value near 0.
chisq.test(M_without_controls)$resi

#The sum of the squares of the Pearson Residuals is the chi-squared statistic 1265.
#Residuals with the largest absolute value point the way to the cells in which the observed and expected counts differed most.


# Statistics to look for difference between number of dEvents among genomic structures
M2 <- as.table(rbind(c(155, 75, 401,1844), c(47,10,12,38))) # création d'une table 2 lignes/3 colonnes
dimnames(M2) <- list(gender=c("Nb species","dEVEs event"), Genome_structure=c("dsDNA","ssDNA", "dsRNA","ssRNA"))#

print(M2)
(test <- chisq.test(M2))
chisq.test(M2)$resi

# Statistics to look for difference between number of dEvents among genomic structures without controls
M2_without_controls <- as.table(rbind(c(155, 75, 401,1844), c(44,10,12,38))) # création d'une table 2 lignes/3 colonnes
dimnames(M2_without_controls) <- list(gender=c("Nb species","dEVEs event"), Genome_structure=c("dsDNA","ssDNA", "dsRNA","ssRNA"))#

print(M2_without_controls)
(test <- chisq.test(M2_without_controls))
chisq.test(M2_without_controls)$resi


# Statistic counting difference between ssRNA(-) and ssRNA(+) Events
M3 <- as.table(rbind(c(1241,603), c(22,117))) # création d'une table 2 lignes/3 colonnes
dimnames(M3) <- list(gender=c("Nb NCBI fam","Nb EVE Viral fam"), Genome_structure=c("+","-"))#
(test <- chisq.test(M3))
chisq.test(M3)$resi
