%This is R code wrapped in LATEX. It uses a R package called knitr to process the code and 
%can output the results in LATEX or a pdf.

\documentclass{article}

\usepackage{times}
\usepackage[cm]{fullpage}
\usepackage{amsmath,amsfonts,amsthm} 

\setlength{\parskip}{6pt}
\setlength{\parindent}{0pt}

\begin{document}

We modify and create the data structures that we will use for regression analysis.
We are given 4 data frames and 1 csv file.
1. Methylation values for the 480,320 sites from 428 patients.
2. Biological annotation from the UCSC database for each site.
	Information includes the region in the genome the site is found, gene name for which the
	site is associated with, etc.
3. The technical annotation for each probe on the array. Each probe is specific to a methylation 
site, hence there are ~480K probes.
4. 169 features for each patient (age, race, BMI at each age, etc.)
5. csv_file: The behavorial scores for each patient.

<<load_data, cache = FALSE, echo = TRUE, include = TRUE>>=

#load the 4 dataframes: mvalues, pvalues, hm450, annot
#hm450: is the biological annotation matrix (480320 * 37), with an entry for each site
#mvalues: contains the methylation values for each patient (480320 * 420)
#pData: 169 patient features for 428 patients (428 * 169)
#annot: technical annotation (483366*13)
#Each dataframe uses the same identifier for the methylation probes
load("/proj/design/il450k/nest13/R/MvalueAnalysisStructures.RData")

dim(mvalues)
#[1] 480320    428
dim(hm450)
#[1] 480320    428
dim(pData)
#[1] 428     37
dim(annot)
#[1] 483366     13

#rename annot to be more descriptive
annot.full <- annot
rm(annot)
#use the nestid as the patient identifier
table(rownames(pData) == pData$nestid)

#load the behavioral scores for each patient
epi<-read.csv("/home/guest/imp9/adhd/NESTSR_IVERSEN_8_22_15.csv", as.is=TRUE)
dim(epi)
#[1] 474  63
length(unique(epi$nestid)) == nrow(epi)
#[1] TRUE
rownames(epi) <-epi$nestid

#delete nestid for which the methylation age doesn't correspond with patient's age
#remove NESTID 198 due to possible age issues
epi<-epi[-grep('198', rownames(epi)),]

#identify patients in pData that are also in epi
table(keep<-rownames(pData) %in% rownames(epi))
#FALSE  TRUE 
#  253   175 
#only keep those patients in pData and mvalues that are also in epi
pData<-pData[keep,]
mvalues<-mvalues[,keep]
#rows will have the same order as pData
#epi will only have nestid which is in pData
epi <- epi[rownames(pData),]
colnames(epi)[colnames(epi) %in% colnames(pData)]
#merge patient info (race, age, etc) with behavioral scores
pData<-cbind(pData,epi[,colnames(epi) != "nest_id"])

#find the variance in the beta values for each site
RowVar <- function(x) {
    rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}
mvar <- RowVar(mvalues)

#only keep those probes where percentage of SNPs are less than half a percent
#remove probes with high percentage of SNPs from hm450 and mvalues
table(hm450$KeepHalfPct)
table(rownames(hm450)==rownames(mvalues))
mvalues<-mvalues[hm450$KeepHalfPct,]
hm450<-hm450[hm450$KeepHalfPct,]
table(rownames(hm450)==rownames(mvalues))

#colnames of mvalues are the barcode
#rownames of pData are the nestid (nest id and sample id are equal)
#for duplicate samples we used the higher quality mvalues
#which leads pData$barcode not corresponding to colnames of mvalues for the duplicate samples
#set rownames of pData to colnames of mvalues which is the barcode
rownames(pData) <- colnames(mvalues)

#create categorical variables for mother's ADHD, education, parity, and pre-pregnancy BMI
pData$education_3cat<- factor(pData$education_3cat, labels = 1:3, levels = 1:3)
pData$parity_3cat <- factor(pData$parity_3cat, labels = 0:2, levels = 0:2)
pData$prePregBMIthreeLev <- rep(NA, nrow(pData))
pData$prePregBMIthreeLev[(!is.na(pData$BMI_LMP_kgm2)) & (pData$BMI_LMP_kgm2 < 30)] <- 0
pData$prePregBMIthreeLev[(!is.na(pData$BMI_LMP_kgm2)) & (pData$BMI_LMP_kgm2 >= 30) & (pData$BMI_LMP_kgm2 < 35)] <- 1
pData$prePregBMIthreeLev[(!is.na(pData$BMI_LMP_kgm2)) & (pData$BMI_LMP_kgm2>=35)] <- 2
pData$prePregBMIthreeLev <-factor(pData$prePregBMIthreeLev, levels = 0:2, labels = c("Lt30", "30toLt35", "Ge35"))

save(pData,mvalues,hm450, file = RegrVars.RData)
@
\end{document}