#This is code from the IMA package that is slightly modified.
#Installing the dependencies for the IMA package can lead to issues, so I copied and 
#modified the code to create a methy450batch object. In addition, I have copied the
#function defn for indexregionfunc(). 

##We create a methy450batch object with 14 slots
##slots: 11 slots for the annnotated regions, remaining slots: mvalues, annotation, pData
#By creating this object we can aggregate statistics for features in the same genomic region
#or chromosomal region.

#Reference: http://ima.r-forge.r-project.org/IMA-manual.pdf

mval.matrix = as.matrix(mvalues)
annotation = as.matrix(hm450)
groupinfo = pData

cat(".......\nSplit the annotation file to 11 
    annotated region categories\n.......\n\n", fill = TRUE)
    annot = annotation
    name = "UCSC_RefGene_Name"
#get all the CpG site names ex. "cg16524862"
    cpGsite = as.character(annot[, 1])
#create a list of lists: each CpG site has a list of gene regions that it is associated with
    genelist = strsplit(as.character(annot[, name]), ";")
#For CpG sites that are not associated with any gene region, put NA as the element in its list.
    genelist[which(genelist == "character(0)")] = "NA"
    name = "UCSC_RefGene_Group"
#create a list of lists: each CpG site has a list of gene features (ex. Body, TSS200, etc.) that it is associated with
    refgene = strsplit(as.character(annot[, name]), ";")
    refgene[which(refgene == "character(0)")] = "NA"
    listlength = lapply(refgene, length)
    listlength[listlength == 0] = 1
##We rep col0, 1, 4, 5 according to how many sites each CpG probe is linked to 
##col0 are the indices, col1 are the Cpg site names
##col2 are the names of the associated genes
##col3 are the names of the gene feature (ex. Body, TSS200, etc.) in which the CpG probe is located 
##col4 are the CpG site's region (i.e. Island, shore, etc.)
##col5 are the chromosomal location of the CpG islands
    col0 = rep(1:nrow(annot), listlength)
    col1 = rep(cpGsite, listlength)
    col2 = unlist(genelist)
    col3 = unlist(refgene)
    col4 = rep(as.character(annot[, "Relation_to_UCSC_CpG_Island"]), 
        listlength)
    col5 = rep(as.character(annot[, "UCSC_CpG_Islands_Name"]), 
        listlength)

    splitToRegionlist = function(grepname = c("TSS1500", "TSS200", 
        "5'UTR", "1stExon", "Gene Body", "3'UTR")) {
        index = col3 == grepname
        col1sub = col1[index]
        col2sub = col2[index]
        temp = split(col1sub, col2sub)
        returnSID = lapply(temp, unique)
        col0sub = col0[index]
        temp = split(col0sub, col2sub)
        returnPID = lapply(temp, unique)
        return(Ind = list(SID = returnSID, PID = returnPID))
    }
##each of the following contains a list such that each element corresponds to a gene.
#each element is 2 lists SID and PID, SID contains the list of indices, and PID contains the corresponding
#Cpgsite ex. 'cg16524862', 'cg16524863'
    TSS1500Ind = splitToRegionlist(grepname = "TSS1500")
    TSS200Ind = splitToRegionlist(grepname = "TSS200")
    UTR5Ind = splitToRegionlist(grepname = "5'UTR")
    EXON1Ind = splitToRegionlist(grepname = "1stExon")
    GENEBODYInd = splitToRegionlist(grepname = "Body")
    UTR3Ind = splitToRegionlist(grepname = "3'UTR")

    cat("TSS1500 region contains:", length(TSS1500Ind$SID), 
        "UCSC REFGENE region \nTSS200 region contains:", 
        length(TSS200Ind$SID), "UCSC REFGENE region\n5'UTR region contains:", 
        length(UTR5Ind$SID), "UCSC REFGENE region\n1st Exon region contains:", 
        length(EXON1Ind$SID), "UCSC REFGENE region\nGene body region contains:", 
        length(GENEBODYInd$SID), "UCSC REFGENE region\n3'UTR region contains:", 
        length(UTR3Ind$SID), "UCSC REFGENE region\n", fill = TRUE)

##Second list of summary stats indexed by chromosomal regions for slots
    splitToRegionlist2 = function(grepname = c("Island", "N_Shore", 
        "S_Shore", "N_Shelf", "S_Shelf")) {
        index = col4 == grepname
        col1sub = col1[index]
        col5sub = col5[index]
        temp = split(col1sub, col5sub)
        returnSID = lapply(temp, unique)
        col0sub = col0[index]
        temp = split(col0sub, col5sub)
        returnPID = lapply(temp, unique)
        return(Ind = list(SID = returnSID, PID = returnPID))
    }
#The following are also list of lists, each element refers to a chromosomal location such as `chr1:533219-534114`.
#Each element of the list is 2 lists (SID and PID) of the CpG sites at that location.
    ISLANDInd = splitToRegionlist2(grepname = "Island")
    NSHOREInd = splitToRegionlist2(grepname = "N_Shore")
    SSHOREInd = splitToRegionlist2(grepname = "S_Shore")
    NSHELFInd = splitToRegionlist2(grepname = "N_Shelf")
    SSHELFInd = splitToRegionlist2(grepname = "S_Shelf")

    cat("Island region contains:", length(ISLANDInd$SID), 
        "UCSC CPG ISLAND region\nN_Shore region contains", 
        length(NSHOREInd$SID), "UCSC CPG ISLAND region\nS_Shore region contains", 
        length(SSHOREInd$SID), "UCSC CPG ISLAND region\nN_Shelf region contains", 
        length(NSHELFInd$SID), "UCSC CPG ISLAND region\nS_Shelf region contains", 
        length(SSHELFInd$SID), "UCSC CPG ISLAND region\n", fill = TRUE)
    setClass("methy450batch", representation(bmatrix = "matrix", 
        annot = "matrix", groupinfo = "data.frame", 
        TSS1500Ind = "list", TSS200Ind = "list", UTR5Ind = "list", 
        EXON1Ind = "list", GENEBODYInd = "list", UTR3Ind = "list", 
        ISLANDInd = "list", NSHOREInd = "list", SSHOREInd = "list", 
        NSHELFInd = "list", SSHELFInd = "list"), where = topenv(parent.frame()))
    x.methy450 = new("methy450batch", bmatrix = as.matrix(mval.matrix), 
        annot = as.matrix(annotation), groupinfo = groupinfo, 
        TSS1500Ind = TSS1500Ind, TSS200Ind = TSS200Ind, 
        UTR5Ind = UTR5Ind, EXON1Ind = EXON1Ind, GENEBODYInd = GENEBODYInd, 
        UTR3Ind = UTR3Ind, ISLANDInd = ISLANDInd, NSHOREInd = NSHOREInd, 
        SSHOREInd = SSHOREInd, NSHELFInd = NSHELFInd, SSHELFInd = SSHELFInd)
    cat("\nA methy450batch class is created and the slotNames are: \n", 
        slotNames(x.methy450), "\n", fill = TRUE)

#function copied from IMA package
#"For each specific region of a gene, IMA will collect the loci within it 
#and derive an index of overall region methylation value. 
#Three different index metrics implemented in IMA:
#mean, median, and Tukey’s Biweight robust average. By default, the mean beta values will be used
#as the region’s methylation index for further analysis"
indexregionfunc <- function (indexlist, beta, indexmethod = c("mean", "median", 
    "tbrm")) 
{
    nr = length(indexlist$PID)
    temp2 = matrix(NA, nrow = nr, ncol = ncol(beta))
    rownames(temp2) = names(indexlist$SID)
    colnames(temp2) = colnames(beta)
    for (i in 1:nr) {
        temp = beta[indexlist$PID[[i]], ]
        if (length(indexlist$PID[[i]]) == 1) {
            temp2[i, ] = temp
        }
        else {
            if (indexmethod == "tbrm") {
                temp2[i, ] = apply(temp, 2, eval(indexmethod))
            }
            else {
                temp2[i, ] = apply(temp, 2, eval(indexmethod), 
                  na.rm = TRUE)
            }
        }
    }
    return(temp2)
}