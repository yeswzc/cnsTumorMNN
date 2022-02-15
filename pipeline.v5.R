##Code use DKFZ classifier mnp.v11b6 to predict samples profiled by EPIC/450k array
##Removed code related to DKFZ classifier
##Requirement:
##R/3.6
##
args = commandArgs(trailingOnly = F)

if(length(args) != 5 + 2){
      stop("Usage: [idat.folder] [output.prefix]\n")
}

idat_path = as.character(args[5+1])
output = as.character(args[5+2])
code.dir <- gsub("--file=", "", args[4])
code.dir <- paste0(dirname(code.dir), "/")

#####################################################################3
#mnpversion = "v11b6"
library(caret) #knn3
library(glmnet)
library(cnsTumorMNN)
package.path <- find.package("cnsTumorMNN")
library(dplyr)
library(plotly)
#library(RANN)
library(knitr)
library(kableExtra)
#suppressMessages(library(mnp.v11b6))
library(IlluminaHumanMethylationEPICmanifest) #CNV
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) #CNV
library(RFpurify)
library(LUMP)

rmd_file = file.path(code.dir, "Generate_HTMLreport.v5.Rmd")
paired.color = RColorBrewer::brewer.pal(12,"Paired")
#####################################################################
#####################################################################

csv_file_name = "Sample_Sheet.csv"
targets <- read.csv(file.path(idat_path, csv_file_name), stringsAsFactors = FALSE, skip = 7)
targets$idat <- paste0(targets$Sentrix_ID,"_",targets$Sentrix_Position);
targets$Basename <- file.path(idat_path, targets$idat);

print(head(targets, 2))
RGset = read.metharray.exp(targets = targets, verbose=T, force=TRUE)
print(colnames(RGset))

###Run DKFZ classifier
if(! file.exists(paste0(output, ".cns_tmp.Rda"))){
    #cat("Run DKFZ classifier...\n");
    cat("Run Tumor purity...\n");
    CNSclassifier = lapply(1:nrow(targets), function(x) classifier(RGset, targets$idat[x], targets$Material_Type[x]));
    CNSclassifier = data.frame(do.call(rbind, CNSclassifier));
    #colnames(CNSclassifier) = c("ID","predFFPE","CNS.MCF","CNS.Subclass", 
    #                            "RFpurify.ABSOLUTE", "RFpurify.ESTIMATE", "LUMP");
    colnames(CNSclassifier) = c("ID", 
                                "RFpurify.ABSOLUTE", "RFpurify.ESTIMATE", "LUMP");
    CNSclassifier = cbind(CNSclassifier, targets[,c(1:(ncol(targets)-1))])
    save(CNSclassifier, targets, file = paste0(output,".cns_tmp.Rda"))
}else{
    load(paste0(output,".cns_tmp.Rda"))
}

#get beta values
if(! file.exists(paste0(output, ".betas32k.rda")) ){
    betas = getBetas32k(RGset, targets);
    save(betas, file = paste0(output, ".betas32k.rda"))
}else{
    load(paste0(output, ".betas32k.rda"))
}
#get pcas
cat("Get PCA transformation...\n")
pca.res = pca_transform(betas)
x.ref = pca.res[[1]]
y.ref = pca.res[[2]]
x.test = pca.res[[3]]
rm(pca.res) 
if(sum(!y.ref %in% names(cancer.colour))> 0 ){ stop("Some cancer type not find matched color:", y.ref[!y.ref %in% names(cancer.colour)]); }
if(sum( !names(ref_distance.subclass) %in% y.ref) > 0){ stop("Some cancer distance density not found: ",names(ref_distance.subclass)[! names(ref_distance.subclass) %in% y.ref]); }

####TSNE
if(! file.exists(paste0(output, ".tsne.rda") ) ){
    cat("TSNE...\n");
    tsne_s = run_tsne(x.ref, x.test, y.ref, n.pc = 300)
    ref.tsne = tsne_s[[1]] 
    test.tsne = tsne_s[[2]]
    ref.tsne.center = tsne_s[[3]]
    rm(tsne_s)
    save(ref.tsne, test.tsne, ref.tsne.center, file = paste0(output, ".tsne.rda"))
}else{
    load(paste0(output, ".tsne.rda"))
}

####UMAP
if(! file.exists(paste0(output, ".umap.rda") ) ){
    cat("UMAP...\n");
    umap_s = predict_umap(x.test[,1:94])
    ref.umap = umap_s[[1]]
    test.umap = umap_s[[2]]
    ref.umap.center = umap_s[[3]]
    rm(umap_s)
    y.test = targets$Sample_Name
    save(ref.umap, test.umap, ref.umap.center, y.ref, y.test, file = paste0(output, ".umap.rda"))
}else{
    load(paste0(output, ".umap.rda"))
}

###
#KNN predict on new test
if(! file.exists(paste0(output, ".KNN_prediciton.rda"))){
    cat("KNN classifer...\n");
    knn.p.raw = predict(cns.knn.v1.1, x.test, type = "prob")
    #calibration
    knn.p.calibrated = predict(knn.score.calibration$glmnet.fit, newx = knn.p.raw,
                               type = "response", s = knn.score.calibration$lambda.1se)[,,1]
    knn.p.calibrated.MCF = get_family_score(knn.p.calibrated) #not very good
    #save(knn.p.raw, knn.p.calibrated, knn.p.calibrated.MCF, file = paste0(output, ".KNN_prediciton.rda"))
}else{
    load(paste0(output, ".KNN_prediciton.rda"))
}

####prepare to generate reports
knn.p1  = get_class_score(knn.p.calibrated.MCF)
knn.p2 = get_class_score(knn.p.calibrated)

CNSclassifier$KNN.pred.MCF = names(knn.p1)
CNSclassifier$KNN.pred.MCF.score = as.vector(knn.p1)
CNSclassifier$KNN.pred.subclass = names(knn.p2)
CNSclassifier$KNN.pred.subclass.score = as.vector(knn.p2)

#CNSclassifier$MCF = sapply( as.character(CNSclassifier$CNS.MCF), function(x) unlist(strsplit(x,";"))[1])
#CNSclassifier$MCF1 = sapply(CNSclassifier$MCF, function(x) unlist(strsplit(x,":"))[1])
#CNSclassifier$MCF1.score = sapply(CNSclassifier$MCF, function(x) unlist(strsplit(x,":"))[2])
#CNSclassifier$Class = sapply(as.character(CNSclassifier$CNS.Subclass), function(x) unlist(strsplit(x,";"))[1])
#CNSclassifier$Class1 = sapply(CNSclassifier$Class, function(x) unlist(strsplit(x,":"))[1])
#CNSclassifier$Class1.score = sapply(CNSclassifier$Class, function(x) unlist(strsplit(x,":"))[2])

CNSclassifier$RFpurify.ABSOLUTE = as.numeric(as.character(CNSclassifier$RFpurify.ABSOLUTE))
CNSclassifier$RFpurify.ESTIMATE = as.numeric(as.character(CNSclassifier$RFpurify.ESTIMATE))
CNSclassifier$LUMP = as.numeric(as.character(CNSclassifier$LUMP))

pwd = getwd()
dir.create(paste0(output,".html"))
out_dir = file.path(pwd, paste0(output,".html"))
targets$idat = paste0(targets$Sentrix_ID,"_",targets$Sentrix_Position)


colnames(targets)[which(colnames(targets) == "Sample_name")] = "Sample_Name";
.. = lapply(1:nrow(targets), function(i){
    sampleIdat = as.character(targets$idat[i])
    sampleName = as.character(targets$Sample_Name[i])
    cat(sampleName,sampleIdat,"...")
    #
    #cns.p1 = CNSclassifier$MCF1[i]
    #cns.p1.score = CNSclassifier$MCF1.score[i]
    #cns.p2 = CNSclassifier$Class1[i]
    #cns.p2.score = CNSclassifier$Class1.score[i]
    cns.p1 = NA
    cns.p1.score = NA
    cns.p2 = NA
    cns.p2.score = NA
    #
    knn.pred.MCFclass = names(knn.p1[i])
    knn.pred.MCFclass.score = round(knn.p1[i],2)
    knn.pred.class = names(knn.p2[i])
    knn.pred.class.score = round(knn.p2[i],2)
    #class.table = data.frame(Classifier = c("CNSv11b6"), MCF = c(cns.p1), Score = c(cns.p1.score), 
    #                         Subclass = c(cns.p2), Score = c(cns.p2.score), check.names = F)
    #if(cns.p2 == cns.p1) class.table = data.frame(Classifier = c("CNSv11b6"), 
    #                                              Prediction = c(cns.p1), Score = c(cns.p1.score));
    class.table = data.frame(Classifier = c("KNN"), MCF = c(knn.pred.MCFclass), Score = c(knn.pred.MCFclass.score), 
                             Subclass = c(knn.pred.class), Score = c(knn.pred.class.score), check.names = F)
    row.names(class.table) = NULL

    histologicDx = paste0("**Diagnosis**: ", targets$Diagnosis[i],". **Age**: ", targets$Age[i], 
                          ", **Gender**:", targets$Gender[i] ,", **Tumor size**: ", 
                          targets$Tumor_site[i],". <br>**Notes**: ", targets$Notes[i],".");
    #prepare data for tsne plot
    plot.col = c(cancer.colour, "#000000");
    names(plot.col)[length(plot.col)] = sampleName;
    x.query = x.test[i, , drop=F]
    x.query.tsne = as.numeric(as.character(test.tsne[i,]))
    tsne = data.frame(rbind(ref.tsne, x.query.tsne));
    tsne$class = c(y.ref, sampleName);
    tsne$size = c(rep(1, length(y.ref)), 4);
    x.query.umap = as.numeric(as.character(test.umap[i,]));
    umap.layout = data.frame(rbind(ref.umap, x.query.umap));
    umap.layout$class = c(y.ref, sampleName);
    umap.layout$size = c(rep(1, length(y.ref)), 4);
    #
    cat("knn.search...")
    res1 = knn.search(x.ref = x.ref, x.query = x.query, y.ref = y.ref, k = 15)
    #
    sex = targets$Gender[i];
    sex = ifelse(sex == "Male" | sex == "M", "Male",ifelse(sex == "Female" | sex == "F", "Female", "unknow"))
    idat_file = targets$Basename[i];
    p.rf.ABSOLUTE = round(CNSclassifier$RFpurify.ABSOLUTE[i], 2);
    p.rf.ESTIMATE = round(CNSclassifier$RFpurify.ESTIMATE[i], 2);
    p.lump = round(CNSclassifier$LUMP[i], 2)
    p.table = data.frame("RFpurify(ABSOLUTE)" = p.rf.ABSOLUTE ,
                         "RFpurify(ESTIMATE)" = p.rf.ESTIMATE, "LUMP" = p.lump, check.names=F);
    cat("read array\n");
    RGset <- minfi::read.metharray(idat_file, verbose=FALSE);  Mset <- minfi::preprocessRaw(RGset); rm(RGset);
    cat("generating report\n");
    rmarkdown::render(rmd_file, output_file = paste0(sampleName,"_MNNv5_",sampleIdat,".html"), 
                      output_dir = out_dir, intermediates_dir = out_dir, knit_root_dir = out_dir,
                      quiet = T);

    nearest.neighbor = NA
    if(sum(res1[[2]]$`N PASS`) > 0 ){
        nearest.neighbor = paste0(as.character(paste0(res1[[2]]$Cancer,"/",
                                  res1[[2]]$`N PASS`)[which(res1[[2]]$`N PASS` > 0)]), collapse = ";");
    }
    nearest.neighbor;
    #break;
})
CNSclassifier$NearestNeighbor = unlist(..)
write.csv(CNSclassifier, file = paste0(output,".final.combined.csv"), row.names = F)
#system(paste0("rm ",output,".cns_tmp.Rda"))


