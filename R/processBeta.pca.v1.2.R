#load("01.DKFZ_2801.pca.rda")
#load("MNPbetas.32kvar.rda")
#anno$`methylation class:ch1`[which( anno$`methylation class:ch1` == "PIN T,  PB A")] = "PIN T, PB A"
#anno$`methylation class:ch1`[which( anno$`methylation class:ch1` == "PIN T,  PB B")] = "PIN T, PB B"
#xx = sapply(rownames(betas.pca$x), function(x) unlist(strsplit(x, "_"))[1])
#pca.y = anno$`methylation class:ch1`[match(xx, rownames(anno))]
#N_non_trivial_pc = 94
#pca.x = betas.pca$x[,1:N_non_trivial_pc]
#pca.center = betas.pca$center
#pca.rotation =  betas.pca$rotation[,1:N_non_trivial_pc]
#save(pca.x, pca.center, pca.rotation, pca.y, file = "DKFZ.2801.pc94.rda")

load("/data/wuz6/project/02.DKFZ.classifier/07.queryNeighbors/R/data/DKFZ.train2801.filter.sd32k.pc300Data.v1.1.rda")
#"pca.center"   "pca.rotation" "x.ref"        "y.ref"
#@ targets: table
getBetas32k = function(RGset, targets, keep_probes = names(pca.center)){
    source("/data/wuz6/software/mnp_training/R/MNPprocessIDAT_functions.R")
    #keep_probes_file = "/data/wuz6/project/02.DKFZ.classifier/07.queryNeighbors/R/data/DKFZ.train2801.filter.sd32k.prbles.gz" 
    #if(! file.exists(keep_probes_file)){ stop("cannot fine the list of 32k probes: ", keep_probes_file,"\n"); }
    batch_coor_data = "/data/wuz6/project/02.DKFZ.classifier/07.queryNeighbors/R/data/data2/ba.coef.RData"
    #message("rgset")
    #RGset <- minfi::read.metharray.exp(targets = targets,force=TRUE, verbose = T)
    message("mset...")
    Mset <- MNPpreprocessIllumina(RGset)
    message("mest done...")
    #keep_probes = read.table(gzfile(keep_probes_file), stringsAsFactors=F, head=F)[,1]
    if(sum(! keep_probes %in% rownames(Mset)) > 0) stop("Error:", sum(!keep_probes %in% rownames(Mset)), " probes not found:\n", keep_probes[!keep_probes %in% rownames(Mset)]);
    Mset = Mset[keep_probes,]
    Mset
    rm(RGset)
    message("keep probes done in Mset")
#
    load(batch_coor_data)

    if( sum(! keep_probes %in% names(methy.coef[[1]]) ) > 0 ){
        stop("Error: probes not found in batch correction data ", batch_coor_data,"\n",
             keep_probes[! keep_probes%in% names(methy.coef[[1]])],"\n")
    }
    message("Batch correction...")
    methy.coef[[1]] = methy.coef[[1]][keep_probes]
    methy.coef[[2]] = methy.coef[[2]][keep_probes]
    unmethy.coef[[1]] = unmethy.coef[[1]][keep_probes]
    unmethy.coef[[2]] = unmethy.coef[[1]][keep_probes]

    methy <- getMeth(Mset)
    unmethy <- getUnmeth(Mset)

    rm(Mset)

    batch <- ifelse(as.character(targets$Material_Type) == "FFPE", "FFPE", "Frozen")
    # perform batch adjustment
    methy.b <- log2(methy +1) + matrix(unlist(methy.coef[match(batch, names(methy.coef))]),ncol=length(batch))
    unmethy.b <- log2(unmethy +1) + matrix(unlist(unmethy.coef[match(batch, names(unmethy.coef))]),ncol=length(batch))
    methy.b[methy.b < 0] <- 0
    unmethy.b[unmethy.b < 0] <- 0
    methy.ba <- 2^methy.b
    unmethy.ba <- 2^unmethy.b
    # illumina-like beta values
    betas<- methy.ba / (methy.ba +unmethy.ba +100)
    betas = t(betas)
    betas
}



#@ betas:
pca_transform = function(betas, n.pc = 300){
    #load("/data/wuz6/project/02.DKFZ.classifier/07.queryNeighbors/R/data/DKFZ.train2801.filter.sd32k.pc300Data.rda")
    #load("/data/wuz6/project/02.DKFZ.classifier/07.queryNeighbors/R/data/DKFZ.train2801.filter.sd32k.pc300Data.v1.1.rda")
    if(! identical(colnames(betas), names(pca.center))){
        stop("Error in pca_trainsform: input betas has different probes from pca transfrom reference.")
    }
    df = t(apply(betas, 1, function(x) x - pca.center))
    rm(betas)

    x.test = df %*% pca.rotation[,1:n.pc]
    #return(list(x.ref, pca.y, x.test))
    return(list(x.ref, y.ref, x.test))
}
#




