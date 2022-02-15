#' getBetas32k
#' @export
getBetas32k = function(RGset, targets, keep_probes = names(pca.center)){
    #batch_coor_data = file.path(package.path <- find.package("cnsTumorMNN"),
     #                           "data/ba.coef.RData")
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
    #load(batch_coor_data)

    if( sum(! keep_probes %in% names(methy.coef[[1]]) ) > 0 ){
        stop("Error: probes not found in batch correction data ", batch_coor_data,"\n",
             keep_probes[! keep_probes%in% names(methy.coef[[1]])],"\n")
    }
    message("Batch correction...")
    methy.coef[[1]] = methy.coef[[1]][keep_probes]
    methy.coef[[2]] = methy.coef[[2]][keep_probes]
    unmethy.coef[[1]] = unmethy.coef[[1]][keep_probes]
    unmethy.coef[[2]] = unmethy.coef[[1]][keep_probes]

    methy <- minfi::getMeth(Mset)
    unmethy <- minfi::getUnmeth(Mset)

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


#' getBetas32k
#' @export
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




