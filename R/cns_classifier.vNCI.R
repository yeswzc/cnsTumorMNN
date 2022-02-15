##Code use DKFZ classifier mnp.v11b4 to predict samples profiled by EPIC/450k array
##Requirement:
##R/3.6
##
#####################################################################3
#@ 
#' @export
get_predicted = function(p, score = 0.3){
    col_names = colnames(p)
    keep = which(p >= score)
    if(length(keep) > 0){
        re_order = order(p[keep], decreasing = T)
        out = paste(col_names[keep][re_order], round(p[keep],3)[re_order], sep=":", collapse =";")
    }else{
        out = "NA"
    }
    out
}
#@ array ID
#@ path
#' @export
#classifier = function(arrayID, path){
classifier = function(RGset, arrayID, ffpe = NULL){
    RGset.subset = RGset[,arrayID]
    Mset = MNPpreprocessIllumina(RGset.subset);
    ### Purity
    p.rf.ABSOLUTE = predict_purity(Mset, method = "ABSOLUTE");
    p.rf.ESTIMATE = predict_purity(Mset, method = "ESTIMATE");
    betas = getBeta(Mset); p.lump = lump(betas);
    ###FFPE/Frozen batch correction
    #if(is.null(ffpe) | (!ffpe %in% c("FFPE", "Frozen")) ){
    #    message("Use FFPE/Frozen prediction.")
    #    p.ffpe = MNPgetFFPE(RGset);
    #    ffpe = p.ffpe
    #}
    #Mset = MNPbatchadjust(Mset,ffpe);
    ###Classifier Prediction
    #calRFscores.cg <- MNPpredict(Mset, type='prob', MCF=TRUE);
    #calRFscores <- MNPpredict(Mset, type='prob');
    #p1 = get_predicted(calRFscores.cg, 0.3);
    #p2 = get_predicted(calRFscores, 0.1);
    #p = c(arrayID, ffpe, p1, p2, p.rf.ABSOLUTE, p.rf.ESTIMATE, p.lump);
    p = c(arrayID, p.rf.ABSOLUTE, p.rf.ESTIMATE, p.lump);
    return(p)
}

