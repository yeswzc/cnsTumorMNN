#' 
#' @export
run_tsne = function(x.ref, x.test, y.ref, n.pc = 94){
    set.seed(123)
    all_tsne = data.frame(Rtsne::Rtsne(rbind(x.ref, x.test)[,1:n.pc], pca=F,max_iter=2500,theta=0,dim=2,verbose=T)$Y)
    #save(all_tsne, file = paste0(output, ".all_tsne.rda"))
    #load(paste0(output, ".all_tsne.rda"))
    ref.tsne = all_tsne[1:nrow(x.ref),]
    test.tsne = all_tsne[-c(1:nrow(x.ref)),]

    #TSNE center of refence tsne
    tsne.center = aggregate(ref.tsne, list(y.ref), median)
    colnames(tsne.center)[1] = "class"
    return(list(ref.tsne, test.tsne, tsne.center))
}
