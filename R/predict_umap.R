#load("01.DKFZ_2801.pca.rda")
#umap = umap(betas.pca$x[,1:94], n_neighbor = 15, n_epochs = 500)

predict_umap = function(pca){
    library(umap)
    #load("/data/wuz6/project/02.DKFZ.classifier/07.queryNeighbors/R/data/data2/DKFZ.2801.pc94umap.rda")
    umap.p = predict(umap, data = pca);
    #
    #all.umap = data.frame(rbind(umap$layout, umap.p));
    #all.umap$Y = c(umap.y, rownames(pca));
    #TSNE center of refence tsne
    ref.umap = data.frame(umap$layout);
    umap.center = aggregate(ref.umap, list(umap.y), median);
    colnames(umap.center)[1] = "class"
    return(list(ref.umap, umap.p, umap.center))
}
