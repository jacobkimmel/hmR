# heteromotility data class

heteromotility <- methods::setClass(
  'heteromotility',
  slots = c(raw.data = 'data.frame', # raw motility statistics data frame
            data = 'data.frame', # scaled scalar data in raw.data
            unscaled.data = 'data.frame', # scalar data in raw.data
            meta.data = 'data.frame', # meta data about specific rows
            feat.data = 'data.frame', # data about specific features
            pcs = 'matrix', # principal components
            tsne = 'matrix', # tsne embedding
            dist.matrix = 'matrix', # distance matrix
            celldataset = 'CellDataSet') # monocle CellDataSet
)

#' Build a `heteromotility` data object from a data.frame
#'
#' @param raw.data : data.frame of raw numeric data with 2 id columns
#' @param meta.data.list : list of character vectors or factors to include in meta.data
#'
#' @return hm : heteromotility data object.
#'
hmMakeObject <- function(raw.data, meta.data.list = list()){

  meta.data = raw.data[,1:2]
  for (i in 1:length(meta.data.list)){
    var.name = names(meta.data.list)[i]
    meta.data[var.name] = as.factor(meta.data.list[[i]])
  }
  hm = heteromotility(raw.data = raw.data,
                      unscaled.data = raw.data[,3:ncol(raw.data)],
                      meta.data = meta.data,
                      dist.matrix = matrix(),
                      feat.data = data.frame(row.names = colnames(raw.data[,3:ncol(raw.data)])))
  return(hm)
}

#' Scales and centers `raw.data` in a heteromotility object `hm`
#'
#' @param hm : `heteromtility` data object.
#'
#' @return hm : `heteromtility` data object.
#'   assigns @data.
#'
hmScaleData <- function(hm){
  scaled = data.frame(scale(hm@unscaled.data))
  hm@data = scaled
  return(hm)
}

#' Gets data from a `heteromotility` object and returns a data.frame
#'
#' @param hm : `heteromtility` data object.
#' @param scalar : string. property for numerical data to be retrieved. {'data', 'raw.data', 'pcs', etc}
#' @param cat : iterable of strings. categorical features to retrieve from `meta.data`
#'
#' @return df : data.frame. contains data in `scalar` and retrieved `meta.data`
#'
hmGetData <- function(hm, scalar='data', cat=c('clust')){

  df = data.frame(attr(hm, scalar))

  for (c in cat){
    df[c] = hm@meta.data[,c]
  }

  return(df)
}

#' Subsets data in heteromotility object `hm`
#'
#' @param hm : `heteromtility` data object.
#' @param cell_index : boolean or numerical index of rows in @data to keep.
#' @param cat : character. category to use when subsetting with `id`.
#' @param id : character. group in `cat` to keep.
#'
#' @return hm : `heteromtility` data object.
#'
hmSubsetData <- function(hm, cell_index=F, cat=F, id=F){
  if (cell_index != F){
    idx = cell_index
  } else if (cat != F & id != F) {
    idx = hm@meta.data[,cat] == id
  } else {
    print('hmSubsetData requires either a `cell_index` or both a category `cat` and `id` to subset.')
    return()
  }

  new = hm
  new@data = new@data[idx,]
  new@unscaled.data = new@unscaled.data[idx,]
  new@raw.data = new@raw.data[idx,]
  new@meta.data = new@meta.data[idx,]

  if (nrow(hm@dist.matrix) == nrow(hm@data)){
    new@dist.matrix = new@dist.matrix[idx,]
  }
  if (nrow(hm@pcs) == nrow(hm@data)){
    new@pcs = new@pcs[idx,]
  }
  if (nrow(hm@tsne) == nrow(hm@data)){
    new@tsne = new@tsne[idx,]
  }

  return(new)
}

#' Performs PCA on `data` in a heteromotility object `hm`
#' adds PCs to `pcs` and loadings to `pc_loadings`
#'
#' @param hm : `heteromtility` data object.
#' @param max_pcs : integer. number of principle components to keep.
#'
#' @return hm : `heteromtility` data object.
#'   assigns @pcs
#'
hmPCA <- function(hm, max_pcs=30){
  pca = prcomp(hm@data)
  hm@pcs = pca$x[,1:max_pcs]

  return(hm)
}

#' Performs hierarchical clustering and sets `meta.data$clust` in `hm`
#'
#' @param hm : `heteromtility` data object.
#' @param k : integer. number of clusters to set with `cutree`.
#' @param linkage : character. method for hierarchical clustering, compatible with `hclust`.
#' @param dist_method : character. method for distance matrix calculation, compatible with `dist`.
#' @param scalar : character. scalar data space to use. {'data', 'unscaled.data', 'pcs'}.
#'
#' @return hm : `heteromtility` data object.
#'   adds 'clust' variable to `hm@meta.data`
#'
hmHClust <- function(hm, k=2, linkage='ward.D2', dist_method='euclidean', scalar='data'){
  dm = dist(attr(hm, scalar), method=dist_method)
  fit = hclust(dm, method=linkage)
  clust = cutree(fit, k=k)
  hm@meta.data$clust = as.factor(clust)

  return(hm)
}

#' Performs tSNE on `pcs` in `hm`
#'
#' @param hm : `heteromtility` data object.
#' @param perplexity : integer. perplexity parameter for `Rtsne`.
#'
#' @return hm : `heteromtility` data object.
#'   assigns `hm@tsne`
#'
hmTSNE <- function(hm, perplexity=50){
  require(Rtsne)
  tsne = Rtsne(hm@pcs, perplexity = perplexity)
  hm@tsne = tsne$Y
  return(hm)
}

#' Perform an ANOVA test for a single feature
anova_feature <- function(df, df.groups, feature_name, plot_path, prefix=''){
  groups <- as.factor(df.groups)
  av <- aov(df[,feature_name] ~ groups, data = df)
  capture.output(summary(av), file = paste(plot_path, prefix, "anova_", feature_name, ".txt", sep =''))
  p <- summary(av)[[1]][['Pr(>F)']][1]
  return(p)
}

#' Perform a t-test for a single feature
ttest_feature <- function(df, df.class, feature_name, class1 = 'Aged', class2 = 'Young', plot_path=FALSE, prefix=''){
  class1_sub <- subset(df, df.class == class1)
  class2_sub <- subset(df, df.class == class2)

  t_result <- t.test(class1_sub[,feature_name], class2_sub[,feature_name], alternative = "two.sided")
  if (is.character(plot_path)){
    capture.output(t_result, file = paste(plot_path, prefix, "ttest_", feature_name, ".txt", sep=''))
  }
  return(t_result$p.value)
}

#' Performs t-tests for each feature in `hm@data` between groups defined by `cat`.
#'
#' Uses a Holm-Bonferroni correction for multiple hypothesis testing.
#'
#' @param hm : `heteromtility` data object.
#' @param cat : character. specifies column name in `hm@meta.data` to use for group definition.
#' @param plot_path : character, optional. path to a directory to save individual t-test summary reports.
#'
#' @return hm : `heteromtility` data object.
#'
hmTTestFeatures <- function(hm, cat='age', plot_path=F){

  classes = levels(hm@meta.data[,cat])

  ttest_results = data.frame(matrix(ncol=3, nrow=ncol(hm@unscaled.data)))
  for (i in 1:ncol(hm@unscaled.data)){
    feature_name <- rownames(hm@feat.data)[i]
    p <- ttest_feature(hm@unscaled.data, hm@meta.data[,cat], feature_name, class1 = classes[1], class2 = classes[2], plot_path = plot_path)
    ttest_results[i,] <- c(feature_name, p, 0)
  }
  colnames(ttest_results) <- c('Feature_Name', 'p_value', 'adj_p_value')
  ttest_results$adj_p_value <- p.adjust(ttest_results$p_value, method = 'holm')
  #ttest_results <- ttest_results[order(ttest_results$adj_p_value),] # order by pvalue, ascending

  hm@feat.data[,paste(cat, '_', 'p_value', sep='')] = ttest_results$p_value
  hm@feat.data[,paste(cat, '_', 'adj_p_value', sep='')] = ttest_results$adj_p_value

  return(hm)
}

#' Performs ANOVAs for each feature in `hm@unscaled.data` between groups defined by `cat`.
#'
#' Uses a Holm-Bonferroni correction for multiple hypothesis testing.
#'
#' @param hm : `heteromtility` data object.
#' @param cat : character. specifies column name in `hm@meta.data` to use for group definition.
#' @param plot_path : character, optional. path to a directory to save individual t-test summary reports.
#'
#' @return hm : `heteromtility` data object.
#'
hmANOVAFeatures <- function(hm, cat='clust', plot_path=F){

  classes = levels(hm@meta.data[,cat])

  anova_results = data.frame(matrix(ncol=3, nrow=ncol(hm@unscaled.data)))
  for (i in 1:ncol(hm@unscaled.data)){
    feature_name <- rownames(hm@feat.data)[i]
    p <- anova_feature(hm@unscaled.data, hm@meta.data[,cat], feature_name, plot_path = plot_path)
    anova_results[i,] <- c(feature_name, p, 0)
  }
  colnames(anova_results) <- c('Feature_Name', 'p_value', 'adj_p_value')
  anova_results$adj_p_value <- p.adjust(anova_results$p_value, method = 'holm')

  hm@feat.data[,paste(cat, '_', 'p_value', sep='')] = anova_results$p_value
  hm@feat.data[,paste(cat, '_', 'adj_p_value', sep='')] = anova_results$adj_p_value

  return(hm)
}

#' Creates a `CellDataSet` object from a `heteromotility` object
#'
#' @param hm : `hetermotility` data object.
#' @param scalar : character. scalar data set to use. {'data', 'unscaled.data', 'pcs'}.
#' @param cat vector of characters. properties from `@meta.data` to include in `phenoData`.
#'
#' @return cds : `CellDataSet` object
#'
hmCreateCellDataset <- function(hm, scalar='data', cat=c('clust')){
  require(monocle)
  # prepare a CellDataSet
  df = attr(hm, scalar)

  fm <- as.matrix( t(df) )

  # Create phenoData object
  pd <- hm@meta.data[,cat]

  # Create featureData object, labeling 'gene_short_name' as the name of the HM feature
  feature_names <- colnames(df)
  fd <- data.frame(feature_names)
  rownames(fd) <- feature_names
  colnames(fd) <- 'gene_short_name'

  exfam=gaussianff()
  # Create CellDataSet object
  pdc <- new("AnnotatedDataFrame", data = pd)
  fdc <- new("AnnotatedDataFrame", data = fd)
  cds <- newCellDataSet(fm, phenoData = pdc, featureData = fdc, expressionFamily = exfam)
  return(cds)
}

#' Performs pseudotiming on motility feature data using `monocle`
#'
#' @param hm : `hetermotility` data object.
#' @param scalar : character. scalar data set to use. {'data', 'unscaled.data', 'pcs'}.
#' @param cat vector of characters. properties from `@meta.data` to include in `phenoData`.
#'
#' @return hm : `heteromtility` data object.
#'   assigns @celldataset.
#'
hmPseudotime <- function(hm, scalar='data', cat=c('clust')){
  df.pc <- princomp(attr(hm, scalar))   # Find ordering features as loadings of PC1
  aload <- abs( with(df.pc, unclass(loadings)) ) # get abs of all loadings
  pload <- sweep(aload, 2, colSums(aload), "/") # proportion per PC
  ordering_features <- tail(sort(pload[,1]), 20) # top 20 loadings in PC1, get names with labels()

  cds = hmCreateCellDataset(hm, scalar=scalar, cat=cat)
  cds <- setOrderingFilter(cds, labels(ordering_features))
  cds <- reduceDimension(cds, max_components = 2, norm_method = 'none', pseudo_expr = 0)
  cds <- orderCells(cds)
  hm@celldataset = cds
  return(hm)
}





