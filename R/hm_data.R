# heteromotility data class
require('cluster')
library(tibble)

attr_exists <- function(object, attr){
  return(attr %in% names(attributes(object)))
}

heteromotility <- methods::setClass(
  'heteromotility',
  slots = c(raw.data = 'data.frame', # raw motility statistics data frame
            data = 'data.frame', # scaled scalar data in raw.data
            unscaled.data = 'data.frame', # scalar data in raw.data
            meta.data = 'data.frame', # meta data about specific rows
            joint.data = 'data.frame', # joint raw and meta data, useful for group summarization
            feat.data = 'data.frame', # data about specific features
            timeseries = 'logical', # boolean if data is a timeseries or not
            pcs = 'matrix', # principal components
            pca = 'ANY', # prcomp object
            tsne = 'matrix', # tsne embedding
            dist.matrix = 'matrix', # distance matrix
            celldataset = 'ANY', # monocle CellDataSet
            state_positions = 'array', # N, T, PCs array of cell positions in PCA space over time
            transitions = 'array', # N, T-1, PCs array of state transition vectors in PCA space over time
            statewise_transitions = 'matrix') # States, PCs array of characteristic state transition vectors
)

#' Build a `heteromotility` data object from a data.frame
#'
#' @param raw.data : data.frame of raw numeric data with 2 id columns
#' @param meta.data.list : list of character vectors or factors to include in meta.data
#' @param timeseries logical. indicates whether a time series is encoded in `raw.data$cell_id`
#'
#' @return hm : heteromotility data object.
#'
hmMakeObject <- function(raw.data, meta.data.list = list(), timeseries = F){

  meta.data = raw.data[,1:2]
  if (length(meta.data.list) > 0){
    for (i in 1:length(meta.data.list)){
      var.name = names(meta.data.list)[i]
      meta.data[var.name] = as.factor(meta.data.list[[i]])
    }
  }
  if (timeseries){
    idents <- data.frame(do.call('rbind', strsplit(as.character(raw.data$cell_id),'-',fixed=TRUE)))
    idents$X1 = as.numeric(as.character(idents$X1))
    idents$X2 = as.numeric(as.character(idents$X2))
    raw.data$cell_id <- idents$X1
    raw.data <- add_column(raw.data, timepoint = as.numeric(idents$X2), .after=2)
    meta.data['Timepoint'] = as.numeric(idents$X2)
    feat_start = 4
  } else{
    feat_start = 3
  }
  
  joint.data = cbind(raw.data, meta.data[,3:ncol(meta.data)])
  
  hm = heteromotility(raw.data = raw.data,
                      unscaled.data = raw.data[,feat_start:ncol(raw.data)],
                      meta.data = meta.data,
                      joint.data = joint.data,
                      dist.matrix = matrix(),
                      feat.data = data.frame(row.names = colnames(raw.data[,feat_start:ncol(raw.data)])),
                      timeseries=timeseries)
  return(hm)
}

#' Collate summary statistics
#' 
#' @param hm : `heteromotility` data object.
#' @param scalar : character. data attribute to use for summary statistics. 
#'                 ('data', 'unscaled.data', 'pcs', 'tsne')
#'
#' @return summary data.frame of summary statistics
#' 
hmSummaryStatistics <- function(hm, scalar='unscaled.data'){
  data = attr(hm, scalar)
  
  count.na <- function(x){sum(is.na(x))}
  central.moment <- function(x, m){mean( (x-mean(x))**m )}
  std.central.moment <- function(x, m){central.moment(x,m)/sd(x)**m}
  
  means = apply(data, 2, mean)
  sds = apply(data, 2, sd)
  nas = apply(data, 2, count.na)
  zeros = apply(data, 2, function(x){sum(x==0)})
  skew = apply(data, 2, function(x){std.central.moment(x, 3)})
  kurtosis = apply(data, 2, function(x){std.central.moment(x, 4)})
  
  summary = data.frame(Mean=means, 
                       SD=sds,
                       CV=sds/means,
                       Skew = skew,
                       Kurt = kurtosis,
                       NACount=nas, 
                       Zeros=zeros, 
                       Feature=colnames(data))
  return(summary)
}

#' Scales and centers `raw.data` in a heteromotility object `hm`
#'
#' Uses `scale()` to set `mean = 0` and `std = 1` for each feature in `@unscaled.data`
#'
#' @param hm : `heteromtility` data object.
#'
#' @return hm : `heteromtility` data object.
#'   assigns @data.
#'
hmScaleData <- function(hm){
  scaled = data.frame(scale(hm@unscaled.data))
  
  # check for nan's
  if (sum(is.na(scaled)) > 0){
    print('NAs encountered in scaling')
    idx = which(is.na(scaled), arr.ind = T)
    print('NAs in : ')
    print(colnames(scaled)[unique(idx[,2])])
    print('Coercing entire feature to 0s...')
    scaled[,unique(idx[,2])] = 0
  }
  hm@data = scaled
  return(hm)
}

#' Gets data from a `heteromotility` object and returns a data.frame
#'
#' @param hm : `heteromtility` data object.
#' @param scalar : character. property for numerical data to be retrieved. {'data', 'raw.data', 'pcs', etc}
#' @param cat : iterable of strings. categorical features to retrieve from `@meta.data`
#'
#' @return df : data.frame. contains data in `scalar` and retrieved `@meta.data`
#'
hmGetData <- function(hm, scalar='data', cat=c('clust')){

  df = data.frame(attr(hm, scalar))

  for (c in cat){
    df[c] = hm@meta.data[,c]
  }

  return(df)
}

#' Retrives a subset of cell data in heteromotility object `hm`
#'
#' Useful for selecting only cells of a certain type or class.
#' `cat` and `id` parameters must be supplied together.
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
  if (attr_exists(hm, 'data')){ new@data = new@data[idx,] }
  if (attr_exists(hm, 'unscaled.data')){ new@unscaled.data = new@unscaled.data[idx,] }
  new@raw.data = new@raw.data[idx,]
  new@meta.data = new@meta.data[idx,]

  if (nrow(hm@dist.matrix) == nrow(hm@raw.data)){
    new@dist.matrix = new@dist.matrix[idx,]
  }
  if (nrow(hm@pcs) == nrow(hm@raw.data)){
    new@pcs = new@pcs[idx,]
  }
  if (nrow(hm@tsne) == nrow(hm@raw.data)){
    new@tsne = new@tsne[idx,]
  }

  return(new)
}

#' Performs PCA on `data` in a heteromotility object `hm`
#'
#' @param hm : `heteromtility` data object.
#' @param max_pcs : integer. number of principle components to keep.
#'
#' @return hm : `heteromtility` data object.
#'   assigns `@pcs`
#'
hmPCA <- function(hm, max_pcs=30){
  pca = princomp(hm@data)
  scores = pca$scores[,1:max_pcs]
  colnames(scores) = paste('PC', seq(1,ncol(scores)), sep='')
  hm@pcs = scores
  hm@pca = pca
  return(hm)
}

#' Performs hierarchical clustering and sets `meta.data$clust` in `hm`
#'
#' @param hm : `heteromtility` data object.
#' @param k : integer. number of clusters to set with `cutree()`.
#' @param linkage : character. method for hierarchical clustering, compatible with `hclust()`.
#' @param dist_method : character. method for distance matrix calculation, compatible with `dist()`.
#' @param scalar : character. scalar data space to use. ['data', 'unscaled.data', 'pcs'].
#'
#' @return hm : `heteromtility` data object.
#'   adds `clust` variable to `@meta.data`
#'
hmHClust <- function(hm, k=2, linkage='ward.D2', dist_method='euclidean', scalar='data'){
  dm = dist(attr(hm, scalar), method=dist_method)
  fit = hclust(dm, method=linkage)
  clust = cutree(fit, k=k)
  hm@meta.data$clust = as.factor(clust)

  return(hm)
}

#' Order clusters by avg_speed
#' 
#' @param hm : `heteromtility` data object.
#' 
hmOrderClusters <- function(hm){
  tmp = hm@data
  tmp$clust = hm@meta.data$clust
  
  mean_speeds = tmp %>% group_by(clust) %>% summarise(mean(avg_speed))
  sorted_by_speed = mean_speeds[order(mean_speeds$`mean(avg_speed)`),]
  
  new_clust = numeric(nrow(hm@meta.data))
  nc = 1
  for (i in sorted_by_speed$clust){
    bidx = hm@meta.data$clust == i
    new_clust[bidx] = nc
    nc = nc + 1
  }
  
  hm@meta.data$clust = as.factor(new_clust)
  return(hm)
}

#' Calculates silhouette values for clusters
#' 
#' @param hm : `heteromotility` data object.
#' @param k_range : numeric, length 2. range of `k` parameter values to try.
#' @param linkage : character. method for hierarchical clustering, compatible with `hclust()`.
#' @param dist_method : character. method for distance matrix calculation, compatible with `dist()`.
#' @param scalar : character. scalar data space to use. ['data', 'unscaled.data', 'pcs'].
#' 
#' @return silhouettes. numeric of silhouette values for each value of `k` in `k_range`.
#' 
hmCalcSilhouettes <- function(hm, k_range=c(2,6), linkage='ward.D2', dist_method='euclidean', scalar='data'){
  
  silhouettes = numeric(k_range[2]-k_range[1]+1)
  dm = dist(attr(hm, scalar), method=dist_method)
  
  fit = hclust(dm, method=linkage)
  i = 1
  for (k in k_range[1]:k_range[2]){
    clust = cutree(fit, k=k)
    s = silhouette(clust, dm)
    silhouettes[i] = summary(s)$avg.width
    i = i + 1
  }
  names(silhouettes) = k_range[1]:k_range[2]
  return(silhouettes)
}

#' Performs tSNE on `pcs` in `hm`
#'
#' @param hm : `heteromtility` data object.
#' @param perplexity : integer. perplexity parameter for `Rtsne`.
#'
#' @return hm : `heteromtility` data object.
#'   assigns `@tsne`
#'
hmTSNE <- function(hm, perplexity=50){
  if (nrow(hm@pcs) < nrow(hm@raw.data)){
    print('PCs must be computed first!')
    stop()
  }
  
  require(Rtsne)
  # use PCA intialization
  tsne = Rtsne(hm@pcs, perplexity = perplexity, Y_init = hm@pcs[,1:2], exaggeration_factor = 20)
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
#' @param method : multiple hypothesis correction method, compatible with `p.adjust`.
#'
#' @return hm : `heteromtility` data object.
#'
hmTTestFeatures <- function(hm, cat='age', plot_path=F, method='holm'){

  classes = levels(hm@meta.data[,cat])

  ttest_results = data.frame(matrix(ncol=3, nrow=ncol(hm@unscaled.data)))
  for (i in 1:ncol(hm@unscaled.data)){
    feature_name <- rownames(hm@feat.data)[i]
    p <- ttest_feature(hm@unscaled.data, hm@meta.data[,cat], feature_name, class1 = classes[1], class2 = classes[2], plot_path = plot_path)
    ttest_results[i,] <- c(feature_name, p, 0)
  }
  colnames(ttest_results) <- c('Feature_Name', 'p_value', 'adj_p_value')
  ttest_results$adj_p_value <- p.adjust(ttest_results$p_value, method = method)
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
#' @param plot_path : character, optional. path to a directory to save individual test summary reports.
#' @param prefix : character, optional. prefix to filenames of individual tests saved in `plot_path`.
#'
#' @return hm : `heteromtility` data object.
#'
hmANOVAFeatures <- function(hm, cat='clust', plot_path=F, prefix=''){

  classes = levels(hm@meta.data[,cat])

  anova_results = data.frame(matrix(ncol=3, nrow=ncol(hm@unscaled.data)))
  for (i in 1:ncol(hm@unscaled.data)){
    feature_name <- rownames(hm@feat.data)[i]
    p <- anova_feature(hm@unscaled.data, hm@meta.data[,cat], feature_name, plot_path = plot_path, prefix=prefix)
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
  hm@meta.data$Pseudotime = hm@celldataset@phenoData@data$Pseudotime
  return(hm)
}

#' Calculate the state positions and transition vectors for each cell
#' 
#' @param hm `heteromotility` data object. must be a timeseries.
#' @param eigenvectors matrix [N, N] of eigenvectors to change basis of the data in `hm`
#'   useful for mapping the state space of variable measurements into the state space defined for 
#'   the whole timelapse experiment.
#' @param n_pcs integer. number of PCs to use to define transition vectors.
#' 
#' @return heteromotility data object with the `state_positions` and `transitions` attributes filled.
#' 
hmCalcTransitions <- function(hm, eigenvectors=NULL, n_pcs=2){
  if (!hm@timeseries){
    stop('Transitions can only be calculated if the object has timeseries data. hm@timeseries is FALSE here.')
  }
  
  if (!is.null(eigenvectors)){
    scores = as.matrix(hm@data) %*% as.matrix(eigenvectors)
  } else{
    scores = hm@pcs
  }
  minT = min(hm@meta.data$Timepoint)
  maxT = max(hm@meta.data$Timepoint)
  
  n_cells = nrow(hm@data) / (maxT + 1 - minT)
  # define the state position matrix [N, T, Dims]
  state_positions = array(NA, dim=c(n_cells, maxT+1-minT, n_pcs))
  # fill the array
  for (t in minT:maxT){
    tpoint = scores[hm@meta.data$Timepoint==t,]
    state_positions[,t+1,] = tpoint[,1:n_pcs]
  }
  # get a transition array
  transitions = state_positions[,2:(maxT+1),] - state_positions[,1:maxT,]
  
  hm@state_positions = state_positions
  hm@transitions = transitions
  return(hm)
}

#' Calculate characteristic state transitions for a predefined categorical label
#' 
#' @param hm `heteromotility` data object with `state_positions` and `transitions` attributes.
#' @param cat character. categorical variable for transition vector calculation, column in `meta.data`.
#' @param func callable. summary function to use, i.e. mean, sum, etc.
#' @param bootstrap integer. number of bootstrap runs to perform to provide error estimates.
#' 
#' @return heteromotility object with [N_States, Dims] matrix of characteristic vectors.
#' 
hmStatewiseTransitions <- function(hm, cat='clust', func=mean, bootstrap=NULL){
  states = as.character(unique(hm@meta.data[cat])[,1])
  
  statewise_transitions = matrix(NA, nrow=length(states), ncol=dim(hm@transitions)[3])
  
  meta.data = hm@meta.data
  transitions = hm@transitions

  meta.data.t0 = meta.data[meta.data$Timepoint=='0',]
  
  for (i in 1:length(states)){
    s = states[i]
    trans = transitions[meta.data.t0[cat]==s,,]
    statewise_transitions[i,] = apply(trans, 3, func) # apply func to each cell across all timepoints
  }
  colnames(statewise_transitions) = paste('PC', 1:ncol(statewise_transitions), sep='')
  rownames(statewise_transitions) = states
  hm@statewise_transitions = statewise_transitions
  return(hm)
}

