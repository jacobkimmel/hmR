# plot from heteromotilty object
library(ggplot2)
library(reshape2)
as_quosure <- function(strs) rlang::parse_quosures(paste(strs, collapse=";"))

sigbra <- function(x.lo, x.hi, y.lo1, y.lo2, y.hi, label = "*", lab.space = .5,
                   text.size = 8, line.size = .3, x.lo.lo = NULL,
                   x.lo.hi = NULL, x.hi.lo = NULL, x.hi.hi = NULL,
                   small.y.len = 1, colour = "black")
{
  out <- list(
    geom_segment(aes(x = x.lo, xend = x.lo, y = y.lo1, yend = y.hi), size = .3,
                 colour = colour),
    geom_segment(aes(x = x.lo, xend = x.hi, y = y.hi, yend = y.hi), size = .3,
                 colour = colour),
    geom_segment(aes(x = x.hi, xend = x.hi, y = y.hi, yend = y.lo2), size = .3,
                 colour = colour),
    annotate("text", x = (x.lo + x.hi) / 2, y = y.hi + 1,
             label = label, size = text.size, colour = colour)
  )
  
  out[[1]]$mapping$x <- x.lo
  out[[1]]$mapping$xend <- x.lo
  out[[1]]$mapping$y <- y.lo1
  out[[1]]$mapping$yend <- y.hi
  out[[1]]$geom_params$size <- line.size
  
  out[[2]]$mapping$x <- x.lo
  out[[2]]$mapping$xend <- x.hi
  out[[2]]$mapping$y <- y.hi
  out[[2]]$mapping$yend <- y.hi
  out[[2]]$geom_params$size <- line.size
  
  out[[3]]$mapping$x <- x.hi
  out[[3]]$mapping$xend <- x.hi
  out[[3]]$mapping$y <- y.hi
  out[[3]]$mapping$yend <- y.lo2
  out[[3]]$geom_params$size <- line.size
  
  out[[4]]$mapping$x <- (x.lo + x.hi) / 2
  out[[4]]$mapping$y <- y.hi + lab.space
  out[[4]]$geom_params$label <- label
  out[[4]]$geom_params$size <- text.size
  
  if (!is.null(x.lo.lo) & !is.null(x.lo.hi))
  {
    i <- length(out) + 1
    
    out[[i]] <- geom_segment(aes(x = x.lo.lo, xend = x.lo.lo, y = y.lo1 - 1,
                                 yend = y.lo1), size = .3, colour = colour)
    out[[i + 1]] <- geom_segment(aes(x = x.lo.lo, xend = x.lo.hi, y = y.lo1,
                                     yend = y.lo1), size = .3, colour = colour)
    out[[i + 2]] <- geom_segment(aes(x = x.lo.hi, xend = x.lo.hi, y = y.lo1,
                                     yend = y.hi - 1), size = .3, colour = colour)
    
    out[[i]]$mapping$x <- x.lo.lo
    out[[i]]$mapping$xend <- x.lo.lo
    out[[i]]$mapping$y <- y.lo1 - small.y.len
    out[[i]]$mapping$yend <- y.lo1
    out[[i]]$geom_params$size <- line.size
    
    out[[i + 1]]$mapping$x <- x.lo.lo
    out[[i + 1]]$mapping$xend <- x.lo.hi
    out[[i + 1]]$mapping$y <- y.lo1
    out[[i + 1]]$mapping$yend <- y.lo1
    out[[i + 1]]$geom_params$size <- line.size
    
    out[[i + 2]]$mapping$x <- x.lo.hi
    out[[i + 2]]$mapping$xend <- x.lo.hi
    out[[i + 2]]$mapping$y <- y.lo1
    out[[i + 2]]$mapping$yend <- y.lo1 - small.y.len
    out[[i + 2]]$geom_params$size <- line.size
  }
  
  if (!is.null(x.hi.lo) & !is.null(x.hi.hi))
  {
    i <- length(out) + 1
    
    out[[i]] <- geom_segment(aes(x = x.hi.lo, xend = x.hi.lo, y = y.lo1 - 1,
                                 yend = y.lo1), size = .3, colour = colour)
    out[[i + 1]] <- geom_segment(aes(x = x.hi.lo, xend = x.hi.hi, y = y.lo1,
                                     yend = y.lo1), size = .3, colour = colour)
    out[[i + 2]] <- geom_segment(aes(x = x.hi.hi, xend = x.hi.hi, y = y.lo1,
                                     yend = y.hi - 1), size = .3, colour = colour)
    
    out[[i]]$mapping$x <- x.hi.lo
    out[[i]]$mapping$xend <- x.hi.lo
    out[[i]]$mapping$y <- y.lo2 - small.y.len
    out[[i]]$mapping$yend <- y.lo2
    out[[i]]$geom_params$size <- line.size
    
    out[[i + 1]]$mapping$x <- x.hi.lo
    out[[i + 1]]$mapping$xend <- x.hi.hi
    out[[i + 1]]$mapping$y <- y.lo2
    out[[i + 1]]$mapping$yend <- y.lo2
    out[[i + 1]]$geom_params$size <- line.size
    
    out[[i + 2]]$mapping$x <- x.hi.hi
    out[[i + 2]]$mapping$xend <- x.hi.hi
    out[[i + 2]]$mapping$y <- y.lo2
    out[[i + 2]]$mapping$yend <- y.lo2 - small.y.len
    out[[i + 2]]$geom_params$size <- line.size
  }
  
  out
}

pl_theme <- theme(axis.line.x = element_line(size=0.75),
                  axis.line.y = element_line(size=0.75),
                  axis.text = element_text(size = rel(1)),
                  axis.ticks.x = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.title = element_text(size=rel(1.3)),
                  legend.position="right",
                  panel.background=element_blank(),
                  panel.border=element_blank(),
                  panel.grid.major=element_blank(),
                  panel.grid.minor=element_blank(),
                  plot.background=element_blank())



#' Plots PCA
#'
#' @param hm : `heteromtility` data object.
#' @param pc1 : integer. first PC to plot.
#' @param pc2 : integer. second PC to plot
#' @param color : string. name of a variable in meta.data to color points.
#' @param save : string. filename to save plots.
#'
#' @return pl : ggplot2 object.
#'
hmPlotPCA <- function(hm, pc1=1, pc2=2, color='clust', save=NA){


  df = hmGetData(hm, scalar='pcs', cat=c(color))

  pl = ggplot(df, aes_string(x=paste('PC', as.character(pc1), sep=''),
                             y=paste('PC', as.character(pc2), sep=''),
                             color=color)) +
        geom_point(size=2) +
        pl_theme

  if (!is.na(save)){
    ggsave(pl, filename=save, width=5, height=4, units='in')
  }

  return(pl)
}

#' Plot statewise transition vectors atop PCA project
#' 
#' @param hm `heteromotility` data object
#' @param color character. column in `meta.data` to use for coloration
#' @param amp float. coefficient to multiply vector length for visualization.
#' @param save character. path to file for saving.
#' @return ggplot2 object
#' 
hmPlotPCATransitions <- function(hm, color, eigenvectors=NULL, amp=1, save=NULL){
  pl = hmPlotPCA(hm, color=color)
  states = as.character(unique(hm@meta.data[color])[,1])
  centroids = matrix(NA, nrow=length(states), ncol=2)
  for (i in 1:length(states)){
    centroids[i,] = apply(hm@pcs[hm@meta.data[color]==states[i],1:2], 2, mean)
  }
  
  arrows = data.frame(State=states, 
             centroid_x = centroids[,1], 
             centroid_y = centroids[,2],
             t_x = hm@statewise_transitions[,1]*amp,
             t_y = hm@statewise_transitions[,2]*amp)
  
  pl = pl + geom_segment(data=arrows, 
                         aes(x=centroid_x, y=centroid_y,
                          xend=centroid_x + t_x, yend=centroid_y + t_y),
                         size=2,
                         arrow=arrow(length = unit(0.5,"cm")))
  return(pl)
}

#' Plot PCA Loadings
#' 
#' @param hm `hetermotility` data object.
#' @param pcs vector of indices for PC loadings to plot
#' @param save character, path to filename for saving loading plots
#' 
#' @return ggplot2 plot
#' 
hmPlotLoadings <- function(hm, pcs, save=NULL){
  
}

#' Plots TSNE
#'
#' @param hm : `heteromtility` data object.
#' @param color : string. name of a variable in meta.data to color points.
#' @param save : string. filename to save plots.
#'
#' @return pl : ggplot2 object.
#'
hmPlotTSNE <- function(hm, color='clust', save=NA){

  df = hmGetData(hm, scalar='tsne', cat=c(color))

  pl = ggplot(df, aes_string(x='X1',
                             y='X2',
                             color=color)) +
    geom_point(size=2) +
    labs(x='tSNE-1', y='tSNE-2') +
    pl_theme

  if (!is.na(save)){
    ggsave(pl, filename=save, width=5, height=4, units='in')
  }

  return(pl)
}


#' Plots feature means for each group defined by a categorical variable `class` in `meta.data`
#'
#' @param hm : `heteromtility` data object.
#' @param class : character. class in `meta.data` to group bars by.
#' @param limited : boolean. use a limited feature set for visualization.
#' @param sig_bars : boolean. plot significance bars. requires `class`_adj_p_value in `@feat.data`.
#' @param width : float. plot width for saving in inches.
#' @param height : float. plot width for saving in inches.
#' @param save : character. filename for saving.
#'
#' @return pl : ggplot2 object.
#'
hmPlotFeatureMeans <- function(hm, class='clust', limited=F, sig_bars=F, width=8, height=5, save=NULL){

  # Set constants
  limited_features <- c("total_distance",
                        "net_distance",
                        "linearity",
                        "spearmanrsq",
                        "progressivity",
                        "avg_speed",
                        "MSD_slope",
                        "hurst_RS",
                        "nongauss",
                        "rw_kurtosis01",
                        "rw_kurtosis05",
                        "avg_moving_speed05",
                        "time_moving05",
                        "autocorr_5")

  feature_theme <- theme(
    panel.background=element_blank(),
    panel.border=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    plot.background=element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1),
    text= element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  )

  # Group data
  df = hm@data

  if (limited){
    df = df[limited_features]
  }
  df$class = as.factor(hm@meta.data[,class])

  sem <- function(x){sd(x)/sqrt(length(x))}
  means <- aggregate(. ~ class, data = df, FUN = mean)
  stderr <- aggregate(. ~ class, data = df, FUN = sem)

  means.m <- melt(means, value.name = "Mean")
  stderr.m <- melt(stderr, value.name = "SE")
  plt_df <- merge(means.m, stderr.m)

  # plot
  mean_limits <- aes(ymax = plt_df$Mean + plt_df$SE, ymin = plt_df$Mean - plt_df$SE)
  p_mean <- ggplot(plt_df, aes(variable, Mean, fill = class)) +
    geom_bar(position = "dodge", stat = "identity") +
    geom_errorbar(mean_limits, position = position_dodge(0.9), width = 0.3) +
    feature_theme +
    labs(title = "Feature Means", x = "Feature", y = "Mean (+/- SE)")

  # add sig bars
  if (sig_bars){
    y_shift = 0.05; alpha = 0.05; text.size=5; lab.space=0.003;

    sig <- data.frame(matrix(nrow = length(levels(plt_df$variable)), ncol=2))
    colnames(sig) <- c('variable', 'adj_p_val')
    for (i in 1:length(levels(plt_df$variable))){
      sig[i,] <- c(levels(plt_df$variable)[i], hm@feat.data[levels(plt_df$variable)[i], paste(class, '_adj_p_value', sep='')])
    }
    sig$adj_p_val <- as.numeric(sig$adj_p_val)

    plt_df$high_pt <- plt_df$Mean + plt_df$SE
    if (sum(colnames(plt_df)=='class') > 0){
      max_arr <- acast(plt_df, variable ~ class, value.var = 'high_pt')
    } else {
      max_arr <- acast(plt_df, variable ~ groups, value.var = 'high_pt')
    }
    y <- apply(max_arr, 1, max) + y_shift

    sig_marks <- data.frame(matrix(nrow=sum(sig$adj_p_val < alpha), ncol=5))
    colnames(sig_marks) <- c('x', 'xend', 'ylow', 'ylow', 'yhigh')
    j <- 1
    for (i in 1:nrow(sig)){
      if (sig[i,]$adj_p_val < alpha){
        sig_marks[j,] <- c(0.8+(i-1), 1.2+(i-1), y[i], y[i], y[i] + 0.02)
        j <- j + 1
      }
    }

    if (nrow(sig_marks) > 0 ){
      for (i in 1:nrow(sig_marks)){
        p_mean <- p_mean + sigbra(sig_marks[i,1], sig_marks[i,2], sig_marks[i,3], sig_marks[i,4], sig_marks[i,5], label = '*', text.size = text.size, lab.space = lab.space)
      }
    }
  }

  # save plot
  if (!is.null(save)){
    ggsave(filename = save, width=width, height=height, units='in')
  }

  return(p_mean)
}

#' Plots a histogram of a single feature
#' 
#' @param hm `heteromotility` data object.
#' @param feature character. feature name in `scalar` data.frame for plotting.
#' @param scalar character. scalar data space to use. ['data', 'unscaled.data', 'pcs'].
#' @param color character. column in `@meta.data` to color subsamples.
#' @param n_bins integer. number of histogram bins.
#' @param height float. height of the plot in inches.
#' @param width float. width of the plot in inches.
#' @param save character. path to output filename.
#' 
#' @return ggplot2 plot object
#' 
hmPlotFeatureHistogram <- function(hm, feature='total_distance', scalar='data', color=NULL, n_bins=50, height=4, width=5, save=NULL){
  
  data = attr(hm, scalar)
  p = ggplot(data=data, aes_string(x=feature))
  if (!is.null(color)){
    p = p + geom_histogram(aes(fill=hm@meta.data[,color]), bins=n_bins)
  } else{
    p = p + geom_histogram(bins=n_bins)
  }
  
  p = p + pl_theme + scale_y_continuous(expand=c(0,0)) + labs(y='Count', title=paste(feature, 'Histogram', sep = ' '))
  
  if (!is.null(save)){
    ggsave(p, filename=save, width=width, height=height, units='in')
  }
  
  return(p)
}

#' Plots a heatmap of Heteromotility feature values
#'
#' @param hm : `heteromotility` object with scaled data.
#' @param scalar : character. scalar data space to use. ['data', 'unscaled.data', 'pcs'].
#' @param plot_path : character. path to save heatmap.
#' @param prefix : character. prefix to append to heatmap file names.
#' @param height : float. height of plot in inches.
#' @param width : float. width of plot in inches.
#' @param ... : all other arguments accepted by `pheatmap()`
#'
#' @exports Saves a heatmap as PREFIXheatmap.png.
#'
hmPlotHeatmap <- function(hm, scalar='data', plot_path=NULL, prefix='', height=7, width=7, ...){
  if (scalar == 'data'){
    mat = hm@data
  } else if (scalar == 'unscaled.data'){
    mat = hm@unscaled.data
  } else if (scalar == 'pcs'){
    mat = hm@pcs
  } else {
    stop('`scalar` must be in (data, unscaled.data, pcs)')
  }

  png(paste(plot_path, prefix, 'heatmap.png', sep=''), width=width, height=height, res=600, units='in')
}

#' Plots monocle pseudotiming
#'
#' @param hm : `heteromotility` object containing @celldataset trajectory data.
#' @param class : character. feature in FeatureData to color by.
#' @param show_state_number : boolean. show state numbers on pseudotiming plot.
#' @param rev_x_axis : boolean. reverse the order of the pseudotiming x-axis.
#' @param save : character. path to file for saving.
#' @param width : integer. Width of plots in inches.
#' @param height : integer. Height of plots in inches.
#'
#' @return pl : ggplot2 object.
#'
hmPlotPseudotime <- function(hm, class='State', show_state_number=F, rev_x_axis=F, save=NULL, width=4, height=4){
  require(monocle)
  cds = hm@celldataset
  pl <- plot_cell_trajectory(cds,
                             color_by=class,
                             show_state_number = show_state_number)

  if (rev_x_axis){
    pl = pl + scale_x_reverse()
  }

  if (!is.null(save)){
    ggsave(pl, filename=save, width=width, height=height, units='in')
  }
  return(pl)
}

sem <- function(x){sd(x)/sqrt(length(x))}

#' Plot a feature over time in time series data
#' 
#' @param hm heteromotility object.
#' @param feature character. feature in joint.data.
#' 
#' @return ggplot2 object
#' 
PlotFeatureTimeseries <- function(hm, feature='avg_speed', group=NULL, save=NULL){
  if (!is.null(group)){
    groupj = c('Timepoint', group)
  } else(
    groupj = c('Timepoint')
  )
  
  feat_v_time = hm@joint.data %>% group_by(!!!as_quosure(groupj))
  S = feat_v_time %>% summarise_at(.vars = feature, .funs = c('mean', 'sem'))
  pl = ggplot(data=S, aes_string(x='Timepoint', y='mean', color=group)) + geom_line() + 
    geom_ribbon(aes_string(ymin='mean-sem', ymax='mean+sem', fill=group), alpha=0.3) +
    labs(y=feature, x='Timepoint')
  
  if (!is.null(save)){
    ggsave(pl, filename=save, width=6, height=4, units='in')
  }
  
  return(pl)
}
