# plot from heteromotilty object
require(ggplot2)
require(reshape2)

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
