# sigbra
# Copyright (C) 2014 Daniel Gromer
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# <http://www.gnu.org/licenses/>
#
#
# A convenient function for displaying significance asterisks together with
# lines in ggplot2 plots.
#
# Parameters:
# x.lo:        numeric; x-position of lower line ending.
# x.hi:        numeric; x-position of upper line ending.
# y.lo1:       numeric; corresponding lower y-position for x.lo.
# y.lo2:       numeric; corresponding lower y-position for x.hi.
# y.hi:        numeric; corresponding upper y-position for x.lo and x.hi
# label:       character; label to display, default is "*".
# lab.space:   numeric; span between y.hi and label, default is .5
# text.size:   numeric; text size for label, default is 8.
# line.size:   numeric; size of line elements, default is .3.
# x.lo.lo:     numeric; if lower line ending should span over multiple values,
#              x.lo.lo defines the lower ending.
# x.lo.hi:     numeric; corresponding upper value for x.lo.lo
# x.hi.lo:     numeric; if upper line ending should span over multiple values,
#              x.hi.lo defines the lower ending.
# x.hi.hi:     nueric; corresponding upper value for x.hi.lo.
# small.y.len: numeric; bracket length for x.lo.lo/hi and/or x.hi.lo/hi.
# colour:      character; colour of label and line elements.
#
# Value:
# A list of ggplot2 geoms
#
# Examples:
# library(ggplot2)
# data <- data.frame(group = c("a", "b"), mean = c(20, 30), se = c(2.5, 3.0))
# ggplot(data, aes(x = group, y = mean)) +
#   geom_bar(stat = "identity") +
#   geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = .5) +
#   sigbra(1, 2, 25, 35, 40)
#
# data <- data.frame(fac.a = c("a", "a", "b", "b"), fac.b = c("c", "d", "c", "d"),
#                    mean = c(10, 12, 20, 35), se = c(2.4, 3, 2.9, 3.1))
# ggplot(data, aes(x = fac.a, y = mean, fill = fac.b)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = .5,
#                 position = position_dodge(.9)) +
#   sigbra(1.775, 2.225, 25, 40, 45) +
#   sigbra(1, 2, 20, 50, 55, x.lo.lo = .55, x.lo.hi = 1.45, x.hi.lo = 1.55,
#          x.hi.hi = 2.45)

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