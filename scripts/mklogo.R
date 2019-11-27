#!/usr/bin/env Rscript
#
# Create the logo (the argument is the filename)

args <- commandArgs(trailingOnly = TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) == 1) {
  stopifnot(is.character(args[1]))
} else if (length(args) == 0) {
  args[1] <- "out.txt"
} else {
  stop("The only argument should be the filename.")
}


# load required packages
library(tidyverse)
library(ggforce)
library(viridis)
library(hexSticker)
library(showtext)

# number of grades in the top tree
d <- 7

# get the grapes
circles <- data.frame(
  x0 = map(1:(d - 1), function(i) seq(
      1 + (i - 1) * 0.5,
      d - (i - 1) * 0.5
    )) %>% unlist(),
  y0 = map(1:(d - 1), function(i) rep(d - i, d - i + 1)) %>% unlist(),
  r = rep(0.4, sum(seq(2, d)))
)

# behold the circles
p <- ggplot() +
  geom_circle(aes(x0 = x0, y0 = y0, r = r, fill = y0), data = circles) +
  scale_fill_viridis(direction = -1, begin = 0, end = 0.25) +
  theme_void() +
  theme_transparent() +
  theme(legend.position = "none")


# Loading Google fonts (http://www.google.com/fonts)
font_add_google("Gochi Hand", "gochi")
# Automatically use showtext to render text for future devices
showtext_auto()

# create the final sticker
cols <- viridis_pal()(20)
sticker(p,
  package = "pyvinecopulib",
  p_size = 22, p_y = 1.3, s_x = 1, s_y = .69, s_width = 1.2, s_height = 0.8,
  p_family = "gochi",
  h_color = cols[15],
  h_fill = cols[10],
  filename = args[1]
)
