###################### Lc50 ssd plots for cu and zn

#### For Cu

# read dataset
# the file argument of read_csv() assumes the file is in your working directory. You may need to change the file path to correctly read your dataset.

data <- read.csv("C:/Users/s443727/OneDrive - University of Canberra/Bd regression analysis/SSD_96h_LC50_Cu_frogs.csv")


# fit distributions
dist <- ssd_fit_dists(data, left = 'Concentration', dists = c('lgumbel'), silent = TRUE, reweight = FALSE, min_pmix = 0, nrow = 6L, computable = TRUE, at_boundary_ok = FALSE, rescale = TRUE)
# plot distributions
ssd_plot_cdf(dist, delta = Inf)
# goodness of fit table
ssd_gof(dist)

# save plot
# width and height are in inches, dpi (dots per inch) sets resolution
ggsave('fit_dist_plot.png', width = 8 , height = 6 , dpi = 300)

# plot model average
# to add confidence intervals set ci = TRUE in predict and ssd_plot
# we recommend using nboot = 10000 in predict, although this may take several minutes to run
pred <- predict(dist, nboot = 10L, ci = TRUE)
new_cu <- ssd_plot(data, pred, left = 'Concentration', label = 'Species', color = NULL, shape = NULL, hc = 5L, ci = TRUE,
                   shift_x = 1.3, xlab = '96 hours LC50 (µg/L)', ylab = 'Species affected (%)') +
  ggtitle('') + ggeasy::easy_all_text_size(8) +ggeasy::easy_all_text_colour("black")+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, colour='black'),
        axis.text = element_text(color = 'black'),
        legend.key = element_rect(fill = NA, colour = NA)) +
  expand_limits(x = 4567.13) +
  scale_color_brewer(palette = 'Dark2', name = '-none-') +
  scale_shape(name = NULL)

new_cu






#### For Zn
### Import the file 

data <- read.csv("C:/Users/s443727/OneDrive - University of Canberra/Bd regression analysis/SSD_96h_LC50_Zn_frogs.csv")

# fit distributions
dist <- ssd_fit_dists(data, left = 'Concentration', dists = c('lgumbel'), silent = TRUE, reweight = FALSE, min_pmix = 0, nrow = 6L, computable = TRUE, at_boundary_ok = FALSE, rescale = TRUE)
# plot distributions
ssd_plot_cdf(dist, delta = Inf)
# goodness of fit table
ssd_gof(dist)

# save plot
# width and height are in inches, dpi (dots per inch) sets resolution
ggsave('fit_dist_plot.png', width = 8 , height = 6 , dpi = 300)

# plot model average
# to add confidence intervals set ci = TRUE in predict and ssd_plot
# we recommend using nboot = 10000 in predict, although this may take several minutes to run
pred <- predict(dist, nboot = 10L, ci = TRUE)
new_zinc <- ssd_plot(data, pred, left = 'Concentration', label = 'Species', color = NULL, shape = NULL, hc = 5L, ci = TRUE,
                     shift_x = 1.3, xlab = '96 hours LC50 (µg/L)', ylab = 'Species affected (%)') +
  ggtitle('') + ggeasy::easy_all_text_size(8) +ggeasy::easy_all_text_colour("black")+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, colour='black'),
        axis.text = element_text(color = 'black'),
        legend.key = element_rect(fill = NA, colour = NA)) +
  expand_limits(x = 28379.99999) +
  scale_color_brewer(palette = 'Dark2', name = '-none-') +
  scale_shape(name = NULL)


new_zinc


##### Add title to each graph 
new_cu <- new_cu +labs(title = "Cu")
new_zinc <- new_zinc +labs(title = "Zn")


SSD_Cu_Zn <- grid.arrange(new_cu, new_zinc, ncol = 2, nrow=1) 

ggsave("SSD_Cu_Zn.png", SSD_Cu_Zn, width = 8, height = 4, units = "in")

























