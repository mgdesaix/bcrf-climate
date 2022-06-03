# Quantifying the amount of climate extrapolation from gradient forest models

Here, I use the Niche Margin Index (Broennimann et al. 2021) to model the niche from the observed sampling locations of our data. This is then used to calculate the distance from the niche margins (whether within or outside the niche space) from all the cells of predicted genetic composition in the gradient forest model. Since the gradient forest model is trained (and tested - inherently with the out-of-bag sampling during the modeling process) on the observed environment from the sampling sites and cells with different environments require extrapolation from the gradient forest turnover functions. By default, this is done with linear extrapolation in the `gradientForest` package (see https://rdrr.io/rforge/gradientForest/man/predict.gradientForest.html). Extrapolation can be problematic. Here, we use the niche margin index to calculate *how* problematic extrapolations may be in different regions based on the climate magnitude from our trained values. Provided are the scripts for calculating the niche margins ([run-nmi-sites.Rmd](./r-scripts/run-nmi-sites.Rmd)) and making pretty maps ([automate-nmi-figures.Rmd](./r-scripts/automate-nmi-figures.Rmd))


## References

Broennimann et al. 2021. Distance to native climatic niche margins explains establishment of alien mammals. Nature Communications. 12, 2353.
