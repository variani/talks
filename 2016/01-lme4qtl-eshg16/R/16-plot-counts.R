### inc
library(devtools)
load_all("~/git/hemostat/lme4qtl")

### GAIT2 model
#m <- relmatGlmer(PTES ~ AGEsc + AGEsc2 + SEXf + (1|ID), phen2, relmat = list(ID = dkin2), family = poisson)

### plot
p <- ggplot(phen2, aes(PTES)) + geom_histogram(binwidth = 20, color = "white", fill = "red")

p <- p + theme_void()
 
