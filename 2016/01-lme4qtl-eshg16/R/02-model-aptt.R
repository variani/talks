### inc
library(devtools)
load_all("~/git/hemostat/lme4qtl")

### model
m <- relmatLmer(APTT ~ AGEc + (1|ID), phen, relmat = list(ID = dkin))
