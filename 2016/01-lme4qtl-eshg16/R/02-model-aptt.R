### inc
library(ggplot2)

library(devtools)
load_all("~/git/hemostat/lme4qtl")

### model
m <- relmatLmer(APTT ~ AGEc + (1|ID), phen, relmat = list(ID = dkin))

### pred
dat <- subset(phen, !is.na(APTT))
dat$APTTpf <- predict(m, re.form = NA)
dat$APTTpfr <- predict(m, re.form = NULL)

p1 <- ggplot(dat, aes(AGE, APTT)) + geom_point() + 
  geom_point(aes(AGE, APTTpfr), color = "red", size = 1) + 
  geom_line(aes(AGE, APTTpf), color = "red", size = 2) + theme_void()

p2 <- ggplot(dat, aes(AGE, APTT)) + geom_point() + 
  geom_line(aes(AGE, APTTpf), color = "red", size = 2) + theme_void()
