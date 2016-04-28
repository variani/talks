### inc
library(ggplot2)

library(devtools)
load_all("~/git/hemostat/lme4qtl")

### par
N <- 20

### model
m <- relmatLmer(APTT ~ AGEc + (1|ID), phen, relmat = list(ID = dkin))

### pred
dat <- subset(phen, !is.na(APTT))

set.seed(1)
ind <- sample(seq(1, nrow(dat)), N)
dat <- dat[ind, ]

dat$APTTpf <- predict(m, dat, re.form = NA)
dat$APTTpfr <- predict(m, dat, re.form = NULL)

p1 <- ggplot(dat) + 
  geom_line(aes(AGE, APTTpf), color = "red", size = 2) + 
  geom_point(aes(AGE, APTT)) + 
  geom_point(aes(AGE, APTTpfr), color = "red") + 
  theme_void()

p1 <- p1 + geom_segment(aes(x = AGE, y = APTTpfr, xend = AGE, yend = APTT), linetype = 3)

p2 <- ggplot(dat, aes(AGE, APTT)) + geom_point() + 
  geom_line(aes(AGE, APTTpf), color = "red", size = 2) + theme_void()
