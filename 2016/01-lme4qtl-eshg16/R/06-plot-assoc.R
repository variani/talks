### inc
library(devtools)
load_all("~/git/hemostat/lme4qtl")

### GAIT1 model

dat <- subset(phen, !is.na(ab0))
#m0 <- relmatLmer(APTT ~ AGEsc + (1|ID), dat, relmat = list(ID = dkin))
#m <- update(m0, . ~. + ab0)

pdat <- subset(dat, select = c("APTT", "AGEsc", "ab0", "ID"))
pdat <- within(pdat, {
  AGEsc <- 0
})

pdat$y <- predict(m, dat, re.form = NA)

### plot
p <- ggplot(pdat, aes(as.factor(ab0), y)) + geom_boxplot(color = "red") 

p <- p + theme_void()
 

