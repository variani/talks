### inc
library(plyr)
library(ggplot2)
library(gridExtra)

library(devtools)
load_all("~/git/ugcd/solarius")
load_all("~/git/hemostat/lme4qtl")

### par
emphcol <- "#d94d3a" # See styles.css
emphcol2 <- "#2a7cdf" # blue #4387fd; blue2 #3c8ef3; blue3 #2a7cdf;

### data
#source("R/01-load-data.R")

#------------------------------
# Figure 1: solarPolygenic
#------------------------------

N <- 20

m <- relmatLmer(APTT ~ AGEc + (1|ID), phen, relmat = list(ID = dkin))

dat <- subset(phen, !is.na(APTT))

set.seed(1)
ind <- sample(seq(1, nrow(dat)), N)
dat <- dat[ind, ]

dat$yf <- predict(m, dat, re.form = NA)
dat$yr <- predict(m, dat, re.form = NULL)

p1 <- ggplot(dat) + 
  geom_line(aes(AGE, yf), color = emphcol, size = 2) + 
  geom_point(aes(AGE, APTT)) + 
  geom_point(aes(AGE, yr), color = emphcol) + 
  geom_segment(aes(x = AGE, y = yr, xend = AGE, yend = APTT), linetype = 3)

p1 <- p1 + labs(title = "solarPolygenic")

p1 <- p1 + theme_void()

#------------------------------
# Figure 2: solarMultipoint
#------------------------------

ids <- c("01101", "01102",
  "01202", "01203", "01205", "01208", "01210",
  "01204", "01209", "01211")

f <- system.file("extdata/gait1/mibd.5.200.gz", package = "gait")
  
mf <- read_mibd_csv_gz(f)
mf <- mf[, 1:3]

mf <- subset(mf, ID1%in% ids & ID2 %in% ids)
mat <- mf2mat(mf)

pf <- melt(mat) # Var1 Var2 value
pf <- mutate(pf,
  Var1 = factor(Var1, levels = as.numeric(ids)),  
  Var2 = factor(Var2, levels = rev(as.numeric(ids))))

p2 <- ggplot(pf, aes(Var1, Var2)) + geom_tile(aes(fill = value), color = "white") + 
  scale_fill_gradient(low = "white", high = emphcol) + guides(fill = "none")

p2 <- p2 + labs(title = "solarMultipoint")

p2 <- p2 + theme_void()

  
#------------------------------
# Figure 3: solarAssoc
#------------------------------
  
dat <- subset(phen, !is.na(ab0))
m0 <- relmatLmer(APTT ~ AGEsc + (1|ID), dat, relmat = list(ID = dkin))
m <- update(m0, . ~. + ab0)

pdat <- subset(dat, select = c("APTT", "AGEsc", "ab0", "ID"))
pdat <- within(pdat, {
  AGEsc <- 0
})

pdat$y <- predict(m, dat, re.form = NA)

p3 <- ggplot(pdat, aes(as.factor(ab0), y)) + geom_boxplot(color = emphcol) 

p3 <- p3 + labs(title = "solarAssoc")

p3 <- p3 + theme_void()

### save figure
p <- marrangeGrob(list(p1, p2, p3), nrow = 1, ncol = 3, top = NULL)
ggsave("solarius-models.pdf", p, width = 12, height = 4)

#pdf("solarius-models.pdf", p, width = 12, height = 4)
#grid.arrange(p1, p2, p3, nrow = 1)
#dev.off()

