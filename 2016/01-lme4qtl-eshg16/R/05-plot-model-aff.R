### inc
library(devtools)
load_all("~/git/hemostat/lme4qtl")

### par
N <- 50

### GAIT1 model
K <- 0.02
#m <- relmatGlmer(aff ~ -1 + AGEsc + SEXfnum + (1|ID), mutate(phen, offset = log(K/(1 - K))), offset = offset, relmat = list(ID = dkin), family = binomial)
#m0 <- relmatGlmer(aff ~ AGEsc + SEXfnum + (1|ID), phen, relmat = list(ID = dkin), family = binomial)

### `dat`
dat <- subset(phen, !is.na(AGEsc))

set.seed(1)
ind <- sample(with(dat, which(SEXfnum == 0)), N)
dat <- dat[ind, ]

### prediction
beta0 <- log(K/(1 - K))
logit <- function(x) 1/(1 + exp(-x))

dat$y <- predict(m, dat, re.form = NULL)
dat <- mutate(dat, prob = logit(y + beta0))

dat$yf <- predict(m, dat, re.form = NA)
dat <- mutate(dat, probf = logit(yf + beta0))

simdat <- data.frame(AGEsc = seq(-2.5, 5, 0.05), SEXfnum = 0)
simdat$y <- predict(m, simdat, re.form = NA)
simdat$y0 <- predict(m0, simdat, re.form = NA)

simdat <- mutate(simdat, prob = logit(y + beta0))
simdat <- mutate(simdat, prob0 = logit(y0))

#simdat <- mutate(simdat, response = predict(m, simdat, re.form = NA, type = "response"))
#simdat <- mutate(simdat, response0 = predict(m0, simdat, re.form = NA, type = "response"))

p1 <- ggplot(simdat) + 
  geom_line(aes(AGEsc, prob)) +
  geom_point(aes(AGEsc, as.numeric(aff), shape = as.factor(aff)), dat, size = 2)

p1 <- p1 + geom_hline(yintercept = 0.5, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2)

p1 <- p1 + scale_shape_manual(values = c(1, 19)) + guides(shape = "none") + theme_void()

p2 <- p1 + geom_line(aes(AGEsc, prob0), linetype = 3)

p3 <- p1 + 
  geom_point(aes(AGEsc, probf, shape = as.factor(aff)), dat) +
  geom_point(aes(AGEsc, prob, shape = as.factor(aff)), dat, color = "red") +
  geom_segment(aes(x = AGEsc, y = probf, xend = AGEsc, yend = prob), dat, linetype = 3)
