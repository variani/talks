### par
N <- 50

K <- 0.05

### GAIT2
dat <- mutate(phen2, offset = log(K/(1 - K)))

m <- relmatGlmer(Throm ~ -1 + AGEsc + SEXfnum + (1|ID), dat, offset = offset, relmat = list(ID = dkin2), family = binomial)
m0 <- relmatGlmer(Throm ~ AGEsc + SEXfnum + (1|ID), dat, relmat = list(ID = dkin2), family = binomial) 
  
dat <- subset(phen2, !is.na(AGEsc))
set.seed(1)
ind <- sample(with(dat, which(SEXfnum == 0)), N)
dat <- dat[ind, ]

beta0 <- log(K/(1 - K))
logit <- function(x) 1/(1 + exp(-x))

dat$y <- predict(m, dat, re.form = NULL)
dat$yf <- predict(m, dat, re.form = NA)
dat <- mutate(dat, prob = logit(y + beta0))
dat <- mutate(dat, probf = logit(yf + beta0))

simdat <- data.frame(AGEsc = seq(-2.5, 5, 0.05), SEXfnum = 0)
simdat$y <- predict(m, simdat, re.form = NA)
simdat$y0 <- predict(m0, simdat, re.form = NA)
simdat <- mutate(simdat, prob = logit(y + beta0))
simdat <- mutate(simdat, prob0 = logit(y0))

p1 <- ggplot(simdat) + 
  geom_line(aes(AGEsc, prob), color = emphcol, size = 2) +
  geom_line(aes(AGEsc, prob0), color = emphcol, size = 2, linetype = 3) +
  geom_point(aes(AGEsc, as.numeric(Throm), shape = as.factor(Throm)), dat, size = 2) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_shape_manual(values = c(1, 19)) + guides(shape = "none") + theme_void()
  
### GAIT1
dat <- mutate(phen, offset = log(K/(1 - K)))

m <- relmatGlmer(aff ~ -1 + AGEsc + SEXfnum + (1|ID), dat, offset = offset, relmat = list(ID = dkin), family = binomial)
m0 <- relmatGlmer(aff ~ AGEsc + SEXfnum + (1|ID), dat, relmat = list(ID = dkin), family = binomial) 
  
dat <- subset(phen, !is.na(AGEsc))
set.seed(1)
ind <- sample(with(dat, which(SEXfnum == 0)), N)
dat <- dat[ind, ]

beta0 <- log(K/(1 - K))
logit <- function(x) 1/(1 + exp(-x))

dat$y <- predict(m, dat, re.form = NULL)
dat$yf <- predict(m, dat, re.form = NA)
dat <- mutate(dat, prob = logit(y + beta0))
dat <- mutate(dat, probf = logit(yf + beta0))

simdat <- data.frame(AGEsc = seq(-2.5, 5, 0.05), SEXfnum = 0)
simdat$y <- predict(m, simdat, re.form = NA)
simdat$y0 <- predict(m0, simdat, re.form = NA)
simdat <- mutate(simdat, prob = logit(y + beta0))
simdat <- mutate(simdat, prob0 = logit(y0))

p2 <- ggplot(simdat) + 
  geom_line(aes(AGEsc, prob), color = emphcol, size = 2) +
  geom_line(aes(AGEsc, prob0), color = emphcol, size = 2, linetype = 3) +
  geom_point(aes(AGEsc, as.numeric(aff), shape = as.factor(aff)), dat, size = 2) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_shape_manual(values = c(1, 19)) + guides(shape = "none") + theme_void()
