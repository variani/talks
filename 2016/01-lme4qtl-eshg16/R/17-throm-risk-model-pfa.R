### inc
library(plyr)

load_all("~/git/ugcd/solarius")
load_all("~/git/ugcd/gait")
load_all("~/git/hemostat/lme4qtl")


### data
load("data/snps.throm.gait2.RData") # -> snps,  gf

K <- 0.05
dat <- mutate(phen2, offset = log(K/(1 - K)))

dat <- within(dat, {
  PFAadp_q10 <- factor(ifelse(PFAadp < quantile(PFAadp[Throm == 1], 0.1, na.rm =  TRUE), "low", "high"))
  PFAepin_q10 <- factor(ifelse(PFAepin < quantile(PFAepin[Throm == 1], 0.1, na.rm =  TRUE), "low", "high"))
  PFAadp_q20 <- factor(ifelse(PFAadp < quantile(PFAadp[Throm == 1], 0.2, na.rm =  TRUE), "low", "high"))
  PFAepin_q20 <- factor(ifelse(PFAepin < quantile(PFAepin[Throm == 1], 0.2, na.rm =  TRUE), "low", "high"))
})

dat <- mutate(dat, 
  PFAadp_q10num = as.numeric(PFAadp_q10) - 1,
  PFAepin_q10num = as.numeric(PFAepin_q10) - 1,
  PFAadp_q20num = as.numeric(PFAadp_q20) - 1,
  PFAepin_q20num = as.numeric(PFAepin_q20) - 1)

dat <- join(dat, gf)

### model
dat1 <- subset(dat, !is.na(PFAadp_q10num))

m10 <- relmatGlmer(Throm ~ -1 + AGEsc + AGEsc2 + SEXfnum + (1|ID), dat1, 
  offset = offset, relmat = list(ID = dkin2), family = binomial)
m1 <- update(m10, . ~ . + PFAadp_q10num)

dat2 <- subset(dat, !is.na(PFAepin_q10num))

m20 <- relmatGlmer(Throm ~ -1 + AGEsc + AGEsc2 + SEXfnum + (1|ID), dat2, 
  offset = offset, relmat = list(ID = dkin2), family = binomial)
m2 <- update(m20, . ~ . + PFAepin_q10num)

### AB0 (rs8176719)
dat1 <- subset(dat, !is.na(PFAadp_q10num) & !is.na(rs8176719))

mod10 <- relmatGlmer(Throm ~ -1 + AGEsc + AGEsc2 + SEXfnum + rs8176719 + (1|ID), dat1, 
  offset = offset, relmat = list(ID = dkin2), family = binomial)
mod1 <- update(mod10, . ~ . + PFAadp_q10num)

dat2 <- subset(dat, !is.na(PFAepin_q10num) & !is.na(rs8176719))

mod20 <- relmatGlmer(Throm ~ -1 + AGEsc + AGEsc2 + SEXfnum + rs8176719 + (1|ID), dat2, 
  offset = offset, relmat = list(ID = dkin2), family = binomial)
mod2 <- update(mod20, . ~ . + PFAepin_q10num)

if(FALSE) 
{

### AB0 
dat1 <- subset(dat, !is.na(PFAadp_q10num) & !is.na(ABOf3num))

mod10 <- relmatGlmer(Throm ~ -1 + AGEsc + AGEsc2 + SEXfnum + ABOf3num + (1|ID), dat1, 
  offset = offset, relmat = list(ID = dkin2), family = binomial)
mod1 <- update(mod10, . ~ . + PFAadp_q10num)

dat2 <- subset(dat, !is.na(PFAepin_q10num) & !is.na(ABOf3num))

mod20 <- relmatGlmer(Throm ~ -1 + AGEsc + AGEsc2 + SEXfnum + ABOf3num + (1|ID), dat2, 
  offset = offset, relmat = list(ID = dkin2), family = binomial)
mod2 <- update(mod20, . ~ . + PFAepin_q10num)


### AB0 (rs8176719) / AGEsc2 skipped
dat1 <- subset(dat, !is.na(PFAadp_q10num) & !is.na(rs8176719))

mod10 <- relmatGlmer(Throm ~ -1 + AGEsc + SEXfnum + rs8176719 + (1|ID), dat1, 
  offset = offset, relmat = list(ID = dkin2), family = binomial)
mod1 <- update(mod10, . ~ . + PFAadp_q10num)

dat2 <- subset(dat, !is.na(PFAepin_q10num) & !is.na(rs8176719))

mod20 <- relmatGlmer(Throm ~ -1 + AGEsc + SEXfnum + rs8176719 + (1|ID), dat2, 
  offset = offset, relmat = list(ID = dkin2), family = binomial)
mod2 <- update(mod20, . ~ . + PFAepin_q10num)

# exp(cbind(fixef(mod2), confint(mod2, method = "Wald", parm = names(fixef(mod2)))))

}
