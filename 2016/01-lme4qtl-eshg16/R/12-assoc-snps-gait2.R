### inc
library(plyr)

load_all("~/git/ugcd/gait")

### data
load("data/snps.throm.gait2.RData") # -> snps,  gf

# update `snpf`
snpf <- within(snpf, {
  maf.imputed <- laply(as.character(Marker), function(x) {
    maf <- NA
    if(x %in% colnames(gf)) {
      maf <- 1 - 0.5 * sum(gf[, x]) / nrow(gf)
      maf <- ifelse(maf > 0.5, 1 - maf, maf)
    }
    
    return(maf)
  })
})

snpf <- mutate(snpf,
  mac.imputed = maf.imputed * nrow(gf))

# create `dat`
dat <- gait2.phen(id.imputed = TRUE, traits = NULL)
dat <- mutate(dat,
  SEXfnum = as.numeric(SEXf) - 1)

dat <- mutate(dat,
  AGEf = as.factor(ifelse(is.na(AGE), NA, 
    ifelse(AGE < 55, 0,
    ifelse(AGE < 64, 1,
    ifelse(AGE < 74, 2,
    ifelse(AGE < 84, 3, 4)))))),
  AGEf2 = as.factor(ifelse(is.na(AGE), NA, 
    ifelse(AGE < 55, 0, 1))),
  AGEf3 = as.factor(ifelse(is.na(AGE), NA, 
    ifelse(AGE < 55, 0,
    ifelse(AGE < 74, 1, 2)))))  

dat <- mutate(dat,
  AGEf3num0 = as.numeric(AGEf3 == 0),
  AGEf3num1 = as.numeric(AGEf3 == 1),
  AGEf3num2 = as.numeric(AGEf3 == 2))

dat <- mutate(dat,
  AGEfnum0 = as.numeric(AGEf == 0),
  AGEfnum1 = as.numeric(AGEf == 1),
  AGEfnum2 = as.numeric(AGEf == 2),
  AGEfnum3 = as.numeric(AGEf == 3),
  AGEfnum4 = as.numeric(AGEf == 4))
  
dat <- join(dat, gf)

K <- 0.05
dat <- mutate(dat, offset = log(K/(1 - K)))

### The null model
#m0 <- relmatGlmer(VT ~ -1 + AGEsc + AGEsc2 + SEXfnum + (1|ID), dat, offset = offset, relmat = list(ID = dkin2), family = binomial)
m0 <- relmatGlmer(VT ~ -1 + AGEf3num1 + AGEf3num2 + SEXfnum + (1|ID), dat, offset = offset, relmat = list(ID = dkin2), family = binomial)

models.snp1 <- dlply(snpf, "Marker", function(x) x$MAF)

### rs8176719 (ABO)
m <- update(m0, . ~ . + rs8176719)

### rs1801020 (f12, RR 5)
m <- update(m0, . ~ . + rs1801020)

### two SNPs
m <- update(m0, . ~ . + rs8176719 + rs1801020)

### rs1799963 (f2, RR 4), MAF ~ 2%
# failed
#m <- update(m0, . ~ . + rs1799963)


