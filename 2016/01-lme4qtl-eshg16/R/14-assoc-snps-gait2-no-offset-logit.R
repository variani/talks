### inc
library(plyr)

library(doMC)

load_all("~/git/ugcd/solarius")
load_all("~/git/ugcd/gait")
load_all("~/git/hemostat/lme4qtl")

### par
cores <- 64
parallel <- (cores > 1)

### parallel
if(parallel) {
  registerDoMC(cores = cores)
}

### data
load("data/snps.throm.gait2.RData") # -> snps,  gf

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


K <- 0.02
dat <- mutate(dat, offset = -qnorm(1 - K))

### Pass 1: one snp
model0 <- relmatGlmer(VT ~ AGEsc + AGEsc2 + SEXfnum + (1|ID), dat, offset = offset, relmat = list(ID = dkin2), family = binomial(probit))

models.snps1 <- dlply(snpf, "Marker", function(x) {
  cat(" * Marker", x$Marker, "\n")
  
  GAIT2.passed <- x$GAIT2
  
  maf.thr <- 0.05
  maf.passed <- ifelse(GAIT2.passed, x$maf.imputed > maf.thr, FALSE)
  
  model <- NULL
  if(maf.passed) {
    model <- try(update(model0, paste0(". ~ . + ", x$Marker)))
  }
  
  anova <- NULL
  if(!is.null(model)) {
    anova <- try(anova(model0, model))
  }
    
  list(Marker = x$Marker, maf.imputed = x$maf.imputed, 
    maf.thr = maf.thr, 
    GAIT2.passed = GAIT2.passed, maf.passed = maf.passed,
    model = model, anova = anova)
}, .parallel = parallel)

### Pass 2: rs8176719 (ABO) + one snp
model0 <- relmatGlmer(VT ~ AGEsc + AGEsc2 + SEXfnum + rs8176719 + (1|ID), dat, offset = offset, relmat = list(ID = dkin2), family = binomial(probit))

models.snps2 <- dlply(snpf, "Marker", function(x) {
  cat(" * Marker", x$Marker, "\n")
  
  null.passed <- !(x$Marker %in% names(fixef(model0)))
  
  GAIT2.passed <- x$GAIT2
  
  maf.thr <- 0.05
  maf.passed <- ifelse(GAIT2.passed, x$maf.imputed > maf.thr, FALSE)
  
  model.passed <- null.passed & maf.passed
  
  model <- NULL
  if(model.passed) {
    model <- try(update(model0, paste0(". ~ . + ", x$Marker)))
  }
  
  anova <- NULL
  if(!is.null(model)) {
    anova <- try(anova(model0, model))
  }
    
  list(Marker = x$Marker, maf.imputed = x$maf.imputed, 
    maf.thr = maf.thr, 
    GAIT2.passed = GAIT2.passed, maf.passed = maf.passed, null.passed = null.passed, model.passed = model.passed,
    model = model, anova = anova)
}, .parallel = parallel)

### Pass 3: rs8176719 (ABO) + rs2036914 (F11) + one snp
model0 <- relmatGlmer(VT ~ AGEsc + AGEsc2 + SEXfnum + rs8176719 + rs2036914 + (1|ID), dat, offset = offset, relmat = list(ID = dkin2), family = binomial(probit))

#m <- relmatGlmer(VT ~ AGEsc + AGEsc2 + SEXfnum + rs8176719 + rs2036914 + rs2288904 + (1|ID), dat, relmat = list(ID = dkin2), family = binomial(probit)(probit))
# m <- relmatGlmer(VT ~ AGEf + SEXfnum + rs8176719 + rs2036914 + rs2288904 + (1|ID), dat, relmat = list(ID = dkin2), family = binomial(probit)) # not converged
#m <- relmatGlmer(VT ~ AGEsc + AGEsc2 + SEXfnum + rs8176719 + rs2036914 + rs2288904 + (1|ID), dat, relmat = list(ID = dkin2), family = binomial(probit)) # not converged

models.snps3 <- dlply(snpf, "Marker", function(x) {
  cat(" * Marker", x$Marker, "\n")
  
  null.passed <- !(x$Marker %in% names(fixef(model0)))
  
  GAIT2.passed <- x$GAIT2
  
  maf.thr <- 0.05
  maf.passed <- ifelse(GAIT2.passed, x$maf.imputed > maf.thr, FALSE)
  
  model.passed <- null.passed & maf.passed
  
  model <- NULL
  if(model.passed) {
    model <- try(update(model0, paste0(". ~ . + ", x$Marker)))
  }
  
  anova <- NULL
  if(!is.null(model)) {
    anova <- try(anova(model0, model))
  }
    
  list(Marker = x$Marker, maf.imputed = x$maf.imputed, 
    maf.thr = maf.thr, 
    GAIT2.passed = GAIT2.passed, maf.passed = maf.passed, null.passed = null.passed, model.passed = model.passed,
    model = model, anova = anova)
}, .parallel = parallel)


### post-process
process_models <- function(models) 
{
  llply(models, function(x) {
    pval <- NA
    if("anova" %in% names(x)) {
      pval <- as.data.frame(x$anova)[2, "Pr(>Chisq)"]
    }
    if(is.null(pval)) {
      pval <- NA
    }
    
    OR <- NA
    OR <- try(as.numeric(fixef(x$model)[x$Marker]), silent = TRUE)
    if(class(OR) == "try-error") {
      OR <- NA
    } 
    
    c(x, list(pval = pval, OR = OR))
  })  
}

# snps1
models.snps1 <- process_models(models.snps1)

tab.snps1 <- ldply(models.snps1, function(x) 
  data.frame(Marker = x$Marker, pval1 = x$pval, OR1 = x$OR, stringsAsFactors = FALSE))
tab.snps1 <- tab.snps1[, -1]

ord <- order(tab.snps1$pval1)
tab.snps1 <- tab.snps1[ord, ]

# snps2
models.snps2 <- process_models(models.snps2)

tab.snps2 <- ldply(models.snps2, function(x) 
  data.frame(Marker = x$Marker, pval2 = x$pval, OR2 = x$OR, stringsAsFactors = FALSE))
tab.snps2 <- tab.snps2[, -1]

ord <- order(tab.snps2$pval2)
tab.snps2 <- tab.snps2[ord, ]

# snps3
models.snps3 <- process_models(models.snps3)

tab.snps3 <- ldply(models.snps3, function(x) 
  data.frame(Marker = x$Marker, pval3 = x$pval, OR3 = x$OR, stringsAsFactors = FALSE))
tab.snps3 <- tab.snps3[, -1]

ord <- order(tab.snps3$pval3)
tab.snps3 <- tab.snps3[ord, ]

### save
save(models.snps1, tab.snps1, models.snps2, tab.snps2, models.snps3, tab.snps3, file = "models.snps.run3.RData")

### rs8176719 (ABO)
#m <- update(m0, . ~ . + rs8176719)

### rs1801020 (f12, RR 5)
#m <- update(m0, . ~ . + rs1801020)

### two SNPs
#m <- update(m0, . ~ . + rs8176719 + rs1801020)

### rs1799963 (f2, RR 4), MAF ~ 2%
# failed
#m <- update(m0, . ~ . + rs1799963)


