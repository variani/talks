### inc
library(plyr)

library(devtools)

load_all("~/git/ugcd/solarius")
load_all("~/git/ugcd/gait")

### GAIT1
phen <- gait1.phen()

phen <- mutate(phen,
  HHID = ifelse(!is.na(HHID), HHID, ID))

# `mtDNA`
f <- "~/git/salambo/projects/14-sex-specific/R/rdata/set.RData"
if(file.exists(f)) {
  load(f) # -> `set`
  dat <- set@sampleInfo
  dat$ID <- rownames(dat)

  phen <- merge(phen, subset(dat, select = c("mtDNA", "ID")), by = "ID")

  phen <- mutate(phen,
    tr_mtDNA = 5 * log(mtDNA + median(mtDNA, na.rm = TRUE)))
}

phen <- mutate(phen,
  SEXfnum = as.numeric(SEXf) - 1)

phen <- mutate(phen,
  AGEsc = AGEc / sd(phen$AGEc, na.rm = TRUE),
  AGEsc2 = AGEsc^2)

phen <- mutate(phen,
  tr1_FXII = log(FXII),
  tr_bmi = log(bmi),
  BMI = log(bmi))

phen <- mutate(phen,
  RID = ID)

phen <- mutate(phen,
  AGEf = as.factor(ifelse(is.na(AGE), NA, 
    ifelse(AGE < 55, 0,
    ifelse(AGE < 64, 1,
    ifelse(AGE < 74, 2,
    ifelse(AGE < 84, 3, 4)))))),
  AGEf2 = as.factor(ifelse(is.na(AGE), NA, 
    ifelse(AGE < 55, 0, 1))),
  AGEf3 = as.factor(ifelse(is.na(AGE), NA, 
    ifelse(AGE < 55, 0,
    ifelse(AGE < 74, 1, 2)))),
  AGEf4 = as.factor(ifelse(is.na(AGE), NA, 
    ifelse(AGE < 55, 0,
    ifelse(AGE < 64, 1,
    ifelse(AGE < 74, 2, 3)))))
)

phen <- mutate(phen,
  AGEf3num0 = as.numeric(AGEf3 == 0),
  AGEf3num1 = as.numeric(AGEf3 == 1),
  AGEf3num2 = as.numeric(AGEf3 == 2))

phen <- mutate(phen,
  AGEf4num0 = as.numeric(AGEf4 == 0),
  AGEf4num1 = as.numeric(AGEf4 == 1),
  AGEf4num2 = as.numeric(AGEf4 == 2),
  AGEf4num3 = as.numeric(AGEf4 == 3))
  
phen <- mutate(phen,
  ab0f2 = as.factor(ifelse(is.na(ab0), NA, 
    ifelse(ab0 < 0, 0, 1))),
  ab0f2num = as.numeric(ab0f2) - 1,
  c46tf2 = as.factor(ifelse(is.na(c46t), NA, 
    ifelse(c46t < 0, 0, 1))))
    
# kinship
dkin <- solarKinship2(phen)

### GAIT2
phen2 <- try(gait2.phen())

if(!class(phen2) == "try-error") {

phen2 <- mutate(phen2, 
  CVI3 = as.factor(ifelse(is.na(CVIlevel), NA,
    ifelse(CVIlevel %in% c(0), 0,
    ifelse(CVIlevel %in% c(1, 2, 3), 1, 
    ifelse(CVIlevel %in% c(4, 5), 2,
    -1))))))
    
phen2 <- mutate(phen2,
  tr_bmi = log(bmi),
  BMI = log(bmi))
  
phen2 <- mutate(phen2,
  RID = ID) 

phen2 <- mutate(phen2,
  SEXfnum = as.numeric(SEXf) - 1)

phen2 <- mutate(phen2,
  ABOf2 = factor(ifelse(is.na(ABO), NA, 
    ifelse(grepl("O", ABO), "O", "Not-O")), levels = c("Not-O", "O")),
  ABOf2num = as.numeric(ABOf2) - 1)
  
phen2 <- within(phen2, {
  ABOf3 <- factor(laply(ABO, function(x)
    ifelse(is.na(x), NA, length(regmatches(x, gregexpr("O", x))[[1]]))))
  ABOf3num <- as.numeric(ABOf3) - 1
})

phen2 <- within(phen2, {
  ABOf3A1 <- factor(laply(ABO, function(x)
    ifelse(is.na(x), NA, length(regmatches(x, gregexpr("A1", x))[[1]]))))
})

phen2 <- mutate(phen2,
  AGEf = as.factor(ifelse(is.na(AGE), NA, 
    ifelse(AGE < 55, 0,
    ifelse(AGE < 64, 1,
    ifelse(AGE < 74, 2,
    ifelse(AGE < 84, 3, 4)))))),
  AGEf2 = as.factor(ifelse(is.na(AGE), NA, 
    ifelse(AGE < 55, 0, 1))),
  AGEf3 = as.factor(ifelse(is.na(AGE), NA, 
    ifelse(AGE < 55, 0,
    ifelse(AGE < 74, 1, 2)))),
  AGEf4 = as.factor(ifelse(is.na(AGE), NA, 
    ifelse(AGE < 55, 0,
    ifelse(AGE < 64, 1,
    ifelse(AGE < 74, 2, 3)))))
)
      
phen2 <- mutate(phen2,
  AGEf3num0 = as.numeric(AGEf3 == 0),
  AGEf3num1 = as.numeric(AGEf3 == 1),
  AGEf3num2 = as.numeric(AGEf3 == 2))

phen2 <- mutate(phen2,
  AGEf4num0 = as.numeric(AGEf4 == 0),
  AGEf4num1 = as.numeric(AGEf4 == 1),
  AGEf4num2 = as.numeric(AGEf4 == 2),
  AGEf4num3 = as.numeric(AGEf4 == 3))
  

phen2 <- mutate(phen2,
  AGEfnum0 = as.numeric(AGEf == 0),
  AGEfnum1 = as.numeric(AGEf == 1),
  AGEfnum2 = as.numeric(AGEf == 2),
  AGEfnum3 = as.numeric(AGEf == 3),
  AGEfnum4 = as.numeric(AGEf == 4))


# kinship    
dkin2 <- solarKinship2(phen2)

}
### end of code for GAIT2
    
