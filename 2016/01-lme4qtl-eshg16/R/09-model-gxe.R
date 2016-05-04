m0 <- relmatLmer(tr_mtDNA ~ SEXf + (1|HHID) + (1|ID), phen, relmat = list(ID = dkin))

m1 <- relmatLmer(tr_mtDNA ~ SEXf + (1|ID) + (0 + SEXf|RID) + (1|HHID), phen, 
  relmat = list(ID = dkin), 
  weights = rep(1e10, nrow(phen)), vcControl = list(rho0 = list(rid = 3)))
  
m2 <- relmatLmer(tr_mtDNA ~ SEXf + (0 + SEXf|ID) + (0 + SEXf|RID) + (1|HHID), phen, 
  relmat = list(ID = dkin), 
  weights = rep(1e10, nrow(phen)), vcControl = list(rho0 = list(rid = 5)))
  
### APTT
  
m0 <- relmatLmer(APTT ~ AGEsc + SEXf + (1|HHID) + (1|ID), phen, relmat = list(ID = dkin))

m1 <- relmatLmer(APTT ~ AGEsc + SEXf + (1|ID) + (0 + SEXf|RID) + (1|HHID), phen, 
  relmat = list(ID = dkin), 
  weights = rep(1e10, nrow(phen)), vcControl = list(rho0 = list(rid = 3)))
  
m2 <- relmatLmer(APTT ~ AGEsc + SEXf + (0 + SEXf|ID) + (0 + SEXf|RID) + (1|HHID), phen, 
  relmat = list(ID = dkin), 
  weights = rep(1e10, nrow(phen)), vcControl = list(rho0 = list(rid = 5)))
  
### FVIII
m0 <- relmatLmer(FVIII ~ AGEsc + AGEsc2 + SEXf + (1|HHID) + (1|ID), phen, relmat = list(ID = dkin))

m1 <- relmatLmer(FVIII ~ AGEsc + AGEsc2 + SEXf + (1|ID) + (0 + SEXf|RID) + (1|HHID), phen, 
  relmat = list(ID = dkin), 
  weights = rep(1e10, nrow(phen)), vcControl = list(rho0 = list(rid = 3)))
  
m2 <- relmatLmer(FVIII ~ AGEsc + AGEsc2 + SEXf + (0 + SEXf|ID) + (0 + SEXf|RID) + (1|HHID), phen, 
  relmat = list(ID = dkin), 
  weights = rep(1e10, nrow(phen)), vcControl = list(rho0 = list(rid = 5)))  
  w
### bmi (GAIT2)
phen2 <- mutate(phen2, BMI = log(bmi))
m0 <- relmatLmer(BMI ~ AGEsc + AGEsc2 + SEXf + (1|HHID) + (1|ID), phen2, relmat = list(ID = dkin2))

m1 <- relmatLmer(BMI ~ AGEsc + AGEsc2 + SEXf + (1|ID) + (0 + SEXf|RID) + (1|HHID), phen2, 
  relmat = list(ID = dkin2), 
  weights = rep(1e10, nrow(phen2)), vcControl = list(rho0 = list(rid = 3)))
  
m2 <- relmatLmer(BMI ~ AGEsc + AGEsc2 + SEXf + (0 + SEXf|ID) + (0 + SEXf|RID) + (1|HHID), phen2, 
  relmat = list(ID = dkin2), 
  weights = rep(1e10, nrow(phen2)), vcControl = list(rho0 = list(rid = 5)))    
  
  
### FVIII
phen2 <- mutate(phen2, FVIII = log(FVIIIc))
m0 <- relmatLmer(FVIII ~ AGEsc + SEXf + (1|HHID) + (1|ID), phen2, relmat = list(ID = dkin2))

m1 <- relmatLmer(FVIII ~ AGEsc + SEXf + (1|ID) + (0 + SEXf|RID) + (1|HHID), phen2, 
  relmat = list(ID = dkin2), 
  weights = rep(1e10, nrow(phen2)), vcControl = list(rho0 = list(rid = 3)))
  
m2 <- relmatLmer(FVIII ~ AGEsc + SEXf + (0 + SEXf|ID) + (0 + SEXf|RID) + (1|HHID), phen2, 
  relmat = list(ID = dkin2), 
  weights = rep(1e10, nrow(phen2)), vcControl = list(rho0 = list(rid = 5)))    
  
