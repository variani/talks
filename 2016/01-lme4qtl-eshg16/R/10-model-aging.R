
### bmi
phen <- mutate(phen, BMI = log(bmi))
m0 <- relmatLmer(BMI ~ AGEsc + AGEsc2 + SEXf + (1|HHID) + (1|ID), phen, relmat = list(ID = dkin))

m1 <- relmatLmer(BMI ~ AGEsc + AGEsc2 + SEXf + (1|ID) + (1 + AGEsc|RID) + (1|HHID), phen, 
  relmat = list(ID = dkin), 
  weights = rep(1e10, nrow(phen)), vcControl = list(rho0 = list(rid = 3)))
  
m2 <- relmatLmer(BMI ~ AGEsc + AGEsc2 + SEXf + (1 + AGEsc|ID) + (1 + AGEsc|RID) + (1|HHID), phen, 
  relmat = list(ID = dkin), 
  weights = rep(1e10, nrow(phen)), vcControl = list(rho0 = list(rid = 5)))
 
  
### bmi (GAIT2)
phen2 <- mutate(phen2, BMI = log(bmi))
m0 <- relmatLmer(BMI ~ AGEsc + AGEsc2 + SEXf + (1|HHID) + (1|ID), phen2, relmat = list(ID = dkin2))

m1 <- relmatLmer(BMI ~ AGEsc + AGEsc2 + SEXf + (1|ID) + (1 + AGEsc|RID) + (1|HHID), phen2, 
  relmat = list(ID = dkin2), 
  weights = rep(1e10, nrow(phen2)), vcControl = list(rho0 = list(rid = 3)))
  
m2 <- relmatLmer(BMI ~ AGEsc + AGEsc2 + SEXf + (1 + AGEsc|ID) + (1 + AGEsc|RID) + (1|HHID), phen2, 
  relmat = list(ID = dkin2), 
  weights = rep(1e10, nrow(phen2)), vcControl = list(rho0 = list(rid = 5)))
 
