### inc
library(devtools)
load_all("~/git/hemostat/lme4qtl")

### ids
#ids <- with(phen, ID[FAMID == "01"])[1:10]
ids <- c("01101", "01102",
  "01202", "01205", "01210", "01203", "01208",
  "01211", "01204", "01209")

### kinsip
ind <- rownames(dkin) %in% ids
mat <- dkin[ind, ind]

pf <- melt(mat) # Var1 Var2 value
pf <- mutate(pf,
  Var1 = factor(Var1, levels = as.numeric(ids)),  
  Var2 = factor(Var2, levels = rev(as.numeric(ids))))

pf <- join(pf, data.frame(Var1 = as.numeric(phen$ID), SEX1 = phen$SEX))
pf <- join(pf, data.frame(Var2 = as.numeric(phen$ID), SEX2 = phen$SEX))

pf <- mutate(pf,
  value = ifelse(SEX1 == 1 & SEX2 == 1, -value, 
    ifelse(SEX1 == 2 & SEX2 == 2, value, 0)))

### plot
p1 <- ggplot(pf, aes(Var1, Var2)) + geom_tile(aes(fill = value), color = "white") + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "grey80")

p1 <- p1 + guides(fill = "none") + theme_void()

### mibd
# "~/git/ugcd/gait/inst/extdata/gait1/mibd.5.200.gz"
f <- system.file("extdata/gait1/mibd.5.200.gz", package = "gait")
  
mf <- read_mibd_csv_gz(f)
mf <- mf[, 1:3]

mf <- subset(mf, ID1%in% ids & ID2 %in% ids)
mat <- mf2mat(mf)

pf <- melt(mat) # Var1 Var2 value
pf <- mutate(pf,
  Var1 = factor(Var1, levels = as.numeric(ids)),  
  Var2 = factor(Var2, levels = rev(as.numeric(ids))))

pf <- join(pf, data.frame(Var1 = as.numeric(phen$ID), SEX1 = phen$SEX))
pf <- join(pf, data.frame(Var2 = as.numeric(phen$ID), SEX2 = phen$SEX))

pf <- mutate(pf,
  value = ifelse(SEX1 == 1 & SEX2 == 1, -value, 
    ifelse(SEX1 == 2 & SEX2 == 2, value, 0)))

p2 <- ggplot(pf, aes(Var1, Var2)) + geom_tile(aes(fill = value), color = "white") + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "grey80")

p2 <- p2 + guides(fill = "none") + theme_void()
 

