### inc
library(plyr)

load_all("~/git/ugcd/gait")

### vars
snps <- c("rs121909548", "rs6025", "rs1801020", "rs1799963", "rs8176719", 
  "rs8176750", "rs2232698", "rs2066865", "rs2036914", "rs2289252", "rs2227589",
  "rs78707713", "rs4981021", "rs1884841", "rs867186",
  "rs4524", "rs169713", "rs2288904", "rs710446", "rs1613662", "rs1063856", 
  "rs5985", "rs1039084", "rs1042579", "rs118203906", "rs118203905")

RR <- c(10, 5, 5, 4, 4, 4, 3.30, 1.47, 1.35, 1.35, 1.30, 1.30, 1.29, 1.22, 1.22, 
  1.21, 1.20, 1.20, 1.19, 1.15, 1.15, 0.85, 0.78, 0.72, NA, NA)
  
snpf <- data.frame(Marker = factor(snps, levels = snps), RR = RR)  

### query NCBI  
ncbif <- NCBI_snp_query(snpf$Marker)
ncbif <- within(ncbif, {
  Gene[Marker == "rs8176719"] <- "ABO"
  Gene[Marker == "rs169713"] <- "HIPEV1"
})

snpf <- join(snpf, ncbif) # subset(ncbif, select = c("Marker", "Chromosome", "BP", "Gene", "MAF")))

### query NCDF files with genotypes
gmat <- gait2.snpsearch.impute.ncdf(snpf$Marker, unique(snpf$Chromosome), cores = 22, verbose = 2)

gf <- data.frame(ID = rownames(gmat), as.data.frame(gmat), stringsAsFactors = FALSE)

### update `snpf`
snpf <- within(snpf,
  GAIT2 <- laply(Marker, function(x) x %in% colnames(gf)))

### save
save(snpf, gf, file = "snps.throm.gait2.RData")

