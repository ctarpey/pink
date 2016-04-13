### Isolation by Distance calculations for populations
###    Distances in Kilometers, from ArcGIS
### Carolyn Tarpey | March 2016
### ---------------------------------------

setwd('G:/Analysis/Pop_analysis/Populations_b3_may/IBD')

#making matrices of 16681 Fst
odd16881 <- read.table ("16681_fst_MATRIX_odd.txt", sep = "\t")
odd16881 <- as.matrix(odd16881)  
even16881 <- read.table ("16681_fst_MATRIX_even.txt", sep = "\t")
even16881 <- as.matrix(even16881)  
all16881 <- read.table ("16681_fst_MATRIX.txt", sep = "\t")
all16881 <- as.matrix(all16881)  

#making matrices of 15996 Fst
odd15996 <- read.table ("15996_fst_MATRIX_odd.txt", sep = "\t")
odd15996 <- as.matrix(odd15996)  
even15996 <- read.table ("15996_fst_MATRIX_even.txt", sep = "\t")
even15996 <- as.matrix(even15996)  
all15996 <- read.table ("15996_fst_MATRIX.txt", sep = "\t")
all15996 <- as.matrix(all15996)  

#making geography matrices
####crosses the ocean
allgeo_ocean <- read.table ("geography_MATRIX_ocean.txt", sep = "\t")
allgeo_ocean  <- as.matrix(allgeo_ocean)  
oddgeo_ocean  <- read.table ("geography_MATRIX_odd_ocean.txt", sep = "\t")
oddgeo_ocean  <- as.matrix(oddgeo_ocean)  
evengeo_ocean  <- read.table ("geography_MATRIX_odd_ocean.txt", sep = "\t")
evengeo_ocean  <- as.matrix(evengeo_ocean)  

####coastline
allgeo_coast <- read.table ("geography_MATRIX_coast.txt", sep = "\t")
allgeo_coast  <- as.matrix(allgeo_coast)  
oddgeo_coast  <- read.table ("geography_MATRIX_odd_coast.txt", sep = "\t")
oddgeo_coast  <- as.matrix(oddgeo_coast)  
evengeo_coast  <- read.table ("geography_MATRIX_odd_coast.txt", sep = "\t")
evengeo_coast  <- as.matrix(evengeo_coast)  

####waterways
allgeo_water <- read.table ("geography_MATRIX_water.txt", sep = "\t")
allgeo_water  <- as.matrix(allgeo_water)  
oddgeo_water  <- read.table ("geography_MATRIX_odd_water.txt", sep = "\t")
oddgeo_water  <- as.matrix(oddgeo_water)  
evengeo_water  <- read.table ("geography_MATRIX_odd_water.txt", sep = "\t")
evengeo_water  <- as.matrix(evengeo_water)

####mantel tests of the fst matrices against the geographic distance
####crosses the ocean
odd16881Mantel <- mantel.test(odd16881, oddgeo_ocean, nperm = 999 )
even16881Mantel <- mantel.test(even16881, evengeo_ocean, nperm = 999)
all16881Mantel <- mantel.test(all16881, allgeo_ocean, nperm = 999)

odd15996Mantel <- mantel.test(odd15996, oddgeo_ocean, nperm = 999 )
even15996Mantel <- mantel.test(even15996, evengeo_ocean, nperm = 999)
all15996Mantel <- mantel.test(all15996, allgeo_ocean, nperm = 999)

# Results
odd16881Mantel 
even16881Mantel
all16881Mantel

odd15996Mantel
even15996Mantel
all15996Mantel

####coastline
odd16881Mantel_coast <- mantel.test(odd16881, oddgeo_coast, nperm = 999 )
even16881Mantel_coast <- mantel.test(even16881, evengeo_coast, nperm = 999)
all16881Mantel_coast <- mantel.test(all16881, allgeo_coast, nperm = 999)

odd15996Mantel_coast <- mantel.test(odd15996, oddgeo_coast, nperm = 999 )
even15996Mantel_coast <- mantel.test(even15996, evengeo_coast, nperm = 999)
all15996Mantel_coast <- mantel.test(all15996, allgeo_coast, nperm = 999)

# Results
odd16881Mantel_coast 
even16881Mantel_coast
all16881Mantel_coast

odd15996Mantel_coast
even15996Mantel_coast
all15996Mantel_coast

####WATERWAy
odd16881Mantel_water <- mantel.test(odd16881, oddgeo_water, nperm = 999 )
even16881Mantel_water <- mantel.test(even16881, evengeo_water, nperm = 999)
all16881Mantel_water <- mantel.test(all16881, allgeo_water, nperm = 999)

odd15996Mantel_water <- mantel.test(odd15996, oddgeo_water, nperm = 999 )
even15996Mantel_water <- mantel.test(even15996, evengeo_water, nperm = 999)
all15996Mantel_water <- mantel.test(all15996, allgeo_water, nperm = 999)

# Results
odd16881Mantel_water
even16881Mantel_water
all16881Mantel_water

odd15996Mantel_water
even15996Mantel_water
all15996Mantel_water
