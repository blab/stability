ddg_1 <- read.delim("~/Desktop/fhcrc/stability/augur/data_sets/0_1HA0_sequences_ddg.txt", header=FALSE)
ddg_2 <- read.delim("~/Desktop/fhcrc/stability/augur/data_sets/0_2YP7_sequences_ddg.txt", header=FALSE)
ddg_1HA0 = ddg_1$V1
ddg_2YP7 = ddg_2$V1
hist(ddg_1HA0, breaks=100, xlim=c(-15,70),main="Distribution of all possible mutations on 1968 structure")
hist(ddg_2YP7, breaks=100, xlim=c(-15,70),main="Distribution of all possible mutations on 2011 structure")
