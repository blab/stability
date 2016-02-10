#MKTIIALSYILCLVFAQKLPGNDNSTATLCLGHHAVPNGTIVKTITNDQIEVTNATELVQSSSTGGICDSPHQILDGENCTLIDALLGDPQCDGFQNKKWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEFNNESFNWTGVTQNGTSSACKRRSDKSFFSRLNWLTQLKYKYPALNVTMPNNEKFDKLYIWGVHHPGTDSDQTSLYVQASGRVTVSTKRSQQTVIPNIGSRPWVRGVSSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGTCNSECITPNGSIPNDKPFQNVNRITYGACPRYVKQNTLKLATGMRNVPEKQTRGIFGAIAGFIENGWEGMVDGWYGFRHQNSEGTGQAADLKSTQAAINQINGKLNRLIGKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFERTKKQLRENAEDMGNGCFKIYHKCDNACIGSIRNGTYDHDVYRDEALNNRFQIKGVELKSGYKDWILWISFAISCFLLCVALLGFIMWACQKGNIRCNICI
#QKLPGNDNSTATLCLGHHAVPNGTIVKTITNDQIEVTNATELVQSSSTGGICDSPHQILDGENCTLIDALLGDPQCDGFQNKKWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEFNNESFNWTGVTQNGTSSACKRRSDKSFFSRLNWLTQLKYKYPALNVTMPNNEKFDKLYIWGVHHPGTDSDQTSLYVQASGRVTVSTKRSQQTVIPNIGSRPWVRGVSSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGTCNSECITPNGSIPNDKPFQNVNRITYGACPRYVKQNTLKLATGMRNVPEKQTRGIFGAIAGFIENGWEGMVDGWYGFRHQNSEGTGQAADLKSTQAAINQINGKLNRLIGKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFERTKKQLRENAEDMGNGCFKIYHKCDNACIGSIRNGTYDHDVYRDEALNNRFQIKGVELKSGYKDWILWISFAISCFLLCVALLGFIMWACQKGNIRCNICI
#signal peptide is 16 sites long
#HA1 starts at "QKL..."
#The sequences I was giving to run_stability cluster included the sp
#In virus_stability.py align_to_outgroup function only compare from site 24 "STAT..."
#So sites in mutator and stability results will be from python index 8, but structure index 9. "STAT"
#So will have to subtract 1 from site position to get python index
#Receptor Binding sites are in the same order as mutator and stability results
#R indexing starts at 1

ddg_1 <- read.delim("~/Desktop/fhcrc/stability/augur/data_sets/0_1HA0_sequences_ddg.txt", header=FALSE)
ddg_2 <- read.delim("~/Desktop/fhcrc/stability/augur/data_sets/0_2YP7_sequences_ddg.txt", header=FALSE)
ddg_1HA0 = ddg_1$V1
ddg_2YP7 = ddg_2$V1
mutation_1HA0 = as.character(ddg_1$V2)
mutation_2YP7 = as.character(ddg_2$V2)

ep = "0000000000000000000000000000000000000000000011111011011001010011000100000001001011110011100110101000001100000100000001000110101011111101011010111110001010011111000101011011111111010010001111101110111001010001110011111111000000111110000000101010101110000000000011100100000001011011100000000000001001011000110111111000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"
ep_ddg_1HA0 = list()
ep_ddg_2YP7 = list()
ep_mutation_1HA0 = list()
ep_mutation_2YP7 = list()

nep_ddg_1HA0 = list()
nep_ddg_2YP7 = list()
nep_mutation_1HA0 = list()
nep_mutation_2YP7 = list()

c1_ddg_1HA0 = list()
c1_ddg_2YP7 = list()
c1_mutation_1HA0 = list()
c1_mutation_2YP7 = list()

c2_ddg_1HA0 = list()
c2_ddg_2YP7 = list()
c2_mutation_1HA0 = list()
c2_mutation_2YP7 = list()

len = length(ddg_1HA0)
for(i in 1:len) {
  mutation = mutation_1HA0[i]
  site = as.numeric(substr(mutation, 2, nchar(mutation)-1))
  if (substr(ep,site,site) == "1") {
    ep_ddg_1HA0 = c(ep_ddg_1HA0, ddg_1HA0[i])
    ep_ddg_2YP7 = c(ep_ddg_2YP7, ddg_2YP7[i])
    ep_mutation_1HA0 = c(ep_mutation_1HA0, mutation_1HA0[i])
    ep_mutation_2YP7 = c(ep_mutation_2YP7, mutation_2YP7[i])
  } else if (substr(ep,site,site) =="0") {
    nep_ddg_1HA0 = c(nep_ddg_1HA0, ddg_1HA0[i])
    nep_ddg_2YP7 = c(nep_ddg_2YP7, ddg_2YP7[i])
    nep_mutation_1HA0 = c(nep_mutation_1HA0, mutation_1HA0[i])
    nep_mutation_2YP7 = c(nep_mutation_2YP7, mutation_2YP7[i])
  }
  if (site<345) {
    c1_ddg_1HA0 = c(c1_ddg_1HA0, ddg_1HA0[i])
    c1_ddg_2YP7 = c(c1_ddg_2YP7, ddg_2YP7[i])
    c1_mutation_1HA0 = c(c1_mutation_1HA0, mutation_1HA0[i])
    c1_mutation_2YP7 = c(c1_mutation_2YP7, mutation_2YP7[i])
  } else if (site>345) {
    c2_ddg_1HA0 = c(c2_ddg_1HA0, ddg_1HA0[i])
    c2_ddg_2YP7 = c(c2_ddg_2YP7, ddg_2YP7[i])
    c2_mutation_1HA0 = c(c2_mutation_1HA0, mutation_1HA0[i])
    c2_mutation_2YP7 = c(c2_mutation_2YP7, mutation_2YP7[i])
  }
}

old.par <- par(mfrow=c(3, 2))

hist(ddg_1HA0, breaks=100, xlim=c(-15,70), main="All possible mutations on 1968 structure")
hist(ddg_2YP7, breaks=100, xlim=c(-15,70), main="All possible mutations on 2011 structure")
hist(as.numeric(c1_ddg_1HA0), breaks=100, xlim=c(-15,70), main="All possible head mutations on 1968 structure")
hist(as.numeric(c1_ddg_2YP7), breaks=100, xlim=c(-15,70), main="All possible head mutations on 2011 structure")
hist(as.numeric(c2_ddg_1HA0), breaks=100, xlim=c(-15,70), main="All possible stalk mutations on 1968 structure")
hist(as.numeric(c2_ddg_2YP7), breaks=100, xlim=c(-15,70), main="All possible stalk mutations on 2011 structure")

hist(ddg_1HA0, breaks=100, xlim=c(-15,70), main="All possible mutations on 1968 structure")
hist(ddg_2YP7, breaks=100, xlim=c(-15,70), main="All possible mutations on 2011 structure")
hist(as.numeric(ep_ddg_1HA0), breaks=100, xlim=c(-15,70), main="All possible epitope mutations on 1968 structure")
hist(as.numeric(ep_ddg_2YP7), breaks=100, xlim=c(-15,70), main="All possible epitope mutations on 2011 structure")
hist(as.numeric(nep_ddg_1HA0), breaks=100, xlim=c(-15,70), main="All possible non-epitope mutations on 1968 structure")
hist(as.numeric(nep_ddg_2YP7), breaks=100, xlim=c(-15,70), main="All possible non-epitope mutations on 2011 structure")
par(old.par)