library(homologene)

updateHomologene(destfile = "../../data/Homologene_240528")

# to convert the name with the new ids use the following
homologene(human_gene,inTax = "9606",outTax = "10090",db = "../../data/Homologene_240528")

# read in the updated reference table
homologeneData2_240528 <- read_tsv("../../data/Homologene_240528") %>%
  as.data.frame()

# once updated pick the annotation of interest
homologene(human_gene, inTax = 9606, outTax = 10090)

homologene(human_gene, inTax = 9606, outTax = 10090,db = homologeneData2)

homologene(human_gene, inTax = 9606, outTax = 10090,db = homologeneData2_240528)