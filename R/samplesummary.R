library(TCGAbiolinks)
library(DT)

tab <-  getSampleFilesSummary(project = c("TCGA-BRCA","TCGA-GBM","TCGA-OV","TCGA-LUAD",
                                          "TCGA-UCEC","TCGA-KIRC","TCGA-HNSC","TCGA-LGG",
                                          "TCGA-THCA","TCGA-LUSC","TCGA-PRAD","TCGA-SKCM",
                                          "TCGA-COAD","TCGA-STAD","TCGA-BLCA","TCGA-LIHC",
                                          "TCGA-CESC","TCGA-KIRP","TCGA-TGCT","TCGA-SARC",
                                          "TCGA-LAML","TCGA-PAAD","TCGA-ESCA","TCGA-PCPG",
                                          "TCGA-READ","TCGA-THYM","TCGA-KICH","TCGA-ACC",
                                          "TCGA-MESO","TCGA-UVM","TCGA-DLBC","TCGA-UCS",
                                          "TCGA-CHOL"))

save(tab, file = "samplesummary.RData")