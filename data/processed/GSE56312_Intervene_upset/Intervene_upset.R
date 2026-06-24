#!/usr/bin/env Rscript
if (suppressMessages(!require("UpSetR"))) suppressMessages(install.packages("UpSetR", repos="http://cran.us.r-project.org"))
library("UpSetR")
pdf("GSE56312_Intervene_upset/Intervene_upset.pdf", width=14, height=8, onefile=FALSE, useDingbats=FALSE)
expressionInput <- c('GSE56312_mES_Sox2_Serum_fixed'=2107,'GSE56312_mES_Nanog_Serum_fixed'=14600,'GSE56312_mES_Nanog_Serum_fixed&GSE56312_mES_Sox2_Serum_fixed'=10935)
upset(fromExpression(expressionInput), nsets=2, nintersects=30, show.numbers="yes", main.bar.color="#ea5d4e", sets.bar.color="#317eab", empty.intersections=NULL, order.by = "freq", number.angles = 0, mainbar.y.label ="No. of Intersections", sets.x.label ="Set size")
invisible(dev.off())
