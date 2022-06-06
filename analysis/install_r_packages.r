if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager",repos='http://cran.us.r-project.org')

#BiocManager::install(version = '3.14') ## for R version 4.1
BiocManager::install(version = '3.15') ## for R version 4.2
BiocManager::install("qvalue")
install.packages("MendelianRandomization")
install.packages("metafor")
install.packages("metacor")