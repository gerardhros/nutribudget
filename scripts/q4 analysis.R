# data analysis RQ4
# Gerard H. ros & Hongzheng
# 13 May 2024

# require packages
require(readxl)
require(data.table)
require(metafor)
require(ggplot2)

# read in the database
d1 <- as.data.table(read_xlsx('data/rq4_database.xlsx',sheet='RQ4_database'))

# make a linear model, likely for a subset for NRECO and RTREAT
m1 <- lm(yi~duration_days * sup_cat,data=d2)
 