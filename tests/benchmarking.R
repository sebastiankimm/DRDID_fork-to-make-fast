## benchmark old estimators vs. C++ implementation

library(DRDID)
# Load data in long format that comes in the DRDID package
data(nsw_long)
# Form the Lalonde sample with CPS comparison group
eval_lalonde_cps <- subset(nsw_long, nsw_long$treated == 0 | nsw_long$sample == 2)

# Implement improved locally efficient DR DID:
out <- drdid(yname = "re", tname = "year", idname = "id", dname = "experimental",
             xformla= ~ age + educ + black + married + nodegree + hisp + re74,
             data = eval_lalonde_cps, panel = TRUE, boot = TRUE, boot.type = "multiplier", nboot = 10)
summary(out)
