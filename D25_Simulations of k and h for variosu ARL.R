# install.packages("spc")   # first-time only
library(spc)

# Values of h resulting in ARL_0 of around 370 observations
xcusum.arl(k = 0.25, h = 8.01, mu = 0, sided = "two")
xcusum.arl(k = 0.50, h = 4.77, mu = 0, sided = "two")
xcusum.arl(k = 0.75, h = 3.34, mu = 0, sided = "two")
xcusum.arl(k = 1.00, h = 2.52, mu = 0, sided = "two")
xcusum.arl(k = 1.25, h = 1.99, mu = 0, sided = "two")
xcusum.arl(k = 1.50, h = 1.61, mu = 0, sided = "two")

# Values of h resulting in ARL_0 of around 40 observations
xcusum.arl(k = 0.25, h = 4.1, mu = 0, sided = "two")
xcusum.arl(k = 0.50, h = 2.7, mu = 0, sided = "two")
xcusum.arl(k = 0.75, h = 1.9, mu = 0, sided = "two")
xcusum.arl(k = 1.00, h = 1.45, mu = 0, sided = "two")
xcusum.arl(k = 1.25, h = 1.1, mu = 0, sided = "two")
xcusum.arl(k = 1.50, h = 0.8, mu = 0, sided = "two")

# Values of h resulting in ARL_0 of around 20 observations
xcusum.arl(k = 0.25, h = 3.1, mu = 0, sided = "two")
xcusum.arl(k = 0.50, h = 2.1, mu = 0, sided = "two")
xcusum.arl(k = 0.75, h = 1.5, mu = 0, sided = "two")
xcusum.arl(k = 1.00, h = 1.1, mu = 0, sided = "two")
xcusum.arl(k = 1.25, h = 0.75, mu = 0, sided = "two")
xcusum.arl(k = 1.50, h = 0.5, mu = 0, sided = "two")
