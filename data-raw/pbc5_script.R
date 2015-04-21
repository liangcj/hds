library(survival)
library(dplyr)

pbc5 <- pbc %>%
  slice(1:312) %>%
  select(time, status, age, edema, bili, albumin, protime) %>%
  mutate(status = (status==2)*1, bili = log(bili), protime = log(protime))
