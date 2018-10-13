library(tidyr)

vemu_wide <- tibble(
  baskets = c("NSCLC", "CRC (vemu)", "CRC (vemu+cetu)", "Bile Duct", 
              "ECD or LCH", "ATC"),
  enrolled = c(20, 10, 27, 8, 18, 7),
  evaluable = c(19, 10, 26, 8, 14, 7),
  responders = c(8, 0, 1, 1, 6, 2),
  one_or_fewer_prior_therapy = c(11, 1, 5, 2, 9, 5),
  two_prior_therapies = c(4, 2, 11, 1, 7, 1),
  three_or_more_therapies = c(5, 7, 11, 5, 2, 1))

vemu <- gather(vemu_wide, key = "prior_therapies", value = "count", 5:7)

dir.create("../data")
save(vemu_wide, file = "../data/vemu_wide.rda")
save(vemu, file = "../data/vemu.rda")

