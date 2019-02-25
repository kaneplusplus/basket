library(basket)
library(igraph)
library(ggplot2)

font_import()
loadfonts(device = "win")

#roxygen2::roxygenise() 

#pack <- "basket"
#path <- find.package(pack)
#system(paste(shQuote(file.path(R.home("bin"), "R")),
#             "CMD", "Rd2pdf", shQuote(path)))

data(vemu_wide)

baskets <- 1:5

vemu_wide1 <- vemu_wide[baskets, ]

# Full Bayes
exact_res <- mem_full_bayes_exact(
  responses = vemu_wide1$responders,
  size = vemu_wide1$evaluable,
  name = vemu_wide1$baskets, p0 = rep(0.15, length(basket))
)



print(exact_res$call)

print(exact_res$basketwise$PEP)
print(exact_res$basketwise$MAP)
print(exact_res$basketwise$post.prob)
print(exact_res$basketwise$HPD)
print(exact_res$basketwise$ESS)
print(exact_res$basketwise$mean_est)
print(exact_res$basketwise$median_est)

print(exact_res$clusterwise$cluster)
print(exact_res$clusterwise$post.prob)
print(exact_res$clusterwise$HPD)
print(exact_res$clusterwise$ESS)
print(exact_res$clusterwise$mean_est)
print(exact_res$clusterwise$median_est)

summary(exact_res$basketwise)
summary(exact_res$clusterwise)


exact_resNew <- update_result(exact_res, 0.25)

summary(exact_res$basketwise)
summary(exact_resNew$basketwise)

summary(exact_res$clusterwise)
summary(exact_resNew$clusterwise)

plot_density(exact_res$basketwise)
plot_posterior_exchangeability(exact_res$basketwise)
plot_all_exchangeability(exact_res$basketwise, c("PRIOR", "MAP", "PEP"), 
                         text_size=3, basket_name_hoffset=-0.2)
plot_all_exchangeability(exact_res$basketwise, c("MAP", "PEP"),
                         text_size=3, basket_name_hoffset=-0.2)

# rm(list=ls())
## Vectors of Observed number of Responses (X) and Patients (N)
## Data <- dget("../inst/code-from-brian/MHAlgorithm/Vemu-data.txt")

baskets <- 1:6
vemu_wide2 <- vemu_wide[baskets, ]

MHResult1 <- mem_full_bayes_mcmc(
  responses = vemu_wide2$responders,
  size = vemu_wide2$evaluable,
  name = c("NSCLC ", "CRC.v ", "CRC.vc", "  BD  ", "ED.LH ", " ATC  "),
  p0 = c(0.15, 0.15, 0.15, 0.2, 0.15, 0.15), Initial = NA
)


print(MHResult1$call)

print(MHResult1$basketwise$PEP)
print(MHResult1$basketwise$MAP)
print(MHResult1$basketwise$post.prob)
print(MHResult1$basketwise$HPD)
print(MHResult1$basketwise$ESS)
print(MHResult1$basketwise$mean_est)
print(MHResult1$basketwise$median_est)


print(MHResult1$clusterwise$cluster)
print(MHResult1$clusterwise$post.prob)
print(MHResult1$clusterwise$HPD)
print(MHResult1$clusterwise$ESS)
print(MHResult1$clusterwise$mean_est)
print(MHResult1$clusterwise$median_est)

summary(MHResult1$basketwise)
summary(MHResult1$clusterwise)

MHResult1New <- update_result(MHResult1, 0.25)

summary(MHResult1$basketwise)
summary(MHResult1New$basketwise)

summary(MHResult1$clusterwise)
summary(MHResult1New$clusterwise)


baskets <- 1:6
vemu_wide2 <- vemu_wide[baskets, ]

MHResult3 <- mem_full_bayes_mcmc(
  responses = c(
    vemu_wide2$responders,
    vemu_wide2$responders
  ),
  size = c(
    vemu_wide2$evaluable,
    vemu_wide2$evaluable
  ),
  name = 1:12,
  p0 = 0.15, niter.MCMC = 4000, Initial = NA
)

print(MHResult3$basketwise$PEP)
print(MHResult3$basketwise$MAP)
print(MHResult3$basketwise$post.prob)
print(MHResult3$basketwise$HPD)
print(MHResult3$basketwise$ESS)
print(MHResult3$basketwise$mean_est)
print(MHResult3$basketwise$median_est)


print(MHResult3$clusterwise$cluster)
print(MHResult3$clusterwise$post.prob)
print(MHResult3$clusterwise$HPD)
print(MHResult3$clusterwise$ESS)
print(MHResult3$clusterwise$mean_est)
print(MHResult3$clusterwise$median_est)
