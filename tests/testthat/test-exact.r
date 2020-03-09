
context("Fit exact models")

data(vemu_wide)

fb_reference <- readRDS("reference_data/mem-exact.rds")

baskets <- c(1, 2, 3)
# baskets <- 1:6

vemu_wide1 <- vemu_wide[baskets, ]

# Full Bayes
exact_res <- mem_exact(
  responses = vemu_wide1$responders,
  size = vemu_wide1$evaluable,
  name = vemu_wide1$baskets,
  cluster_analysis = TRUE,
  p0 = 0.25
)
# plot_mem(exact_res, text_size = 1.75, expand = c(0.5, 0.5), basket_name_hjust = 0.85)
# plot_pep(exact_res$basket)

#
expect_true(class(summary(exact_res)) == "mem_summary")
expect_true(class(print(summary(exact_res))) == "mem_summary")

expect_true(class(summary(exact_res$basket)) == "mem_basket_summary")
expect_true(class(print(summary(exact_res$basket))) == "mem_basket_summary")

expect_true(class(summary(exact_res$cluster)) == "mem_cluster_summary")
expect_true(class(print(summary(exact_res$cluster))) == "mem_cluster_summary")


fb <- exact_res$basket

expect_equal(fb$PEP, fb_reference$basketwise$PEP)
expect_equal(fb$HPD, fb_reference$basketwise$HPD, tolerance = 5e-2)
expect_equal(fb$post.prob, fb_reference$basketwise$post.prob)
expect_equal(fb$ESS, fb_reference$basketwise$ESS, tolerance = 5e-4)
# summary(exact_res)
# print(exact_res)

res1 <- update(exact_res, p0 = 0.15)
expect_true(inherits(res1, "mem_exact"))

expect_equal(basket_name(res1), c("NSCLC", "CRC (vemu)", "CRC (vemu+cetu)"))

expect_true(inherits(
  exact_res_up <- update_p0(exact_res, alternative = "less"),
  "exchangeability_model"
))

time_taken <- system.time({
  exact_res1 <- mem_exact(
    responses = c(0),
    size = c(1),
    name = "Basket 1",
    p0 = 0.15
  )
})
