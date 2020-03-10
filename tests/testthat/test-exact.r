
context("Fit exact models")

data(vemu_wide)

fb_reference <- readRDS("reference_data/mem-exact.rds")

baskets <- c(1, 2, 3)

vemu_wide1 <- vemu_wide[baskets, ]

# Full Bayes
exact_res <- mem_exact(
  responses = vemu_wide1$responders,
  size = vemu_wide1$evaluable,
  name = vemu_wide1$baskets,
  cluster_analysis = TRUE,
  p0 = 0.25
)

exact_res_basket <- basket(
  responses = vemu_wide1$responders,
  size = vemu_wide1$evaluable,
  name = vemu_wide1$baskets,
  cluster_analysis = TRUE,
  p0 = 0.25,
  method = "exact")

exact_res2 <- mem_exact(
  responses = c(4, 3, 0),
  size = c(10, 3, 0),
  name = letters[1:3],
  cluster_analysis = TRUE,
  p0 = 0.25)

exact_res3 <- mem_exact(
  responses = c(4, 3),
  size = c(10, 3),
  name = letters[1:2],
  cluster_analysis = TRUE,
  p0 = 0.25)

exact_res3 <- mem_exact(
  responses = c(4, 3),
  size = c(10, 3),
  name = letters[1:2],
  cluster_analysis = FALSE,
  p0 = 0.25)

# Remove the call and check for equality.
expect_equal(exact_res$basket[-18], 
             exact_res_basket$basket[-18], 
             tolerance = 0.5)


expect_true(class(summary(exact_res)) == "mem_summary")
expect_true(class(print(summary(exact_res))) == "mem_summary")

expect_true(class(summary(exact_res$basket)) == "mem_basket_summary")
expect_true(class(print(summary(exact_res$basket))) == "mem_basket_summary")

expect_true(class(summary(exact_res$cluster)) == "mem_cluster_summary")
expect_true(class(print(summary(exact_res$cluster))) == "mem_cluster_summary")


fb <- exact_res$basket

expect_equal(fb$pep, fb_reference$basketwise$PEP)
expect_equal(fb$hpd, fb_reference$basketwise$HPD, tolerance = 5e-2)
expect_equal(fb$post_prob, fb_reference$basketwise$post.prob)
expect_equal(fb$ess, fb_reference$basketwise$ESS, tolerance = 5e-4)

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
