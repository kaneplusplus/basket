# Full Bayes can be calculated correctly.

    Code
      unclass(summary(exact_res))
    Output
      $call
      mem_exact(responses = vemu_wide1$responders, size = vemu_wide1$evaluable, 
          name = vemu_wide1$baskets, p0 = 0.25, cluster_analysis = TRUE)
      
      $basket
      
      The Null Response Rates (alternative is greater):
                     NSCLC CRC (vemu) CRC (vemu+cetu)
      Null           0.250      0.250            0.25
      Posterior Prob 0.951      0.002            0.00
      
      Posterior Mean and Median Response Rates:
             NSCLC CRC (vemu) CRC (vemu+cetu)
      Mean   0.425      0.041           0.042
      Median 0.422      0.032           0.033
      
      Highest Posterior Density Interval with Coverage Probability 0.95:
            NSCLC CRC (vemu) CRC (vemu+cetu)
      lower 0.208      0.000           0.000
      upper 0.628      0.105           0.109
      
      Posterior Effective Sample Size:
       NSCLC CRC (vemu) CRC (vemu+cetu)
       19.93     55.647          50.867
      
      
      $cluster
      
      Cluster 1        
       "NSCLC"
      Cluster 2                               
       "CRC (vemu)" "CRC (vemu+cetu)"
      
      The Null Response Rates (alternative is greater):
                                 Cluster 1 Cluster 2
      Posterior for null of 0.25     0.951     0.001
      
      Posterior Mean and Median Response Rates:
             Cluster 1 Cluster 2
      Mean       0.425     0.041
      Median     0.422     0.033
      
      Highest Posterior Density Interval with Coverage Probability 0.95:
                  Cluster 1 Cluster 2
      Lower Bound     0.214     0.000
      Upper Bound     0.631     0.107
      
      Posterior Effective Sample Size:
       Cluster 1 Cluster 2
           20.28     52.91
      
      

# Exact corner case models

    Code
      summary(basket(responses = vemu_wide1$responders, size = vemu_wide1$evaluable,
      name = vemu_wide1$baskets, cluster_analysis = TRUE, p0 = 0.25, method = "exact"))
    Output
      
      -- The MEM Model Call ----------------------------------------------------------
      
      mem_exact(responses = responses, size = size, name = name, p0 = p0, 
          shape1 = shape1, shape2 = shape2, prior = prior, hpd_alpha = hpd_alpha, 
          alternative = alternative, seed = seed, cluster_analysis = cluster_analysis, 
          call = call, cluster_function = cluster_function)
      
      -- The Basket Summary ----------------------------------------------------------
      
      The Null Response Rates (alternative is greater):
                     NSCLC CRC (vemu) CRC (vemu+cetu)
      Null           0.250      0.250            0.25
      Posterior Prob 0.951      0.002            0.00
      
      Posterior Mean and Median Response Rates:
             NSCLC CRC (vemu) CRC (vemu+cetu)
      Mean   0.424      0.041           0.042
      Median 0.422      0.032           0.033
      
      Highest Posterior Density Interval with Coverage Probability 0.95:
            NSCLC CRC (vemu) CRC (vemu+cetu)
      lower 0.215      0.000           0.000
      upper 0.625      0.107           0.105
      
      Posterior Effective Sample Size:
        NSCLC CRC (vemu) CRC (vemu+cetu)
       20.933     52.153          55.786
      
      -- The Cluster Summary ---------------------------------------------------------
      
      Cluster 1        
       "NSCLC"
      Cluster 2                               
       "CRC (vemu)" "CRC (vemu+cetu)"
      
      The Null Response Rates (alternative is greater):
                                 Cluster 1 Cluster 2
      Posterior for null of 0.25     0.951     0.001
      
      Posterior Mean and Median Response Rates:
             Cluster 1 Cluster 2
      Mean       0.424     0.041
      Median     0.422     0.033
      
      Highest Posterior Density Interval with Coverage Probability 0.95:
                  Cluster 1 Cluster 2
      Lower Bound     0.219     0.000
      Upper Bound     0.637     0.107
      
      Posterior Effective Sample Size:
       Cluster 1 Cluster 2
          20.091    52.918
      

---

    Code
      summary(mem_exact(responses = c(4, 3, 0), size = c(10, 3, 0), name = letters[1:
      3], cluster_analysis = TRUE, p0 = 0.25))
    Output
      
      -- The MEM Model Call ----------------------------------------------------------
      
      mem_exact(responses = c(4, 3, 0), size = c(10, 3, 0), name = letters[1:3], 
          p0 = 0.25, cluster_analysis = TRUE)
      
      -- The Basket Summary ----------------------------------------------------------
      
      The Null Response Rates (alternative is greater):
                         a     b     c
      Null           0.250 0.250 0.250
      Posterior Prob 0.873 0.997 0.855
      
      Posterior Mean and Median Response Rates:
                 a     b     c
      Mean   0.420 0.848 0.590
      Median 0.415 0.919 0.574
      
      Highest Posterior Density Interval with Coverage Probability 0.95:
                a     b     c
      lower 0.146 0.473 0.067
      upper 0.693 1.000 1.000
      
      Posterior Effective Sample Size:
            a     b     c
       11.073 5.199 2.187
      
      -- The Cluster Summary ---------------------------------------------------------
      
      Cluster 1    
       "a"
      Cluster 2    
       "b"
      Cluster 3    
       "c"
      
      The Null Response Rates (alternative is greater):
                                 Cluster 1 Cluster 2 Cluster 3
      Posterior for null of 0.25     0.873     0.997     0.855
      
      Posterior Mean and Median Response Rates:
             Cluster 1 Cluster 2 Cluster 3
      Mean       0.420     0.848     0.590
      Median     0.415     0.919     0.574
      
      Highest Posterior Density Interval with Coverage Probability 0.95:
                  Cluster 1 Cluster 2 Cluster 3
      Lower Bound     0.147     0.472     0.063
      Upper Bound     0.699     1.000     1.000
      
      Posterior Effective Sample Size:
       Cluster 1 Cluster 2 Cluster 3
          10.817     5.171      2.12
      

---

    Code
      summary(mem_exact(responses = c(4, 3), size = c(10, 3), name = letters[1:2],
      cluster_analysis = TRUE, p0 = 0.25))
    Output
      
      -- The MEM Model Call ----------------------------------------------------------
      
      mem_exact(responses = c(4, 3), size = c(10, 3), name = letters[1:2], 
          p0 = 0.25, cluster_analysis = TRUE)
      
      -- The Basket Summary ----------------------------------------------------------
      
      The Null Response Rates (alternative is greater):
                         a     b
      Null           0.250 0.250
      Posterior Prob 0.873 0.997
      
      Posterior Mean and Median Response Rates:
                 a     b
      Mean   0.419 0.847
      Median 0.414 0.918
      
      Highest Posterior Density Interval with Coverage Probability 0.95:
                a     b
      lower 0.148 0.474
      upper 0.699 1.000
      
      Posterior Effective Sample Size:
            a     b
       10.816 5.251
      
      -- The Cluster Summary ---------------------------------------------------------
      
      Cluster 1    
       "a"
      Cluster 2    
       "b"
      
      The Null Response Rates (alternative is greater):
                                 Cluster 1 Cluster 2
      Posterior for null of 0.25     0.872     0.997
      
      Posterior Mean and Median Response Rates:
             Cluster 1 Cluster 2
      Mean       0.419     0.847
      Median     0.414     0.918
      
      Highest Posterior Density Interval with Coverage Probability 0.95:
                  Cluster 1 Cluster 2
      Lower Bound     0.148     0.471
      Upper Bound     0.699     1.000
      
      Posterior Effective Sample Size:
       Cluster 1 Cluster 2
          10.841     5.188
      

