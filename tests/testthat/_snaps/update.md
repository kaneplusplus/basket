# Update works.

    Code
      summary(update_p0(mh1, 0.25))
    Output
      
      -- The MEM Model Call ----------------------------------------------------------
      
      mem_mcmc(responses = vemu_wide$responders, size = vemu_wide$evaluable, 
          name = vemu_wide$baskets, p0 = c(0.15, 0.15, 0.15, 0.2, 0.15, 
              0.15), mcmc_iter = 10, mcmc_burnin = 10)
      
      -- The Basket Summary ----------------------------------------------------------
      
      The Null Response Rates (alternative is greater):
                     NSCLC CRC (vemu) CRC (vemu+cetu) Bile Duct ECD or LCH   ATC
      Null           0.250      0.250           0.250     0.250      0.250 0.250
      Posterior Prob 0.873      0.068           0.002     0.758      0.793 0.746
      
      Posterior Mean and Median Response Rates:
             NSCLC CRC (vemu) CRC (vemu+cetu) Bile Duct ECD or LCH   ATC
      Mean   0.359      0.113           0.053     0.319      0.330 0.311
      Median 0.349      0.096           0.042     0.332      0.322 0.308
      
      Highest Posterior Density Interval with Coverage Probability 0.95:
                  NSCLC CRC (vemu) CRC (vemu+cetu) Bile Duct ECD or LCH   ATC
      Lower Bound 0.177      0.000           0.000     0.088      0.150 0.084
      Upper Bound 0.562      0.268           0.135     0.516      0.523 0.512
      
      Posterior Effective Sample Size:
        NSCLC CRC (vemu) CRC (vemu+cetu) Bile Duct ECD or LCH    ATC
       22.412     21.181          41.655    17.932     23.163 17.608
      

