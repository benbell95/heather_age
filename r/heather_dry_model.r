##########################################################################
### Function to simulate age distribution of plants on dry moorland, with change of mortality over time + optional managed burns
### Copyright 2023 - 2025, Benjamin Bell.
### Code and updates: https://github.com/benbell95/heather_age
### Blog: https://www.benjaminbell.co.uk
### Paper: ... Ritson et al. 2025

### This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

### This function has been tested with R v4.4.3. Use at your own risk.

### See README.md file for usage. 

### Mortality data for heather (as used in the paper) is provided in file "exp_mortality.csv".
### Additional functions available for analysing and plotting model output.

### Additional version of the model, which uses a random mortality chance also available for experimentation.

### Arguments (options)
# plants = How many plants
# m = Mortality % - should be simple numeric vector, e.g. (1:10, for 1 % chance at plant age 1, and 10% chance at plant age 10)
# r = Reset age of plant to this value on death (can be negative value for lag time)
# years = How many years to simulate
# sa = Start age of plants - the age of the plants for the first run (default 1). This can also be a vector of ages, for example, to start the scenario with a community of plants with a typical age distribution.
# managed (optional) = if plants are managed, set to true, and set following arguments:
# ms = Management start year (i.e. start at year 10)
# me = Management end year. Optional, leave blank for management to carry on indefinitely.
# mf = Management frequency (How often does management occur)
# mk = Percentage of plants to kill (reset age)
# write = logical - whether to write the output to a csv
# fn = file name if write=TRUE - include system path

heather_dry <- function(plants=1000, m, r=0, years=100, sa=1, managed=FALSE, ms=10, me, mf=10, mk=10, write=FALSE, fn) {
    # Check max mortality value - if lower than 100, repeat max value so length of m = years + 1 (otherwise generates NA values)
    # This does NOT increase mortality chance - change input values if want to do this
    if(max(m) < 100) {
        m <- c(m, rep(max(m), times=(years + 1) - length(m)))
    }
    # Convert % to decimals (for random number since using runif() for probability)
    m <- m / 100
    # Create matrix (with additional row [year 1]) to hold results
    mat <- matrix(NA, ncol=plants, nrow=years+1)
    # First row needs to be start age (e.g all plants start age 1 - default)
    mat[1,] <- sa
    ### Management
    if(managed==TRUE) {
        # Create a sequence for when to do management
        m_s1 <- seq(ms, years, mf)
        # sequence needs to match length of years
        # fill gaps with 0s
        m_s <- rep(0, times=years)
        m_s[m_s1] <- m_s1
        # Convert mk (kill percentage) into number of groups
        mkg <- plants / (mk / 100 * plants)
        # Create groups for plants to kill (reset age)
        m_g <- split(1:plants, cut(1:plants, mkg, labels=FALSE))
        # sequence to select group
        m_sgr <- rep(as.numeric(names(m_g)), times=ceiling(length(m_s1) / length(m_g))) |> {\(x) x[1:length(m_s1)]}()
        m_sg <- m_s 
        m_sg[m_s1] <- m_sgr
        # Management end
        if(hasArg(me)) {
            m_s[me:years] <- 0
        }
    }
    ### Loop
    for(i in 1:years) {
        ### Management (must occur first)
        if(managed==TRUE) {
            if(m_s[i]==i) {
                # Burn % of plants (reset age) - select row and columns 
                mat[i, c(m_g[[m_sg[i]]]) ] <- r              
            }   
        }
        # Increment age / mortality test
        # Additional check if age is negative (if so, always increment)
        mat[i + 1,] <- sapply(mat[i,], function(x) ifelse(ifelse(x < 1, TRUE, runif(1) >= m[x]), x + 1, r))
    }
    # Remove final row from matrix (extra year is run)
    mat <- mat[-nrow(mat),]
    # Write data
    if(write==TRUE) {
        write.csv(mat, fn)
    }
    # Return matrix
    return(mat)
}
