##########################################################################
### Function to simulate age distribution of plants on dry moorland, with change of mortality over time + optional managed burns + optional random component to mortality.
### Copyright 2023 - 2025, Benjamin Bell.
### Code and updates: https://github.com/benbell95/heather_age
### Blog: https://www.benjaminbell.co.uk

### This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

### This function has been tested with R v4.4.3. Use at your own risk.

### This is a modified version of the "heather_dry" model, with the following additional arguments which introduce an element of random chance at mortality. Refer to original model for details on other arguments.

# random = logical. Random chance for random mortality
# rc = % chance of random mortality occurring
# ri = logical. TRUE = Always increase mortality. FALSE = completely random mortality.
# rm = random multiplier. "random" = random, else numeric for set increase. (this only works when ri=TRUE)

### increased mortality is a multiplier - either random, or increase

heather_dry_ran <- function(plants=1000, m, r=0, years=100, sa=1, random=TRUE, rc=10, ri=TRUE, rm="random", managed=FALSE, ms=10, me, mf=10, mk=10) {
    # Check max mortality value - if lower than 100, repeat max value so length of m = years + 1 (otherwise generates NA values)
    # This does NOT increase mortality chance - change input values if want to do this
    if(max(m) < 100) {
        m <- c(m, rep(max(m), times=(years + 1) - length(m)))
    }
    # Convert % to decimals (for random number since using runif() for probability)
    m <- m / 100
    rc <- rc / 100
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
    # Copy original mortality (important, as otherwise loop overwrites original when random=TRUE)
    if(random==TRUE) {mr <- m}
    ### Loop
    for(i in 1:years) {
        ### Random chance of different mortality 
        if(random==TRUE) {
            if(runif(1) >= rc) {   
                m <- mr 
            } else {
                # Always increase mortality
                if(ri==TRUE && rm=="random") { m <- mr * (runif(1, 1, 10) / 10 + 1) }
                if(ri==TRUE && rm!="random") { m <- mr * rm }              
                # random mortality change
                if(ri==FALSE) {m <- mr * runif(1)}
            }
        }
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
    # Return matrix
    return(mat)
}
