##########################################################################
### Helper functions to use with heather_dry model for stats, plots and analysis
### Copyright 2023 - 2025, Benjamin Bell.
### Code and updates: https://github.com/benbell95/heather_age
### Blog: https://www.benjaminbell.co.uk
### Paper: ... Ritson et al. 2025


### These functions have been tested with R v4.4.3. Use at your own risk.

### See README.md file for usage and examples. 


### Analysis
# Function to reformat data for ANOVA/non-parametric testing (Wilcoxon rank sum test)
simdf <- function(x, stack=TRUE) {
    # Get no. years
    y <- nrow(x)
    # Transpose
    x <- transpose(x)
    # rename columns
    setnames(x, as.character(1:y))
    # convert values below 1 to NA (dead plants)
    for(col in names(x)) set(x, i=which(x[[col]] < 1), j=col, value=NA)
    # Stack
    if(stack==TRUE) {
        x <- melt(x, measure.vars=1:y, variable.name="year", value.name="plant_age")
    }
    return(x)
}

### Function to calculate adjusted p value using Wilcox test on a stacked data.table
# cny = Distribution is compared to how many "final" years (default = 10, i.e. last 10 years)
# pa = P adjustment method, see ?p.adjust for options
# future = use parallel processing or not
# write = logical. write data, or return as object
# fn = to write data.
# Data is technically "independent" - since the values are only dependent on the previous year
# If comparing year to year, use a paired test

wctp <- function(x, cny=10, pair=FALSE, pa="BY", future=FALSE, write=FALSE, fn) {
    # number of years (get number of levels)
    l <- x[, nlevels(year)]
    # Years to compare distribution to
    cy <- (l-cny):l
    # Calculate p
    if(future==TRUE) {
        r <- future_lapply(1:l, \(y) wilcox.test(x[year %in% cy, plant_age], x[year==y, plant_age], paired=pair)$p.value)
    } else {
        # store
        r <- vector("numeric", l)
        # loop - get p values
        # loop is just as fast as lapply :)
        for(i in 1:l) {
            r[i] <- wilcox.test(x[year %in% cy, plant_age], x[year==i, plant_age], paired=pair)$p.value
        } 
    }
    # Adjust p
    r <- p.adjust(r, pa)
    # Write data?
    if(write==TRUE) {
        write.csv(r, fn)
    } else {
        return(r)
    }
}


### Function to find stable periods
# x should be vector of p values
# t = threshold for p value to define similarity (use high value since test is of variance, not similarity)
# sp = stable period length (i.e. number of years in a given period to be considered stable)
stable <- function(x, t=0.90, sp=10) {
    # Get values above threshold
    a <- ifelse(x > t, TRUE, FALSE)
    # When do periods occur
    a_cs <- cumsum(rle(a)$lengths)
    a_t <- as.numeric(a_cs[which(rle(a)$values==FALSE)] + 1)
    # No FALSE means all stable, thus starts at first year
    if(length(a_t)==0) {
        a_t <- 1
    }
    # Get lengths
    a_l <- as.numeric(rle(a)$lengths[which(rle(a)$values==TRUE)])
    # Results
    # Check lengths match
    if(length(a_l)!=length(a_t)) {
        a_t <- a_t[-length(a_t)]
    }
    mat <- data.frame(year=a_t, length=a_l)
    # Remove years where length is less than sp (stable period)
    # But check to see if any match first!
    if(length(which(a_l < sp))>0) {
        mat <- mat[-c(which(a_l < sp)),]
    }
    return(mat)
}


### Function for stats table
sim_stats <- function(m, m_fy, med, med_fy, d, sp, sp1, sp1len) {   
     # Mean age entire distribution
     m_mean <- mean(m, na.rm=TRUE)
     m_sd <- sd(m, na.rm=TRUE)
     m_r <- round(range(m), 2)
     # Mean age final year
     m_fy_mean <- mean(m_fy, na.rm=TRUE)
     m_fy_sd <- sd(m_fy, na.rm=TRUE)
     m_fy_r <- round(range(m_fy), 2)
     # Median age entire distribution
     med_mean <- mean(med, na.rm=TRUE)
     med_sd <- sd(med, na.rm=TRUE)
     med_r <- round(range(med), 2)
     # Median age final year
     med_fy_mean <- mean(med_fy, na.rm=TRUE)
     med_fy_sd <- sd(med_fy, na.rm=TRUE)
     med_fy_r <- round(range(med_fy), 2)
     # Dead plants
     d_mean <- mean(d, na.rm=TRUE)
     d_sd <- sd(d, na.rm=TRUE)
     d_r <- round(range(d), 2)
     # Stable periods
     sp_mean <- mean(sp, na.rm=TRUE)
     sp_sd <- sd(sp, na.rm=TRUE)
     sp_r <- round(range(sp, na.rm=TRUE), 0) 
     # First stable period
     sp1_mean <- mean(sp1, na.rm=TRUE)
     sp1_sd <- sd(sp1, na.rm=TRUE)
     sp1_r <- round(range(sp1, na.rm=TRUE), 0) 
     # First stable period length
     sp1len_mean <- mean(sp1len, na.rm=TRUE)
     sp1len_sd <- sd(sp1len, na.rm=TRUE)
     sp1len_r <- round(range(sp1len, na.rm=TRUE), 0) 
     ### Results table
     # Mean, SD, Range (min/max)
     result <- data.frame(mean=round(c(m_mean, m_fy_mean, med_mean, med_fy_mean, d_mean, sp_mean, sp1_mean, sp1len_mean), 2), sd=round(c(m_sd, m_fy_sd, med_sd, med_fy_sd, d_sd, sp_sd, sp1_sd, sp1len_sd), 5), mean_min=c(m_r[1], m_fy_r[1], med_r[1], med_fy_r[1], d_r[1], sp_r[1], sp1_r[1], sp1len_r[1]), mean_max=c(m_r[2], m_fy_r[2], med_r[2], med_fy_r[2], d_r[2], sp_r[2], sp1_r[2], sp1len_r[2]))
     rownames(result) <- c("Mean age (all years)", "Mean age (final year)", "Median age (all years)", "Median age (final year)", "Dead Plants", "Stable periods", "First stable period", "First stable period length")
     return(result)
}


#######################################
### Plot - compare boxplot and p values
# Works with stacked data.table
plot_stable <- function(ad, p, sp, t=0.9, main="") {
    # number of years (get number of levels)
    years <- ad[, nlevels(year)]
    # age range
    ar <- range(ad[, plant_age], na.rm=TRUE)
    # x labels
    xl <- seq(0, years, 10)
    # dead plants
    dp <- lapply(1:years, \(x)  sum(is.na(ad[year==x, plant_age]))) |> unlist()
    dpm <- mean(dp)
    # Mean plant age
    mp <- lapply(1:years, \(x)  mean(ad[year==x, plant_age], na.rm=TRUE)) |> unlist()
    mpm <- mean(mp)
    # Median age at final year (final year is first year)
    ym <- median(ad[year==years, plant_age], na.rm=TRUE)
    #########################################################
    ### Plot layout
    layout(matrix(c(4,4,2,3,1,1), ncol=1))
    par(oma=c(7,5,7,5), bg="white")
    #########################################################
    ### PLOT 1 (p values) [appears at bottom]
    # Dummy plot
    par(mar=c(1,3,1,3))
    plot(p, xlim=c(0, years+1), type="n", xaxs="i", axes=FALSE, ann=FALSE)
    # Stable periods
    #if(lengths(sims_sp[[66]])[1]!=0) {
    if(lengths(sp)[1]!=0) {
        rect(xleft=sp[,1], ybottom=-1, xright=sp[,1] + (sp[,2] - 1), ytop=ar[2], col=adjustcolor("green", 0.1), border=NA, xpd=NA)
        text(x=sp[,1]+1, y=t - 0.05, pos=4, labels=paste("Length:", sp[,2], "years"), col="green4", cex=2)
    }
    # Plot p values
    lines(p)
    # Axis
    axis(2, cex.axis=2, las=1)
    axis(1, at=1:years, labels=NA, cex.axis=2)
    axis(1, at=xl, lwd.ticks=1.5, tck=-0.02, cex.axis=2, mgp=c(3,2,0))
    axis(4, cex.axis=2, las=1)
    mtext("Year", side=1, line=4.5, cex=1.5)
    mtext("Adjusted p value", side=2, line=5, cex=1.5)
    mtext("Adjusted p value", side=4, line=5, cex=1.5)
    text(x=3, y=1, labels="D", cex=4, font=2, pos=1)
    # Threshold
    abline(h=t, col="green4", lwd=2, lty=2)
    text(x=years, y=t - 0.05, pos=2, labels=paste("Threshold:", t), col="green4", font=2, cex=2)
    #########################################################
    ### PLOT 2 (mean plant age) [appears in middle]
    par(mar=c(1,3,1,3))
    plot(x=mp, xlim=c(0, years+1), col="black", type="l", xaxs="i", lwd=2, axes=FALSE, ann=FALSE)
    # Mean
    abline(h=mpm, col=adjustcolor("grey50", 0.75), lwd=3, lty=2) 
    text(x=years-15, y=min(mp), labels=paste0("Mean plant age: ", round(mpm,1)), col="black", cex=2, xpd=NA)
    # Axis
    axis(2, cex.axis=2, las=1)
    axis(4, cex.axis=2, las=1)
    mtext("Mean plant age", side=2, line=5, cex=1.5)
    mtext("Mean plant age", side=4, line=5, cex=1.5)
    text(x=3, y=max(mp), labels="B", cex=4, font=2, pos=1)
    #########################################################
    ### PLOT 3 (dead plants) [appears in middle]
    par(mar=c(1,3,1,3))
    plot(x=dp, xlim=c(0, years+1), col="red", type="l", xaxs="i", lwd=2, axes=FALSE, ann=FALSE)
    # Mean
    abline(h=dpm, col=adjustcolor("red3", 0.75), lwd=3, lty=2) 
    text(x=years-15, y=min(dp), labels=paste0("Mean dead plants: ", round(dpm,1)), col="red3", cex=2, xpd=NA)
    # Axis
    axis(2, cex.axis=2, las=1)
    axis(4, cex.axis=2, las=1)
    mtext("Dead plants", side=2, line=5, cex=1.5)
    mtext("Dead plants", side=4, line=5, cex=1.5)
    text(x=3, y=max(dp), labels="C", cex=4, font=2, pos=1)
    #########################################################
    ### PLOT 4 (Age distribution boxplot) [appears at top]
    # Dummy plot
    par(mar=c(1,3,1,3))
    plot(1:years, xlim=c(0, years+1), ylim=ar, type="n", xaxs="i", axes=FALSE, ann=FALSE)
    # Boxplot
    # Add boxplot
    boxplot(ad[, plant_age ~ year], at=1:years, col=adjustcolor("grey50", 0.5), add=TRUE, axes=FALSE, ann=FALSE) 
    # Median final year
    abline(h=ym, col=adjustcolor("blue", 0.75), lwd=3, lty=2) 
    text(x=years-15, y=-0.5, labels=paste0("Median age at year ", years, ": ", ym), col="blue", cex=2, xpd=NA)
    # Axis
    axis(2, cex.axis=2, las=1)
    axis(3, at=1:years, labels=NA, cex.axis=2)
    axis(3, at=xl, lwd.ticks=1.5, tck=-0.02, cex.axis=2)
    axis(4, cex.axis=2, las=1) 
    mtext("Plant age", side=2, line=5, cex=1.5)
    mtext("Plant age", side=4, line=5, cex=1.5) 
    text(x=3, y=ar[2], labels="A", cex=4, font=2, pos=1)
    ### Plot title 
    title(main=main, outer=TRUE, line=3.5, cex.main=3)
}



### Generate stats
# e = Exclude values below x
# r = age mean/range after x number of years (default 50)
heather_stats <- function(x, e, r=50, conf=0.05) {
    if(!missing(e)) {
        x[which(x < e)] <- NA
    }
    # Number of cols
    ns <- ncol(x)
    # Mean
    x_me <- rowMeans(x, na.rm=TRUE)
    # Median
    x_med <- apply(x, 1, median, na.rm=TRUE)
    # Std Dev
    x_sd <- apply(x, 1, sd, na.rm=TRUE)
    # std error
    x_se <- x_sd / sqrt(ns)
    # t score
    x_ts <- qt(p=conf/2, df=ns - 1, lower.tail=FALSE)
    # Margin of error
    x_mar <- x_ts * x_se
    # Confidence
    c_up <- x_me + x_mar
    c_lo <- x_me - x_mar
    # Mean range after x years
    x_r <- range(x_me[-c(0:r)])
    # Combine
    xl <- list(mean_age_year=x_me, median_age_year=x_med, sd_age_year=x_sd, std_error=x_se, t_score=x_ts, margin_of_error=x_mar, upper_confidence=c_up, lower_confidence=c_lo, max_age=max(x, na.rm=TRUE), mean_age=mean(x, na.rm=TRUE), plants=ncol(x), mean_age_range_after=as.numeric(c(x_r, paste(r))), mean_age_after=mean(x[-c(0:r)], na.rm=TRUE))
    return(xl)
}

# Plot mean age/stats
# x = stats
heather_plot_mean <- function(x,  type="o", col="black", pch=1, lwd=2, t="", mean=TRUE, error=TRUE, conf=TRUE, conf_col="blue", show_stats=TRUE, ...) {
    # Main plot
    plot(x[[1]], type=type, col=col, pch=pch, lwd=lwd, xlab="Year", ylab="Mean Age", main="", ...)
    # x axis position
    xa <- 1:length(x[[1]])
    # Error bars
    if(error==TRUE) {
        arrows(x0=xa, y0=x[[1]] - x[[3]], x1=xa, y1=x[[1]] + x[[3]], code=3, angle=90, length=0.1)

    }
    # Confidence interval
    if(conf==TRUE) {
        polygon(x=c(xa, rev(xa)), y=c(x[[7]], rev(x[[8]])), col=adjustcolor(conf_col, 0.25), border=NA)
        points(xa, x[[1]], type=type, col=col, pch=pch, lwd=lwd)
    }
    # Get coordinates
    co <- par("usr")
    # Mean
    if(mean==TRUE) {
        rect(xleft=co[1], ybottom=x[[7]][1], xright=co[2], ytop=x[[7]][2], col=adjustcolor("blue4", 0.25), border=NA)
        abline(h=x[[5]], col="blue")
    }
    # Additional tick marks
    axis(1, at=seq(0, length(x[[1]]), 10), labels=FALSE, tcl=-0.35)
    # Labels
    if(show_stats==TRUE) {
        text(length(x[[1]]), co[3] + co[4] * 0.2, paste("Plants:", x[[11]], "\nYears:", length(x[[1]]), "\nMax age (entire population):", x[[9]], "\nMean age (entire population):", round(x[[10]], 1), "\nMean age range after", x[[12]][3], "years:", paste(round(x[[12]][1], 1)), "to", paste(round(x[[12]][2], 1)), "\nMean age after", x[[12]][3], "years:", paste(round(x[[13]], 1))), adj=1)
    }
    title(main=t)
}
