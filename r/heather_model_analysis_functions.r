##########################################################################
### Helper functions to use with heather_dry model for stats, plots and analysis
### Copyright 2023 - 2025, Benjamin Bell.
### Code and updates: https://github.com/benbell95/heather_age
### Blog: https://www.benjaminbell.co.uk
### Paper: ... Ritson et al. 2025

### These functions have been tested with R v4.5.1. Use at your own risk.

### See README.md file for usage and examples. 

### Analysis
# Function to reformat data as a stack
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


### Function to find stable periods
# x should be vector of p or d values
# pd = Specify whether p or d values ("p" / "d").
# t = threshold to define similarity. For "p" values, it looks at values greater than the threshold. For D values, it looks at values less than.
# sp = stable period length (i.e. number of years in a given period to be considered stable)
stable <- function(x, pd="d", t=0.90, sp=10) {
    # Get values against threshold
    if(pd=="p") a <- ifelse(x > t, TRUE, FALSE)
    if(pd=="d") a <- ifelse(x < t, TRUE, FALSE)
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
        # If first year = TRUE, it is not included in a_t, need to check and add it in
        if(a[1]==TRUE) {
            a_t <- c(1, a_t)
        } else {
            a_t <- a_t[-length(a_t)]
        }
    }
    mat <- data.frame(year=a_t, length=a_l)
    # Remove years where length is less than sp (stable period)
    # But check to see if any match first!
    if(length(which(a_l < sp))>0) {
        mat <- mat[-c(which(a_l < sp)),]
    }
    return(mat)
}


### Stable period plot
# d = d values (from ks test)
# sp = stable period (from stable function)
# sp_sd = stable period std dev, plotted as diagonal shading. (optional)
# m = management period (start, end) (optional)
# y = number of years
# xa = logical, plot x axis?
# xlabs = logical, plot x axis labels
# ya1 = logical, plot y axis on the left (side = 2)
# ya2 = logical, plot y axis on the right (side = 4)
# cols = colour for stable stage
# colm = colour for managment stage
# cole = colour for establishing stage
# ... = additional arguments to plot() e.g. use xlim to control x axis etc.

# Colorblind space colours = cols="#44AA99", colm="#d88a97", cole="#ddcc77"

plot_stable <- function(d, sp, sp_sd=0, m, y=300, xa=TRUE, xlabs=TRUE, ya1=TRUE, ya2=TRUE, cols="#44AA99", colm="#d88a97", cole="#ddcc77", ...) {
    # Empty plot
    plot(d, type="n", ylim=c(0, 1), xaxs="i", ann=FALSE, axes=FALSE, ...)
    pu <- par("usr")
    # Plot time periods
    if(hasArg(m)) { e <- m[2] } else { e <- 1}
    # Stable
    if(sp >= e) {
        rect(xleft=sp, ybottom=pu[3], xright=y, ytop=pu[4], col=cols, border=NA)
    }  
    if(sp > 2) {
        # Developing/establishing
        rect(xleft=e, ybottom=pu[3], xright=sp+sp_sd, ytop=pu[4], col=cole, border=NA)
        # stable/establish cross over
        rect(xleft=sp-sp_sd, ybottom=pu[3], xright=sp+sp_sd, ytop=pu[4], col=cols, density=5, lwd=10, border=NA, lend=2)
    }   
    # Management
    if(hasArg(m)) {
        rect(xleft=m[1], ybottom=pu[3], xright=m[2], ytop=pu[4], col=colm, border=NA)
        e <- m[2]
        #if(m[2] > sp) {
        #    rect(xleft=sp, ybottom=pu[3], xright=m[2], ytop=pu[4], col=cols, density=2, lwd=20, border=NA, lend=1, angle=90)
        #} 
    }
    # Plot range
    polygon(x=c(1:y, y:1), y=c(d[,5], rev(d[,4])), col="grey50", border=NA)
    # Plot mean
    lines(x=1:y, y=d[,2], lwd=3, col="black")
    # Axes
    if(xlabs==TRUE) {
        xlab1 <- 1
        xlabs <- seq(0, y, 50)
    } else {
        xlab1 <- NA
        xlabs <- NA
    }
    if(xa==TRUE) {
    axis(1, at=1, lab=xlab1)
    axis(1, at=seq(0, y, 50), lab=xlabs, lwd=0, lwd.ticks=1, tcl=-0.85)
    axis(1, at=seq(10, y-10, 10), lwd=0, lwd.ticks=1, lab=NA)
    }
    if(ya1==TRUE) { axis(2, las=1, lwd=0, lwd.ticks=1) }
    if(ya2==TRUE) { axis(4, las=1, lwd=0, lwd.ticks=1) }
    box()
}
