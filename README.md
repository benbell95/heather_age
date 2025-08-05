# Heather Age Modelling

This script models heather age distribution for a dry heath/moorland environment, with or without management (e.g. managed burns). Heather plants age with time, with each plant having a chance at death, which increases with age. 

This model can be used to determine the age distribution of a heather plant community at different time periods, and determine how long the plant community takes to "stabilise". Please refer to the paper (...) for scientific background. 

There are two models available, heather_dry_model, and heather_dry_ran_model. The first is the main model (used in the paper), that uses a set mortality chance, which increases with age (100% mortality at plant age 44). The second model adds the ability to apply random mortality modifiers (e.g. random catastrophic events), and is available for experimentation.

Additionally, helper functions are available to generate stats, analysis and plots from the model run.

This guide provides example usage based on the heather_dry_model.

## Getting the code

You can load the model code directly into R from Github using the following:

```
# Load model code
source("https://raw.githubusercontent.com/benbell95/heather_age/r/heather_dry_model.r")
```

Heather plant mortality data is also available in the .csv file "exp_mortality.csv" which you can download and add to your project folder.

Alternatively, download the code for a local copy.

## Example

Example scenario 1: Simulate the age of 100 heather plants for 100 years, where there is a 5-year lag time after a plant dies until it regrows.

```
# Load mortality data (this should be a vector)
rmo <- read.csv("exp_mortality.csv")[,2]

# Simulate age of 100 heather plants for 100 years, with 5 year lag time.
sim1 <- heather_dry(plants=100, m=rmo, r=-4, years=100)
```

This will produce a 100 x 100 matrix of numeric data, representing the age of each "plant" (columns) for each year simulated (rows). Numbers below 1 represent "dead" plants. For a 5 year lag time, the reset argument (r) is set to -4 to account for 0.

Example scenario 2: Simulate the age of 100 heather plants for 100 years, where there is a 5-year lag time after a plant dies until it regrows, with managed burns, burning 10% of plants every 10 years starting at year 10.

```
# Simulate age of 100 heather plants for 100 years, with 5 year lag time, + managed burned (10% every 10 years, starting at year 10).
sim2 <- heather_dry(plants=100, m=rmo, r=-4, years=100, managed=TRUE, ms=10, mf=10, mk=10)
```

Details of each argument (and additional available arguments) are available in the code file.

The results matrix by itself is not that useful, it is the "raw" model output, but you can use this data to generate useful stats.

## Model analysis and plots

### Quick descriptive stats

Several "helper" functions are available which allow you to generate meaningful insights into the model output, and plot the results.

```
# Load model code
source("https://raw.githubusercontent.com/benbell95/heather_age/r/heather_model_analysis_functions.r")
```

For example, to generate some quick descriptive stats and plot the results, you can use the following code:

```
# Generate quick stats
sim1_stats <- heather_stats(sim1, e=1, r=50)
sim2_stats <- heather_stats(sim2, e=1, r=50)

# Plot results
heather_plot_mean(sim1_stats, t="sim 1", ylim=c(0,25))
heather_plot_mean(sim2_stats, t="sim 2", ylim=c(0,25))
```

Refer to the script for available options. After running heather_stats(), it will produce a list of descriptive statistics for each year simulated, as well as averages for the entire period.

For example, in sim1, the max and mean age of the heather plants is 41 and 11.6 respectively. In sim2, the max and mean ages are 36 and 10.5. (Your results will vary).

The plot shows mean age of the heather plant community and how this changes changes for each year simulated. 

### "Stability" period analysis

While the basic stats show changes over time, we also want to know when the age distribution reaches "stability" - that is, when there is no longer large variance in the plant age distribution from year to year.

For this, we use the Wilcoxon rank sum test, and compare the "final" age distribution (calculated as the mean of the final 10 years of the model run) to each year in the model run, looking at the p values for significant results. Since this is actually a test of variability, a significant p value (typically < 0.05) would suggest that the age distribution of the compared years (e.g. "final" years and year 1) varies and are not "stable". Therefore, we are actually interested in non-significant p values, and in this example, we set that threshold at < 0.9. Additionally, a "stable period" is defined when 10 or more consecutive years have p values that are not significant. 

To run the analysis, the raw model output data needs to be reformatted to work with the base r wilcox.test() function. 

The helper functions are designed to work with data.tables (since typically you might simulate thousands of plants or years), therefore it is necessary to convert the model output data first, before using the simdf() function to reformat the data.

```
# Load data.table library 
library("data.table")

# Convert to data.table
sim1 <- as.data.table(sim1)
sim2 <- as.data.table(sim2)

# Reformat data
sim1r <- simdf(sim1, stack=FALSE)
sim2r <- simdf(sim2, stack=FALSE)
```

Next, we'll run the Wilcoxon rank sum test to compare "final" year age distribution to each year, to check for variability, extract the p values, and adjust them. Adjusting the p values is important to reduce false positives since we are running many tests.

```
# Wilcoxon rank sum test
sim1_wt <- lapply(sim1r[,-1], function(x) wilcox.test(sim1r[[1]], x))
sim2_wt <- lapply(sim2r[,-1], function(x) wilcox.test(sim2r[[1]], x))

# Get p values
sim1_wt_p <- lapply(sim1_wt, function(x) x$p.value) |> unlist()
sim2_wt_p <- lapply(sim2_wt, function(x) x$p.value) |> unlist()

# Adjust p values
sim1_wt_pa <- p.adjust(sim1_wt_p, "BY")
sim2_wt_pa <- p.adjust(sim2_wt_p, "BY")
```

Next, we'll calculate the stable periods using the helper function.

```
# Stable period
sim1_sp <- stable(sim1_wt_pa)
sim2_sp <- stable(sim2_wt_pa)
```

Looking at the results, a single stable period was found in "sim1", which started at year 37, and lasted for 64 years. In "sim2", there were three stable periods, the first starting at year 9, lasting 26 years, the second in year 36 lasting 19 years, and the third starting in year 58 for 43 years. (Your results will vary)

You can then plot the simulation using the helper functions.

```
# Plot simulations/stable period
pdf("sim1_stable_plot.pdf", width=24, height=18)
plot_stable(ad=sim1r, p=sim1_wt_pa, sp=sim1_sp, t=0.9)
dev.off()

pdf("sim2_stable_plot.pdf", width=24, height=18)
plot_stable(ad=sim2r, p=sim2_wt_pa, sp=sim2_sp, t=0.9)
dev.off()
```

## Parallel processing

The two examples run very fast since they are only simulating a small number of plants and years.

However, you may wish to simulate thousands of plants over a longer time period. You might also want to run the simulation multiple times.

Whilst the function by itself does not use parallel processing, you can use the Future framework and wrap the function in future_lapply() to significantly speed up processing.

### Install, load and plan parallel processing in R

```
# Install package (Do this once only)
install.packages("future.apply")
# Load package and setup
library(future.apply)
plan(multisession)
```

### Example

Lets run the first simulation 1000 times.

```
# Run sim 1, 1000 times
sr <- 1:1000

sim1_m <- future_lapply(sr, function(x) heather_dry(plants=100, m=rmo, r=-4, years=100), future.seed=TRUE)
```

Since the function code uses a loop() and random number generator, you must include the argument future.seed=TRUE to avoid errors.

This should run fairly fast on most systems, and the result is a list of the simulation raw data. However, since R keeps results in memory, you may run into problems if you have low system resources, or you are running very large simulations.

Instead, you can write the results to disk as .csv files, and import them back into R later. This is also useful for saving the raw model output for future analysis.

```
# Run model and save output directly to disk
future_lapply(sr, function(x) heather_dry(plants=100, m=rmo, r=-4, years=100, write=TRUE, fn=paste0("~/model_run/sim1/", x, ".csv")), future.seed=TRUE)
```

Notice that the data was not saved as an R object. Each model run will save output in a separate file - in this example, 1000 .csv files are created.

To import this data back into R for analysis, use the following code (uses data.table).

```
# Load raw data
ds <- "~/model_run/"
run <- "sim1/"

fl <- list.files(path=paste0(ds, run), ".csv")
sims <- lapply(fl, function(x) fread(paste0(ds, run, x), drop=1, skip=1))
# rename
names(sims) <- fl
```

The result is a list of raw model output data as data.table objects. You can now run the analysis (as previous) on this data, by using lapply() and future_lapply().

```
# Reformat data
sims_r <- lapply(sims, function(x) simdf(x, stack=FALSE))
# Wilcoxon rank sum test
sims_wt <- future_lapply(sims_r, function(x) lapply(x[,-1], function(x2) wilcox.test(x[,1], x2)))
# Get p values
sims_wt_p <- lapply(sims_wt, function(x) unlist(lapply(x, function(x2) unlist(x2$p.value))))
# Adjust p values
sims_wt_pa <- lapply(sims_wt_p, function(x) p.adjust(x, "BY"))
# Stable period
sims_sp <- lapply(sims_wt_pa, function(x) stable(x))
```

You can also plot the results of individual simulations as per previous example.

### Statistics of all simulations

Having run the simulation multiple times, you will likely want to know some average stats across all of these simulations. For example, the mean age of the heather plants in the final year of each simulation run.

```
# Mean age at final year (column 1 in formatted data)
mean_fy <- lapply(sims_r, function(x) mean(x[[1]], na.rm=TRUE)) |> unlist() |> unname()

# Overall mean
mean(mean_fy, na.rm=TRUE)
# Standard deviation
sd(mean_fy, na.rm=TRUE)
```

Or you might want to know the mean age of the heather plants for every year of each simulation run.

```
# Every year, every sim
mean_ev <- future_lapply(sims_r, function(x) lapply(x, function(x) mean(x, na.rm=TRUE))) 

# Each sim average
mean_a <- lapply(mean_ev, function(x) mean(unlist(x), na.rm=TRUE)) |> unlist() |> unname() 
```

If you want information about the stable periods for each simulation

```
# Stable periods
stable_a <- lapply(sims_sp, function(x) nrow(x)) |> unlist() |> unname() 

# Which year does the first stable period occur
stable1_a <- lapply(sims_sp, function(x) x$year[1]) |> unlist() |> unname() 

# How long is the first stable period
stable1len_a <- lapply(sims_sp, function(x) x$length[1]) |> unlist() |> unname() 

# Overall mean
mean(stable_a, na.rm=TRUE)
mean(stable1_a, na.rm=TRUE)
mean(stable1len_a, na.rm=TRUE)
```
