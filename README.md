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
source("https://raw.githubusercontent.com/benbell95/heather_age/refs/heads/main/r/heather_dry_model.r")
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
source("https://raw.githubusercontent.com/benbell95/heather_age/refs/heads/main/r/heather_model_analysis_functions.r")
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

### "Stable distribution" analysis

While the basic stats show changes over time, we also want to know when the age distribution reaches "stability" - that is, when the plant age distribution has very high similarity from year to year.

For this, we use the [Kolmogorov-Smirnov test](https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test) to assess the age distribution of each year in the model run, which is compared to the "final age distribution", calculated as the mean of the final 10 years of the model run. This produces the Kolmogorov's D statistic, a value between 0 and 1, with values close to 0 meaning high similarity, and values close to 1 meaning low similarity in the age distribution. The test also produces a p value, with significant p values (typically < 0.05) indicating significant differences in the distribution. 

The analysis is run within the scenario - therefore, to produce a "natural" age distribution, the model should be run for at least 100 years without any management occurring. The final years age distribution could then be considered as the "stable natural" age distribution. If management occurs for the duration of the model run, the age distribution will also become stable (resulting in low D statistic like the non-management scenario), but this should be considered as a "stable management" age distribution, rather than a natural one. 

Alternatively, you could compare the age distribution of any model scenario to a pre-defined age distribution of your choosing.

To run the analysis, we will use the base r ks.test() function. 

The helper functions are designed to work with data.tables (since typically you might simulate thousands of plants or years), therefore it is necessary to convert the model output data first.

```
# Load data.table library 
library("data.table")

# Convert to data.table
sim1 <- as.data.table(sim1)
sim2 <- as.data.table(sim2)
```

Next we'll run the test on each and every year of the model run, comparing to the final 10 year's age distribution.

```
# Run ks tests
years <- 100
sim1_ks <- lapply(1:years, \(y) ks.test(unlist(sim1[91:100,]), unlist(sim1[y,])))
sim2_ks <- lapply(1:years, \(y) ks.test(unlist(sim2[91:100,]), unlist(sim2[y,])))
```

Then, we'll extract the D and P values, which we can use to assess the similarity of the age distributions.

```
# Extract D values
sim1_ks_d <- lapply(sim1_ks, "[", 1) |> unlist() |> unname()
sim2_ks_d <- lapply(sim2_ks, "[", 1) |> unlist() |> unname()

# Extract P values
sim1_ks_p <- lapply(sim1_ks, "[", 2) |> unlist() |> unname()
sim2_ks_p <- lapply(sim2_ks, "[", 2) |> unlist() |> unname()
```

For the p values, since we have run multiple tests, you should adjust these values to avoid false positives.

```
# Adjust p values
sim1_ks_pa <- p.adjust(sim1_ks_p, method="BY")
sim2_ks_pa <- p.adjust(sim2_ks_p, method="BY")
```

You can create a quick plot to take a look at the results.

```
plot(sim1_ks_d, type="o", ylim=c(0, 1), xlab="Year")
lines(sim1_ks_pa, col="red")
```

On the plot, you will see that near the beginning of the model run, D values are high, while p values are low. As time progresses, this relationship is reversed, with lower D values indicating more similar age distributions.

Since in this example we are only looking at 100 plants, and we've only run the model once, the values are not particularly smooth - ideally, you would simulate thousands of plant ages and repeat the model several times.

You can compare these results to the second scenario.

You could also take a look at the actual distributions by plotting a histogram or density plot of the ages from the raw model output. For example, to compare the first 10 years and the last 10 years:

```
# Compare age distributions
layout(matrix(1:2, ncol=2))
hist(sim1[1:10,] |> unlist(), freq=FALSE)
hist(sim1[91:100,] |> unlist(), freq=FALSE)
```

The analysis so far gives a good indication that the age distributions differ, but we can define the "stable" period more clearly using one of the helper functions. The function uses the Kolmogorov D statistic to define the stable period, based on a threshold value that you set.

We'll use the mean D statistic for the last 10 years of the model run, plus 0.1. For a model run which uses thousands of plants and is run multiple times, you'll want to use a lower value (typically 0.015).

First, calculate the mean D statistic, then run the stable function.

```
# Stable period
sim1_sp <- stable(sim1_wt_pa)
sim2_sp <- stable(sim2_wt_pa)
```

The results will show the year which the "stable" period starts, and how long it lasts. There may be multiple stable periods.

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

The result is a list of raw model output data as data.table objects. You can now run the analysis as per previous examples on this data, by using lapply() and future_lapply().
