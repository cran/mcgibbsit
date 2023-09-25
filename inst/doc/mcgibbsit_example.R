## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

nformat <- function(x)
  format(x, big.mark = ',')

## ----setup--------------------------------------------------------------------
library(mcgibbsit)

set.seed(42)        # for reproducibility
tmpdir <- tempdir()

## -----------------------------------------------------------------------------
# Define a function to generate the output of our imaginary MCMC sampler
gen_samples <- function(run_id, nsamples=200)
{
  x <- matrix(nrow = nsamples+1, ncol=4)
  colnames(x) <- c("alpha","beta","gamma", "nu")
  
  x[,"alpha"] <- exp(rnorm (nsamples+1, mean=0.025, sd=0.025))
  x[,"beta"]  <- rnorm (nsamples+1, mean=53,    sd=14)
  x[,"gamma"] <- rbinom(nsamples+1, 20,         p=0.15) + 1
  x[,"nu"]    <- rnorm (nsamples+1, mean=x[,"alpha"] * x[,"beta"], sd=1/x[,"gamma"])
#'
  # induce serial correlation of 0.25
  x <- 0.75 * x[2:(nsamples+1),] + 0.25 * x[1:nsamples,]

  # induce ~50% acceptance rate
  accept <- runif(nsamples) > 0.50
  for(i in 2:nsamples)
    if(!accept[i]) x[i,] <- x[i-1,]

  write.table(
    x,
    file = file.path(
      tmpdir,
      paste("mcmc", run_id, "csv", sep=".")
      ),
    sep = ",",
    row.names = FALSE
  )
}

## -----------------------------------------------------------------------------
# Generate and load 3 runs 
for(i in 1:3)
  gen_samples(i, 200)
  
mcmc.3 <- read.mcmc(
  3, 
  file.path(tmpdir, "mcmc.#.csv"), 
  sep=",",
  col.names=c("alpha","beta","gamma", "nu")
  )

## ----trace_density_plot_3,fig.width=8,fig.height=8----------------------------
# Trace and Density Plots
plot(mcmc.3)

## -----------------------------------------------------------------------------
# And check the necessary run length 
mcg.3 <- mcgibbsit(mcmc.3)
print(mcg.3)

## -----------------------------------------------------------------------------
# Generate and load 7 more runs 
for(i in 3 + (1:7))
  gen_samples(i, 200)
  
mcmc.10 <- read.mcmc(
  10, 
  file.path(tmpdir, "mcmc.#.csv"), 
  sep=",",
  col.names=c("alpha","beta","gamma", "nu")
  )

## ----trace_density_plot_10,fig.width=8,fig.height=8---------------------------
# Trace and Density Plots
plot(mcmc.10)

## -----------------------------------------------------------------------------
# And check the necessary run length 
mcg.10 <- mcgibbsit(mcmc.10)
print(mcg.10)

## -----------------------------------------------------------------------------
summary(mcmc.10)

