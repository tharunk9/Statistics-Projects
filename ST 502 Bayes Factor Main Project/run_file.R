args = commandArgs(trailingOnly=TRUE)


                      ### FUNCTIONS ###



# References:
#   Verhagen, J., & Wagenmakers, E.-J. (2014). Bayesian tests to quantify the
#   result of a replication attempt. Journal of Experimental Psychology:
#   General, 143(4), 1457–1475. http://doi.org/10.1037/a0036731

## to download the ReplicationBF package from the paper install devtools and download git code
# library("devtools")
# devtools::install_github("neurotroph/ReplicationBF")

suppressMessages(library("plyr"))
suppressMessages(library("MCMCpack"))
suppressMessages(library("R.utils"))
suppressMessages(library("ggplot2"))
suppressMessages(library("MBESS"))
suppressMessages(library("tidyverse"))
suppressMessages(library("ReplicationBF"))
suppressMessages(library("ggthemes"))
suppressMessages(library("BayesFactor"))


# REPLICATION BAYES FACTOR FOR F-TESTS -----------------------------------------

BFrep.ANOVA <- function( ... ) { BFrep.ANOVA.MCMC(...) }

BFrep.ANOVA.MCMC <- function(F.orig, df.effect.orig, df.error.orig, N.orig,
                             F.rep,  df.effect.rep,  df.error.rep, N.rep,
                             plot.posterior = F, M = 1e5, burnin = 500)
{
  # Loglikelihood of Original study
  loglik <- function(theta, XF, Xdf1, Xdf2, XN) {
    ifelse(theta < 0,
           log(0),
           df(XF, Xdf1, Xdf2, ncp = theta * XN, log = T))
  }
  
  captureOutput(expr={sample.prior <- MCMCmetrop1R(loglik, theta.init = c(0),
                                                   logfun = T, burnin=burnin, mcmc=M,
                                                   thin=1, verbose=0,
                                                   #optim.method="L-BFGS-B", optim.lower=0, 
                                                   V = matrix(ncol=1, nrow=1, c(1)),
                                                   XF= F.orig, Xdf1 = df.effect.orig,
                                                   Xdf2 = df.error.orig, XN = N.orig)})
  if (plot.posterior)
    plot(sample.prior)
  
  # H0: Likelihood for Skeptic
  lik.h0 = df(F.rep, df.effect.rep, df.error.rep)
  
  # Hr: Likelihood for Proponent
  lik.hr = mean(df(F.rep, df.effect.rep, df.error.rep, sample.prior * N.rep))
  
  return (list( bf = lik.hr / lik.h0, orig.posterior = sample.prior ))
}

# REPLICATION BAYES FACTOR FOR T-TESTS -----------------------------------------

BFrep.t <- function( ... ) { BFrep.t.MCMC(...) }

BFrep.t.MCMC <- function(t.orig, n1.orig, n2.orig,
                         t.rep, n1.rep, n2.rep,
                         plot.posterior = F, M = 1e5, burnin = 500)
{
  df.orig = (n1.orig + n2.orig - 2)
  sqrt.n.orig = sqrt( 1 / (1/n1.orig + 1/n2.orig))
  
  df.rep = (n1.rep + n2.rep - 2)
  sqrt.n.rep = sqrt( 1 / (1/n1.rep + 1/n2.rep))
  
  # Step 1: approximate posterior of original study
  loglik <- function(theta, Xt, Xdf, XN) { dt(Xt, Xdf, ncp=theta*XN, log=T) }
  captureOutput(expr={sample.prior <- MCMCmetrop1R(loglik, theta.init=c(0),
                                                   logfun=T, mcmc=M,
                                                   burnin=burnin,
                                                   verbose=0,
                                                   Xt=t.orig, Xdf=df.orig,
                                                   XN=sqrt.n.orig)})
  
  # Step 2: Likelihood for Skeptic (H0)
  likelihood.h0 = dt(t.rep, df.rep)
  
  # Step 3: Likelihood for Proponent (Hr)
  if (plot.posterior)
    plot(sample.prior)
  likelihood.hr = mean(dt(t.rep, df.rep, ncp=sample.prior * sqrt.n.rep))
  
  return (list(bf = likelihood.hr / likelihood.h0,
               orig.posterior = sample.prior ))
}

BFrep.t.normalApprox <- function(t.orig, n1.orig, n2.orig,
                                 t.rep, n1.rep, n2.rep,
                                 plot.posterior = F, M = 1e6) {
  df.orig = (n1.orig + n2.orig - 2)
  sqrt.n.orig = sqrt( 1 / (1/n1.orig + 1/n2.orig))
  
  df.rep = (n1.rep + n2.rep - 2)
  sqrt.n.rep = sqrt( 1 / (1/n1.rep + 1/n2.rep))
  
  # Step 1: Calculate CI for NCP of t-Distribution under H1 of original study
  #delta.lower <- confint.ncp.t(t.orig, df.orig)$lower
  delta.lower <- conf.limits.nct(t.orig, df.orig)$Lower.Limit
  
  # Step 2: Calculate parameters for normal approximation to posterior
  mu.delta = t.orig / sqrt.n.orig
  sd.delta = abs(((t.orig - delta.lower) / qnorm(.025)) / sqrt.n.orig)
  
  # Step 3: Likelihood for Skeptic (H0)
  likelihood.h0 = dt(t.rep, df.rep)
  
  # Step 4: Likelihood for Proponent (Hr)
  sample.prior = rnorm(M, mu.delta, sd.delta)
  likelihood.hr = mean(dt(t.rep, df.rep, sample.prior * sqrt.n.rep))
  
  return( list(bf = likelihood.hr / likelihood.h0,
               orig.posterior = sample.prior) )
}

# HELPER FUNCTION --------------------------------------------------------------

check.equality.all <- function(t.orig, n1.orig, n2.orig, t.rep, n1.rep, n2.rep)
{
  # check.equality.all(...) is a helper function, that calcules three variants
  # of the Replication Bayes factor and shows a plot with the densities of the
  # original study's posterior distribution. This plot was used in a previous
  # version of the manuscript.
  
  ret <- list(
    t.normal = BFrep.t.normalApprox(t.orig, n1.orig, n2.orig,
                                    t.rep, n1.rep, n2.rep),
    t.MCMC = BFrep.t.MCMC(t.orig, n1.orig, n2.orig,
                          t.rep, n1.rep, n2.rep),
    F.MCMC = BFrep.ANOVA.MCMC(t.orig^2, 1, (n1.orig+n2.orig-2),
                              (n1.orig+n2.orig),
                              t.rep^2, 1, (n1.rep+n2.rep-2),
                              (n1.rep+n2.rep))
  )
  
  # Create plot visualizing different posteriors
  p <- ggplot() +
    geom_density(aes(ret$t.normal$orig.posterior, linetype="1", color="3")) +
    geom_density(aes(ret$t.MCMC$orig.posterior, linetype="1", color="1")) +
    geom_density(aes(ret$F.MCMC$orig.posterior, linetype="2", color="1")) +
    scale_x_continuous(name="Effect Size") +
    scale_y_continuous(name="Posterior Density") +
    scale_color_discrete(labels=c("MCMC", "Grid Approximation",
                                  "Normal Approximation"),
                         breaks=c("1", "2", "3")) +
    scale_linetype_discrete(labels=c("t-Test", "F-Test"), breaks=c("1", "2")) +
    guides(color=guide_legend(title=""), linetype=guide_legend(title=""))
  print(p)
  
  return(lapply(ret, '[[', 'bf'))
}
                      ### FUNCIONS (END) ###




##########################################################################
#------------------------------------------------------------------------#
##########################################################################



                ###  BFrep Simulation 1 ### 



# This simulation shows the Replication Bayes factor in different combinations
# of parameters, i.e. in different scenarios.

# Simulation parameters and set up ---------------------------------------------
n.set  <- c(10,25,50)
k.set  <- c(2,3,4,5)
f2.set <- c(0.00001, 0.15, 0.35)

combinations <- expand.grid(n.orig = n.set, n.rep = n.set, k = k.set,
                            f2.orig = f2.set, f2.rep = f2.set)
# Perform simulations only for cases, where we use a larger sample in the
# replication than in the original study to be more realistic (and save time
# when running the simulations).
combinations <- combinations[combinations$n.rep >= combinations$n.orig,]

# Run simulations --------------------------------------------------------------
results <- apply(combinations, 1, function(row) {
  # Calculated observed F-values
  
  F.obs.orig <- (row['f2.orig'] * row['n.orig'] * row['k']) / (row['k'] - 1)
  F.obs.rep <- (row['f2.rep'] * row['n.rep'] * row['k']) / (row['k'] - 1)
  
  # Calculate Replication Bayes factor for F-tests
  BFrep <- RBF_Ftest(# Original study .........................
    F.obs.orig, c(row['k'] - 1,
                  row['n.orig'] - row['k']),
    row['n.orig'] * row['k'],
    # Replication study ......................
    F.obs.rep, c(row['k'] - 1,
                 row['n.rep'] - row['k']),
    row['n.rep'] * row['k'],
    # Settings ...............................
    M = 1e4,
    store.samples = FALSE)
  
  # Collect results
  c(row,
    # Observed F-values
    F.orig = F.obs.orig, F.rep = F.obs.rep,
    # Replication Bayes factor
    BFrep = BFrep$bayesFactor)
})
# Prepare data frae with results
results <- as.data.frame(t(results))
colnames(results) <- c("n.orig", "n.rep", "k", "f2.orig", "f2.rep",
                       "F.orig", "F.rep", "BFrep")

# Save and/or load simulation results to disk ----------------------------------
temp <- paste0(args[1], "\\BFrep_Sim1.rds")
saveRDS(results, file=temp)#temp


               ###  BFrep Simulation 1 (END) ### 





##########################################################################
#------------------------------------------------------------------------#
##########################################################################


                 
                 ###  BFrep Simulation 3  ###



# The simulations in this script show the behavior of the proposed Replication
# Bayes Factor in comparison to the originally proposed approach and in some
# additional scenarios.


# Parameter configuration ------------------------------------------------------
n.orig.set = c(15, 50)
n.rep.set = c(50, 100)
d.set = c(0.2, 0.6, 1)
sim.count = 2

sim.combinations <- expand.grid(n.orig = n.orig.set, n.rep = n.rep.set,
                                d.orig = d.set, d.rep = d.set,
                                i = 1:sim.count)
# Perform simulations only for cases, where we use a larger sample in the
# replication than in the original study to be more realistic (and save time
# when running the simulations).
sim.combinations <- sim.combinations[sim.combinations$n.rep >=
                                       sim.combinations$n.orig,]
sim.combinations$sim.no <- seq.int(nrow(sim.combinations))

# Fixed seed for reproducible simulation results
set.seed(20032018)

# Setup ------------------------------------------------------------------------
.pb = txtProgressBar(min = 0, max = nrow(sim.combinations), style = 3)

# Function for the simulation --------------------------------------------------
simulate.studypair <- function(r) {
  # Original study .............................................................
  if ((sim.count) > 1) {
    # If `sim.count` is set to a value larger than 1, d.orig is true population
    # effect size and random samples are drawn. This introduces sampling
    # variation in the results of the simulation.
    sample1 = rnorm(r['n.orig'], mean = 0, sd = 1)
    sample2 = rnorm(r['n.orig'], mean = r['d.orig'], sd = 1)
    
    t.orig = as.numeric(t.test(sample1, sample2)$statistic)
    d.obs.orig = mean(sample2 - sample1)
  } else {
    # If `sim.count` is 1, d.orig is observed effect size and the t-value is
    # calculated from d.orig for use in the Replication Bayes factor.
    t.orig <- r['d.orig'] * r['n.orig'] / sqrt(2*r['n.orig'])
    d.obs.orig <- NA
  }
  
  # Replication study ..........................................................
  if ((sim.count) > 1) {
    # If `sim.count` is set to a value larger than 1, d.rep is true population
    # effect size and random samples are drawn. This introduces sampling
    # variation in the results of the simulation.
    sample1.rep = rnorm(r['n.rep'], mean = 0, sd = 1)
    sample2.rep = rnorm(r['n.rep'], mean = r['d.rep'], sd = 1)
    
    t.rep = as.numeric(t.test(sample1.rep, sample2.rep)$statistic)
    d.obs.rep = mean(sample2.rep - sample1.rep)
  } else {
    # If `sim.count` is 1, d.rep is observed effect size and the t-value is
    # calculated from d.rep for use in the Replication Bayes factor.
    t.rep <- r['d.rep'] * r['n.rep'] / sqrt(2*r['n.rep'])
    d.obs.rep <- NA
  }
  
  # Combined study
  sample1.combined = c(sample1, sample1.rep)
  sample2.combined = c(sample2, sample2.rep)
  
  t.combined = as.numeric(t.test(sample1.combined,
                                 sample2.combined)$statistic)
  
  # Update progress bar
  setTxtProgressBar(.pb, r['sim.no'])
  
  # Collect results
  c(
    r,
    # observed t-values .......................................................
    t.orig = t.orig,
    t.rep = t.rep,
    # Observed effect sizes ...................................................
    d.obs.orig = d.obs.orig,
    d.obs.rep = d.obs.rep,
    # RBF for t-Test (Monte Carlo estimate) ...................................
    BFrep.t.MCe = BFrep.t.normalApprox(t.orig, r['n.orig'], r['n.orig'],
                                       t.rep, r['n.rep'], r['n.rep'])$bf,
    # RBF for t-Test (Importance sampling estimate) ...........................
    BFrep.t.ISe = RBF_ttest(t.orig, c(r['n.orig'], r['n.orig']),
                            t.rep, c(r['n.rep'], r['n.rep']))$bayesFactor,
    # RBF for F-tests (with F = t^2) ..........................................
    BFrep.F.ISe = RBF_Ftest(t.orig^2, c(1, (2 * r['n.orig'] - 2)),
                            2 * r['n.orig'],
                            t.rep^2, c(1, (2 * r['n.rep'] - 2)),
                            2 * r['n.rep'],
                            M = 1e5)$bayesFactor
    # Bayes Factor using Evidence Updating (see Ly et al., 2017) ..............
    #BFrep.t.EU = NA
  ) 
}

# Run simulation ---------------------------------------------------------------
sim.results <- apply(sim.combinations, 1, simulate.studypair)
data <- as.data.frame(t(sim.results))

# Calculate some variables -----------------------------------------------------
data$diff.MCvIS  <- data$BFrep.t.MCe - data$BFrep.t.ISe
data$ratio.MCvIS <- data$BFrep.t.MCe / data$BFrep.t.ISe

data$diff.tF  <- data$BFrep.t.ISe - data$BFrep.F.ISe
data$ratio.tF <- data$BFrep.t.ISe / data$BFrep.F.ISe

#data$diff.ISvEU  <- data$BFrep.t.ISe - data$BFrep.t.EU
#data$ratio.ISvEU <- data$BFrep.t.ISe / data$BFrep.t.EU

# Save and/or load simulation results to disk ----------------------------------
temp2 <- paste0(args[1], "\\BFrep_Sim3.rds")
saveRDS(data, file=temp2)#temp
                
                   ###  BFrep Simulation 3 (END)  ###





##########################################################################
#------------------------------------------------------------------------#
##########################################################################


lapply(names(sessionInfo()$basePkgs), require, character.only = TRUE)
invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE, force=TRUE))

suppressMessages(library("plyr"))
suppressMessages(library("tidyverse"))
suppressMessages(library("ReplicationBF"))


                 ###  BFrep Example 3 ###



set.seed(26032018)
n.orig <- 15
n.rep <- 30

# Set up original study
orig.group.1 <- rnorm(n.orig, mean = 1.5, sd = 1)
orig.group.2 <- rnorm(n.orig, mean = 2.2, sd = 1)
orig.group.3 <- rnorm(n.orig, mean = 2.9, sd = 1)
orig.study <- data.frame(y = c(orig.group.1, orig.group.2, orig.group.3),
                         x = c(rep(1, n.orig),
                               rep(2, n.orig),
                               rep(3, n.orig)))
orig.study$x <- as.factor(orig.study$x)
# ANOVA on original study
anova.orig <- unlist(summary(aov(y ~ x, data = orig.study)))
# Extract ANOVA statistics
F.orig <- anova.orig["F value1"]
df.effect.orig <- anova.orig["Df1"]
df.error.orig <- anova.orig["Df2"]
N.tot.orig <- n.orig * 3
f2.orig <- F.orig * df.effect.orig / df.error.orig

# Set up replication study - with reversed pattern of effects
rep.group.1 <- rnorm(n.rep, mean = 2.9, sd = 1)
rep.group.2 <- rnorm(n.rep, mean = 2.2, sd = 1)
rep.group.3 <- rnorm(n.rep, mean = 1.5, sd = 1)
rep.study <- data.frame(y = c(rep.group.1, rep.group.2, rep.group.3),
                        x = c(rep(1, n.rep),
                              rep(2, n.rep),
                              rep(3, n.rep)))
rep.study$x <- as.factor(rep.study$x)
# ANOVA on replication study
anova.rep <- unlist(summary(aov(y ~ x, data = rep.study)))
# Extract ANOVA statistics
F.rep <- anova.rep["F value1"]
df.effect.rep <- anova.rep["Df1"]
df.error.rep <- anova.rep["Df2"]
N.tot.rep <- n.rep * 3
f2.rep <- F.rep * df.effect.rep / df.error.rep

# Calculate RBF on F-test
rbf.ex3 <- RBF_Ftest(F.orig, c(df.effect.orig, df.error.orig), N.tot.orig,
                     F.rep, c(df.effect.rep, df.error.rep), N.tot.rep,
                     store.samples = T)
#rbf.ex3$bayesFactor

# Generate compile data.frame for interaction plot
studies.ex3 <- data.frame(y = c(orig.study$y, rep.study$y),
                          x = as.factor(c(orig.study$x, rep.study$x)),
                          study = c(rep("Original", 3*n.orig),
                                    rep("Replication", 3*n.rep)))
studies.ex3.agg <- studies.ex3 %>%
  as_tibble() %>%
  group_by(study, x) %>%
  summarise(mean = mean(y),
            sd = sd(y),
            se = sd(y)/n())


# Perform post-hoc analysis
ph.t.orig.12 <- t.test(orig.group.1, orig.group.2)$statistic
ph.t.rep.12 <- t.test(rep.group.1, rep.group.2)$statistic
value1 <- RBF_ttest(ph.t.orig.12, c(n.orig, n.orig),
          ph.t.rep.12, c(n.rep, n.rep))$bayesFactor

ph.t.orig.13 <- t.test(orig.group.1, orig.group.3)$statistic
ph.t.rep.13 <- t.test(rep.group.1, rep.group.3)$statistic
value2 <- RBF_ttest(ph.t.orig.13, c(n.orig, n.orig),
          ph.t.rep.13, c(n.rep, n.rep))$bayesFactor

ph.t.orig.23 <- t.test(orig.group.2, orig.group.3)$statistic
ph.t.rep.23 <- t.test(rep.group.2, rep.group.3)$statistic
value3 <- RBF_ttest(ph.t.orig.23, c(n.orig, n.orig),
          ph.t.rep.23, c(n.rep, n.rep))$bayesFactor

RBF_data <- data.frame(matrix(NA, nrow = 3, ncol = 3))

RBF_data[1,1] <- ph.t.orig.12
RBF_data[1,2] <- ph.t.rep.12
RBF_data[1,3] <- value1
RBF_data[2,1] <- ph.t.orig.13
RBF_data[2,2] <- ph.t.rep.13
RBF_data[2,3] <- value2
RBF_data[3,1] <- ph.t.orig.23
RBF_data[3,2] <- ph.t.rep.23
RBF_data[3,3] <- value3

row.names(RBF_data) <- c("Group1vs2","Group1vs3","Group2vs3")
colnames(RBF_data) <- c("t.original","t.replication","bayesFactor")

# Save and/or load simulation results to disk ----------------------------------
temp3 <- paste0(args[1], "\\BFrep_Ex3.rds")
saveRDS(studies.ex3.agg, file=temp3)#temp
temp4 <- paste0(args[1], "\\BFrep_Ex3data.rds")
saveRDS(RBF_data, file=temp4)#temp
          

        ###  BFrep Example 3 (END) ###



##########################################################################
#------------------------------------------------------------------------#
##########################################################################