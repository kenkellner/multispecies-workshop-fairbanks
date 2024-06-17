#' ---
#' title: An R Markdown document converted from "nimble_nosolutions.ipynb"
#' output: html_document
#' ---
#' 
#' # Introduction to NIMBLE
#' 
#' Ken Kellner
#' 
#' **What is NIMBLE/`nimble`?**:
#' 
#' * Software for fitting statistical models (de Valpine et al., 2017)
#' * Best-known for Bayesian analysis
#' * Modeling language (BUGS-based)
#' * Algorithm language
#' 
#' # Outline
#' 
#' 1. Simulate dataset
#' 2. The BUGS language
#' 3. Simple `nimble` workflow
#' 4. Customized `nimble` workflow
#' 5. Other algorithms
#' 6. Macros
#' 
#' # Simulate dataset
#' 
#' * Probability of wolf occurrence at a site
#' * Probability depends on forest cover
#' * Ignore imperfect detection
#' 
#' ![](wolf.jpeg)
#' 
#' ## Simulation
#' 
#' First some design and covariate data:
#' 
## -----------------------------------------------------------------------------
set.seed(1)
M <- 100 # number of sites
forest <- runif(M, 0, 100)
forest_scale <- as.numeric(scale(forest))

#' 
#' Calculate occurrence probability $\psi$:
#' 
#' $$\mu_i = \beta_0 + \beta_1 \cdot forest_i $$
#' 
#' $$ \psi_i = \mathrm{ilogit}(\mu_i) $$
#' 
## -----------------------------------------------------------------------------
beta0 <- qlogis(0.3) # average occurrence is 0.3
beta1 <- 0.5 # positive effect of forest
truth <- c(beta0=beta0, beta1=beta1)

#' 
## -----------------------------------------------------------------------------
mu <- psi <- numeric(M)
for (i in 1:M){
  mu[i] <- beta0 + beta1 * forest_scale[i]
  psi[i] <- plogis(mu[i])
}

#' 
#' Simulate observed data $y$:
#' 
#' $$ y_i \sim \mathrm{Binomial}(\psi_i, 1) $$
#' 
#' 
## -----------------------------------------------------------------------------
y <- numeric(M)
for (i in 1:M){
  y[i] <- rbinom(1, 1, psi[i]) # simulate random Bernoulli
}
y

#' 
#' ## Visualize simulated dataset
#' 
#' Start by breaking up the response data into categories by forest cover
#' 
## -----------------------------------------------------------------------------
breaks <- quantile(forest, c(0, 0.2, 0.4, 0.6, 0.8, 1))
forest_cat <- cut(forest, breaks = breaks)
tab <- aggregate(y ~ forest_cat, FUN=mean)
tab$n <- table(forest_cat)
tab$se <- sqrt(tab[,2] * (1-tab[,2]) / tab$n)

#' 
#' Plot the forest - occurrence relationship
#' 
## -----------------------------------------------------------------------------
library(ggplot2)
options(repr.plot.width=10, repr.plot.height=10)

ggplot(data = tab) +
  geom_point(aes(x=forest_cat, y = y)) +
  geom_errorbar(aes(x = forest_cat, ymin= y - se, ymax=y + se), width=0.2) +
  labs(y = "Wolf occurrence", x = "Forest cover") +
  theme_bw(base_size = 24) +
  theme(panel.grid = element_blank())

#' 
#' ## Fit the model with `glm`
#' 
## -----------------------------------------------------------------------------
fit_glm <- glm(y ~ forest_scale, family=binomial)
summary(fit_glm)

est_glm <- coef(fit_glm)

#' 
#' # The BUGS language
#' 
#' * Used mainly (but not exclusively) for Bayesian analysis
#' * Defines relationships between *nodes*
#' * Nodes may be data or parameters to estimate
#' 
#' * Model consists of *declarations* about nodes
#' * Nodes declared deterministic using `<-`
#' * or coming from a probability distribution (using `~`)
#' 
#' ```r
#' x <- 2 + 2
#' y <- alpha + beta
#' z ~ dnorm(0, 1)
#' 
#' ```
#' 
#' * Order of declarations generally doesn't matter (resulting in a directed acyclic graph or DAG)
#' 
#' ## BUGS Code for the model
#' 
#' Data nodes: `M`, `y`, `forest_scale`
#' 
#' Parameter nodes: `beta0`, `beta1`, `mu`, `psi`
#' 
#' ```r
#' # Likelihood
#' for (i in 1:M){
#'   # Calculate linear predictor
#'   mu[i] <- beta0 + beta1 * forest_scale[i]
#' 
#'   # Transform it to a probability
#'   logit(psi[i]) <- mu[i]
#'   #psi[i] <- ilogit(mu[i]) # equivalent
#' 
#'   # Data model
#'   y[i] ~ dbern(psi[i])
#'   #y[i] ~ dbinom(psi[i], 1) # equivalent
#' }
#' 
#' # Priors (uninformative/vague)
#' beta0 ~ dnorm(0, sd = 10)
#' beta1 ~ dnorm(0, sd = 10)
#' 
#' ```
#' 
#' # Simple `nimble` workflow
#' 
#' 1. Create list of data
#' 2. Write BUGS code with `nimbleCode`
#' 3. Set parameters to monitor
#' 4. Run `nimbleMCMC`
#' 5. Do things with the output
#' 
## -----------------------------------------------------------------------------
library(nimble)

#' 
#' 1. Create list of data
#' 
## -----------------------------------------------------------------------------
nimble_data <- list(M = M, forest_scale = forest_scale, y = y)

#' 
#' 2. Write BUGS code
#' 
## -----------------------------------------------------------------------------
code <- nimbleCode({
  # Likelihood
  for (i in 1:M){
    mu[i] <- beta0 + beta1 * forest_scale[i]   
    logit(psi[i]) <- mu[i]
    
    y[i] ~ dbern(psi[i])
  }

  # Priors
  beta0 ~ dnorm(0, sd = 10)
  beta1 ~ dnorm(0, sd = 10)
})

#' 
## -----------------------------------------------------------------------------
code

#' 
#' 3. Set parameters to monitor
#' 
## -----------------------------------------------------------------------------
pars <- c("beta0", "beta1")

#' 
#' 4. Run `nimbleMCMC`
#' 
## -----------------------------------------------------------------------------
fit_nimble <- nimbleMCMC(code = code,
                         constants = nimble_data,
                         monitors = pars,
                         niter = 2000,
                         nburnin = 1000,
                         nchains = 3,
                         samplesAsCodaMCMC = TRUE)

#' 
#' ## 5. Do things with the output
#' 
#' Look at the basic structure of the output object
#' 
## -----------------------------------------------------------------------------
class(fit_nimble)
lapply(fit_nimble, head)

#' 
#' ### Model diagnostics
#' 
#' Using the `coda` package (`MCMCvis` is also a good option that Chris will show tomorrow)
#' 
## -----------------------------------------------------------------------------
par(mfrow=c(2,1))
coda::traceplot(fit_nimble)

#' 
## -----------------------------------------------------------------------------
# R-hat
coda::gelman.diag(fit_nimble)

#' 
## -----------------------------------------------------------------------------
# Effective sample size
coda::effectiveSize(fit_nimble)

#' 
#' ### Inference
#' 
## -----------------------------------------------------------------------------
summary(fit_nimble)

#' 
#' Look at the posterior distributions
#' 
## -----------------------------------------------------------------------------
plot(fit_nimble)

#' 
#' Compare estimates
#' 
## -----------------------------------------------------------------------------
est_nimble <- apply(as.matrix(fit_nimble), 2, mean)
comp <- cbind(truth = truth, glm = est_glm, nimble = est_nimble)
round(comp, 3)

#' 
#' ## Create a figure from the output
#' 
#' We'll make a figure that shows the relationship between forest and occurrence probability as a smooth curve.
#' 
#' First, convert the posterior samples into a matrix:
#' 
## -----------------------------------------------------------------------------
post <- as.matrix(fit_nimble)
dim(post)
head(post)

#' 
#' Make the x-axis of the plot, a sequence of forest values:
#' 
## -----------------------------------------------------------------------------
n <- 100
forest_rng <- range(forest)
forest_seq <- seq(forest_rng[1], forest_rng[2], length.out = n)
head(forest_seq)
forest_seq_scale <- (forest_seq - mean(forest)) / sd(forest)

#' 
#' For each value of `forest_seq_scale`, calculate a posterior distribution for the corresponding value of `psi`:
#' 
## -----------------------------------------------------------------------------
new_psi <- matrix(NA, nrow = nrow(post), ncol = n)

# Iterate over MCMC samples
for (i in 1:nrow(post)){
  # Iterate over forest values
  for (j in 1:n){
    # Calculate psi
    new_psi[i, j] <- plogis(post[i, "beta0"] + post[i, "beta1"] * forest_seq_scale[j])
  }
}

# posterior for 50th value of forest_seq
hist(new_psi[,50])

#' 
#' Calculate posterior summary stats for each value of `forest_seq`
#' 
## -----------------------------------------------------------------------------
plot_dat <- data.frame(
  forest = forest_seq,
  psi = apply(new_psi, 2, mean),
  low = apply(new_psi, 2, quantile, 0.025),
  high = apply(new_psi, 2, quantile, 0.975)
)

head(plot_dat)

#' 
#' Use the data frame of summary stats to make a plot:
#' 
## -----------------------------------------------------------------------------
ggplot(data = plot_dat) +
  geom_ribbon(aes(x = forest, ymin = low, ymax = high), alpha=0.2) +
  geom_line(aes(x = forest, y = psi)) +
  labs(x = "Forest cover", y = "Predicted psi and 95% CRI") +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank())

#' 
#' Now also show the "true" relationship between forest and $\psi$
#' 
## -----------------------------------------------------------------------------
true_line <- data.frame(x = forest, y = psi)
ggplot(data = plot_dat) +
  geom_ribbon(aes(x = forest, ymin = low, ymax = high), alpha=0.2) +
  geom_line(aes(x = forest, y = psi)) +
  geom_line(data = true_line, aes(x=forest, y = psi), col='red') +
  labs(x = "Forest cover", y = "Predicted psi and 95% CRI") +
  theme_bw(base_size = 18) +
  theme(panel.grid = element_blank())

#' 
#' # Customized `nimble` workflow
#' 
#' * Done in several steps instead of just one step (`nimbleMCMC`)
#' * Allows for more control
#' 
#' Workflow:
#' 
#' 1. Create `nimbleModel` object
#' 2. Configure MCMC
#' 3. Create MCMC object
#' 4. Compile model and MCMC
#' 5. Run analysis
#' 
#' ## 1. Create `nimbleModel` object
#' 
#' We will explicitly initialize our two parameters in this example
#' 
## -----------------------------------------------------------------------------
inits <- list(beta0 = 0, beta1 = 0)

#' 
## -----------------------------------------------------------------------------
# Create a NIMBLE model from BUGS code
r_model <- nimbleModel(code = code, constants = nimble_data,
                       inits = inits) 

#' 
#' What's contained in the model object?
#' 
## -----------------------------------------------------------------------------
r_model$getCode()

#' 
## -----------------------------------------------------------------------------
r_model$y
r_model$beta0

#' 
## -----------------------------------------------------------------------------
# The DAG, which is not very helpful in this case
plot(r_model$modelDef$graph)

#' 
#' ## 2. Configure MCMC
#' 
#' Create an MCMC configuration using `nimble`'s defaults
#' 
#' The default is Metropolis-Hastings random walk samplers for our parameters
#' 
## -----------------------------------------------------------------------------
mcmc_config <- configureMCMC(r_model, monitors = pars)

#' 
#' ## 3. Build the MCMC
#' 
## -----------------------------------------------------------------------------
r_mcmc <- buildMCMC(mcmc_config)

#' 
#' Technically we can now sample using the MCMC object, but it's very slow.
#' 
## -----------------------------------------------------------------------------
system.time(r_mcmc$run(10))
as.matrix(r_mcmc$mvSamples)

#' 
#' ## 4. Compile the model and MCMC sampler
#' 
#' `nimble` converts our R code to C++, which is then compiled to fast machine code
#' 
## -----------------------------------------------------------------------------
c_model <- compileNimble(r_model)
c_mcmc <- compileNimble(r_mcmc, project = r_model) # note project argument

#' 
#' Now we can sample much faster:
#' 
## -----------------------------------------------------------------------------
system.time(c_mcmc$run(10))

#' 
#' ## 5. Run the analysis
#' 
#' We'll run multiple chains, and we want to initialize each chain differently.
#' 
## -----------------------------------------------------------------------------
inits <- function(){
  list(beta0 = rnorm(1), beta1 = rnorm(1))
}
inits()

#' 
## -----------------------------------------------------------------------------
inits()

#' 
#' The `runMCMC` function generates the final samples using options similar to `nimbleMCMC`.
#' 
## -----------------------------------------------------------------------------
samples <- runMCMC(c_mcmc, inits = inits, nchains = 3, niter = 2000, nburnin = 1000, 
                   samplesAsCodaMCMC = TRUE)

summary(samples)

#' 
#' # Other algorithms
#' 
#' Examples:
#' 
#' * Hamiltonian Monte Carlo (HMC)
#' * Laplace approximation
#' * Pólya-Gamma sampler
#' 
#' ## HMC
#' 
#' * Used most famously by `Stan`. 
#' * Reduces correlation between successive samples (thus increasing effective sample size)
#' * Each iteration generally slower `nimble`'s default samplers
#' * Requires calculating gradients (derivatives) and thus, practically, automatic differentation (AD)
#' * Certain kinds of models are more complicated to fit (such as models with latent variables like occupancy)
#' 
## -----------------------------------------------------------------------------
library(nimbleHMC)

#' 
#' Start by building a new model object, and specifying that we'll need derivative capabilities
#' 
## -----------------------------------------------------------------------------
r_model <- nimbleModel(code = code, constants = nimble_data, buildDerivs = TRUE) 

#' 
#' Configure an HMC sampler and build it:
#' 
## -----------------------------------------------------------------------------
hmc_config <- configureHMC(r_model, monitors = pars)
r_hmc <- buildMCMC(hmc_config)

#' 
#' NUTS stands for "No U-Turn Sampler" which is a type of HMC (that Stan also uses)
#' 
#' Compile the model and sampler, same as before:
#' 
## -----------------------------------------------------------------------------
c_model <- compileNimble(r_model)
c_hmc <- compileNimble(r_hmc, project = r_model)

#' 
#' Take samples same as before.
#' 
## -----------------------------------------------------------------------------
samples_hmc <- runMCMC(c_hmc, inits = inits, nchains = 3, niter = 2000, nburnin = 1000,
                      samplesAsCodaMCMC = TRUE)

summary(samples_hmc)

#' 
#' As expected, our effective sample size will be higher than with the default sampler:
#' 
## -----------------------------------------------------------------------------
round(coda::effectiveSize(samples_hmc))
round(coda::effectiveSize(samples)) # MH-RW has fewer effective samples

#' 
#' ## Pólya-Gamma sampler
#' 
#' A sampler for logistic regression that uses Pólya-Gamma data augmentation.
#' 
#' * Just released in `nimble` 1.2
#' * Sampler goes on the parameters in the linear predictor
#' 
## -----------------------------------------------------------------------------
pg_config <- configureMCMC(r_model, monitors = pars, print=FALSE)
pg_config$removeSamplers(pars)
pg_config$addSampler(type = "sampler_polyagamma", target = pars, 
                     control = list(fixedDesignColumns = TRUE))
pg_config
r_pg <- buildMCMC(pg_config)

#' 
## -----------------------------------------------------------------------------
c_model <- compileNimble(r_model)
c_pg <- compileNimble(r_pg, project = r_model)

#' 
## -----------------------------------------------------------------------------
system.time(out_mcmc <- runMCMC(c_mcmc, niter = 10000, nburnin = 8000, samplesAsCodaMCMC = TRUE))
round(coda::effectiveSize(out_mcmc))

#' 
## -----------------------------------------------------------------------------
system.time(out_pg <- runMCMC(c_pg, niter = 10000, nburnin = 8000, samplesAsCodaMCMC = TRUE))
round(coda::effectiveSize(out_pg))

#' 
#' # Exercise: Adding a random intercept
#' 
#' Suppose individual sites are nested within a larger spatial unit called a "group"
#' 
#' We'll assign each site to one of 10 groups randomly
#' 
## -----------------------------------------------------------------------------
ngroups <- 10
group <- sample(1:ngroups, M, replace=TRUE)
group

#' 
## -----------------------------------------------------------------------------
# Add to our data
nimble_data_rand <- nimble_data
nimble_data_rand$ngroups <- ngroups
nimble_data_rand$group <- group

#' 
#' A logical way specify the model would be to have each `group` have its own intercept:
#' 
#' $$ \mu_i = \beta_{0_{group_i}} + \beta_1 \cdot forest_i $$
#' $$ \psi_i = \mathrm{ilogit}(\mu_i) $$
#' 
#' **Hint:** This requires "nested" indexing:
#' 
#' $$ \mu_i = \beta_0[group_i] + \beta_1 \cdot forest_i $$
#' 
#' Or in BUGS:
#' 
#' ```r
#' beta0[group[i]]
#' ```
#' 
#' The group-specific intercepts can then be specified as coming from a common normal distribution with hyperparameters $\beta_{0_{mean}}$ and $\beta_{0_{sd}}$:
#' 
#' $$ \beta_{0_{group}} \sim \mathrm{Normal}(\beta_{0_{mean}}, \beta_{0_{sd}}) $$
#' 
#' **Hint:** You'll need to loop over `ngroups` and specify the distribution for each $\beta_0$
#' 
#' ```r
#' for (g in 1:ngroups){
#'  ...   
#' }
#' ```
#' 
#' **Final hint:** You'll also need to specify priors on the hyperparameters $\beta_{0_{mean}}$ and $\beta_{0_{sd}}$
#' 
#' Here's the final code:
#' 
## -----------------------------------------------------------------------------
code_random <- nimbleCode({
 # insert here
})

#' 
## -----------------------------------------------------------------------------
# use nimbleMCMC here

#' 
#' ## Laplace Approximation
#' 
#' * Used by, e.g., `TMB`
#' * Based on maximum likelihood, not Bayesian
#' * Much faster
#' * Much easier to estimate random effects vs. standard max likelihood approaches
#' * Also requires gradient/derivatives
#' 
#' Create the model object (and specify we'll need derivatives):
#' 
## -----------------------------------------------------------------------------
r_model <- nimbleModel(code = code_random, constants = nimble_data_rand, buildDerivs = TRUE) 

#' 
#' Create a Laplace approximation algorithm object:
#' 
## -----------------------------------------------------------------------------
alg_laplace <- buildLaplace(r_model, paramNodes = c("beta0_mean", "beta0_sd", "beta1"), 
                            randomEffectsNodes=c("beta0"))

#' 
#' Compile everything:
#' 
## -----------------------------------------------------------------------------
c_model <- compileNimble(r_model)
c_alg <- compileNimble(alg_laplace, project = r_model)

#' 
#' Use the algorithm to find the maximum likelihood estimates:
#' 
## -----------------------------------------------------------------------------
mle_inits <- c(0, 1, 0)
mle <- c_alg$findMLE(mle_inits)
mle$par

#' 
#' Standard errors:
#' 
## -----------------------------------------------------------------------------
sqrt(diag(solve(-mle$hessian)))

#' 
#' You can see `nimble` struggled to estimate the SD.
#' 
#' Estimates of `beta0`:
#' 
## -----------------------------------------------------------------------------
c_alg$summary(mle)$randomEffects$estimates

#' 
#' # Macros
#' 
#' * Special code chunks that, when processed, turn into new code
#' * Users can create their own macros
#' * Under development
#' 
#' For example:
#' 
#' ```r
#' mu[1:M] <- myDistribution()
#' ```
#' 
#' Could "expand" into valid BUGS code:
#' 
#' ```r
#' for (i in 1:M){
#'   mu[i] ~ dnorm(0, 1)
#' }
#' ```
#' 
#' **A more realistic example**
#' 
## -----------------------------------------------------------------------------
library(nimbleMacros) # not available on CRAN yet

#' 
#' Create a new block of model code with a macro, `linPred`.
#' 
#' The `linPred` macro generates a linear predictor (optionally transformed) using R's formula interface
#' 
## -----------------------------------------------------------------------------
code <- nimbleCode({
  psi[1:M] <- linPred(~forest_scale, link = logit)
  
  for (i in 1:M){
    y[i] ~ dbern(psi[i])    
  }
})

#' 
#' Creating the `nimbleModel` will process the macro code, generating the final valid BUGS code.
#' 
## -----------------------------------------------------------------------------
mod <- nimbleModel(code = code, constants = nimble_data)
mod$getCode()

#' 
