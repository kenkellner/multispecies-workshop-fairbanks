#' ---
#' title: An R Markdown document converted from "occupancy_nosolutions.ipynb"
#' output: html_document
#' ---
#' 
#' # Single-species occupancy models
#' 
#' Ken Kellner
#' 
#' # Outline
#' 
#' 1. Introduction
#' 2. Simulate dataset
#' 3. Fit in `unmarked`
#' 4. Fit in `nimble`
#' 5. Other methods
#' 
#' # Introduction
#' 
#' The probability a species occupies (is present at) a site ($\psi$)
#' 
#' A key consideration of species distribution models (SDMs)
#' 
#' Problem: detection is imperfect (*p* < 1)
#' 
#' Just because we don't detect a species at a site doesn't mean it isn't there
#' 
#' Solution: **Repeated samples**
#' 
#' Visiting the site multiple times allows estimation of detection probability *p* and thus correction of $\psi$
#' 
#' ![](occupancy.png)
#' 
#' ## State (occupancy) process
#' 
#' **Parameters**
#' 
#' $z_i$: Latent occupancy at site $i$
#' 
#' $\psi_i$: Occupancy probability at site $i$
#' 
#' **Math**
#' 
#' $$z_i \sim \mathrm{Bernoulli}(\psi_i)$$
#' 
#' Equivalent to:
#' 
#' $$z_i \sim \mathrm{Binomial}(1, \psi_i)$$
#' 
#' ## Detection process
#' 
#' **Parameters/data**
#' 
#' $y_{ij}$: Observed detection at site $i$ for repeated sample $j$
#' 
#' $p_{ij}$: Probability of detecting an individual at site $i$ in sample $j$
#' 
#' **Math**
#' 
#' $$y_{ij} \sim \mathrm{Bernoulli}(p_{ij} \cdot z_i)$$
#' 
#' ## Key assumptions
#' 
#' * Population is closed during the repeated samples
#' * Repeated samples are independent
#' * No false positive detections (only false negative)
#' 
#' # Simulate a dataset
#' 
#' Example adapted from Kéry and Kellner, *Applied Statistical Modelling for Ecologists* (2024)
#' 
#' ![](book.jpg)
#' 
#' 
#' Occupancy of plant (gentian) at a series of sites
#' 
#' Both occupancy and detection depend on humidity
#' 
#' ![](gentern.png)
#' 
#' **Study design and covariates**
#' 
## -----------------------------------------------------------------------------
set.seed(123)
nsites <- 150
nvisits <- 3
humidity <- runif(nsites, -1, 1)

#' 
#' **Occupancy parameters**
#' 
#' $occ_{int}$: The occupancy intercept
#' 
#' $occ_{hum}$: The effect of humidity on occupancy
#' 
## -----------------------------------------------------------------------------
occ_int <- 0
occ_hum <- 2

#' 
#' **Simulate occupancy at each site**
#' 
#' $$\mu_{occ_i} = \mathrm{occ}_{int} + \mathrm{occ}_{hum} \cdot \mathrm{humidity_i} $$
#' $$ \psi_i = \mathrm{ilogit}(\mu_{occ_i}) $$
#' 
## -----------------------------------------------------------------------------
mu_occ <- occ_int + occ_hum * humidity
psi <- plogis(mu_occ)

#' 
#' Plot the relationship:
#' 
## -----------------------------------------------------------------------------
plot_psi <- data.frame(humidity=humidity, est=psi)
plot_psi <- plot_psi[order(humidity),]

library(ggplot2)
options(repr.plot.width=10, repr.plot.height=10)

ggplot(data=plot_psi, aes(x=humidity, y=est)) +
  theme_bw(base_size=24) +
  theme(panel.grid = element_blank()) +
  geom_line() +
  labs(y = "psi")

#' 
#' Simulate the latent occupancy state:
#' 
#' $$ z_i \sim \mathrm{Bernoulli}(\psi_i) $$
#' 
## -----------------------------------------------------------------------------
z <- rbinom(nsites, 1, p = psi)
z

#' 
#' Plot the raw latent occupancy data:
#' 
## -----------------------------------------------------------------------------
breaks <- quantile(humidity, c(0, 0.2, 0.4, 0.6, 0.8, 1))
hum_cat <- cut(humidity, breaks = breaks)
tab <- aggregate(z ~ hum_cat, FUN=mean)
tab$n <- table(hum_cat)
tab$se <- sqrt(tab[,2] * (1-tab[,2]) / tab$n)

#' 
## -----------------------------------------------------------------------------
ggplot(data = tab) +
  geom_point(aes(x=hum_cat, y = z)) +
  geom_errorbar(aes(x = hum_cat, ymin= z - se, ymax=z + se), width=0.2) +
  labs(y = "True plant occupancy", x = "Humidity") +
  theme_bw(base_size = 24) +
  theme(panel.grid = element_blank())

#' 
#' **Detection parameters**
#' 
#' $det_{int}$: The detection intercept
#' 
#' $det_{hum}$: The effect of humidity on detection
#' 
## -----------------------------------------------------------------------------
det_int <- 0
det_hum <- -3
truth <- c(occ_int = occ_int, occ_hum = occ_hum, det_int = det_int, det_hum = det_hum)

#' 
#' Note that humidity has oppposite effects on occupancy (2) and detection (-3)! A "pathological" example.
#' 
#' **Simulate detection process**
#' 
#' $$\mu_{det_i} = \mathrm{det}_{int} + \mathrm{det}_{hum} \cdot \mathrm{humidity_i} $$
#' $$ p_i = \mathrm{ilogit}(\mu_{det_i}) $$
#' 
## -----------------------------------------------------------------------------
mu_det <- det_int + det_hum * humidity
p <- plogis(mu_det)

#' 
## -----------------------------------------------------------------------------
plot_p <- data.frame(humidity=humidity, est=p)
plot_p <- plot_p[order(humidity),]

plot_psi$param <- "psi"
plot_p$param <- "p"
plot_dat <- rbind(plot_psi, plot_p)

ggplot(data=plot_dat) +
  theme_bw(base_size=24) +
  theme(panel.grid = element_blank()) +
  geom_line(aes(x=humidity, y=est, col=param)) +
  labs(y="Parameter value", col="Parameter")

#' 
#' Simulate observed detection/non-detection:
#' 
#' $$ y_{ij} \sim \mathrm{Bernoulli}(p_{i} \cdot z_i) $$
#' 
## -----------------------------------------------------------------------------
y <- matrix(NA, nsites, nvisits)

for (i in 1:nsites){
  for (j in 1:nvisits){
    y[i,j] <- rbinom(1, 1, p = p[i] * z[i])
  }
}

#' 
## -----------------------------------------------------------------------------
head(y)

#' 
#' **Plot observed detection/non-detection as a function of humidity**
#' 
#' First determine if gentian was observed at least once at each site:
#' 
## -----------------------------------------------------------------------------
y_naive <- apply(y, 1, max) # ever observed

#' 
#' Next bin `y_naive` and plot the same way we did with `z` earlier:
#' 
## -----------------------------------------------------------------------------
breaks <- quantile(humidity, c(0, 0.2, 0.4, 0.6, 0.8, 1))
hum_cat <- cut(humidity, breaks = breaks)
tab <- aggregate(y_naive ~ hum_cat, FUN=mean)
tab$n <- table(hum_cat)
tab$se <- sqrt(tab[,2] * (1-tab[,2]) / tab$n)

#' 
## -----------------------------------------------------------------------------
ggplot(data = tab) +
  geom_point(aes(x=hum_cat, y = y_naive)) +
  geom_errorbar(aes(x = hum_cat, ymin= y_naive - se, ymax=y_naive + se), width=0.2) +
  labs(y = "Naive occupancy", x = "Humidity") +
  theme_bw(base_size = 24) +
  theme(panel.grid = element_blank())

#' 
#' ### Naive analysis
#' 
#' Fit a standard logistic regression, ignoring imperfect detection
#' 
## -----------------------------------------------------------------------------
fit_glm <- glm(y_naive~humidity, family=binomial)
glm_est <- c(coef(fit_glm), NA, NA)
summary(fit_glm)

#' 
## -----------------------------------------------------------------------------
nd <- data.frame(humidity = sort(humidity))
mod_est <- predict(fit_glm, newdata=nd, se.fit=TRUE)

nd$est <- plogis(mod_est$fit)
nd$lower <- plogis(mod_est$fit - 1.96*mod_est$se.fit)
nd$upper <- plogis(mod_est$fit + 1.96*mod_est$se.fit)

ggplot(data=nd, aes(x=humidity, y=est)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
  geom_line() +
  theme_bw(base_size=24) +
  theme(panel.grid=element_blank()) +
  ylim(0, 1) +
  labs(y="naive occupancy")

#' 
## -----------------------------------------------------------------------------
fit_glm2 <- glm(y_naive~humidity+I(humidity^2), family=binomial)
summary(fit_glm2)

#' 
## -----------------------------------------------------------------------------
nd <- data.frame(humidity = sort(humidity))
mod_est <- predict(fit_glm2, newdata=nd, se.fit=TRUE)

nd$est <- plogis(mod_est$fit)
nd$lower <- plogis(mod_est$fit - 1.96*mod_est$se.fit)
nd$upper <- plogis(mod_est$fit + 1.96*mod_est$se.fit)

ggplot(data=nd, aes(x=humidity, y=est)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
  geom_line() +
  theme_bw(base_size=18) +
  theme(panel.grid=element_blank()) +
  ylim(0, 1) +
  labs(y="naive occupancy")

#' 
#' Both analyses capture the wrong pattern (actually measuring $\psi \cdot p$)
#' 
#' # Fit in R: `unmarked`
#' 
#' * A package for fitting a wide variety of occupancy and abundance models, accounting for imperfect detection
#' * Uses maximum likelihood
#' 
## -----------------------------------------------------------------------------
library(unmarked)

#' 
#' ## Organize data
#' 
#' Using an `unmarkedFrame`, analogous to an R data frame
#' 
#' Components:
#' 
#' * observed data `y`
#' * site covariates `siteCovs`
#' * obs covariates `obsCovs`
#' 
## -----------------------------------------------------------------------------
umf <- unmarkedFrameOccu(y = y, siteCovs = data.frame(humidity = humidity))
head(umf)

#' 
## -----------------------------------------------------------------------------
plot(umf)

#' 
#' ## Fit the model
#' 
#' Using the `occu` function
#' 
#' Need to specify a "double" formula: first detection, then occupancy
#' 
#' `~humidity ~humidity`
#' 
## -----------------------------------------------------------------------------
fit_unm <- occu(~humidity ~humidity, data = umf)

#' 
## -----------------------------------------------------------------------------
summary(fit_unm)

#' 
#' ## Model fit
#' 
#' Residuals:
#' 
#' $$ r_{ij} = (y_{ij} - \hat{y}_{ij}) $$
#' 
#' $$\hat{y}_{ij} = \psi_i p_{ij} $$
#' 
## -----------------------------------------------------------------------------
plot(fit_unm)

#' 
#' ### Check model fit with parametric bootstrap
#' 
#' 1. Calculate test statistic using fitted model and real dataset
#' 2. Simulate new datasets based on model
#' 3. Fit models to new datasets
#' 4. Calculate test statistics using new datasets/models
#' 
#' 1. Basic test statistic: sum of squared errors (SSE)
#' 
#' $$ \mathrm{SSE} = \sum{r_{ij}^2} $$
#' 
## -----------------------------------------------------------------------------
SSE(fit_unm)

#' 
#' 2. We can simulate new datasets using `simulate`
#' 
## -----------------------------------------------------------------------------
new_datasets <- simulate(fit_unm, 2)
lapply(new_datasets, head)

#' 
#' 3-4: We could fit the model to each new dataset manually, but `parboot` does it for us:
#' 
## -----------------------------------------------------------------------------
pb <- parboot(fit_unm, nsim = 30)

#' 
## -----------------------------------------------------------------------------
plot(pb)
pb

#' 
#' ### An alternative fit statistic
#' 
#' The MacKenzie-Bailey chi-square (MacKenzie & Bailey, 2004)
#' 
## -----------------------------------------------------------------------------
library(AICcmodavg)

mb.chisq(fit_unm)

#' 
## -----------------------------------------------------------------------------
mb_test <- mb.gof.test(fit_unm, nsim=30)

#' 
#' ## Inference
#' 
## -----------------------------------------------------------------------------
summary(fit_unm)

#' 
#' ### Compare with truth
#' 
## -----------------------------------------------------------------------------
unm_est  <- coef(fit_unm)

comp <- cbind(truth = truth, glm = glm_est, unm = unm_est)
round(comp, 3)

#' 
## -----------------------------------------------------------------------------
plogis(0)
plogis(0.457)

#' 
#' ### Get parameter estimates for each site
#' 
#' Use `predict`!
#' 
#' For $\psi$:
#' 
## -----------------------------------------------------------------------------
psi_site <- predict(fit_unm, type = "state")
head(psi_site)

#' 
#' Or detection ($p$, results are in site-major order):
#' 
## -----------------------------------------------------------------------------
p_all <- predict(fit_unm, type = "det")
dim(p_all)
head(p_all)

#' 
#' ### Plot covariate effects
#' 
#' 1. Manually
#' 2. Using `plotEffects`
#' 
#' **Manually**
#' 
#' Using `predict` and an approach similar to the last module
#' 
#' 1. Make sequence of humidity values
#' 2. Predict $\psi$ for each value of humidity
#' 3. Plot line and error bars
#' 
## -----------------------------------------------------------------------------
hum_rng <- range(humidity)
hum_seq <- seq(hum_rng[1], hum_rng[2], length.out = 100)
nd <- data.frame(humidity = hum_seq)

#' 
## -----------------------------------------------------------------------------
pr <- predict(fit_unm, type = "state", newdata = nd, appendData = TRUE)
head(pr)

#' 
## -----------------------------------------------------------------------------
ggplot(data = pr, aes(x=humidity)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha = 0.2) +
  geom_line(aes(y=Predicted)) +
  theme_bw(base_size=24) + 
  theme(panel.grid = element_blank()) +
  labs(y = "Occupancy and 95% CI")

#' 
#' **Using `plotEffects`**
#' 
## -----------------------------------------------------------------------------
plotEffects(fit_unm, type='state', covariate='humidity')

#' 
## -----------------------------------------------------------------------------
plotEffects(fit_unm, type='det', covariate='humidity')

#' 
#' ### Estimate total number of occupied sites
#' 
#' The `ranef` function: uses empirical Bayes to to estimate posterior distribution of $z$
#' 
## -----------------------------------------------------------------------------
r <- ranef(fit_unm)
round(sum(bup(r)), 2)

#' 
#' ## Model selection
#' 
#' First fit another (null) model
#' 
## -----------------------------------------------------------------------------
fit_null <- occu(~1~1, umf)

#' 
#' Combine the models into a `fitList`
#' 
## -----------------------------------------------------------------------------
fits <- fitList(hum=fit_unm, null=fit_null)

#' 
#' Get AIC model selection table
#' 
## -----------------------------------------------------------------------------
modSel(fits)

#' 
#' ## NIMBLE
#' 
#' Workflow:
#' 
#' 1. Bundle the data
#' 2. Write the BUGS code
#' 3. Set initial values
#' 4. Run with `nimbleMCMC`
#' 5. Diagnostics
#' 5. Inference
#' 
## -----------------------------------------------------------------------------
library(nimble)

#' 
#' ## Bundle the data
#' 
## -----------------------------------------------------------------------------
data_list <- list(nsites = nsites, nvisits = nvisits, y = y, humidity = humidity)
str(data_list)

#' 
#' ## Write the BUGS code
#' 
#' ### State model
#' 
#' **Math**
#' 
#' 
#' $$\mu_{occ_i} = \mathrm{occ}_{int} + \mathrm{occ}_{hum} \cdot \mathrm{humidity_i} $$
#' 
#' $$\psi_i = \mathrm{expit}(\mu_{occ_i})$$
#' 
#' $$z_i \sim \mathrm{Bernoulli}(\psi_i)$$
#' 
#' ```r
#' for (i in 1:nsites){
#'   mu_occ[i] <- occ_int + occ_hum * humidity[i]
#'   psi[i] <- expit(mu_occ[i])
#'   z[i] ~ dbern(psi[i])
#' }
#' ```
#' 
#' ### Detection model
#' 
#' **Math**
#' 
#' $$ \mu_{det_i} = \mathrm{det}_{int} + \mathrm{det}_{hum} \cdot \mathrm{humidity_i} $$
#' 
#' $$ p_i = \mathrm{expit}(\mu_{det_i}) $$
#' 
#' $$y_{ij} \sim \mathrm{Bernoulli}(p_i \cdot z_i)$$
#' 
#' **Code**
#' 
#' ```r
#' for (i in 1:nsites){
#'   mu_det[i] <- det_int + det_hum * humidity[i]
#'   p[i] <- expit(mu_det[i])
#'   
#'    for (j in 1:nvisits){
#'     y[i,j] ~ dbern(p[i] * z[i])
#'   }
#' }
#' ```
#' 
#' ### Number of occupied sites
#' 
#' **Math**
#' 
#' $$z_{sum} = \sum{z}$$
#' 
#' **Code**
#' 
#' ```r
#' zsum <- sum(z[1:nsites])
#' ```
#' 
#' ### Goodness-of-fit
#' 
#' We'll calculate SSE for the real dataset and for simulated datasets, same as with `parboot`.
#' 
#' **Math**
#' 
#' $$ \mathrm{SSE} = \sum\left(y_{ij} - \hat{y}_{ij}\right)^2 $$
#' $$ \hat{y}_{ij} = \psi_i \cdot p_{ij} $$
#' 
#' **Code**
#' 
#' ```r
#' for (i in 1:nsites){
#'   for (j in 1:nvisits){
#'     yhat[i,j] <- psi[i] * p[i]
#'     r2[i,j] <- (y[i,j] - yhat[i,j])^2
#'       
#'     # New dataset
#'     y_new[i,j] ~ dbern(p[i] * z[i])
#'     r2_new[i,j] <- (y_new[i,j] - yhat[i,j])^2
#'   }    
#' }
#' SSE <- sum(r2[1:nsites, 1:nvisits])
#' SSE_new <- sum(r2_new[1:nsites, 1:nvisits])
#' ```
#' 
## -----------------------------------------------------------------------------
code_occu <- nimbleCode({
  # State model
  for (i in 1:nsites){
    mu_occ[i] <- occ_int + occ_hum * humidity[i]
    psi[i] <- expit(mu_occ[i])
    z[i] ~ dbern(psi[i])
  }
  
  # Detection model
  for (i in 1:nsites){
    mu_det[i] <- det_int + det_hum * humidity[i]
    p[i] <- expit(mu_det[i])
    
    for (j in 1:nvisits){
      y[i,j] ~ dbern(p[i] * z[i])
    }
  }

  # Number of occupied sites
  zsum <- sum(z[1:nsites])

  # Priors
  occ_int ~ dunif(-5, 5)
  occ_hum ~ dnorm(0, sd = 10)
  
  det_int ~ dunif(-5, 5)
  det_hum ~ dnorm(0, sd = 10)
    
  # GOF
  for (i in 1:nsites){
    for (j in 1:nvisits){
      yhat[i,j] <- psi[i] * p[i]
      r2[i,j] <- (y[i,j] - yhat[i,j])^2
      
      # New dataset
      y_new[i,j] ~ dbern(p[i] * z[i])
      r2_new[i,j] <- (y_new[i,j] - yhat[i,j])^2
    }    
  }
  SSE <- sum(r2[1:nsites, 1:nvisits])
  SSE_new <- sum(r2_new[1:nsites, 1:nvisits])
})

#' 
#' Quick test run:
#' 
## -----------------------------------------------------------------------------
fit_nim <- nimbleMCMC(code_occu, constants = data_list,
                       nchain = 4, niter = 300, nburnin = 100, samplesAsCodaMCMC = TRUE)

#' 
#' ### Initial values 
#' 
#' **Issue**: need to initialize `z` at some reasonable values
#' 
#' Initialize at the maximum observation (i.e., `z = 1` if species was observed at least once)
#' 
## -----------------------------------------------------------------------------
inits <- function() list(z = apply(y, 1, max)) # same as y_naive

#' 
#' ### Fit the model
#' 
## -----------------------------------------------------------------------------
fit_nim <- nimbleMCMC(code_occu, 
                      constants = data_list, 
                      inits = inits,
                      monitors = c("occ_int","occ_hum","det_int","det_hum","zsum", "SSE", "SSE_new"),
                      nchain = 4, 
                      niter = 2000, 
                      nburnin = 1000, 
                      samplesAsCodaMCMC = TRUE)

#' 
## -----------------------------------------------------------------------------
lapply(fit_nim, head)

#' 
#' ### Diagnostics
#' 
#' **Traceplots:**
#' 
## -----------------------------------------------------------------------------
par(mfrow=c(2,2))
coda::traceplot(fit_nim)

#' 
#' **Rhat:**
#' 
## -----------------------------------------------------------------------------
coda::gelman.diag(fit_nim)

#' 
#' **Effective samples:**
#' 
## -----------------------------------------------------------------------------
round(coda::effectiveSize(fit_nim))

#' 
#' **Posterior predictive check**
#' 
#' Compare SSE for real and simulated datasets
#' 
#' Similar to what `parboot` did
#' 
## -----------------------------------------------------------------------------
samples <- as.data.frame(as.matrix(fit_nim))
head(samples)

#' 
## -----------------------------------------------------------------------------
ggplot(data = samples) +
  geom_point(aes(x=SSE, y=SSE_new)) +
  geom_abline(col = "red") +
  theme_bw(base_size = 24) +
  theme(panel.grid = element_blank()) +
  labs(x = "SSE real", y = "SSE simulated")

#' 
#' **Bayesian p-value**
#' 
#' * Essentially the proportion of points above/below the line
#' * Ideally should be ~0.5
#' * Values close to 0 or 1 indicate poor fit
#' 
## -----------------------------------------------------------------------------
mean(samples$SSE > samples$SSE_new)

#' 
#' ## Inference
#' 
## -----------------------------------------------------------------------------
summary(fit_nim)

#' 
#' Compare estimates:
#' 
## -----------------------------------------------------------------------------
samps <- as.matrix(fit_nim)[,c("occ_int","occ_hum","det_int","det_hum")]
nim_est <- apply(samps, 2, mean)

#' 
## -----------------------------------------------------------------------------
comp <- cbind(comp, nimble=nim_est)
round(comp, 3)

#' 
#' **Make a plot of humidity effect on occupancy**
#' 
#' **Math**
#' 
#' 
#' $$\mu_{occ_i} = \mathrm{occ}_{int} + \mathrm{occ}_{hum} \cdot \mathrm{humidity_i} $$
#' 
#' $$\psi_i = \mathrm{ilogit}(\mu_{occ_i})$$
#' 
#' **Code**
#' 
#' ```r
#' for (i in 1:n_hum_values){
#'   mu_occ[i] <- occ_int + occ_hum * hum_values[i]
#'   psi[i] <- plogis(mu_occ[i])
#' }
#' ```
#' 
#' **Values of humidity**
#' 
## -----------------------------------------------------------------------------
hum_seq <- seq(min(humidity), max(humidity), length.out=100)

#' 
#' **Posteriors of `occ_int` and `occ_hum`**
#' 
## -----------------------------------------------------------------------------
samps <- as.matrix(fit_nim)
head(samps)

#' 
#' **Calculate posterior of `mu` for each value of humidity**
#' 
## -----------------------------------------------------------------------------
mu_post <- matrix(NA, nrow=nrow(samps), ncol = 100)

for (i in 1:100){
  mu_post[,i] <- samps[,"occ_int"] + samps[,"occ_hum"] * hum_seq[i] 
}

#' 
## -----------------------------------------------------------------------------
psi_post <- plogis(mu_post)

#' 
#' **Calculate mean psi and 95% CI at each value of humidity**
#' 
## -----------------------------------------------------------------------------
psi_mean <- apply(psi_post, 2, mean)
psi_lower <- apply(psi_post, 2, quantile, 0.025)
psi_upper <- apply(psi_post, 2, quantile, 0.975)

#' 
#' **Make plot**
#' 
## -----------------------------------------------------------------------------
plot_dat <- data.frame(humidity = hum_seq, est=psi_mean, lower=psi_lower, upper=psi_upper)

ggplot(data=plot_dat, aes(x=humidity)) +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
  geom_line(aes(y=est)) +
  theme_bw(base_size=24) +
  theme(panel.grid=element_blank()) +
  labs(y = "Occupancy and 95% CRI")

#' 
#' ### Exercise: modifying the model
#' 
#' * Simulate a new site-level covariate and add it to the model
#' * Plot its effect if you have time
#' 
#' # Bonus content: Occupancy macro in `nimble`
#' 
#' Macros are an in-development feature in `nimble`
#' 
#' They are "shortcut" code chunks which are expanded into larger blocks of model code.
#' 
## -----------------------------------------------------------------------------
library(nimbleEcology) # development version
nimbleOptions(enableModelMacros = TRUE)

#' 
## -----------------------------------------------------------------------------
code <- nimbleCode({
  y[1:nsites, 1:nvisits] ~ occupancy(~humidity[1:nsites], ~humidity[1:nsites])
})

#' 
## -----------------------------------------------------------------------------
mod <- nimbleModel(code, constants = data_list)
mod$getCode()

#' 
#' # Bonus content: `ubms`
#' 
#' Fit Bayesian occupancy models using samed syntax as `unmarked`.
#' 
#' Just need to change `occu` to `stan_occu`
#' 
## -----------------------------------------------------------------------------
library(ubms)
fit_ubms <- stan_occu(~humidity ~humidity, data = umf)

#' 
## -----------------------------------------------------------------------------
fit_ubms

#' 
## -----------------------------------------------------------------------------
traceplot(fit_ubms)

#' 
#' **Residuals**
#' 
#' Separated by occupancy and detection process
#' 
## -----------------------------------------------------------------------------
plot(fit_ubms)

#' 
#' **Goodness of fit with posterior predictive check**
#' 
#' * Uses the specialized MacKenzie-Bailey chi-square test for occupancy models
#' * This test compares the expected frequencies of encounter histories given the model ([001] [110] and so on) with what we actually observed
#' * We calculate MB statistic for real dataset and for simulated datasets and compare
#' 
## -----------------------------------------------------------------------------
g <- gof(fit_ubms, draws = 1000)
plot(g)

#' 
#' Extract full posterior distributions:
#' 
## -----------------------------------------------------------------------------
post <- extract(fit_ubms)
lapply(post, head)

#' 
## -----------------------------------------------------------------------------
plot_posteriors(fit_ubms)

#' 
#' Plot effects:
#' 
## -----------------------------------------------------------------------------
plot_effects(fit_ubms, 'state')

#' 
## -----------------------------------------------------------------------------
plot_effects(fit_ubms, 'det')

#' 
#' **Estimate total occupied sites**
#' 
## -----------------------------------------------------------------------------
zpost <- posterior_predict(fit_ubms, "z")

#' 
## -----------------------------------------------------------------------------
zpost_sum <- apply(zpost, 1, sum)
mean(zpost_sum)
hist(zpost_sum)

#' 
#' **Information criterion**
#' 
## -----------------------------------------------------------------------------
waic(fit_ubms)

#' 
## -----------------------------------------------------------------------------
loo(fit_ubms) # approximation based on single fit

#' 
## -----------------------------------------------------------------------------
kfold(fit_ubms, K=5) # actual real cross-validation

#' 
#' **Model selection**
#' 
## -----------------------------------------------------------------------------
ubms_null <- stan_occu(~1~1, umf)

#' 
## -----------------------------------------------------------------------------
fl <- fitList(null=ubms_null, hum=fit_ubms)

#' 
#' Use LOO-CV (approximation) for model selection
#' 
## -----------------------------------------------------------------------------
modSel(fl)

#' 
#' `ubms` supports spatial random effects via RSR (restricted spatial regression) - see vignette.
#' 
#' Another very useful package, particularly for multispecies/spatial analyses: `spOccupancy`
#' 
#' # Bonus content: "Do it yourself" Maximum Likelihood
#' 
#' We need to calculate the marginal or integrated likelihood (i.e., "integrating out" the latent occupancy states at each site $z$)
#' 
#' **Two possible situations we observe:**
#' 
#' 1. The species was detected at least once at site $i$
#' 2. The species was never detected at site $i$
#' 
#' **If the species was detected at least once at site $i$:**
#' 
#' We know the species is present at the site ($z = 1$).
#' 
#' The likelihood is the product of the probability of occupancy * likelihood of detection process
#' 
#' $$L_i = \psi \cdot \mathrm{dbinom}(n_i, J, p)$$
#' 
#' Here $n_i$ is the number of times the species was detected in J samples
#' 
#' For example an encounter history [01010] would have $J = 5$ and $n = 2$
#' 
#' Or if $p$ varies with detection period $j$, we split the binomial into a series of Bernoullis and multiply them:
#' 
#' $$L_i = \psi \cdot \prod_{j=1}^J\mathrm{dbinom}(y_{ij}, 1, p_i)$$
#' 
#' **If the species was not detected at site $i$:**
#' 
#' Either 
#' 
#' 1. $z = 1$ (species present but undetected) 
#' 2. $z = 0$ (species not present)
#' 
#' We calculate the likelhood for each possibility and add them together.
#' 
#' 1. $\psi \cdot \mathrm{dbinom}(n_i, J, p)$
#' 2. $(1 - \psi)$
#' 
#' $$ L_i = \psi \cdot \mathrm{dbinom}(n_i, J, p) + (1-\psi)$$
#' 
#' **In summary:**
#' 
#' $$
#' L_i = \begin{cases}
#' \psi \cdot \mathrm{dbinom}(n_i, J, p), & any\ detection \\
#' \psi \cdot \mathrm{dbinom}(n_i, J, p) + (1-\psi), & no detection \\
#' \end{cases}
#' $$
#' 
#' or
#' 
#' $$L_i = \psi \cdot \mathrm{dbinom}(n_i, J, p) + (1-\psi) \cdot (\mathrm{no\_detects} == 1)$$
#' 
#' Need to add new data ($n$ and $\mathrm{no\_detects}$) to our list
#' 
## -----------------------------------------------------------------------------
data_list$n <- apply(y, 1, sum)
data_list$no_detects <- 1 - apply(y, 1, max)
data_list$y <- NULL
str(data_list)

#' 
#' **Math:**
#' 
#' $$L_i = \psi \cdot \mathrm{dbinom}(n_i, J, p) + (1-\psi) \cdot (\mathrm{no\_detects} == 1)$$
#' 
#' **Code:**
#' 
#' ```r
#' lik[i] <- psi[i] * dbinom(data$n[i], data$nvisits, p[i]) +
#'           (1 - psi[i]) * data$no_detects[i]
#' ```
#' 
## -----------------------------------------------------------------------------
NLL <- function(pars, data){
  occ_int <- pars[1]
  occ_hum <- pars[2]
  det_int <- pars[3]
  det_hum <- pars[4]
  
  # Calculate occupancy probability
  mu_occ <- occ_int + occ_hum * data$humidity
  psi <- plogis(mu_occ)
  
  # Calculate detection probability
  mu_det <- det_int + det_hum * data$humidity
  p <- plogis(mu_det)
    
  lik <- rep(0, data$nsites)
  for (i in 1:data$nsites){
    lik[i] <- psi[i] * dbinom(data$n[i], data$nvisits, p[i]) +
      (1 - psi[i]) * data$no_detects[i]
  }
  ll <- log(lik)
  return(-sum(ll))
}

#' 
## -----------------------------------------------------------------------------
par_start <- c(occ_int = 0, occ_hum = 0, det_int = 0, det_hum = 0)
occ_diy <- optim(par_start, fn = NLL, hessian = TRUE, data = data_list)

#' 
## -----------------------------------------------------------------------------
occ_diy$par

#' 
