#' ---
#' title: An R Markdown document converted from "multiocc-unmarked_nosolutions.ipynb"
#' output: html_document
#' ---
#' 
#' # Multispecies occupancy models in `unmarked`
#' 
#' ### Ken Kellner
#' 
#' # Outline
#' 
#' 1. Model review
#' 2. Fit model with `unmarked`
#' 3. Marginal and conditional occupancy
#' 4. Model fitting challenges
#' 
#' # Model review
#' 
#' * Suppose you sample for $S = 2$ species at N sites, $i=1...N$
#' * Then the latent state $\textbf{Z}_i$ at site $i$ is a vector of length 2, and there are four possible states:
#' 
#' $$ \textbf{Z}_i = \{11\}, \{10\},\{01\}, or \{00\} $$
#' 
#' $$ \textbf{Z}_i \sim \mathrm{MultivariateBernoulli}(\boldsymbol\psi) $$
#' 
#' $$ pr(Z_i) = \{11\} \Longrightarrow \psi_{11} $$
#' $$ pr(Z_i) = \{10\} \Longrightarrow \psi_{10} $$
#' $$ pr(Z_i) = \{01\} \Longrightarrow \psi_{01} $$
#' $$ pr(Z_i) = \{00\} \Longrightarrow \psi_{00} $$
#' 
#' $$ \sum{\boldsymbol\psi} = 1$$
#' 
#' ### How do we get the values of $\psi$?
#' 
#' ### Start with the "natural parameters", $f$
#' 
#' $f_1$: Term associated with occupancy of species 1 by itself
#' 
#' $f_2$: Term associated with occupancy of species 2 by itself
#' 
#' $f_{12}$: Interaction term between the two species
#' 
#' If $f_{12}$ is positive, presence of one species increases occupancy probability of the other
#' 
#' If $f_{12}$ is negative, presence of one species decreases occupancy probability of the other
#' 
#' ### Calculating $\psi$: the multinomial logit link
#' 
#' Each element of $\mathbf{\psi}$ is a function of the natural parameters $f$.
#' 
#' $$ 
#' \begin{align}
#' \psi_{11} & \propto exp(f_{1} + f_{2} + f_{12}) \\
#' \psi_{10} & \propto exp(f_{1}) \\
#' \psi_{01} & \propto exp(f_{2}) \\
#' \psi_{00} & \propto exp(0) = 1
#' \end{align}
#' $$
#' 
#' Making sure $\mathbf{\psi}$ sums to 1:
#' 
#' $$ \psi_{11} = \frac{exp(f_{1} + f_{2} + f_{12})}{1 + exp(f_{1}) + exp(f_{2}) + exp(f_{1} + f_{2} + f_{12})} $$
#' 
#' $$ \psi_{10} = \frac{exp(f_{1})}{1 + exp(f_{1}) + exp(f_{2}) + exp(f_{1} + f_{2} + f_{12})} $$
#' 
#' $$ \psi_{01} = \frac{exp(f_{2})}{1 + exp(f_{1}) + exp(f_{2}) + exp(f_{1} + f_{2} + f_{12})} $$
#' 
#' $$ \psi_{00} = \frac{1}{1 + exp(f_{1}) + exp(f_{2}) + exp(f_{1} + f_{2} + f_{12})} $$
#' 
#' ### A note on removing interaction terms
#' 
#' Original math:
#' 
#' $$ 
#' \begin{align}
#' \psi_{11} & \propto exp(f_{1} + f_{2} + f_{12}) \\
#' \psi_{10} & \propto exp(f_{1}) \\
#' \psi_{01} & \propto exp(f_{2}) \\
#' \psi_{00} & \propto 1
#' \end{align}
#' $$
#' 
#' What if we force $f_{12} = 0$?
#' 
#' $$ 
#' \begin{align}
#' \psi_{11} & \propto exp(f_{1} + f_{2}) \\
#' \psi_{10} & \propto exp(f_{1}) \\
#' \psi_{01} & \propto exp(f_{2}) \\
#' \psi_{00} & \propto 1
#' \end{align}
#' $$
#' 
#' This is equivalent to fitting two separate occupancy models. I'll show this later.
#' 
#' ### Modeling detection
#' 
#' * It's the same as single-species occupancy
#' * Each species gets its own separate detection model
#' 
#' For example, detection/nondetection of species 1 at site $i$ and occasion $j$ can be modeled as
#' 
#' $$ y_{1,ij} \sim \mathrm{Bernoulli}(p_{1} \cdot z_{1,i}) $$
#' 
#' where $z_{1,i}$ is the element 1 of $Z_i$ (i.e., latent presence/absence of species 1 at site $i$) and $p_{1}$ is the species 1 detection probability
#' 
#' # Multispecies occupancy model in `unmarked`
#' 
## -----------------------------------------------------------------------------
library(unmarked)

#' 
#' We will use the function `occuMulti`
#' 
#' ![](occuMulti_setup.png)
#' 
#' ## Set up the input data
#' 
#' ### Example data
#' 
## -----------------------------------------------------------------------------
data(MesoCarnivores)
lapply(MesoCarnivores, head) # look at raw data

#' 
#' ### Formatting the data
#' 
#' We will create an `unmarkedFrameOccuMulti`. 
#' 
#' The main difference from `unmarkedFrameOccu` is that you have >1 matrix of observations: we now have one `y` per species.
#' 
#' ![](umf_om_construction.png)
#' 
#' Get some more detail by looking at the help file for `unmarkedFrameOccuMulti`:
#' 
## -----------------------------------------------------------------------------
?unmarkedFrameOccuMulti

#' 
#' ### y list
#' 
#' The first three elements of `MesoCarnivores` are the y-matrices. To keep things simple, we'll just use coyote and red fox.
#' 
#' Notice that the resulting list is **named**: this is important!
#' 
## -----------------------------------------------------------------------------
ylist <- MesoCarnivores[c(2:3)]
lapply(ylist, head) # look at first few rows

#' 
#' ### Site covariates
#' 
#' The last element of the `MesoCarnivores` list is a data frame of site covariates:
#' 
## -----------------------------------------------------------------------------
site_covs <- MesoCarnivores$sitecovs
head(site_covs)

#' 
#' ### Construct the unmarked frame
#' 
## -----------------------------------------------------------------------------
umf <- unmarkedFrameOccuMulti(y=ylist, siteCovs=site_covs)
head(umf)

#' 
## -----------------------------------------------------------------------------
summary(umf)

#' 
## -----------------------------------------------------------------------------
plot(umf)

#' 
#' ## Set up the formulas
#' 
#' ![](occuMulti_setup.png)
#' 
#' ### `stateformulas`
#' 
#' While `occu` had a single formula for occupancy, `occuMulti` requires one formula per natural parameter $f$, organized into a `stateformulas` vector.
#' 
#' It can be hard to keep track of how many natural parameters there are and what each one represents.
#' 
#' It can be helpful to look at the $f$-design matrix, which is generated by `unmarkedFrameOccuMulti`.
#' 
## -----------------------------------------------------------------------------
umf@fDesign

#' 
#' The number and order of the formulas in the `stateformulas` vector should match the column names of this matrix.
#' Therefore we'll need 3 formulas total.
#' For this model we'll set the 1st-order $f$ parameters to be a function of (standardized) housing density. The interaction term will be intercept-only.
#' 
#' Our `stateformulas` should look like this:
#' 
## -----------------------------------------------------------------------------
stateformulas <- c("~scale(HDens_5km)","~scale(HDens_5km)","~1")
stateformulas

#' 
#' Notice that the formulas are **strings** (each wrapped in `""`). This is required and is just working around some R limitations.
#' 
#' ### `detformulas`
#' 
#' This one is easier, there is just one formula per species, so there should be 2 total.
#' 
#' Intercept-only models for both species:
#' 
## -----------------------------------------------------------------------------
detformulas <- c("~1","~1")
detformulas

#' 
#' ## Run `occuMulti`
#' 
#' We now have all the pieces we need (`unmarkedFrameOccuMulti`, `stateformulas`, `detformulas`) needed to run a model.
#' 
## -----------------------------------------------------------------------------
?occuMulti

#' 
## -----------------------------------------------------------------------------
mod_hdens <- occuMulti(detformulas=detformulas, stateformulas=stateformulas, data=umf)
mod_hdens

#' 
#' ## Model goodness-of-fit
#' 
#' Fit statistic: sum of squared errors (SSE) across all species
#' 
#' $$ SSE = \sum_{s=1}^S \sum_{i=1}^{M\cdot J} \left(y_{i} - \hat{y}_{i}\right)^2 $$
#' 
#' $$ \hat{y}_i = \psi_{marg_i} \cdot p_i $$
#' 
#' $\psi_{marg_i}$ is the marginal occupancy for the species (I'll talk about it more in a minute)
#' 
#' Run a parametric bootstrap (takes a while):
#' 
## -----------------------------------------------------------------------------
pb <- parboot(mod_hdens, nsim=30)
plot(pb)

#' 
#' Dotted line (real dataset) should fall within the histogram (simulated datasets).
#' 
#' There isn't a function available for the MB test for multispecies models (yet?)
#' 
#' ## Inference
#' 
## -----------------------------------------------------------------------------
mod_hdens

#' 
#' * Coyote occupancy is higher than red fox
#' * No relationship between coyote occupancy and housing density
#' * Positive relationship between red fox and housing density
#' * Positive interaction between coyote and fox
#' 
#' ### Exercise
#' 
#' * Add `Dist_5km` (disturbance) as a site-level covariate on one or more occupancy parameters
#' * Add `Trail` as a covariate on detection probability for both species
#' * Name the model object `mod_dist`
#' * If you have time, check goodness-of-fit with `parboot`
#' 
## -----------------------------------------------------------------------------
head(umf)

#' 
#' ### Occupancy probabilities
#' 
#' To get the expected probability $\boldsymbol\psi$ for each occupancy state at each site, use `predict`.
#' This gives you the probabilities along with standard errors and a 95% CI.
#' 
## -----------------------------------------------------------------------------
occ_prob <- predict(mod_hdens, type="state")
lapply(occ_prob, head)

#' 
#' Rows should sum to 1:
#' 
## -----------------------------------------------------------------------------
head(apply(occ_prob$Predicted, 1, sum))

#' 
#' ### Marginal occupancy
#' 
#' You might want to know the *marginal* occupancy of species 2 (red fox) at a site, which is the probability of red fox occupancy regardless of the occupancy states of other species.
#' 
#' This is calculated by summing all $\psi$ states in which red fox is present:
#' 
#' $$ \psi_{marg,i2} = \psi_{11}+\psi_{01} $$
#' 
#' Calculate manually:
#' 
## -----------------------------------------------------------------------------
head(occ_prob$Predicted)

#' 
## -----------------------------------------------------------------------------
marg_manual <- occ_prob$Predicted[,1] + occ_prob$Predicted[,3]
head(marg_manual)

#' 
#' You can also do this by specifying the `species` argument in `predict`.
#' 
## -----------------------------------------------------------------------------
redfox_marginal <- predict(mod_hdens, type="state", species="redfox")
head(redfox_marginal)

#' 
#' #### Plotting marginal occupancy
#' 
#' Compare across species (at site 1 only) with a plot:
#' 
## -----------------------------------------------------------------------------
coy_marginal <- predict(mod_hdens, type='state', species="coyote") # get coyote

marg_plot_dat <- rbind(redfox_marginal[1,], coy_marginal[1,])
marg_plot_dat$Species <- c("Red fox", "Coyote")
marg_plot_dat

#' 
## -----------------------------------------------------------------------------
library(ggplot2)
options(repr.plot.width=10, repr.plot.height=10)

ggplot(data = marg_plot_dat, aes(x=Species, col=Species)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width = 0.2) +
  geom_point(aes(y=Predicted), size=3) +
  theme_bw(base_size=24) +
  theme(panel.grid = element_blank(), legend.position='none') +
  labs(y = "Marginal occupancy and 95% CI", x = "Species")

#' 
#' ### Conditional occupancy
#' 
#' Alternatively, you might want to know the probability of occupancy of one species, conditional on the presence of another.
#' 
#' For example, the probability of red fox occupancy, conditional on coyote presence:
#' 
#' $$ \psi_{i, fox, coy} = \frac{\psi_{i{11}}}{\psi_{i{10}} + \psi_{i{11}}}$$
#' 
## -----------------------------------------------------------------------------
lapply(occ_prob, head)

#' 
## -----------------------------------------------------------------------------
cond_manual <- occ_prob$Predicted[,1] / (occ_prob$Predicted[,2] + occ_prob$Predicted[,1])
head(cond_manual)

#' 
#' 
#' And on coyote absence:
#' 
#' $$ \psi_{i, fox, -coy} = \frac{\psi_{i{01}}}{\psi_{i{01}} + \psi_{i{00}}}$$
#' 
## -----------------------------------------------------------------------------
cond_manual <- occ_prob$Predicted[,3] / (occ_prob$Predicted[,3] + occ_prob$Predicted[,4])
head(cond_manual)

#' 
#' With `predict`, use the `species` and `cond` arguments together for this.
#' 
## -----------------------------------------------------------------------------
redfox_coy <- predict(mod_hdens, type="state", species="redfox", cond="coyote")
head(redfox_coy)

#' 
#' What about conditional on coyote *absence*?
#' 
#' Use `-`
#' 
## -----------------------------------------------------------------------------
redfox_nocoy <- predict(mod_hdens, type="state", species="redfox", cond="-coyote")
head(redfox_nocoy)

#' 
#' #### Plotting conditional occupancy
#' 
#' You can use this output from `predict` to generate comparison plots.
#' 
#' Here we'll compare conditional occupancy for just site 1.
#' 
## -----------------------------------------------------------------------------
plot_data <- rbind(redfox_coy[1,], redfox_nocoy[1,])
plot_data$Coyote_status <- c("Present","Absent")
head(plot_data)

#' 
## -----------------------------------------------------------------------------
ggplot(data = plot_data, aes(x = Coyote_status, col=Coyote_status)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  geom_point(aes(y = Predicted), size=3) +
  theme_bw(base_size = 24) +
  theme(panel.grid = element_blank(), legend.position='none') +
  labs(y = "Fox conditional occupancy and 95% CI", x = "Coyote status")

#' 
#' ### Co-occurence probability
#' 
#' Probability that *both* species 1 and 2 are at a site (or some other arbitrary combo of species)
#' 
#' $$ \psi_{i, fox+coy} = \frac{\psi_{i{11}}}{\psi_{i{11}}+\psi_{i{10}} + \psi_{i{01}} + \psi_{i{00}}}$$
#' 
#' (Trivial for only 2 species)
#' 
## -----------------------------------------------------------------------------
colnames(occ_prob$Predicted)

#' 
## -----------------------------------------------------------------------------
cooccur_manual <- occ_prob$Predicted[,1] / rowSums(occ_prob$Predicted)
head(cooccur_manual)

#' 
#' Or with predict:
#' 
## -----------------------------------------------------------------------------
psi_coccur <- predict(mod_hdens, type="state", species=c("coyote","redfox"))
head(psi_coccur)

#' 
#' ### Plotting covariate effects
#' 
#' To plot the effect of housing density on **marginal** occupancy, we again use `predict`.
#' 
#' 1. Generate sequence of housing density values
#' 2. Predict for each value
#' 3. Plot 95% CI and estimate line
#' 
#' 1. Generate sequence of possible `Hdens_5km` values for X-axis and put in data frame
#' 
#' NOTE: Even though we scaled this covariate in the analysis, we don't need to manually scale here. 
#' 
#' `unmarked` will do it for us because the scaling was specified *in the formula*.
#' 
## -----------------------------------------------------------------------------
hdens_range <- range(siteCovs(umf)$HDens_5km)
hdens_range
hdens_seq <- seq(hdens_range[1], hdens_range[2], length.out=100)
nd <- data.frame(HDens_5km = hdens_seq)

#' 
#' 2. `predict` occupancy at each value of `Hdens_5km` along our sequence. Supply `nd` to `newdata` argument
#' 
#' We do this for both coyote and red fox separately:
#' 
## -----------------------------------------------------------------------------
occ_hdens_coy <- predict(mod_hdens, type="state", species="coyote", newdata=nd)
occ_hdens_coy$Species <- "Coyote"

occ_hdens_fox <- predict(mod_hdens, type="state", species="redfox", newdata=nd)
occ_hdens_fox$Species <- "Red fox"

#' 
#' Then combine the two pieces:
#' 
## -----------------------------------------------------------------------------
occ_hdens <- rbind(occ_hdens_coy, occ_hdens_fox)
occ_hdens$Hdens_5km = hdens_seq
head(occ_hdens)

#' 
#' 3. Make the plot
#' 
## -----------------------------------------------------------------------------
ggplot(data = occ_hdens, aes(x = Hdens_5km)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Species), alpha = 0.2) +
  geom_line(aes(y = Predicted, col = Species)) +
  theme_bw(base_size = 24) +
  theme(panel.grid = element_blank()) +
  labs(y = "Marginal occupancy and 95% CI", x = "Housing density")

#' 
#' ### Exercise: Plot conditional occupancy for a covariate
#' 
#' * Plot conditional occupancy of red fox as a function of housing density
#' * Do this for coyote present, coyote absent, or both (if you have time)
#' * Should be able to adapt code from previous section on marginal occupancy, just add the `cond` argument.
#' 
#' # Removing interaction terms
#' 
#' (or "forcing interactions to 0")
#' 
#' This is pretty common. The two main reasons:
#' 
#' * 3-way and higher interactions are tough to interpret and rarely part of hypotheses
#' * Higher-order interactions are hard to estimate
#' 
#' Two options to remove interactions:
#' 
#' * Set corresponding formula to `"0"`
#' * Set `maxOrder` argument
#' 
#' Setting a formula to `"0"` removes just that particular interaction term.
#' 
#' Review the $f$ design matrix:
#' 
## -----------------------------------------------------------------------------
umf@fDesign

#' 
## -----------------------------------------------------------------------------
stateformulas_0 <- c("~scale(HDens_5km)", "~scale(HDens_5km)", "0")
mod_noint <- occuMulti(detformulas=detformulas, stateformulas=stateformulas_0, data=umf)
mod_noint

#' 
#' Setting `maxOrder` gets rid of all interactions of above that level.
#' 
#' For example, `maxOrder = 1` means no interactions at all.
#' 
## -----------------------------------------------------------------------------
stateformulas <- c("~scale(HDens_5km)", "~scale(HDens_5km)")
mod_noint <- occuMulti(detformulas=detformulas, stateformulas=stateformulas, data=umf, maxOrder=1)
mod_noint

#' 
#' This is equivalent to fitting two separate occupancy models.
#' 
#' For example, the coyote single-species model:
#' 
## -----------------------------------------------------------------------------
umf_coy <- unmarkedFrameOccu(y=ylist$coyote, siteCovs=site_covs)
occu(~1~scale(HDens_5km), umf_coy)

#' 
#' ### Exercise
#' 
#' * Set up a 3-species model with the `MesoCarnivores` dataset (i.e., add bobcat)
#' * Add more covariates if you want
#' * Don't include a 3-species interaction term (i.e., force it to be 0)
#' 
#' # Model fitting challenges
#' 
#' 1. Amount of data required
#' 2. Poor estimates
#' 
#' ## Amount of data required
#' 
#' * You need a surprising amount of data to have high power to detect interaction effects
#' * We can explore this with simulation + power analysis
#' 
#' ### Quick and very dirty review of power analysis
#' 
#' * Type I error: We reject the null hypothesis when in fact it is true.
#' * Type II error: We fail to reject the null hypothesis when in fact it is false.
#' 
#' * Power: 1 - Type II error (rule-of-thumb we want > 0.8)
#' * Power analysis: Determining power for a particular experimental design, effect size, and Type I error ($\alpha$)
#' 
#' We need to use simulation for doing power analysis with this model.
#' 
#' ### Power analysis in `unmarked`
#' 
#' Information we need:
#' 
#' * A set of hypothetical effect sizes (= parameter values)
#' * A "template" model that contains all the information about the experimental design
#' * Value of $\alpha$ (we'll leave at 0.05)
#' 
#' `unmarked` then automatically:
#' 
#' * Plugs effect sizes (= parameter values) into the template model
#' * Simulates new datasets and re-fits the models (many times)
#' * Checks parameter estimates for statistical significance at $\alpha$
#' * Summarizes the results as power
#' 
#' #### Effect sizes
#' 
#' We define a positive interaction between two species: when species 2 is present, species 1 is twice as likely to occupy a site.
#' 
#' First define $f$ values:
#' 
## -----------------------------------------------------------------------------
set.seed(123)

f1 <- 0.2
f2 <- 0.1
f12 <- 0.5 # <- EFFECT SIZE OF INTERACTION

#' 
#' For illustration, calculate the corresponding occupancy probabilities:
#' 
## -----------------------------------------------------------------------------
psi <- numeric(4)
psi[1] <- exp(f1 + f2 + f12)
psi[2] <- exp(f1)
psi[3] <- exp(f2)
psi[4] <- 1
psi <- psi / sum(psi)
names(psi) <- c("[11]","[10]","[01]","[00]")
round(psi, 3)

#' 
#' Note probability of both species present [11] is twice as high as either species by itself.
#' 
#' This is a large interaction effect!
#' 
#' #### Simulate a template dataset
#' 
#' Basically, a placeholder dataset with the correct design (number of sites, species, etc)
#' 
#' We can use a more flexible version of the `unmarked` `simulate` method which can simulate datasets from scratch and requires three specific inputs.
#' 
#' 1. Model structure / covariates (with formulas)
#' 
#' Two species and no covariates in this example
#' 
## -----------------------------------------------------------------------------
formulas <- list(
  state = c("~1", "~1", "~1"),
  det = c("~1", "~1")
)

#' 
#' 2. Experimental design
#' 
#' Number of sites (`M`) and occasions (`J`)
#' 
## -----------------------------------------------------------------------------
design <- list(M = 700, J = 5)

#' 
#' 3. Coefficient values for state and det models
#' 
#' The coefficient values don't actually matter for this "template" dataset (as they are later replaced with the effect sizes we choose), but we'll use the ones we specified above
#' 
## -----------------------------------------------------------------------------
# Coefficient values
coef <- list(state = c(f1, f2, f12), det = c(0, 0))

#' 
#' Finally, simulate the dataset:
#' 
## -----------------------------------------------------------------------------
umf_temp <- simulate("occuMulti", formula = formulas, design = design, coef = coef)
head(umf_temp)

#' 
#' ### Create template model
#' 
#' Fit a model using our template dataset.
#' 
#' Remember, the actual estimates from this model don't matter, just the structure of the model (sample size, covariates included, etc.)
#' 
## -----------------------------------------------------------------------------
fit_temp <- occuMulti(detformulas = formulas$det, stateformulas = formulas$state, umf_temp)

#' 
#' ## Run power analysis
#' 
#' Required info:
#' 
#' 1. The template model (to tell `unmarked` about the study design etc.)
#' 2. Effect sizes to get power for (called `coefs`)
#' 3. Value of $\alpha$ if we want something other than 0.05
#' 
#' Here's the 'true' effect sizes we specified earlier (we're focused on the interaction effect size, 0.5)
#' 
## -----------------------------------------------------------------------------
effects <- list(state = c(f1, f2, f12), det = c(0, 0))
effects

#' 
#' Now run the power analysis with the `powerAnalysis` function.
#' 
#' This is the power at the sample size we specified (700 sites)
#' 
## -----------------------------------------------------------------------------
pa700 <- powerAnalysis(fit_temp, coefs = effects) # coef = effect sizes to test
pa700

#' 
#' We can also test smaller sample sizes using the `design` argument, which randomly sub-samples the existing sites:
#' 
## -----------------------------------------------------------------------------
pa500 <- powerAnalysis(fit_temp, coefs = coef, design = list(M=500, J=5))
pa300 <- powerAnalysis(fit_temp, coefs = coef, design = list(M=300, J=5))
pa100 <- powerAnalysis(fit_temp, coefs = coef, design = list(M=100, J=5))

#' 
#' Finally, combine all the power analyses into a list for comparison and plot power for the interaction term.
#' 
## -----------------------------------------------------------------------------
pl <- unmarkedPowerList(list(pa100, pa300, pa500, pa700))

plot(pl, param = "[sp1:sp2] (Intercept)", power=0.8)

#' 
#' We need quite a few sites (~700) to have adequate power to detect a significant positive interaction, even when the effect size is relatively large!
#' 
#' I think many studies using this model are significantly under-powered.
#' 
#' ## Poor estimates
#' 
#' * You might get poor estimates under certain conditions:
#'     - Sparse data (many 0s)
#'     - Boundary estimates (occupancy close to 0 or 1)
#'     - Few observations where multiple species are detected
#'     - Separation (perfect correlation with covariate)
#' * What do I mean by poor estimates? Very large absolute values and SEs
#' 
## -----------------------------------------------------------------------------
state_complex <- c(rep("~scale(Dist_5km)+scale(HDens_5km)", 6), 0)
det_complex <- rep("~Trail",3)

mod_complex <- occuMulti(stateformulas=state_complex, detformulas=det_complex, umf)
mod_complex

#' 
#' ### Potential Solutions
#' 
#' * Fewer covariates
#' * Fewer species
#' * Adjust observation period length if possible
#' * Set higher-order interaction terms to 0
#' * Penalized likelihood with `optimizePenalty` function (Clipp et al. 2021)
#' * Bayesian framework
#' 
#' ### Penalized likelihood
#' 
#' "Bayes-inspired" penalty (Hutchinson et al. 2015):
#' 
#' $$
#' -\lambda\frac{1}{2}\sum_i{}\theta_i^2
#' $$
#' 
#' Where $\theta$ is the vector of parameter estimates.
#' 
#' Size of penalty increases as values of $\theta$ increase.
#' 
#' Tradeoff:
#' * Introduce small bias in parameter estimates
#' * Potentially large reduction in variance
#' 
#' ```r
#' mod_penalty <- occuMulti(stateformulas=state_complex, detformulas=det_complex, umf, penalty=1)
#' mod_penalty
#' ```
#' 
#' ```r
#' Bootstraping covariance matrix
#' 
#' 
#' Call:
#' occuMulti(detformulas = det_complex, stateformulas = state_complex, 
#'     data = umf, penalty = 1, maxOrder = 3L)
#' 
#' Occupancy:
#'                                  Estimate    SE      z  P(>|z|)
#' [bobcat] (Intercept)              -1.7810 0.269 -6.627 3.43e-11
#' [bobcat] scale(Dist_5km)          -1.3143 0.293 -4.491 7.08e-06
#' [bobcat] scale(HDens_5km)         -2.8200 0.539 -5.230 1.70e-07
#' [coyote] (Intercept)              -0.6049 0.214 -2.830 4.66e-03
#' [coyote] scale(Dist_5km)           0.0285 0.106  0.270 7.87e-01
#' [coyote] scale(HDens_5km)         -1.0908 0.378 -2.887 3.88e-03
#' [redfox] (Intercept)              -1.5659 0.305 -5.140 2.75e-07
#' [redfox] scale(Dist_5km)          -0.3068 0.157 -1.958 5.02e-02
#' [redfox] scale(HDens_5km)          0.4730 0.546  0.867 3.86e-01
#' [bobcat:coyote] (Intercept)        1.1871 0.434  2.736 6.22e-03
#' [bobcat:coyote] scale(Dist_5km)    0.9347 0.346  2.705 6.84e-03
#' [bobcat:coyote] scale(HDens_5km)  -0.3218 1.008 -0.319 7.50e-01
#' [bobcat:redfox] (Intercept)       -0.8831 0.347 -2.545 1.09e-02
#' [bobcat:redfox] scale(Dist_5km)    0.0364 0.305  0.119 9.05e-01
#' [bobcat:redfox] scale(HDens_5km)   2.5609 1.074  2.384 1.71e-02
#' [coyote:redfox] (Intercept)        1.0001 0.310  3.227 1.25e-03
#' [coyote:redfox] scale(Dist_5km)    0.0236 0.167  0.141 8.87e-01
#' [coyote:redfox] scale(HDens_5km)   1.3920 0.345  4.030 5.57e-05
#' 
#' Detection:
#'                      Estimate     SE      z  P(>|z|)
#' [bobcat] (Intercept)    -2.44 0.1755 -13.89 7.63e-44
#' [bobcat] Trail           1.74 0.1966   8.87 7.25e-19
#' [coyote] (Intercept)    -1.89 0.0926 -20.45 5.64e-93
#' [coyote] Trail           2.10 0.1179  17.84 3.31e-71
#' [redfox] (Intercept)    -1.49 0.1782  -8.35 6.64e-17
#' [redfox] Trail           1.72 0.1714  10.04 9.73e-24
#' 
#' AIC: 6135.555 
#' ```
#' 
#' How to choose penalty value? `optimizePenalty()` function.
#' 
#' ```r
#' mod_penalty <- optimizePenalty(mod_complex, penalties=c(0.5,1))
#' ```
#' 
#' ### Literature Cited
#' 
#' Clipp, H.L., Evans, A.L., Kessinger, B.E., Kellner, K. and Rota, C.T. (2021). A penalized likelihood for multispecies occupancy models improves predictions of species interactions. Ecology, p.e03520.
#' 
#' Rota, C.T., Ferreira, M.A., Kays, R.W., Forrester, T.D., Kalies, E.L., McShea, W.J., Parsons, A.W. and Millspaugh, J.J. (2016). A multispecies occupancy model for two or more interacting species. Methods in Ecology and Evolution, 7(10), pp.1164-1173.
#' 
