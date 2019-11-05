## ----include = FALSE-----------------------------------------------------
options(knitr.kable.NA = "")


## ----ex1a, message = FALSE, echo = FALSE, results = "hide"---------------
## Load data and convert binary outcomes to factors
library(tidyverse)
perulung <- read.csv("perulung_ems.csv")
perulung$sex <- factor(perulung$sex, c(0, 1), c("female", "male"))
perulung$respsymptoms = factor(perulung$respsymptoms,
                               c(0, 1), c("no symptoms", "symptoms"))

## Calculate regression coefficient estimates
y <- perulung$fev1
x <- perulung$height
ybar <- mean(y)
xbar <- mean(x)

beta1_hat <- sum((x - xbar) * (y - ybar)) / sum((x - xbar)^2)
beta0_hat <- ybar - beta1_hat * xbar

## Calculate residual standard deviation
ypred <- beta0_hat + beta1_hat * x
rss <- sum((y - ypred)^2)
df <- length(x) - 2
sigma_hat <- sqrt(rss / df)

## Standard error for regression parameters
se_beta1 <- sigma_hat / sqrt(sum((x - xbar)^2))
se_beta0 <- sigma_hat * sqrt(1/length(x) + xbar^2 / sum((x - xbar)^2))

## t-statistic, 95% CI, and p-value
beta0_tt <- beta0_hat / se_beta0
beta0_ci <- beta0_hat + c(-1,1) * qt(0.975, df) * se_beta0
beta0_p <- 2 * pt(abs(beta0_tt), df, lower.tail = FALSE)

beta1_tt <- beta1_hat / se_beta1
beta1_ci <- beta1_hat + c(-1,1) * qt(0.975, df) * se_beta1
beta1_p <- 2 * pt(abs(beta1_tt), df, lower.tail = FALSE)

## Compare with lm(...)
fit1a <- lm(fev1 ~ height, data = perulung)


## ------------------------------------------------------------------------
data.frame(parameter = c("beta0 (Intercept)", "beta1 (height)", 
                         "sigma (residual std. dev.)"),
           estimate = c(beta0_hat, beta1_hat, sigma_hat),
           std_error = c(se_beta0, se_beta1, NA),
           t_stat    = c(beta0_tt, beta1_tt, NA),
           ci_lower  = c(beta0_ci[1], beta1_ci[1], NA),
           ci_upper  = c(beta0_ci[2], beta1_ci[2], NA),
           p_value   = c(beta0_p, beta1_p, NA)) %>%
  knitr::kable(digits = 3)


## ---- fig.height = 2.7, fig.width = 3, fig.align = "center", fig.show = "hold"----
## Using base R graphics
par(mar = c(3, 3, 0.5, 0.5), tcl = -0.25, mgp = c(2, 0.5, 0))
plot(fev1 ~ height, data = perulung, 
     las = 1, xlab = "Height (cm)", ylab = "FEV1 (litres/second)")
abline(a = beta0_hat, b = beta1_hat, col = "red", lwd = 2)

## Using ggplot
ggplot(perulung, aes(height, fev1)) + 
  geom_point() + 
  geom_abline(slope = beta1_hat, intercept = beta0_hat, color = "red") +
  labs(x = "Height (cm)", y = "FEV1 (litres/second")


## ------------------------------------------------------------------------
fit1a <- lm(fev1 ~ height, data = perulung)
summary(fit1a)
confint(fit1a)


## ----ex1a, message = FALSE, echo = TRUE, results = "hide"----------------
## Load data and convert binary outcomes to factors
library(tidyverse)
perulung <- read.csv("perulung_ems.csv")
perulung$sex <- factor(perulung$sex, c(0, 1), c("female", "male"))
perulung$respsymptoms = factor(perulung$respsymptoms,
                               c(0, 1), c("no symptoms", "symptoms"))

## Calculate regression coefficient estimates
y <- perulung$fev1
x <- perulung$height
ybar <- mean(y)
xbar <- mean(x)

beta1_hat <- sum((x - xbar) * (y - ybar)) / sum((x - xbar)^2)
beta0_hat <- ybar - beta1_hat * xbar

## Calculate residual standard deviation
ypred <- beta0_hat + beta1_hat * x
rss <- sum((y - ypred)^2)
df <- length(x) - 2
sigma_hat <- sqrt(rss / df)

## Standard error for regression parameters
se_beta1 <- sigma_hat / sqrt(sum((x - xbar)^2))
se_beta0 <- sigma_hat * sqrt(1/length(x) + xbar^2 / sum((x - xbar)^2))

## t-statistic, 95% CI, and p-value
beta0_tt <- beta0_hat / se_beta0
beta0_ci <- beta0_hat + c(-1,1) * qt(0.975, df) * se_beta0
beta0_p <- 2 * pt(abs(beta0_tt), df, lower.tail = FALSE)

beta1_tt <- beta1_hat / se_beta1
beta1_ci <- beta1_hat + c(-1,1) * qt(0.975, df) * se_beta1
beta1_p <- 2 * pt(abs(beta1_tt), df, lower.tail = FALSE)

## Compare with lm(...)
fit1a <- lm(fev1 ~ height, data = perulung)


## ---- fig.height = 4, fig.width = 7--------------------------------------
par(mfrow = c(1,2))
plot(fit1a, 1:2)


## ------------------------------------------------------------------------
perulung$height_cat <- cut(perulung$height, c(-Inf, 120, 130, Inf), 
                           labels = c("<=120cm", "120-130cm", ">130cm"))


## ------------------------------------------------------------------------
table(perulung$height_cat)


## ------------------------------------------------------------------------
fit1d <- lm(fev1 ~ height_cat, data = perulung)
summary(fit1d)


## ------------------------------------------------------------------------
aggregate(fev1 ~ height_cat, perulung, mean)


## ------------------------------------------------------------------------
levels(perulung$height_cat)


## ------------------------------------------------------------------------
perulung$height_cat <- relevel(perulung$height_cat, "120-130cm")
levels(perulung$height_cat)


## ------------------------------------------------------------------------
fit1e <- lm(fev1 ~ height_cat, perulung)
summary(fit1e)


## ------------------------------------------------------------------------
fit1g <- lm(fev1 ~ sex, perulung)
summary(fit1g)
confint(fit1g)

t.test(perulung$fev1[perulung$sex == "male"], 
       perulung$fev1[perulung$sex == "female"], var.equal = TRUE)



## ----results = "hide"----------------------------------------------------
library(NHANES)
data(NHANES)
nhanes_child <- subset(NHANES, AgeMonths < 120)


## ----results="hide", fig.height = 3, fig.width = 3, fig.show = "hold", fig.align = "center", out.width = "40%"----
ex3a <- subset(nhanes_child, !is.na(Height) & !is.na(Length))
nrow(ex3a)
summary(ex3a[c("AgeMonths", "Height", "Length")])
sd(ex3a$Height)
sd(ex3a$Length)
cor(ex3a$Height, ex3a$Length)

par(mar = c(3, 3, 0.5, 0.5), tcl = -0.25, mgp = c(2, 0.5, 0))
plot(ex3a$Height, ex3a$Length)
abline(a = 0, b = 1)

t.test(ex3a$Length, ex3a$Height, paired = TRUE)


## ----ex3b----------------------------------------------------------------
nhanes_child$height_all <- ifelse(nhanes_child$AgeMonths < 24, nhanes_child$Length, nhanes_child$Height)
fit2b <- lm(height_all ~ AgeMonths, nhanes_child)
summary(fit2b)


## ---- fig.height = 3, fig.width = 3, fig.align = "center"----------------
par(mar = c(3, 3, 0.5, 0.5), tcl = -0.25, mgp = c(2, 0.5, 0))
plot(height_all ~ AgeMonths, nhanes_child)
abline(fit2b, col = "red", lwd = 2)


## ----ex3c, fig.height = 4, fig.width = 7, fig.align = "center"-----------
par(mfrow = c(1,2), tcl = -0.25, mgp = c(2, 0.5, 0))
plot(fit2b, 1, lwd = 3); plot(fit2b, 2)


## ----ex3d, fig.height = 4, fig.width = 7, fig.align = "center"-----------
## Square-root transform AgeMonths to improve linearity
fit2d1 <- lm(height_all ~ sqrt(AgeMonths), nhanes_child)
summary(fit2d1)

## Plot predicted height versus observed data on transformed and
## and original scales
par(mfrow = c(1,2), tcl = -0.25, mgp = c(2, 0.5, 0))
plot(height_all ~ sqrt(AgeMonths), nhanes_child,
     main = "Transformed scale")
abline(fit2d1, col = "red", lwd = 3)

fit2d1_pred <- predict(fit2d1, newdata = data.frame(AgeMonths = 0:120))
plot(height_all ~ AgeMonths, nhanes_child,
     main = "Original scale")
lines(0:120, fit2d1_pred, col = "red", lwd = 3)


## Residual analysis
plot(fit2d1, 1:2, lwd = 3)


## ----fig.height = 4, fig.width = 7, fig.align = "center"-----------------
## Log-transform height_all to address heteroskedasticity
fit2d2 <- lm(log(height_all) ~ sqrt(AgeMonths), nhanes_child)

## Plot predicted height versus observed data on transformed and
## and original scales
par(mfrow = c(1,2), tcl = -0.25, mgp = c(2, 0.5, 0))
plot(log(height_all) ~ sqrt(AgeMonths), nhanes_child,
     main = "Transformed scale")
abline(fit2d2, col = "red", lwd = 3)

fit2d2_pred <- predict(fit2d2, newdata = data.frame(AgeMonths = 0:120))
plot(height_all ~ AgeMonths, nhanes_child, main = "Original scale")
lines(0:120, exp(fit2d2_pred), col = "red", lwd = 3)

## Residual analysis
plot(fit2d2, 1:2, lwd = 3)


## ----fig.height = 4, fig.width = 3.5, fig.align = "center"---------------
## Add both AgeMonths and sqrt(AgeMonths) as predictors to improve linear trend
fit2d3 <- lm(log(height_all) ~ AgeMonths + sqrt(AgeMonths), nhanes_child)

## Plot predicted height versus observed data on original scale
par(mfrow = c(1,1), tcl = -0.25, mgp = c(2, 0.5, 0))   
fit2d3_pred <- predict(fit2d3, newdata = data.frame(AgeMonths = 0:120))
plot(height_all ~ AgeMonths, nhanes_child, main = "Original scale")
lines(0:120, exp(fit2d3_pred), col = "red", lwd = 3)


## ----fig.height = 4, fig.width = 7, fig.align = "center"-----------------
## Residual analysis
par(mfrow = c(1,2), tcl = -0.25, mgp = c(2, 0.5, 0))
plot(fit2d3, 1:2, lwd = 3)


## ------------------------------------------------------------------------
summary(fit2d3)


## ------------------------------------------------------------------------
## install.packages("sandwich")
library(sandwich)

confint_robust <- function(object, parm, level = 0.95, type = "HC3") {

  cf <- coef(object)
  vcov <- sandwich::vcovHC(object, type = type)

   pnames <- names(cf)
   if (missing(parm))
       parm <- pnames
   else if (is.numeric(parm))
       parm <- pnames[parm]
   a <- (1 - level)/2
   a <- c(a, 1 - a)
   pct <- stats:::format.perc(a, 3)
   fac <- qnorm(a)
   ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, pct))
   ses <- sqrt(diag(vcov))[parm]
   ci[] <- cf[parm] + ses %o% fac
   ci
}


## ------------------------------------------------------------------------
confint(fit2b)
confint_robust(fit2b)


## ----ex3a, cache = TRUE--------------------------------------------------
#' @param n simulated sample size.
#' @param data a data frame from which to subsample with replacement.
#' @param formula a formula with a numeric outcome and covariate (y ~ x).
#' @param beta1_true the true value for the slope.
#'
#' @return
#' A numeric vector consisting of three elements:
#'   1. The estimated regression slope from the simulated dataset.
#'   2. A binary outcome (0/1) whether the least squares 95% CI contains
#'      the true slope.
#'   3. A binary outcome (0/1) whether the robust 95% CI contains the true
#'      slope.
#'
sim_lm <- function(n, data, formula, beta1_true) {
  
  df_sim <- data[sample.int(nrow(data), n, replace = TRUE), ]
  fit <- lm(formula, df_sim)
  
  ci_ls <- confint(fit)[2,]
  ci_robust <- confint_robust(fit)[2,]
  
  covers <- c(coef(fit)[2],
              ci_ls[1] < beta1_true & beta1_true < ci_ls[2],
              ci_robust[1] < beta1_true & beta1_true < ci_robust[2])
  names(covers) <- c("estimate", "least_squares", "robust")
  covers
}

nsim <- 3000
data <- nhanes_child
formula <- height_all ~ AgeMonths
b1_true <- coef(lm(height_all ~ AgeMonths, nhanes_child))[2]

sim_n10 <- replicate(nsim, sim_lm(10, data, formula, b1_true))
sim_n25 <- replicate(nsim, sim_lm(25, data, formula, b1_true))
sim_n50 <- replicate(nsim, sim_lm(50, data, formula, b1_true))
sim_n100 <- replicate(nsim, sim_lm(100, data, formula, b1_true))
sim_n500 <- replicate(nsim, sim_lm(500, data, formula, b1_true))

ex3a <- bind_rows(
  data.frame(n = 10, t(sim_n10)),
  data.frame(n = 25, t(sim_n25)),
  data.frame(n = 50, t(sim_n50)),
  data.frame(n = 100, t(sim_n100)),
  data.frame(n = 500, t(sim_n500))
)


## ---- fig.height = 2, fig.width = 7--------------------------------------
ex3a %>%
  ggplot(aes(estimate)) +
  geom_histogram(aes(y = ..density..), bins = 100) +
  geom_vline(xintercept = b1_true, color = "red") +
  facet_wrap(~n, scales = "free_y", nrow = 1) +
  theme_light() +
  ggtitle("beta1: slope for change in height by AgeMonths")


## ---- fig.height = 3, fig.width = 5, fig.align = "center"----------------
ex1a %>%
  mutate(n = as.factor(n)) %>%
  group_by(n) %>%
  summarise(lsq_coverage = mean(least_squares)) %>%
  ggplot(aes(n, lsq_coverage)) +
  geom_col() +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  coord_cartesian(ylim = c(0.8, 1.0)) +
  theme_light() +
  labs(title = "Coverage of least-squares 95% CI",
       y = "95% CI coverage probability")


## ---- fig.height = 3, fig.width = 5, fig.align = "center"----------------
ex1a %>%
  mutate(n = as.factor(n)) %>%
  gather(CI_type, covers, least_squares, robust) %>%
  group_by(n, CI_type) %>%
  summarise(coverage = mean(covers)) %>%
  ggplot(aes(n, coverage, fill = CI_type)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  coord_cartesian(ylim = c(0.8, 1.0)) +
  theme_light() +
  labs(title = "Coverage of least-squares 95% CI",
       y = "95% CI coverage probability")


## ---- cache = TRUE, warning = FALSE--------------------------------------
set.seed(6841871)
  
nhanes20pl <- NHANES[NHANES$Age >= 20, ]

## Height ~ Age 
data <- nhanes20pl[!is.na(nhanes20pl$Height), ]
formula <- Height ~ Age
b1_true <- coef(lm(formula, data))[2]

height_n10 <- replicate(nsim, sim_lm(10, data, formula, b1_true))
height_n25 <- replicate(nsim, sim_lm(25, data, formula, b1_true))
height_n50 <- replicate(nsim, sim_lm(50, data, formula, b1_true))
height_n100 <- replicate(nsim, sim_lm(100, data, formula, b1_true))
height_n500 <- replicate(nsim, sim_lm(500, data, formula, b1_true))

## BMI ~ Age 
data <- nhanes20pl[!is.na(nhanes20pl$BMI), ]
formula <- BMI ~ Age
b1_true <- coef(lm(formula, data))[2]

bmi_n10 <- replicate(nsim, sim_lm(10, data, formula, b1_true))
bmi_n25 <- replicate(nsim, sim_lm(25, data, formula, b1_true))
bmi_n50 <- replicate(nsim, sim_lm(50, data, formula, b1_true))
bmi_n100 <- replicate(nsim, sim_lm(100, data, formula, b1_true))
bmi_n500 <- replicate(nsim, sim_lm(500, data, formula, b1_true))

## AlcoholYear ~ Age 
data <- nhanes20pl[!is.na(nhanes20pl$AlcoholYear), ]
formula <- AlcoholYear ~ Age
b1_true <- coef(lm(formula, data))[2]

alc_n10 <- replicate(nsim, sim_lm(10, data, formula, b1_true))
alc_n25 <- replicate(nsim, sim_lm(25, data, formula, b1_true))
alc_n50 <- replicate(nsim, sim_lm(50, data, formula, b1_true))
alc_n100 <- replicate(nsim, sim_lm(100, data, formula, b1_true))
alc_n500 <- replicate(nsim, sim_lm(500, data, formula, b1_true))

ex3e <- bind_rows(
  data.frame(outcome = "Height", n = 10, t(height_n10)),
  data.frame(outcome = "Height", n = 25, t(height_n25)),
  data.frame(outcome = "Height", n = 50, t(height_n50)),
  data.frame(outcome = "Height", n = 100, t(height_n100)),
  data.frame(outcome = "Height", n = 500, t(height_n500)),
  data.frame(outcome = "BMI", n = 10, t(bmi_n10)),
  data.frame(outcome = "BMI", n = 25, t(bmi_n25)),
  data.frame(outcome = "BMI", n = 50, t(bmi_n50)),
  data.frame(outcome = "BMI", n = 100, t(bmi_n100)),
  data.frame(outcome = "BMI", n = 500, t(bmi_n500)),
  data.frame(outcome = "Alcohol days per year", n = 10, t(alc_n10)),
  data.frame(outcome = "Alcohol days per year", n = 25, t(alc_n25)),
  data.frame(outcome = "Alcohol days per year", n = 50, t(alc_n50)),
  data.frame(outcome = "Alcohol days per year", n = 100, t(alc_n100)),
  data.frame(outcome = "Alcohol days per year", n = 500, t(alc_n500))  
)



## ---- fig.height = 7, fig.width = 5, fig.align = "center"----------------

ex3e %>%
  mutate(n = as.factor(n)) %>%
  gather(CI_type, covers, least_squares, robust) %>%
  group_by(outcome, n, CI_type) %>%
  summarise(coverage = mean(covers)) %>%
  ungroup() %>%
  mutate(outcome = fct_relevel(outcome, "Height", "BMI")) %>%
  ggplot(aes(n, coverage, fill = CI_type)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  facet_wrap(~outcome, nrow = 3) +
  coord_cartesian(ylim = c(0.8, 1.0)) +
  theme_light() +
  labs(title = "Coverage of least-squares 95% CI",
       y = "95% CI coverage probability")

