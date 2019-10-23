## ----setup, echo = FALSE-------------------------------------------------
my_plot_hook <- function(x, options)
  paste("\n", knitr::hook_plot_tex(x, options), "\n")
knitr::knit_hooks$set(plot = my_plot_hook)


## ---- message = FALSE----------------------------------------------------
library(tidyverse)
perulung <- read.csv("perulung_ems.csv")
perulung <- perulung %>%
  mutate(sex = factor(sex, c(0, 1), c("female", "male")),
         respsymptoms = factor(respsymptoms, c(0, 1), c("no symptoms", "symptoms")))


## ---- results = "hold"---------------------------------------------------
x <- perulung$height

## A function to return requested outputs as a vector
calc_ex1a <- function(x, digits = 3) {
  mu <- mean(x)
  se <- sd(x) / sqrt(length(x))

  print_ci <- function(ci) paste0("(", round(ci[1], digits), ", ", 
                                  round(ci[2], digits), ")")
 
  c(mean = round(mu, digits),  
    ci_norm = print_ci(mu + c(-1, 1) * qnorm(0.975) * se),
    ci_tdist = print_ci(mu + c(-1, 1) * qt(0.975, length(x)-1) * se),
    ci_t.test = print_ci(t.test(x)$conf.int))
}

rbind(
  height = calc_ex1a(perulung$height),
  fev1 = calc_ex1a(perulung$fev1),
  height_female  = calc_ex1a(perulung$height[perulung$sex == "female"]),
  height_male  = calc_ex1a(perulung$height[perulung$sex == "male"]),
  fev1_nosymp  = calc_ex1a(perulung$fev1[perulung$respsymptoms == "no symptoms"]),
  fev1_symptom  = calc_ex1a(perulung$fev1[perulung$respsymptoms == "symptoms"])
) %>%
knitr::kable(align = "lcccc")


## ---- echo = FALSE-------------------------------------------------------
print_ci <- function(ci, digits = 3) paste0("(", round(ci[1], digits), ", ", 
                                             round(ci[2], digits), ")")


## ------------------------------------------------------------------------
mean(perulung$height)
t.test(perulung$height)$conf.int
t.test(perulung$height, mu = 124, alternative = "greater")


## ------------------------------------------------------------------------
x_female <- perulung$height[perulung$sex == "female"]
length(x_female)
mean(x_female)
t.test(x_female, mu = 123.5)



## ------------------------------------------------------------------------
table(perulung$sex)
x_female <- perulung$height[perulung$sex == "female"]
x_male <- perulung$height[perulung$sex == "male"]

mean(x_male) - mean(x_female)
t.test(x_male, x_female, var.equal = TRUE)


## ------------------------------------------------------------------------
fit <- lm(height ~ age, data = perulung)
summary(fit)
confint(fit)


## ------------------------------------------------------------------------
fev1_symp <- perulung$fev1[perulung$respsymptoms == "symptoms"]
fev1_nosymp <- perulung$fev1[perulung$respsymptoms == "no symptoms"]
mean(fev1_symp) - mean(fev1_nosymp)
t.test(fev1_symp, fev1_nosymp)
t.test(fev1_symp, fev1_nosymp, alternative = "less")


## ------------------------------------------------------------------------
fit6 <- lm(fev1 ~ height + age, data = perulung)
summary(fit6)
confint(fit6)



## ------------------------------------------------------------------------
library(NHANES)
data(NHANES)
nhanes20pl <- NHANES[NHANES$Age >= 20, ]


## ---- fig.height = 2.5, fig.width = 3, warning = FALSE-------------------
ggplot(nhanes20pl, aes(Height)) +
geom_histogram(bins = 50)

ggplot(nhanes20pl, aes(BMI)) +
geom_histogram(bins = 50)

ggplot(nhanes20pl, aes(AlcoholYear)) +
geom_histogram(bins = 50)


## ------------------------------------------------------------------------
height_mean <- mean(nhanes20pl$Height, na.rm = TRUE)
height_sd <- sd(nhanes20pl$Height, na.rm = TRUE)
bmi_mean <- mean(nhanes20pl$BMI, na.rm = TRUE)
bmi_sd <- sd(nhanes20pl$BMI, na.rm = TRUE)
alc_mean <- mean(nhanes20pl$AlcoholYear, na.rm = TRUE)
alc_sd <- sd(nhanes20pl$AlcoholYear, na.rm = TRUE)

data.frame(mean = c(height_mean, bmi_mean, alc_mean),
           sd = c(height_sd, bmi_sd, alc_sd),
           row.names = c("height", "bmi", "alcohol_year")) %>%
knitr::kable(digits = 1, align = "cc")


## ------------------------------------------------------------------------
ans2c <- c(1-pnorm(165, height_mean, height_sd),
            pnorm(160, height_mean, height_sd) - pnorm(153, height_mean, height_sd),
            qnorm(0.9, height_mean, height_sd),
            1-pnorm(30, bmi_mean, bmi_sd),
            pnorm(30, bmi_mean, bmi_sd) - pnorm(25, bmi_mean, bmi_sd),
            qnorm(0.25, bmi_mean, bmi_sd),
            1-pnorm(100, alc_mean, alc_sd),
            pnorm(10, alc_mean, alc_sd),
            qnorm(0.75, alc_mean, alc_sd) - qnorm(0.25, alc_mean, alc_sd))


## ------------------------------------------------------------------------
ans2d <- c(mean(nhanes20pl$Height > 165, na.rm=TRUE),
           mean(nhanes20pl$Height > 153 & nhanes20pl$Height < 160, na.rm=TRUE),
           quantile(nhanes20pl$Height, 0.9, na.rm=TRUE),
           mean(nhanes20pl$BMI > 30, na.rm=TRUE),
           mean(nhanes20pl$BMI > 25 & nhanes20pl$BMI <= 30, na.rm=TRUE),
           quantile(nhanes20pl$BMI, 0.25, na.rm=TRUE),
           mean(nhanes20pl$AlcoholYear > 100, na.rm=TRUE),
           mean(nhanes20pl$AlcoholYear <= 10, na.rm=TRUE),
           diff(quantile(nhanes20pl$AlcoholYear, c(0.25, 0.75), na.rm=TRUE)))
    
labels <- c("Height above 165cm", "Height between 153-160cm", "Height 90th percentile",
             "BMI above 30", "BMI between 25-30", "BMI 25th percentile",
             "Alcohol days above 100", "Alcohol days below 10", "Alcohol days IQR")
data.frame(normal_approx = ans2c,
           empirical_dist = ans2d,
           row.names = labels) %>%
           knitr::kable(digits = 2)


## ------------------------------------------------------------------------
height <- nhanes20pl$Height[!is.na(nhanes20pl$Height)]
bmi <- nhanes20pl$BMI[!is.na(nhanes20pl$BMI)]
alcohol <- nhanes20pl$AlcoholYear[!is.na(nhanes20pl$AlcoholYear)]


## ---- cache = TRUE-------------------------------------------------------
set.seed(77316870)

nsim <- 10000 # number of simulated datasets
nsamp <- c(5, 10, 25, 50, 100, 500)   # sample size

sim_ex3 <- function(x, n) {

  mu <- mean(x)  # true population mean
  samp <- replicate(nsim, sample(x, n, replace = TRUE))

  est <- colMeans(samp)
  se <- apply(samp, 2, sd) / sqrt(n)

  ci_l_z <- est - qnorm(0.975) * se
  ci_u_z <- est + qnorm(0.975) * se

  ci_l_t <- est - qt(0.975, n-1) * se
  ci_u_t <- est + qt(0.975, n-1) * se
  
  list(est = est,
       cover_z = mean(ci_l_z < mu & mu < ci_u_z),
       cover_t = mean(ci_l_t < mu & mu < ci_u_t))
}

height_res <- list()
bmi_res <- list()
alc_res <- list()

for(n in nsamp){
  height_res[[as.character(n)]] <- sim_ex3(height, n)
  bmi_res[[as.character(n)]] <- sim_ex3(bmi, n)
  alc_res[[as.character(n)]] <- sim_ex3(alcohol, n)
}


## ---- fig.height = 3.5, fig.width=6, fig.show = "hold", cache = TRUE-----
## sample mean summary
sample_mean_summary <-
  bind_rows(
    data.frame(outcome = "height", 
               sapply(height_res, "[[", "est"), 
               check.names = FALSE, stringsAsFactors = FALSE),
    data.frame(outcome = "bmi", 
               sapply(bmi_res, "[[", "est"), 
               check.names = FALSE, stringsAsFactors = FALSE),
    data.frame(outcome = "alcohol", 
               sapply(alc_res, "[[", "est"), 
               check.names = FALSE, stringsAsFactors = FALSE)
  ) %>%
  gather(n, sample_mean, `5`:`500`) %>%
  mutate(outcome = fct_relevel(outcome, "height", "bmi", "alcohol"),
         n = as.integer(n))
    
sample_mean_summary %>%
  filter(outcome == "height") %>%
  ggplot(aes(sample_mean)) +
  geom_histogram(aes(y = ..density..), bins = 100) +
  geom_vline(xintercept = mean(height), color = "red") +
  facet_wrap(~n, scales = "free_y") +
  theme_light() +
  ggtitle("Height: sample mean")

sample_mean_summary %>%
  filter(outcome == "bmi") %>%
  ggplot(aes(sample_mean)) +
  geom_histogram(aes(y = ..density..), bins = 100) +
  geom_vline(xintercept = mean(bmi), color = "red") +
  facet_wrap(~n, scales = "free_y") +
  theme_light() +
  ggtitle("BMI: sample mean")

sample_mean_summary %>%
  filter(outcome == "alcohol") %>%
  ggplot(aes(sample_mean)) +
  geom_histogram(aes(y = ..density..), bins = 100) +
  geom_vline(xintercept = mean(alcohol), color = "red") +
  facet_wrap(~n, scales = "free_y") +
  theme_light() +
  ggtitle("Alcohol days per year: sample mean")


## ---- fig.height = 2.5, fig.width=7, cache = TRUE------------------------
## Compile coverage estimates
cover_summary <- 
  bind_rows(
    tibble(outcome = "height", 
           n = nsamp, 
           ci_type = "z", 
           coverage = sapply(height_res, "[[", "cover_z")),
    tibble(outcome = "bmi", 
           n = nsamp, 
           ci_type = "z",
           coverage = sapply(bmi_res, "[[", "cover_z")),
    tibble(outcome = "alcohol", 
           n = nsamp, 
           ci_type = "z",
           coverage = sapply(alc_res, "[[", "cover_z")),
    tibble(outcome = "height", 
           n = nsamp, 
           ci_type = "t", 
           coverage = sapply(height_res, "[[", "cover_t")),
    tibble(outcome = "bmi", 
           n = nsamp, 
           ci_type = "t",
           coverage = sapply(bmi_res, "[[", "cover_t")),
    tibble(outcome = "alcohol", 
           n = nsamp, 
           ci_type = "t",
           coverage = sapply(alc_res, "[[", "cover_t")),
    
  ) %>%
  mutate(outcome = fct_relevel(outcome, "height", "bmi", "alcohol"))

## summary table
cover_summary %>% 
  unite(key, outcome, ci_type) %>%
  spread(key, coverage) %>%
  select(n, height_z, height_t, bmi_z, bmi_t, alcohol_z, alcohol_t) %>%
  knitr::kable(caption = "Coverage of theoretical 95% confidence intervals")

## summary figure
cover_summary %>%
  mutate(n = as.factor(n)) %>%
  ggplot(aes(n, coverage, color = ci_type, group = ci_type)) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  geom_point() +
  geom_line() +
  scale_color_brewer(palette = "Set1") +
  facet_wrap(~outcome) +
  theme_light()


