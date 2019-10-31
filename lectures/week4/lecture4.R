library(tidyverse)



#' ## Odds versus probabilty


data.frame(Probability = seq(0, 0.99, 0.01)) %>%
  mutate(Odds = Probability / (1 - Probability)) %>%
  ggplot(aes(Probability, Odds)) +
  geom_abline(intercept = 0.0, slope = 1.0, linetype = "dashed") +
  geom_line() +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 2.)) +
  theme_bw()

ggsave("figures/odds-probability-comparison.png",
       h = 4, w = 3)

#' 

p <- seq(0, 1, 0.001)

tibble(Probability = p,
       `n = 5` = sqrt(p*(1-p)/5),
       `n = 10` = sqrt(p*(1-p)/10),
       `n = 25` = sqrt(p*(1-p)/25),
       `n = 50` = sqrt(p*(1-p)/50)) %>%
  gather(n, SE, -1) %>%
  mutate(n = fct_relevel(n, "n = 5")) %>%
  ggplot(aes(Probability, SE, color = n)) +
  geom_line(size = 1.5) +
  scale_color_brewer(element_blank(), palette = "Set1") +
  xlab("Sample proportion") +
  theme_bw(16) +
  theme(panel.grid = element_blank(),
        legend.position = "bottom")

ggsave("figures/proportion_se.png", height = 5, width = 6)
  
plot(p, p * (1-p)/sqrt(5), type = "l")


#' ## Binomial

size <- c(5, 10, 25, 50)
p <- 0.4

df <- tibble(
  n = unlist(Map(rep, size, size)),
  x = unlist(lapply(size, seq_len)),
  p = p) %>%
  mutate(Probability = dbinom(x, n, p),
         group = paste0("n = ", n, ", p = ", p, "; np = ", n*p) %>% fct_inorder())

ggplot(df, aes(x, Probability)) +
  geom_col(width=1.0, color = "grey50") +
  facet_wrap(~group, scale = "free", nrow = 1) +
  theme_bw() +
  theme(panel.grid = element_blank())
  


#' ## Proportion example

prop.test(123, 1000, 0.13, correct = FALSE)
prop.test(123, 1000, 0.13)


#' ## Continuity correction

x <- 7
n <- 12
p <- 0.4
e_x <- n*p
sd_x <- sqrt(n*p*(1-p))

df <- tibble(x = 0:12) %>%
  mutate(Probability = dbinom(x, 12, 0.4),
         tail = x >= 7)

sum(dbinom(x:n, n, p))
pnorm(7, e_x, sd_x, FALSE)
pnorm(6.5, e_x, sd_x, FALSE)

## Binomial 
fig1 <- df %>%
  ggplot(aes(x, Probability, fill = tail)) +
  geom_col(width = 1.0, color = "black", alpha = 0.8) +
  scale_fill_manual(values = c("grey90", "steelblue")) +
  scale_x_continuous(breaks = 0:12) + 
  theme_bw(16) +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  stat_function(fun = dnorm,
                geom = "blank",
                xlim = c(-0.5, 12.5),
                args = list(mean = e_x, sd = sd_x))

ggsave("figures/binom_continuity_correction_1.png", fig1, h=4, w=6)

## Binomial + normal approximation
fig2 <- fig1 +
  stat_function(fun = dnorm,
                xlim = c(-0.5, 12.5),
                args = list(mean = e_x, sd = sd_x),
                color = "red3", size = 1.5) +
  stat_function(fun = dnorm,
                geom = "area",
                xlim = c(x, 12),
                args = list(mean = e_x, sd = sd_x),
                fill = "red3", size = 1.5, alpha = 0.5)

ggsave("figures/binom_continuity_correction_2.png", fig2, h=4, w=6)

## Binomial + normal approximation w/ continuity correction
fig3 <- fig2 +
  stat_function(fun = dnorm,
                xlim = c(-0.5, 12.5),
                args = list(mean = e_x, sd = sd_x),
                color = "red3", size = 1.5) +
  stat_function(fun = dnorm,
                geom = "area",
                xlim = c(x-0.5, 12),
                args = list(mean = e_x, sd = sd_x),
                fill = "forestgreen", size = 1.5, alpha = 0.5)

  ggsave("figures/binom_continuity_correction_3.png", fig3, h=4, w=6)
