## Given the properties we can read off from the
## PH diagram in Eavor's patent, basic properties
## of ethane, ammonia, and the 60:40 mixing 
## proportion, we want to complete the property list
## for the vertices in the power cycle - e.g., 
## temperature, enthalpy, and others if needed.

## First, we source the aforementioned values we
## could read off, which are loaded in 
## parameter_values.R

library(tidyverse)
library(here)

source("parameter_values.R")

## Next, we download isobaric properties for ethane
## and ammonia for p_r values matching those of the 
## 5 power cycle vertices. There are 3 distinct
## p_r values: 0.778 (A), 0.322 (B,C), 0.939 (D,E).

p_r_vec <- c(0.322, 0.778, 0.939)
ammonia_p_vec <- nh3$p_c * p_r_vec
ethane_p_vec  <- c2h6$p_c * p_r_vec

## We haven't yet been able to script the downloading
## of data from NIST, so the data files were gathered
## manually and collected in the data/ subdirectory.

## Next, we proceed through the vertices one at a time
## and work through the fitting procedure. This is
## essentially the same code block for vertices
## A, C, D, E, and B; these could be combined at 
## some point. We did vertex B out of sequence because
## that was the only one where we thought the vertex
## might wind up in the biphasic region and thus
## require special handling. As it happens, it wound
## up in the vapor region, so no special handling was
## required. 


##############################################################
## fit vertex A
##############################################################

## load isobar data
nh3_isobar <- 
  read_delim(here::here("data", "nh3_0p78.tsv"), delim = "\t")
c2h6_isobar <- 
  read_delim(here::here("data", "c2h6_0p78.tsv"), delim = "\t")
## compute reduced temperatures
nh3_isobar$t_r  <- nh3_isobar$`Temperature (K)` / nh3$t_c
c2h6_isobar$t_r <- c2h6_isobar$`Temperature (K)` / c2h6$t_c
## identify range of common t_r values
range(nh3_isobar$t_r)
range(c2h6_isobar$t_r)
## here, this runs from 0.49 to 1.72
## interpolate enthalpies to this common scale
common_t_r <- seq(from = 0.49, to = 1.72, by = 0.01)
mix_h_mat <- matrix(NA, nrow = length(common_t_r), ncol = 4)
colnames(mix_h_mat) <- c("t_r", "nh3_h", "c2h6_h", "mix")
mix_h_mat[, "t_r"] <- common_t_r
mix_h_mat[, "nh3_h"] <- 
  approx(x = nh3_isobar$t_r, y = nh3_isobar$`Enthalpy (kJ/mol)`,
         xout = mix_h_mat[, "t_r"])$y
mix_h_mat[, "c2h6_h"] <- 
  approx(x = c2h6_isobar$t_r, y = c2h6_isobar$`Enthalpy (kJ/mol)`,
         xout = mix_h_mat[, "t_r"])$y
mix_h_mat[, "mix"] <-
  0.4 * (mix_h_mat[, "nh3_h"] - nh3$h_c) + 
  0.6 * (mix_h_mat[, "c2h6_h"] - c2h6$h_c)
new_t_r <-
  approx(x = mix_h_mat[, "mix"], y = mix_h_mat[, "t_r"],
         xout = (power_cycle$A$h - power_cycle$critical$h))$y
## using this new t_r, compute t_nh3, t_c2h6, and t_mix
new_t_r * nh3$t_c
new_t_r * c2h6$t_c
new_t_r * (0.4 * nh3$t_c + 0.6 * c2h6$t_c)
## over at NIST, look up enthalpies for NH3 and C2H6
## at their respective P,T values: 27.383 and 17.203, resp
## compute h_mix as weighted sum
0.6 * 17.203 + 0.4 * 27.383
## compute h_c_mix
0.6 * c2h6$h_c + 0.4 * nh3$h_c
## compare h_mix - h_c_mix with target value (here 4.4)
21.275 - 16.8738

power_cycle$A$t_r <- new_t_r
power_cycle$A$t <- new_t_r * (0.4 * nh3$t_c + 0.6 * c2h6$t_c)




##############################################################
## fit vertex C
##############################################################

## load isobar data
nh3_isobar <- 
  read_delim(here::here("data", "nh3_0p32.tsv"), delim = "\t")
c2h6_isobar <- 
  read_delim(here::here("data", "c2h6_0p32.tsv"), delim = "\t")
## compute reduced temperatures
nh3_isobar$t_r  <- nh3_isobar$`Temperature (K)` / nh3$t_c
c2h6_isobar$t_r <- c2h6_isobar$`Temperature (K)` / c2h6$t_c
## identify range of common t_r values
range(nh3_isobar$t_r)
range(c2h6_isobar$t_r)
## here, this runs from 0.49 to 1.72
## interpolate enthalpies to this common scale
common_t_r <- seq(from = 0.49, to = 1.72, by = 0.01)
mix_h_mat <- matrix(NA, nrow = length(common_t_r), ncol = 4)
colnames(mix_h_mat) <- c("t_r", "nh3_h", "c2h6_h", "mix")
mix_h_mat[, "t_r"] <- common_t_r
mix_h_mat[, "nh3_h"] <- 
  approx(x = nh3_isobar$t_r, y = nh3_isobar$`Enthalpy (kJ/mol)`,
         xout = mix_h_mat[, "t_r"])$y
mix_h_mat[, "c2h6_h"] <- 
  approx(x = c2h6_isobar$t_r, y = c2h6_isobar$`Enthalpy (kJ/mol)`,
         xout = mix_h_mat[, "t_r"])$y
mix_h_mat[, "mix"] <-
  0.4 * (mix_h_mat[, "nh3_h"] - nh3$h_c) + 
  0.6 * (mix_h_mat[, "c2h6_h"] - c2h6$h_c)
new_t_r <-
  approx(x = mix_h_mat[, "mix"], y = mix_h_mat[, "t_r"],
         xout = (power_cycle$C$h - power_cycle$critical$h))$y ## update!
## using this new t_r, compute t_nh3, t_c2h6, and t_mix
new_t_r * nh3$t_c
new_t_r * c2h6$t_c
new_t_r * (0.4 * nh3$t_c + 0.6 * c2h6$t_c)
## over at NIST, look up enthalpies for NH3 and C2H6
## at their respective P,T values: 10.54 and 5.0956, resp
## compute h_mix as weighted sum
0.6 * 5.0956 + 0.4 * 10.54
## compute h_c_mix
0.6 * c2h6$h_c + 0.4 * nh3$h_c
## compare h_mix - h_c_mix with target value (here -9.6)
7.273 - 16.8738

power_cycle$C$t_r <- new_t_r                                  ## update!
power_cycle$C$t <- new_t_r * (0.4 * nh3$t_c + 0.6 * c2h6$t_c) ## update!



##############################################################
## fit vertex D
##############################################################

## load isobar data
nh3_isobar <- 
  read_delim(here::here("data", "nh3_0p94.tsv"), delim = "\t")  ## update!
c2h6_isobar <- 
  read_delim(here::here("data", "c2h6_0p94.tsv"), delim = "\t") ## update!
## compute reduced temperatures
nh3_isobar$t_r  <- nh3_isobar$`Temperature (K)` / nh3$t_c
c2h6_isobar$t_r <- c2h6_isobar$`Temperature (K)` / c2h6$t_c
## identify range of common t_r values
range(nh3_isobar$t_r)
range(c2h6_isobar$t_r)
## here, this runs from 0.49 to 1.72
## interpolate enthalpies to this common scale
common_t_r <- seq(from = 0.49, to = 1.72, by = 0.01)
mix_h_mat <- matrix(NA, nrow = length(common_t_r), ncol = 4)
colnames(mix_h_mat) <- c("t_r", "nh3_h", "c2h6_h", "mix")
mix_h_mat[, "t_r"] <- common_t_r
mix_h_mat[, "nh3_h"] <- 
  approx(x = nh3_isobar$t_r, y = nh3_isobar$`Enthalpy (kJ/mol)`,
         xout = mix_h_mat[, "t_r"])$y
mix_h_mat[, "c2h6_h"] <- 
  approx(x = c2h6_isobar$t_r, y = c2h6_isobar$`Enthalpy (kJ/mol)`,
         xout = mix_h_mat[, "t_r"])$y
mix_h_mat[, "mix"] <-
  0.4 * (mix_h_mat[, "nh3_h"] - nh3$h_c) + 
  0.6 * (mix_h_mat[, "c2h6_h"] - c2h6$h_c)
new_t_r <-
  approx(x = mix_h_mat[, "mix"], y = mix_h_mat[, "t_r"],
         xout = (power_cycle$D$h - power_cycle$critical$h))$y ## update!
## using this new t_r, compute t_nh3, t_c2h6, and t_mix
new_t_r * nh3$t_c
new_t_r * c2h6$t_c
new_t_r * (0.4 * nh3$t_c + 0.6 * c2h6$t_c)
## over at NIST, look up enthalpies for NH3 and C2H6
## at their respective P,T values: 10.647 and 5.191, resp
## compute h_mix as weighted sum
0.6 * 5.191 + 0.4 * 10.647
## compute h_c_mix
0.6 * c2h6$h_c + 0.4 * nh3$h_c
## compare h_mix - h_c_mix with target value (here -9.5)
7.3734 - 16.8738

power_cycle$D$t_r <- new_t_r                                  ## update!
power_cycle$D$t <- new_t_r * (0.4 * nh3$t_c + 0.6 * c2h6$t_c) ## update!




##############################################################
## fit vertex E
##############################################################

## load isobar data
nh3_isobar <- 
  read_delim(here::here("data", "nh3_0p94.tsv"), delim = "\t")  ## update!
c2h6_isobar <- 
  read_delim(here::here("data", "c2h6_0p94.tsv"), delim = "\t") ## update!
## compute reduced temperatures
nh3_isobar$t_r  <- nh3_isobar$`Temperature (K)` / nh3$t_c
c2h6_isobar$t_r <- c2h6_isobar$`Temperature (K)` / c2h6$t_c
## identify range of common t_r values
range(nh3_isobar$t_r)
range(c2h6_isobar$t_r)
## here, this runs from 0.49 to 1.72
## interpolate enthalpies to this common scale
common_t_r <- seq(from = 0.49, to = 1.72, by = 0.01)
mix_h_mat <- matrix(NA, nrow = length(common_t_r), ncol = 4)
colnames(mix_h_mat) <- c("t_r", "nh3_h", "c2h6_h", "mix")
mix_h_mat[, "t_r"] <- common_t_r
mix_h_mat[, "nh3_h"] <- 
  approx(x = nh3_isobar$t_r, y = nh3_isobar$`Enthalpy (kJ/mol)`,
         xout = mix_h_mat[, "t_r"])$y
mix_h_mat[, "c2h6_h"] <- 
  approx(x = c2h6_isobar$t_r, y = c2h6_isobar$`Enthalpy (kJ/mol)`,
         xout = mix_h_mat[, "t_r"])$y
mix_h_mat[, "mix"] <-
  0.4 * (mix_h_mat[, "nh3_h"] - nh3$h_c) + 
  0.6 * (mix_h_mat[, "c2h6_h"] - c2h6$h_c)
new_t_r <-
  approx(x = mix_h_mat[, "mix"], y = mix_h_mat[, "t_r"],
         xout = (power_cycle$E$h - power_cycle$critical$h))$y ## update!
## using this new t_r, compute t_nh3, t_c2h6, and t_mix
new_t_r * nh3$t_c
new_t_r * c2h6$t_c
new_t_r * (0.4 * nh3$t_c + 0.6 * c2h6$t_c)
## over at NIST, look up enthalpies for NH3 and C2H6
## at their respective P,T values: 27.831 and 17.579, resp
## compute h_mix as weighted sum
0.6 * 17.579 + 0.4 * 27.831
## compute h_c_mix
0.6 * c2h6$h_c + 0.4 * nh3$h_c
## compare h_mix - h_c_mix with target value (here 4.8)
21.6798 - 16.8738

power_cycle$E$t_r <- new_t_r                                  ## update!
power_cycle$E$t <- new_t_r * (0.4 * nh3$t_c + 0.6 * c2h6$t_c) ## update!





##############################################################
## fit vertex B
##############################################################

## load isobar data
nh3_isobar <- 
  read_delim(here::here("data", "nh3_0p32.tsv"), delim = "\t")  ## update!
c2h6_isobar <- 
  read_delim(here::here("data", "c2h6_0p32.tsv"), delim = "\t") ## update!
## compute reduced temperatures
nh3_isobar$t_r  <- nh3_isobar$`Temperature (K)` / nh3$t_c
c2h6_isobar$t_r <- c2h6_isobar$`Temperature (K)` / c2h6$t_c
## identify range of common t_r values
range(nh3_isobar$t_r)
range(c2h6_isobar$t_r)
## here, this runs from 0.49 to 1.72
## interpolate enthalpies to this common scale
common_t_r <- seq(from = 0.49, to = 1.72, by = 0.01)
mix_h_mat <- matrix(NA, nrow = length(common_t_r), ncol = 4)
colnames(mix_h_mat) <- c("t_r", "nh3_h", "c2h6_h", "mix")
mix_h_mat[, "t_r"] <- common_t_r
mix_h_mat[, "nh3_h"] <- 
  approx(x = nh3_isobar$t_r, y = nh3_isobar$`Enthalpy (kJ/mol)`,
         xout = mix_h_mat[, "t_r"])$y
mix_h_mat[, "c2h6_h"] <- 
  approx(x = c2h6_isobar$t_r, y = c2h6_isobar$`Enthalpy (kJ/mol)`,
         xout = mix_h_mat[, "t_r"])$y
mix_h_mat[, "mix"] <-
  0.4 * (mix_h_mat[, "nh3_h"] - nh3$h_c) + 
  0.6 * (mix_h_mat[, "c2h6_h"] - c2h6$h_c)
new_t_r <-
  approx(x = mix_h_mat[, "mix"], y = mix_h_mat[, "t_r"],
         xout = (power_cycle$B$h - power_cycle$critical$h))$y ## update!
## using this new t_r, compute t_nh3, t_c2h6, and t_mix
new_t_r * nh3$t_c
new_t_r * c2h6$t_c
new_t_r * (0.4 * nh3$t_c + 0.6 * c2h6$t_c)
## over at NIST, look up enthalpies for NH3 and C2H6
## at their respective P,T values: 27.678 and 16.817, resp
## compute h_mix as weighted sum
0.6 * 16.817 + 0.4 * 27.678
## compute h_c_mix
0.6 * c2h6$h_c + 0.4 * nh3$h_c
## compare h_mix - h_c_mix with target value (here 3.9)
21.1614 - 16.8738
## We're notably off here; we get 4.3 vs 3.9. 
## Their shift into the biphasic region hurts our modeling. 

power_cycle$B$t_r <- new_t_r                                  ## update!
power_cycle$B$t <- new_t_r * (0.4 * nh3$t_c + 0.6 * c2h6$t_c) ## update!

