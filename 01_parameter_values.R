## parameter values used in various analyses. 
## Most of these come from Eavor's public patents
## and press releases.
##
## Fluid properties are drawn 
## from the National Institute of Standards Webbook
## (NIST).

## ammonia properties

nh3 <- list(
  t_c = 405.40,  ## K
  p_c = 113.330, ## bar 
  h_c = 22.416,  ## kJ / mol
  s_c = 70.823,  ## J / (mol K)
  m_w = 17.0305  ## g / mol 
)

## ethane properties

c2h6 <- list(
  t_c = 305.33, ## K
  p_c = 48.718, ## bar
  h_c = 13.179, ## kJ / mol
  s_c = 50.722, ## J / (mol K)
  m_w = 30.069  ## g / mol
)

## fluid mix properties

fluid_mix <- list(
  components = c("ammonia", "ethane"),
  component_weights = c(0.4, 0.6),
  t_c = 345.36, ## K
  p_c = 74.563, ## bar
  h_c = 16.874, ## kJ / mol
  s_c = 58.762, ## J / (mol K)
  m_w = 24.854  ## g / mol
)

## water properties

## power_cycle

power_cycle <- list(
  A = list(
    p = 58,  ## bar
    h = 9.5, ## kJ / mol
    t = NA,  ## K
    s = NA,  ## J / (mol K)
    x_v = 1, ## unitless
    p_r = NA, ## unitless
    t_r = NA, ## unitless
    location = "entering turbine"
  ),
  B = list(
    p = 24,   ## bar
    h = 8.0,  ## kJ / mol
    t = NA,   ## K
    s = NA,   ## J / (mol K)
    x_v = NA, ## unitless
    p_r = NA, ## unitless
    t_r = NA, ## unitless
    location = "exiting turbine"
  ),
  C = list(
    p = 24,   ## bar
    h = -4.5, ## kJ / mol
    t = NA,   ## K
    s = NA,   ## J / (mol K)
    x_v = 0,  ## unitless
    p_r = NA, ## unitless
    t_r = NA, ## unitless
    location = "descending vertical bore"
  ),
  D = list(
    p = 70,   ## bar
    h = -4.4, ## kJ / mol
    t = NA,   ## K
    s = NA,   ## J / (mol K)
    x_v = 0,  ## unitless
    p_r = NA, ## unitless
    t_r = NA, ## unitless
    location = "entering lateral"
  ),
  E = list(
    p = 70,   ## bar
    h = 9.9,  ## kJ / mol
    t = NA,   ## K
    s = NA,   ## J / (mol K)
    x_v = 1,  ## unitless
    p_r = NA, ## unitless
    t_r = NA, ## unitless
    location = "exiting lateral"
  ),
  critical = list(
    p = 74.563, ## bar
    t = 345.36, ## K
    h = 5.1     ## kJ / mol
  )
)

for(i1 in 1:5){
  power_cycle[[i1]]$p_r <- power_cycle[[i1]]$p / power_cycle$critical$p
}

## loop

eavor_loop <- list(
  n_laterals = 12,
  depth = 2500, ## m
  width = 5000, ## m
  lateral_radius = 0.078, ## m
  vertical_radius = 0.28  ## m
)
