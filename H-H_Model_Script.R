# Hodgkin-Huxley Nerve Cell Model 
# See https://en.wikipedia.org/wiki/Hodgkin%E2%80%93Huxley_model
#################################
# Clear environment
cat("\f") # clears console
rm(list = ls()) # clears environment

# Load packages
pkgNames <- c("simecol", "zoo")
for (pkgName in pkgNames)
{
  if(!require(pkgName, character.only = TRUE))
  {
    install.packages(pkgName)
  }
  try(library(pkgName, character.only = TRUE))
}

# Model 
i <- 0
V <- -70
while ( V < 0 ) 
{
  simTime <- 100 # in ms
  delta_T <- 0.01 # time interval in ms
  
  t <- seq(from = 0, to = simTime, by = delta_T) # set 0.01ms time intervals from 0-100ms
  
  # External current 
  changeTimes <- 0 # in ms
  currentLevels <- i # adjust for different current levels
  
  I <- NA # initialize empty vector for storing to-be applied current levels
  I[1:500] <- currentLevels 
  I[501:2000] <- 0
  I[2001:length(t)] <- currentLevels
  
  # Constant parameters
  # refer to table 3 (Hodkin & Huxley, 1952)
  Gbar_Na <- 120; Gbar_K <- 36; Gbar_l <- 0.3
  C_M <- 1; V_Na <- 115; V_K <- -12; V_l <- 10.6
  
  ls_HH <- rep(list(list(list())), 4) 
  names(ls_HH) <- c("V", "Coefficients", "Currents", "Derivates")
  
  # Initial states (at t = 1)
  ls_HH$V <- rep(0, 10001)
  ls_HH$V[1] <- 0 # voltage at baseline
  
  ls_HH[[2]]$alpha_n <- rep(0, 10001); ls_HH[[2]]$alpha_n[1] <- .01 * ( (10-ls_HH$V[1]) / (exp((10-ls_HH$V[1])/10)-1) )
  ls_HH[[2]]$beta_n <- rep(0, 10001); ls_HH[[2]]$beta_n[1] <- .125*exp(-ls_HH$V[1]/80)
  ls_HH[[2]]$alpha_m <- rep(0, 10001); ls_HH[[2]]$alpha_m[1] <- .1*( (25-ls_HH$V[1]) / (exp((25-ls_HH$V[1])/10)-1) ) 
  ls_HH[[2]]$beta_m <- rep(0, 10001); ls_HH[[2]]$beta_m[1] <- 4*exp(-ls_HH$V[1]/18)
  ls_HH[[2]]$alpha_h <- rep(0, 10001); ls_HH[[2]]$alpha_h[1] <- .07*exp(-ls_HH$V[1]/20)
  ls_HH[[2]]$beta_h <- rep(0, 10001); ls_HH[[2]]$beta_h[1] <- 1/(exp((30-ls_HH$V[1])/10)+1) 
  
  ls_HH[[2]]$n_inf <- rep(0, 10001)
  ls_HH[[2]]$m_inf <- rep(0, 10001) 
  ls_HH[[2]]$h_inf <- rep(0, 10001) 
  ls_HH[[2]]$n_inf[1] <- ls_HH[[2]]$alpha_n[1] / (ls_HH[[2]]$alpha_n[1] + ls_HH[[2]]$beta_n[1]) 
  ls_HH[[2]]$m_inf[1] <- ls_HH[[2]]$alpha_m[1] / (ls_HH[[2]]$alpha_m[1] + ls_HH[[2]]$beta_m[1]) 
  ls_HH[[2]]$h_inf[1] <- ls_HH[[2]]$alpha_h[1] / (ls_HH[[2]]$alpha_h[1] + ls_HH[[2]]$beta_h[1]) 
  
  # Loop over coefficients, currents, and derivatives at each time interval
  for ( j in 1:(length(t) - 1) ) 
  {
    # Calculate coefficients
    ls_HH[[2]]$alpha_n[j] <- .01 * ( (10-ls_HH$V[j]) / (exp((10-ls_HH$V[j])/10)-1) ) 
    ls_HH[[2]]$beta_n[j] <- .125*exp(-ls_HH$V[j]/80) 
    ls_HH[[2]]$alpha_m[j] <- .1*( (25-ls_HH$V[j]) / (exp((25-ls_HH$V[j])/10)-1) ) 
    ls_HH[[2]]$beta_m[j] <- 4*exp(-ls_HH$V[j]/18) 
    ls_HH[[2]]$alpha_h[j] <- .07*exp(-ls_HH$V[j]/20) 
    ls_HH[[2]]$beta_h[j] <- 1/(exp((30-ls_HH$V[j])/10)+1) 
    
    # Calculate currents
    ls_HH[[3]]$I_Na[j] <- (ls_HH[[2]]$m_inf[j]^3) * Gbar_Na * ls_HH[[2]]$h_inf[j] * (ls_HH$V[j]-V_Na) 
    ls_HH[[3]]$I_K[j] <- (ls_HH[[2]]$n_inf[j]^4) * Gbar_K * (ls_HH$V[j]-V_K) 
    ls_HH[[3]]$I_L[j] <- Gbar_l *(ls_HH$V[j]-V_l) 
    ls_HH[[3]]$I_ion[j] <- I[j] - ls_HH[[3]]$I_K[j] - ls_HH[[3]]$I_Na[j] - ls_HH[[3]]$I_L[j] 
    
    # Calculate derivatives 
    ls_HH$V[j  +1] <- ls_HH$V[j] + delta_T * (ls_HH[[3]]$I_ion[j] / C_M)
    ls_HH[[2]]$n_inf[j + 1] <- ls_HH[[2]]$n_inf[j] + delta_T*(ls_HH[[2]]$alpha_n[j] *(1-ls_HH[[2]]$n_inf[j]) - ls_HH[[2]]$beta_n[j] * ls_HH[[2]]$n_inf[j])
    ls_HH[[2]]$m_inf[j + 1] <- ls_HH[[2]]$m_inf[j] + delta_T*(ls_HH[[2]]$alpha_m[j] *(1-ls_HH[[2]]$m_inf[j]) - ls_HH[[2]]$beta_m[j] * ls_HH[[2]]$m_inf[j])
    ls_HH[[2]]$h_inf[j + 1] <- ls_HH[[2]]$h_inf[j] + delta_T*(ls_HH[[2]]$alpha_h[j] *(1-ls_HH[[2]]$h_inf[j]) - ls_HH[[2]]$beta_h[j] * ls_HH[[2]]$h_inf[j])
  
  }
  ls_HH$V <- ls_HH$V - 70 # set resting potential to -70mV
  
  # Plots
  par(mfrow = c(1, 2))
  
  # jpeg("Voltage Plot.jpg")
  # Plot voltage
  plot(t, ls_HH$V, 
       type = "l",
       col = "blue",
       ylim = c(-90, ifelse(max(ls_HH$V) <= 30, 30, max(ls_HH$V))),
       xlab = "Time (ms)",
       ylab = "Voltage (mV)")
  par(lwd = 1)
  abline(h = max(ls_HH$V), lty = 2)
  text(80, (max(ls_HH$V) + 2), "Peak action potential")
  abline(h = 0, lty = 2)
  abline(h = -61, lty = 2)
  text(80, -59, "Threshold of excitation")
  text(15, ifelse(max(ls_HH$V) <= 30, 30, max(ls_HH$V)), paste0("Current level: ", currentLevels))
  
  # Plot conductance for sodium and potassium
  # jpeg("Conductance for Sodium and Potassium Plot.jpg")
  g_K <- Gbar_K *ls_HH[[2]]$n_inf^4
  g_Na <- Gbar_Na * (ls_HH[[2]]$m_inf^3) * ls_HH[[2]]$h_inf
  
  if ( max(g_K) > max(g_Na))
  {
    y_max <- max(g_K)
  } else {
    y_max <- max(g_Na)
  }
  plot(t, g_K, 
       type = "l",
       col = "blue",
       ylim = c(0, ifelse(max(y_max) <= 30, 30, max(y_max))), 
       xlab = "Time (ms)",
       ylab = "Conductance (g)")
  par(lwd = 1)
  lines(t, g_Na, 
        type = "l",
        col = "red",
        ylim = c(0, ifelse(max(y_max) <= 30, 30, max(y_max))), 
        xlab = "Time (ms)",
        ylab = "Conductance (g)")
  par(lwd = 1)
  legend(60, ifelse(max(y_max) <= 30, 30, max(y_max)), 
         c("Potassium", "Sodium"),
         lty = c(1, 1),
         lwd = c(1, 1),
         col = c("blue", "red"))
  text(15, ifelse(max(y_max) <= 30, 30, max(y_max)), 
       paste0("Current level: ", currentLevels))
  
  V <- max(ls_HH$V)
  i <- i + 0.1
  print(paste0("Current level: ", currentLevels, "; Max Voltage: ", V))
}
# Find local maxima
localMaxima <- function(x)
{
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if ( x[[1]] == x[[2]]) 
  {
    y <- y[-1]
  }
  y
}

V_Smooth <- smooth(ls_HH$V) # smooth voltage
