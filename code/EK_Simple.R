rm(list = ls())

library(R.matlab)

# Citation: Artuc, Erhan and Ortega, Johan (2026). "International Trade Policy and Quantitative Models:
# A Practitioner’s Guide," World Bank Policy Research Working Paper Series.

#X: trade matrix --> row: n(importer) x column: m(exporter) 

# SET WORKING DIRECTORY HERE
#setwd('/home/ea/Dropbox/Research_new/WTM/CODE/Artuc_Ortega_2026/model_code')
#load data from matlab format
d=readMat('data_tiva25_simple.mat')

#setup
d$Niter <- 5000 #max iterations
d$err_tol <- 1e-6 #tolerance
d$rho <- 0.2 #speed
d$theta <- 4.0 #trade eleasticity


# Shocks
s <- list()
s$deltahat <- array(1, dim = c(d$DC, d$DC))

###############################################
# Prepare for 10% trade cost reduction simulation

INDX <- 1:d$DC
for (i in INDX) {
  i_ <- setdiff(1:d$DC, i)  
  s$deltahat[i, i_] <- 1.1
}

#----------------------------------------------------------------------
#DEFINE FUNCTION FOR SOLUTION

Solve_EK_simple <- function(s, d) {
  
  #number of countries
  DC <- d$DC
  
  #shock
  deltahat <- s$deltahat
  
  #parameters and shares
  Niter <- d$Niter
  err_tol <- d$err_tol
  theta <- d$theta
  
  Xbar <- d$X
  pibar <- d$pibar
  
  #speed
  rho1 <- d$rho 
  
  #initialize
  Xhat <- array(1, dim = c(DC, DC))
  wagehat <- array(1, dim = c(DC))
  
  
  iter <- 0
  err_big <- TRUE
  
  print("current error=")
 
  
  while (iter < Niter && err_big) {
    iter <- iter + 1
    
    #SOLUTION STEPS BEGIN
  
    #Step 1: price including iceberg trade costs (Eq.1)
    costhat_nm <- aperm( replicate(DC, wagehat), c(2, 1) )
    pnmhat <- deltahat * costhat_nm
    
    #Step 2: price index (Eq.2)
    Phihat <- (apply( pibar * (pnmhat^(-theta)), 1, sum))^(-1/theta)
    Phihat_nm <- replicate(DC, Phihat)
    
    #Step 3: trade flows (Eq.3)
    pihat <- (pnmhat / Phihat_nm)^(-theta)
    
    #Step 4: market clearing and demand (Eq.4)
    incomehat <- replicate(DC, wagehat)
    Xhat <- incomehat * pihat
    
    #Step 5: new wages (Eq.5)
    wagehat_star <- apply(Xbar*Xhat, 2, sum)/apply(Xbar, 2, sum)
    
    #Step 6: scale - need to scale guessed nominal variable
    scale <- sum(Xhat*Xbar)/sum(Xbar)
    wagehat_star <- wagehat_star / scale
    
    #Step 7: error computation then update variables
    err <- sum(abs(wagehat_star - wagehat))   
    wagehat <- rho1 * wagehat_star + (1 - rho1) * wagehat
    
    #Step 8: check convergence
    if (iter %% 20 == 0) {
      err_big <- (err > err_tol)
      print(err)
    }
    
    #SOLUTION STEPS END

  } #end of iteration
  
  # GENERATE OUTPUTS
  rGDPhat=wagehat/Phihat

  # percent changes
  cGDP    <- 100 * (rGDPhat - 1.0)
  
  list(
    Xhat = Xhat,  #trade flows
    wagehat = wagehat,
    pihat = pihat, #import share
    Phihat = Phihat, #price index
    err = err, #error in the end
    iter = iter, #how many iterations it took
    scale = scale, #normalization otherwise price can go to infinity or zero
    CPIhat = Phihat, #same as price index
    rGDPhat=rGDPhat, #real gdp in hat
    cGDP=cGDP #gdp change in percent
  )
}

#END FUNCTION FOR SOLUTION
#----------------------------------------------------------------------

# RUN SIMULATION


ptm <- proc.time()
r <-Solve_EK_simple(s,d)
elapsed <- proc.time() - ptm
print(elapsed)
print(r$err)

# SAVE RESULTS
results<-data.frame( Country=unlist(d$countries), GDP=r$cGDP)
write.csv(results, "results_ek.csv", row.names=FALSE)
