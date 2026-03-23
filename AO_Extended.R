rm(list = ls())

library(R.matlab)

# Citation: Artuc, Erhan and Ortega, Johan (2026). "International Trade Policy and Quantitative Models:
# A Practitioner’s Guide," World Bank Policy Research Working Paper Series.

#betaM row country column sector producing using inputs from third dimension
#X: trade matrix --> row: n(importer) x column: m(exporter)  x 3rd dim: j(sector)

# SET WORKING DIRECTORY HERE
#setwd('/home/ea/Dropbox/Research_new/WTM/CODE/Artuc_Ortega_2026/model_code')
#load data from matlab format
d=readMat('data_tiva25_tariff.mat')

#setup
d$Niter <- 5000 #max iter
d$err_tol <- 1e-6 #tolerance 
d$rho <- 0.5 #speed
d$theta <- 4.0 #trade elasticity
d$nu <- 2.0 # labor elasticity


###############################################
# Shocks

DC <- d$DC
DJ <- d$DJ

#------------------------------------------------------------------
#SHOCK: High Income vs Non-High Income Tariff Hike
s <- list(
  dhat = array(1, dim = c(DC, DC, DJ)),
  taunew = d$taubar
)

INDX <- 1:DC

for (i in INDX) {
  for (j in INDX){
    
    taunew_ <- d$taubar[i, j, 1:DJ-1]
    if ( d$high[i] == 1 & d$high[j] != 1) {
      taunew_=d$taubar[i, j, 1:d$DJ-1 ] + 0.25
    } 
    if ( d$high[i] != 1 & d$high[j] == 1) {
      taunew_=d$taubar[i, j, 1:d$DJ-1 ] + 0.25
    } 
    
    taunew_[taunew_ < 0] <- 0
    s$dhat[i, j, 1:DJ-1] <- (1 + taunew_) / (1 + d$taubar[i, j, 1:DJ-1])
    s$taunew[i, j, 1:DJ-1] <- taunew_
  }
}

# END SHOCK SETUP
#-------------------------------------------------------------------



#----------------------------------------------------------------------
#DEFINE FUNCTION FOR SOLUTION

Solve_Model_AO <- function(s, d) {
  
  # number of countries and sectors
  DC <- d$DC
  DJ <- d$DJ
  
  #shock
  taunew <- s$taunew
  dhat <- s$dhat
  
  #parameters and shares
  Niter <- d$Niter
  err_tol <- d$err_tol
  theta <- d$theta
  nu <- d$nu
  
  Xibar <- d$Xibar
  Xbar <- d$X
  Ymj_alphaL=d$Ymj.alphaL
  DbarL <- d$DbarL
  Dbar <- d$Dbar
  alphabarL <- d$alphabarL
  alphabarM <- d$alphabarM
  betaM <- d$betaM
  gammaL <- d$gammaL
  
  savings <- c(d$savings)
  taurevshare <- c(d$taurevshare)
  incomebar <- c(d$income)
  
  pibar <- d$pibar
  Lbar <- d$Lbar
  
  #speed
  rho1 <- d$rho * 0.5
  rho2 <- d$rho
  rho3 <- d$rho
  
  # proeduction function
  alphabarL_ <- alphabarL / (alphabarL + alphabarM)
  alphabarM_ <- alphabarM / (alphabarL + alphabarM)
  
  #initialize
  Xhat <- array(1, dim = c(DC, DC, DJ))
  costMhat <- array(1, dim = c(DC, DJ))
  wagehat <- array(1, dim = c(DC, DJ))
  
  iter <- 0
  err_big <- TRUE
  
  print("current error=")
 
  while (iter < Niter && err_big) {
    iter <- iter + 1
    
    #SOLUTION STEPS BEGIN 
  
    #Step 1: production cost (Eq.6) - in different shapes
    costhat <- (wagehat^alphabarL_) * (costMhat^alphabarM_)
    costhat_nmj <- aperm( replicate(DC, costhat), c(3, 1, 2) )

    #Step 2: price including tariffs (Eq.7)
    pnmjhat <- dhat * costhat_nmj
    
    #Step 3: price index (Eq.8) - in different shapes
    Phihat <- (apply( pibar * (pnmjhat^(-theta)), c(1,3), sum))^(-1/theta)
    Phihat_nmj <- aperm( replicate(DC, Phihat), c(1, 3, 2))
    
    #Step 4: new intermediate input costs (Eq.9)
    costMhat_star <- apply((aperm(replicate(DJ, Phihat),c(1,3,2)))^betaM,c(1,2),prod)
    
    #Step 5: trade flows (Eq.10)
    pihat <- (pnmjhat / Phihat_nmj)^(-theta)
    
    #Step 6: labor income (Eq.11)
    What <- (rowSums(Lbar * (wagehat^nu)))^(1/nu)
    
    #Step 7: labor allocation (Eq.12)
    Lhat <- (wagehat / replicate(DJ, What))^(nu)
    
    #step 8: sectoral output (Eq.13)
    Ymj_hat <- apply(Xibar * Xhat, c(2, 3), sum)

    #step 9: new wages (Eq.14)
    wagehat_star <- Ymj_hat * (Lhat^(-(nu - 1)/nu))

    #step 10: market clearing and demand (Eq.15,16,17)
    incomehat <- What*(1.0-taurevshare)+(apply(Xhat*Xbar*taunew,1,sum)/incomebar)
    expenditurehat <- (incomehat - savings) / (1 - savings)
    
    pihat_ <- aperm(replicate(DJ,pihat),c(1, 2, 4, 3) ) #alt shape
    expenditurehat_ <- replicate(DJ,replicate(DC,expenditurehat)) #alt shape
    
    cont_C <- DbarL * expenditurehat_ * pihat #demand by consumers
    cont_Y <- aperm(replicate(DC, Ymj_hat ),c(1, 3, 2)) #alternative shape
    cont_Y <- replicate(DJ, cont_Y) #final needed shape
    cont_I <- apply(Dbar * pihat_ * cont_Y, c(1, 2, 4), sum) #demand by producers
    
    Xhat_star <- (cont_C + cont_I)/dhat #demand corridor by corridor
       
    #Step 11: scale - all guessed nominal vars
    scale <- sum(Xhat*Xbar)/sum(Xbar)
    Xhat_star <- Xhat_star/scale
    costMhat_star <- costMhat_star/scale
    wagehat_star <- wagehat_star/scale
    
    #step 12: computer errors then update variables
    err1 <- sum(abs(Xhat_star - Xhat) )
    err2 <- sum(abs(costMhat_star - costMhat))
    err3 <- sum(abs(wagehat_star - wagehat))
    err <- c(err1, err2, err3)
    
    Xhat <- rho1 * Xhat_star + (1 - rho1) * Xhat
    costMhat <- rho2 * costMhat_star + (1 - rho2) * costMhat
    wagehat <- rho3 * wagehat_star + (1 - rho3) * wagehat
    
    #step 13: check convergence
    if (iter %% 20 == 0) {
      err_big <- (err1 > err_tol || err2 > err_tol || err3 > err_tol)
      print(err)
    }
    
    #SOLUTION STEPS END

  } #end of iteration
  
  # GENERATE OUTPUTS
  CPIhat <- apply(Phihat^gammaL, 1, prod)
  rGDPhat <- incomehat / CPIhat       
  rWhat      <- What / CPIhat

  trade_mask <- array(1, dim = c(DC, DC, DJ))
  for (i in 1:DC) {
    trade_mask[i, i, ] <- 0.0
  }
  
  tmp1 <- Xhat * d$X * trade_mask
  tmp2 <- d$X * trade_mask
  
  sum_tmp1 <- apply(tmp1, c(2), sum)
  sum_tmp2 <- apply(tmp2, c(2), sum)
  Exporthat   <- matrix(sum_tmp1 / sum_tmp2, ncol = 1)
  
  sum_tmp1 <- apply(tmp1, c(1), sum)
  sum_tmp2 <- apply(tmp2, c(1), sum)
  Importhat   <- matrix(sum_tmp1 / sum_tmp2, ncol = 1)
  
  #percent change
  cGDP    <- 100 * (rGDPhat - 1.0)
  cIncome <- 100 * (rWhat      - 1.0)
  cExport <- 100 * (Exporthat/CPIhat  - 1.0) 
  cImport <- 100 * (Importhat/CPIhat  - 1.0) 
  
  list(
    Xhat = Xhat, #trade
    Ymj_hat = Ymj_hat, #sectoral output
    wagehat = wagehat, #wage (per effective unit)
    Lhat = Lhat, #labor allocation
    What = What, #expected wage
    costMhat = costMhat, #material cost
    pihat = pihat, #import share
    Phihat = Phihat, #price index
    CPIhat = CPIhat, #cpi
    err = err,
    iter = iter,
    scale = scale,
    rGDPhat=rGDPhat, #real gdp in hat
    cIncome=cIncome, #in percent
    cGDP=cGDP, #in percent
    cExport=cExport, #in percent
    cImport=cImport #in percent
  )
}

#END FUNCTION FOR SOLUTION
#----------------------------------------------------------------------

# RUN SIMULATION
ptm <- proc.time()
r <-Solve_Model_AO(s,d)
elapsed <- proc.time() - ptm
print(elapsed)
print(r$err)

# SAVE RESULTS
results<-data.frame( Country=unlist(d$countries), DIncome=r$cIncome, GDP=r$cGDP, Export=r$cExport, Import=r$cImport)
write.csv(results, "results_ao.csv", row.names=FALSE)
