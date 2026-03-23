import numpy as np
import time
import pandas as pd
import hdf5storage

# Citation: Artuc, Erhan and Ortega, Johan (2026). "International Trade Policy and Quantitative Models:
# A Practitioner’s Guide," World Bank Policy Research Working Paper Series.

#betaM row country column sector producing using inputs from third dimension
#X: trade matrix --> row: n(importer) x column: m(exporter) x 3rd dim: j(sector)

# Load data from matlab format
d = hdf5storage.loadmat('data_tiva25_tariff.mat')

#Setup
d["Niter"] = 5000 #max iterations
d["err_tol"] = 1e-6 #tolerance
d["rho"] = 0.5 #speed
d["theta"] = 4.0 #trade elasticity
d["nu"] = 2.0 #labor elasticity

#scalarize
d["DC"]=d["DC"].item()
d["DJ"]=d["DJ"].item()


#------------------------------------------------------------------
#SHOCK: High Income vs Non-High Income Tariff Hike
s = {}
s["dhat"]=np.ones((d["DC"],d["DC"],d["DJ"]))
s["taunew"]=d["taubar"]
INDX = np.arange(1, d["DC"]+1)
for i in INDX:
    for j in INDX:

        taunew_=d["taubar"][i-1,j-1, 0:d["DJ"]-1]
        
        if d["high"][i-1]==1.0 and d["high"][j-1] != 1.0:
            taunew_=d["taubar"][i-1, j-1, 0:d["DJ"]-1 ] + 0.25        
        
        if d["high"][i-1]!=1.0 and d["high"][j-1]==1.0:
            taunew_=d["taubar"][i-1, j-1, 0:d["DJ"]-1 ] + 0.25        
        
        taunew_[taunew_ < 0.0] = 0.0
    
        s["dhat"][i-1, j-1, 0:d["DJ"]-1] = (1.0 + taunew_) / (1.0 + d["taubar"][i-1,j-1, 0:d["DJ"]-1] )
        s["taunew"][i-1,  j-1, 0:d["DJ"]-1]=taunew_

# END SHOCK SETUP
#-------------------------------------------------------------------



#----------------------------------------------------------------------
#DEFINE FUNCTION FOR SOLUTION

def Solve_Model_AO(s, d):

    #number of countries and sectors
    DC = d["DC"]
    DJ = d["DJ"]
    
    #shocks
    taunew = s["taunew"]
    dhat = s["dhat"]
    
    #params and shares
    Niter = d["Niter"]
    err_tol = d["err_tol"]
    theta = d["theta"]
    nu = d["nu"]
    
    Xibar = d["Xibar"]
    Xbar = d["X"]
    DbarL = d["DbarL"]
    Dbar = d["Dbar"]

    alphabarL = d["alphabarL"]
    alphabarM = d["alphabarM"]
    betaM = d["betaM"]
    gammaL = d["gammaL"]
    
    pibar = d["pibar"]
    Lbar = d["Lbar"]

    savings=d["savings"].reshape(DC,)
    taurevshare=d["taurevshare"].reshape(d["DC"],)
    incomebar=d["income"].reshape(DC,)
    
    #speed
    rho1 = d["rho"] * 0.5
    rho2 = d["rho"]
    rho3 = d["rho"]
    
    #production function
    alphabarL_ = alphabarL / (alphabarL + alphabarM)
    alphabarM_ = alphabarM / (alphabarL + alphabarM)
    
    #initialize
    Xhat = np.ones((DC, DC, DJ))
    costMhat = np.ones((DC, DJ))
    wagehat = np.ones((DC, DJ))
     
    iter_count = 0
    err_big = True
    err=100
    
    while iter_count <  Niter and err_big:
        iter_count += 1
        
        #SOLUTION STEPS BEGIN 

        #Step 1: production cost (Eq.6) - in different shapes
        costhat = (wagehat ** alphabarL_) * (costMhat ** alphabarM_)
        costhat_ = costhat[np.newaxis, :, :]
 
        #Step 2: price including tariffs (Eq.7)
        pnmjhat = dhat * costhat_
 
        #Step 3: price index (Eq.8) - in different shapes
        Phihat = np.sum(pibar * (pnmjhat ** (-theta)), axis=1) ** (-1/theta)
        Phihat_nmj = Phihat[:, np.newaxis, :]

        #Step 4: new intermediate input costs (Eq.9)
        costMhat_star = np.prod( Phihat_nmj**betaM , axis=2)

        #Step 5: trade flows (Eq.10)
        pihat = (pnmjhat / Phihat_nmj) ** (-theta)
        
        #Step 6: labor income (Eq.11)
        What = np.sum(Lbar * (wagehat ** nu), axis=1) ** (1/nu)

        #Step 7: labor allocation (Eq.12)
        Lhat = (wagehat / What[:,np.newaxis]) ** nu

        #Step 8: sectoral output (Eq.13)
        Ymj_hat = np.sum(Xibar * Xhat, axis=0)

        #Step 9: new wages (Eq.14)
        wagehat_star = Ymj_hat * (Lhat ** (-(nu - 1)/nu))

        #Step 10: market clearing and demand (Eqs.15,16,17)
        incomehat = What*(1.0-taurevshare)+(np.sum( Xhat*Xbar*taunew,axis=(1,2))/incomebar)
        expenditurehat = (incomehat - savings) / (1 - savings)

        pihat_ = pihat[:, :, np.newaxis,:] #alternative shape
        expenditurehat_ = expenditurehat[:,np.newaxis,np.newaxis] #alternative shape
        
        cont_C = DbarL * expenditurehat_ * pihat #demand by consumers 
        cont_Y = Ymj_hat[:, np.newaxis, :, np.newaxis] #alternative shape
        cont_I = np.sum(Dbar * pihat_ * cont_Y, axis=2) #demand by producers 

        Xhat_star = (cont_C + cont_I )/dhat #change in demand corridor by corridor
        
        #Step 11: scale - all guessed nominal variables
        scale=np.sum(Xhat*Xbar, axis=(0,1,2))/np.sum(Xbar, axis=(0,1,2))
        Xhat_star=Xhat_star/scale
        costMhat_star=costMhat_star/scale
        wagehat_star=wagehat_star/scale

        # Step 12: error computation then update variables
        err1 = np.sum(np.abs(Xhat_star - Xhat))
        err2 = np.sum(np.abs(costMhat_star - costMhat))
        err3 = np.sum(np.abs(wagehat_star - wagehat))
        err = [err1, err2, err3]

        Xhat = rho1 * Xhat_star + (1 - rho1) * Xhat
        costMhat = rho2 * costMhat_star + (1 - rho2) * costMhat
        wagehat = rho3 * wagehat_star + (1 - rho3) * wagehat

        # Step 13: check convergence
        if iter_count % 20 == 0:
            err_big = err1 > err_tol or err2 > err_tol or err3 > err_tol
            print(f"Iteration {iter_count}, Errors: {err}")

        #SOLUTION STEPS END

    # GENERATE OUTPUTS
    CPIhat = np.prod(Phihat ** gammaL, axis=1)
    rGDPhat=incomehat/CPIhat 
    rWhat = What/CPIhat

    trade_mask = np.ones((DC, DC, DJ))
    i = np.arange(DC)
    trade_mask[i, i, :] = 0.0
    
    tmp1=Xhat*d["X"]*trade_mask
    tmp2=d["X"]*trade_mask
    Exporthat=tmp1.sum(axis=(0, 2))/tmp2.sum(axis=(0, 2)) 
    Importhat=tmp1.sum(axis=(1, 2))/tmp2.sum(axis=(1, 2)) 

    #percent change   
    cIncome = 100*(rWhat- 1.0)
    cGDP = 100*(rGDPhat- 1.0)
    cExport = 100*(Exporthat/CPIhat- 1.0)
    cImport = 100*(Importhat/CPIhat- 1.0) 

    return {
        'Xhat': Xhat, #trade flows
        'Ymj_hat': Ymj_hat, #sectoral output
        'wagehat': wagehat, #wage per effective unit
        'Lhat': Lhat, #labor alloc
        'What': What, #expected wage
        'costMhat': costMhat, #material cost
        'pihat': pihat, #import share
        'Phihat': Phihat, #price index
        'CPIhat': CPIhat, #cpi
        'err': err, #final error
        'iter': iter_count, #iter count
        'scale': scale,         
        'rGDPhat': rGDPhat, #real gdp in hat
        'cImport': cImport, #in percent
        'cExport': cExport, #in percent
        'cIncome': cIncome, #in percent
        'cGDP': cGDP #in percent
    }


#END FUNCTION FOR SOLUTION
#----------------------------------------------------------------------


# RUN SIMULATION
start_time = time.time()
r = Solve_Model_AO(s, d)
elapsed_time = time.time() - start_time
print(f"Elapsed time: {elapsed_time:.2f} seconds")
print(r["err"])

# SAVE RESULTS 
iso=[]
for i in range(0,d["DC"]):
    iso.append(d["countries"][i][0][0])
results = pd.DataFrame({ 'Country': iso})
results["DIncome"]=r["cIncome"]
results["GDP"]=r["cGDP"]
results["Export"]=r["cExport"]
results["Import"]=r["cImport"]

results.to_csv("results_ao.csv", index=False)
