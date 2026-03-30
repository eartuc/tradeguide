import numpy as np
import time
import pandas as pd
import hdf5storage 

# Citation: Artuc, Erhan and Ortega, Johan (2026). "International Trade Policy and Quantitative Models:
# A Practitioner’s Guide," World Bank Policy Research Working Paper Series.

#X: trade matrix --> row: n(importer) x column: m(exporter)

# Load data from matlab format
d = hdf5storage.loadmat('data_tiva25_simple.mat')

#Setup
d["Niter"] = 5000 #max iterations
d["err_tol"] = 1e-6 #tolerance
d["rho"] = 0.2 #speed
d["theta"] = 4.0 #trade elasticity

#scalarize
d["DC"]=d["DC"].item()

# shocks
s = {}
s["deltahat"] = np.ones((d["DC"], d["DC"]))


###############################################
# Prepare for 10% trade cost reduction simulation

INDX = np.arange(1, d["DC"] + 1)

for i in INDX:
    i_ = np.setdiff1d(np.arange(1, d["DC"] + 1), i)
    s["deltahat"][i-1, i_-1] =  1.1  



#----------------------------------------------------------------------
#DEFINE FUNCTION FOR SOLUTION

def Solve_EK_simple(s, d):
    
    #number of countries
    DC = d["DC"]

    #shock
    deltahat = s["deltahat"]
    
    #paramaters and shares
    Niter = d["Niter"]
    err_tol = d["err_tol"]
    theta = d["theta"]

    Xbar = d["X"]
    pibar = d["pibar"]
    
    #speed
    rho1 = d["rho"] 
        
    #initialize    
    Xhat = np.ones((DC, DC))
    wagehat = np.ones(DC)
    
    iter_count = 0
    err_big = True
    err=100
    
    while iter_count <  Niter and err_big:
        
        iter_count += 1
        
        #SOLUTION STEPS BEGIN

        #Step 1: price including iceberg trade costs (Eq.1)
        costhat_nm = wagehat[np.newaxis, :]
        pnmhat = deltahat * costhat_nm
 
        #Step 2: price index (Eq.2)
        Phihat = np.sum(pibar * (pnmhat ** (-theta)), axis=1, keepdims=True) ** (-1/theta)
        
        #Step 3: trade flows (Eq.3)
        pihat = (pnmhat / Phihat) ** (-theta)
        
        #Step 4: market clearing and demand (Eq.4)                
        incomehat = np.tile(wagehat[:,np.newaxis], (1, DC))         
        Xhat = incomehat * pihat     
   
        #Step 5: new wages (Eq.5)
        wagehat_star=np.sum(Xbar*Xhat,axis=0)/np.sum(Xbar,axis=0)

        #Step 6: scale - need to scale guessed nominal variable
        scale=np.sum(Xhat*Xbar, axis=(0,1))/np.sum(Xbar, axis=(0,1))
        wagehat_star =  wagehat_star / scale
     
        #Step 7: error computation then update variables
        err = np.sum(np.abs(wagehat_star - wagehat))
        wagehat = rho1 * wagehat_star + (1 - rho1) * wagehat

        #Step 8: check convergence
        if iter_count % 20 == 0:        
            err_big = err > err_tol  
            print(f"Iteration {iter_count}, Errors: {err}")

        #SOLUTION STEPS END

    # GENERATE OUTPUTS
    rGDPhat=wagehat/Phihat[:,0]

    #percent changes   
    cGDP = 100*(rGDPhat- 1.0)

    return {
        'Xhat': Xhat,#trade flow
        'wagehat': wagehat, #wage
        'pihat': pihat, #import share
        'Phihat': Phihat, #price index
        'err': err, #final error
        'iter': iter_count,
        'scale': scale, #normalization otherwise price can go to infinity or zero
        'CPIhat': Phihat, #same as price index here as there are no sectors
        'rGDPhat': rGDPhat, #real gdp in hats
        'cGDP': cGDP #gdp change in percent
    }

#END FUNCTION FOR SOLUTION
#----------------------------------------------------------------------

# RUN SIMULATION
start_time = time.time()
r = Solve_EK_simple(s, d)
elapsed_time = time.time() - start_time
print(f"Elapsed time: {elapsed_time:.2f} seconds")
print(r["err"])

# SAVE RESULTS 
iso=[]
for i in range(0,d["DC"]):
    iso.append(d["countries"][i][0][0])
results = pd.DataFrame({ 'Country': iso})
results["GDP"]=r["cGDP"]

results.to_csv("results_ek.csv", index=False)