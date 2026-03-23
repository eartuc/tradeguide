using MAT
using LinearAlgebra
using Printf
using DataFrames
using CSV

# Citation: Artuc, Erhan and Ortega, Johan (2026). "International Trade Policy and Quantitative Models:
# A Practitioner’s Guide," World Bank Policy Research Working Paper Series.

#X: trade matrix --> row: n(importer) x column: m(exporter)

# Load data from matlab format
d= matread("data_tiva25_simple.mat")

#Setup
d["Niter"] = 5000 #iterations
d["err_tol"] = 1.0e-6 #tolerance
d["rho"] = 0.2 #speed
d["theta"] = 4.0 #trade eleasticity

#integerize
d["DC"]=Int(d["DC"])

# Shocks
s = Dict()
s["deltahat"] = ones(d["DC"], d["DC"])

###############################################
# Prepare for 10% trade cost reduction simulation
INDX = vec(1:d["DC"])

for i in INDX
    i_ = setdiff(1:d["DC"], i)
    s["deltahat"][i, i_, :] .= 1.1  # Julia uses 1-based indexing
end


#----------------------------------------------------------------------
#DEFINE FUNCTION FOR SOLUTION

function Solve_EK_simple(s, d)

    #number of countries
    DC = d["DC"]
    
    #shock
    deltahat = s["deltahat"]
    
    #parameters and shares
    Niter = d["Niter"]
    err_tol = d["err_tol"]
    theta = d["theta"]
     
    Xbar = d["X"]
    pibar = d["pibar"]
    
    #speed
    rho1 = d["rho"]

    #initialize
    Xhat = ones(DC, DC)
    wagehat = ones(DC, 1)
    
    #for scope purposes
    scale=1
    Phihat = ones(DC, 1)
    pihat = ones(DC,DC)
    
    iter_count = 0
    err_big = true
    err = 100.0
    
    while iter_count < Niter && err_big
        
        iter_count += 1

        #SOLUTION STEPS BEGIN

        #Step 1: price including iceberg trade costs (Eq.1)
        costhat_nm = permutedims(wagehat, [2, 1]) 
        pnmhat = deltahat .* costhat_nm

        #Step 2: price index (Eq.2)
        Phihat = sum(pibar .* (pnmhat .^ (-theta)), dims=2) .^ (-1/theta)
        
        #Step 3: trade flows (Eq.3)
        pihat = (pnmhat ./ Phihat) .^ (-theta)

        #Step 4: market clearing and demand (Eq.4)
        incomehat = repeat(wagehat, 1, DC)
        Xhat = incomehat .* pihat

        #Step 5: new wages (Eq.5)
        wagehat_star=reshape(sum(Xbar.*Xhat,dims=1)./sum(Xbar,dims=1),(DC,1))
        
        #Step 6: scale - need to scale guessed nominal variable
        scale=sum(Xhat[:].*Xbar[:])/sum(Xbar[:])
        wagehat_star = wagehat_star/scale

        #Step 7: error computation then update variables
        err = sum(abs.(wagehat_star - wagehat))        
        wagehat = rho1 * wagehat_star + (1 - rho1) * wagehat
        
        #Step 8: check convergence
        if iter_count % 20 == 0
            err_big = err > err_tol 
            @printf("Iteration %d, Errors: %s\n", iter_count, err)
        end

        #SOLUTION STEPS END

    end

    #generate outputs
    rGDPhat=wagehat./Phihat
    
    #percent changes   
    cGDP = 100*(rGDPhat .- 1.0)


    return Dict(
        "Xhat" => Xhat, #trade flows
        "wagehat" => wagehat, #wages
        "pihat" => pihat, #import shares
        "Phihat" => Phihat, #price index
        "err" => err, #final error
        "iter" => iter_count,
        "scale" => scale, #normalization otherwise price can go to infinity or zero
        "CPIhat" => Phihat, #cpi change same as price index
        "rGDPhat"=>rGDPhat, #real gdp in hat
        "cGDP"=>cGDP #in percent
    )
end

#END FUNCTION FOR SOLUTION
#----------------------------------------------------------------------

# RUN SIMULATION

start_time = time()
r = Solve_EK_simple(s, d)
elapsed_time = time() - start_time
@printf("Elapsed time: %.2f seconds\n", elapsed_time)
print(r["err"])

#SAVE RESULTS
results = DataFrame(d["countries"] ,["Country"])
results.GDP = r["cGDP"][:,1]

CSV.write("results_ek.csv", results)

