using MAT
using LinearAlgebra
using Printf
using DataFrames
using CSV

# Citation: Artuc, Erhan and Ortega, Johan (2026). "International Trade Policy and Quantitative Models:
# A Practitioner’s Guide," World Bank Policy Research Working Paper Series.

#betaM row country column sector producing using inputs from third dimension
#X: trade matrix --> row: n(importer) x column: m(exporter)  x 3rd dim: j(sector)

# Load data from matlab format
d = matread("data_tiva25_tariff.mat")

#Setup
d["Niter"] = 5000 #max iter
d["err_tol"] = 1.0e-6 #tolerance
d["rho"] = 0.5 #speed
d["theta"] = 4.0 #trade elas
d["nu"] = 2.0 #labor elas

#integerize
d["DC"]=Int(d["DC"])
d["DJ"]=Int(d["DJ"])


#------------------------------------------------------------
#SHOCK: High Income vs Non-High Income Tariff
s = Dict()
s["dhat"]=ones(d["DC"],d["DC"],d["DJ"])
s["taunew"]=d["taubar"]

INDX=vec(1:d["DC"])

for i=INDX
    for j=INDX
    taunew_=d["taubar"][i,j,1:d["DJ"]-1];
    if d["high"][i]==1 && d["high"][j]!=1
        taunew_ = d["taubar"][i,j,1:d["DJ"]-1] .+ 0.25    
    end
    if d["high"][i] != 1 && d["high"][j]==1
        taunew_ = d["taubar"][i,j,1:d["DJ"]-1] .+ 0.25    
    end

    taunew_[taunew_ .< 0.0].= 0.0
    s["dhat"][i,j,1:d["DJ"]-1] .= (1.0 .+ taunew_)./(1.0 .+ d["taubar"][i,j,1:d["DJ"]-1])
    s["taunew"][i,j,1:d["DJ"]-1].= taunew_

    end
end
# END SHOCK SETUP
#-------------------------------------------------------------------



#----------------------------------------------------------------------
#DEFINE FUNCTION FOR SOLUTION

function Solve_Model_AO(s, d)
    
    #number of countries and sectors
    DC = d["DC"]
    DJ = d["DJ"]
    
    #shock
    taunew = s["taunew"]
    dhat = s["dhat"]
    
    #parameters and shares
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

    savings=d["savings"]
    taurevshare=d["taurevshare"]
    incomebar=d["income"]
    
    #speed
    rho1 = d["rho"] * 0.5
    rho2 = d["rho"]
    rho3 = d["rho"]
    
    #production function
    alphabarL_ = alphabarL ./ (alphabarL + alphabarM)
    alphabarM_ = alphabarM ./ (alphabarL + alphabarM)
    
    #initialize
    Xhat = ones(DC, DC, DJ)
    Ymj_hat = ones(DC, DJ)
    costMhat = ones(DC, DJ)
    wagehat = ones(DC, DJ)
    incomehat=ones(DC,1)
    
    #for scope purposes
    scale=1.0
    Phihat = ones(DC, DJ)
    Lhat = ones(DC, DJ)
    What = ones(DC, 1)
    pihat = ones(DC,DC,DJ)

    iter_count = 0
    err_big = true
    err = 100.0
    
    while iter_count < Niter && err_big
        iter_count += 1
        
        #Step 1: production cost (Eq.6) - in different shapes
        costhat = (wagehat .^ alphabarL_) .* (costMhat .^ alphabarM_)
        costhat_nmj = permutedims(repeat(costhat, 1, 1, DC), [3, 1, 2]) 

        #Step 2: price including tariffs (Eq.7)
        pnmjhat = dhat .* costhat_nmj

        #Step 3: price index (Eq.8) - in different shapes
        Phihat_nmj = sum(pibar .* (pnmjhat .^ (-theta)), dims=2) .^ (-1/theta) 
        Phihat = reshape(Phihat_nmj, (DC, DJ))
        
        #Step 4: new intemediate input costs (Eq.9)
        costMhat_star = dropdims(prod( Phihat_nmj.^ betaM, dims=3), dims=3)


        #Step 5: trade flows (Eq.10)
        pihat = (pnmjhat ./ Phihat_nmj) .^ (-theta)

        #Step 6: labor income (Eq.11)
        What = sum(Lbar .* (wagehat .^ nu), dims=2) .^ (1/nu)

        #Step 7: labor allocation (Eq.12)
        Lhat = (wagehat ./ What) .^ nu
        
        #Step 8: sectoral output (Eq.13)
        Ymj_hat = dropdims( sum(Xibar .* Xhat, dims=1), dims=1)

        #Step 9: new wages (Eq.14)
        wagehat_star = Ymj_hat .* (Lhat .^ (-(nu - 1)/nu))
       
        #Step 10: market clearing and demand (Eqs.15,16,17)
        incomehat = What.*(1.0 .- taurevshare)+(sum(Xhat.*Xbar.*taunew,dims=(2,3))./incomebar)
        expenditurehat = (incomehat - savings) ./ (1.0 .- savings)

        pihat_ = permutedims(repeat(pihat, 1, 1, 1, DJ), (1, 2, 4, 3)) #alternative shape      

        cont_C = DbarL .* expenditurehat .* pihat #demand by consumers
        cont_Y = permutedims(reshape(Ymj_hat,(DC,DJ,1)), [1,3,2] ) #alternative shape
        cont_I = dropdims(sum(Dbar.*pihat_.*cont_Y,dims=3),dims=3) #demand by producers

        Xhat_star = (cont_C + cont_I)./dhat #demand corridor by corridor
        
        #Step 11: scale
        scale=sum(Xhat[:].*Xbar[:])/sum(Xbar[:])
        Xhat_star=Xhat_star./scale
        costMhat_star=costMhat_star./scale
        wagehat_star=wagehat_star./scale

        #Step 12: compute errors then update variables
        err1 = sum(abs.(Xhat_star - Xhat) )
        err2 = sum(abs.(costMhat_star - costMhat))
        err3 = sum(abs.(wagehat_star - wagehat))
        err = [err1, err2, err3]

        Xhat = rho1 * Xhat_star + (1 - rho1) * Xhat
        costMhat = rho2 * costMhat_star + (1 - rho2) * costMhat
        wagehat = rho3 * wagehat_star + (1 - rho3) * wagehat

        #Step 13: check convergence
        if iter_count % 20 == 0            
            err_big = err1 > err_tol || err2 > err_tol || err3 > err_tol
            @printf("Iteration %d, Errors: %s\n", iter_count, err)
        end

        #SOLUTION STEPS END

    end

    # GENERATE OUTPUTS
    CPIhat = prod(Phihat .^ gammaL, dims=2)
    rGDPhat=incomehat./CPIhat;
    rWhat=What./CPIhat;

    trade_mask = ones(DC, DC, DJ)
    for i = 1:DC
        trade_mask[i, i, :] .= 0.0
    end
    
    tmp1=Xhat.*d["X"].*trade_mask
    tmp2=d["X"].*trade_mask
    Exporthat=reshape(sum(tmp1,dims=(1, 3))./sum(tmp2, dims=(1, 3)), (DC,1) ) 
    Importhat=reshape(sum(tmp1, dims=(2, 3))./sum(tmp2, dims=(2, 3)), (DC,1) ) 

    #convert to percent   
    cGDP = 100*(rGDPhat .- 1.0)
    cIncome = 100*(rWhat .- 1.0)
    cExport = 100*(Exporthat./CPIhat .- 1.0)
    cImport = 100*(Importhat./CPIhat .- 1.0)

    return Dict(
        "Xhat" => Xhat, #trade flow change
        "Ymj_hat" => Ymj_hat, #sectoral output
        "wagehat" => wagehat, #wage per effective units
        "Lhat" => Lhat, #labor allocation
        "What" => What, #expected wage
        "costMhat" => costMhat, #material cost
        "pihat" => pihat, #import share change
        "Phihat" => Phihat, #price index
        "CPIhat" => CPIhat, #cpi
        "err" => err, #final error
        "iter" => iter_count, #number of iterations
        "scale"=>scale,
        "rGDPhat"=> rGDPhat, #real gdp change in hat
        "cImport"=> cImport, #in percent
        "cExport"=> cExport, #in percent
        "cIncome"=> cIncome, #in percent
        "cGDP"=> cGDP #in percent
    )
end

#END FUNCTION FOR SOLUTION
#----------------------------------------------------------------------

# RUN SIMULATION
start_time = time()
r = Solve_Model_AO(s, d)
elapsed_time = time() - start_time
@printf("Elapsed time: %.2f seconds\n", elapsed_time)
print(r["err"])

#SAVE RESULTS
results = DataFrame(d["countries"] ,["Country"])
results.DIncome = r["cIncome"][:,1]
results.GDP = r["cGDP"][:,1]
results.Export = r["cExport"][:,1]
results.Import = r["cImport"][:,1]

CSV.write("results_ao.csv", results)
