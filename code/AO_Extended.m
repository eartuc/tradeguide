%MAIN SCRIPT
clear

% Citation: Artuc, Erhan and Ortega, Johan (2026). "International Trade Policy and Quantitative Models:
% A Practitioner’s Guide," World Bank Policy Research Working Paper Series.

% Updated March 21, 2026

%betaM row country column sector producing using inputs from third dimension
%X: trade matrix --> row: n(importer) x column: m(exporter)  x 3rd dim: j(sector)

%load data
d=load("data_tiva25_tariff.mat");

%setup
d.Niter=5000; %iterations
d.err_tol=1e-6; %tolerance
d.rho=0.5; %speed
d.theta=4.0; %trade elas
d.nu=2.0; %labor elas

%------------------------------------------------------------------
%SHOCK: High Income vs Non-High Income Tariff
s.dhat=ones(d.DC,d.DC,d.DJ);
s.taunew=d.taubar;
INDX(1,:)=1:d.DC;
for i=INDX
    for j=INDX
    taunew_=d.taubar(i,j,1:d.DJ-1);
    if d.high(i)==1 && d.high(j)~=1
        taunew_=d.taubar(i,j,1:d.DJ-1) + 0.25;    
    end
    if d.high(i)~=1 && d.high(j)==1
        taunew_=d.taubar(i,j,1:d.DJ-1) + 0.25;    
    end

    taunew_(taunew_<0)=0;
    s.dhat(i,j,1:d.DJ-1) =(1+taunew_)./(1+d.taubar(i,j,1:d.DJ-1));
    s.taunew(i,j,1:d.DJ-1)=taunew_;

    end
end
% END SHOCK SETUP
%-------------------------------------------------------------------




% RUN SIMULATION
tic
r=Solve_Model_AO(s,d);
toc
r.err

% SAVE RESULTS
resultTable = table(d.countries, r.cIncome, r.cGDP, r.cExport, r.cImport ,'VariableNames', {'Country','DIncome', 'GDP', 'Export', 'Import'}); 
writetable(resultTable,'results_ao.csv');  

%----------------------------------------------------------------------
%DEFINE FUNCTION FOR SOLUTION 
%WARNING: SOME VERSIONS OF MATLAB REQUIRE THE DEFINITION TO BE IN THE END

function r = Solve_Model_AO(s,d)

%number of countries and sectors
DC=d.DC;
DJ=d.DJ;

%SHOCK
taunew=s.taunew;
dhat=s.dhat;

%PARAMETERS AND SHARES
Niter=d.Niter;
err_tol=d.err_tol;
theta=d.theta;
nu=d.nu;

Xibar=d.Xibar;
Xbar=d.X;
DbarL=d.DbarL;
Dbar=d.Dbar;

alphabarL=d.alphabarL; 
alphabarM=d.alphabarM;
betaM=d.betaM;
gammaL=d.gammaL;

savings=d.savings;
taurevshare=d.taurevshare;
incomebar=d.income;

pibar=d.pibar;
Lbar=d.Lbar;

%speed variables
rho1=d.rho*0.5;
rho2=d.rho;
rho3=d.rho;

%prod function
alphabarL_=(alphabarL)./(alphabarL+alphabarM);
alphabarM_=(alphabarM)./(alphabarL+alphabarM);

%allocate space
Xhat=(ones(DC,DC,DJ)); 
costMhat=(ones(DC,DJ));
wagehat=(ones(DC,DJ));

iter=0;
err_big=true;
err=100;

while (iter <  Niter && err_big )

    iter=iter+1;

    %SOLUTION STEPS BEGIN 

    %Step 1: production cost (Eq.6) - in different shapes
    costhat = ( wagehat.^(alphabarL_) ).*(costMhat.^(alphabarM_)); 
    costhat_nmj = permute(costhat ,[3,1,2]);

    %Step 2: price including tariffs (Eq.7)
    pnmjhat= dhat.*costhat_nmj; 
    
    %Step 3: price index (Eq.8) - in different shapes
    Phihat_nmj = ( sum( pibar.*(pnmjhat.^(-theta)) , 2 ) ).^(-(1/theta)); 
    Phihat = permute(Phihat_nmj,[1,3,2]); 

    %Step 4: new intermediate input costs (Eq.9)
    costMhat_star = prod( Phihat_nmj.^betaM,3); 
    
    %Step 5: trade flows (Eq.10)
    pihat = ( pnmjhat./Phihat_nmj ).^(-theta);

    %Step 6: labor income (Eq.11)
    What = ( sum( Lbar.*(wagehat.^nu),2 ) ).^(1/nu); 

    %Step 7: labor allocation (Eq.12)
    Lhat= (wagehat./What).^(nu);
        
    %Step 8: sectoral output (Eq.13)
    Ymj_hat=permute(sum(Xibar.*Xhat,1),[2,3,1]);

    %Step 9: new wages (Eq.14)
    wagehat_star=Ymj_hat.*( Lhat.^(-(nu-1)./nu));
    
    %Step 10: market clearing and demand (Eqs.15,16,17)
    incomehat = What.*(1-taurevshare)+(sum(Xhat.*Xbar.*taunew ,[2,3])./incomebar);
    expenditurehat = (incomehat - savings)./(1-savings);    

    pihat_ = permute(repmat(pihat,1,1,1,DJ),[1,2,4,3]); %alternative shape

    cont_C = DbarL.*expenditurehat.*pihat; %demand by consumers
    cont_Y = permute(Ymj_hat,[1,3,2]); %need alternative shape
    cont_I = permute( sum( Dbar.*pihat_.*cont_Y , 3), [1, 2, 4, 3]); %by producers

    Xhat_star = (cont_C + cont_I)./dhat; %demand corridor by corridor
    
    %Step 11: scale - guessed variables in nominal terms
    scale=sum(Xhat(:).*Xbar(:))/sum(Xbar(:));   
    Xhat_star=Xhat_star./scale;
    costMhat_star=costMhat_star./scale;
    wagehat_star=wagehat_star./scale;

    %Step 12: calculate errors then update variables
    err1= sum( abs(Xhat_star(:) - Xhat(:)));
    err2= sum( abs(costMhat_star(:) - costMhat(:)) );
    err3= sum( abs(wagehat_star(:) - wagehat(:)) );
    err=[err1 err2 err3 ];

    Xhat = rho1*Xhat_star + (1-rho1)*Xhat;
    costMhat = rho2*costMhat_star + (1-rho2)*costMhat;
    wagehat = rho3*wagehat_star + (1-rho3)*wagehat;

    %Step 13: check convergence
    if mod(iter,20)==0            
        err_big=(err1 > err_tol || err2 > err_tol || err3 > err_tol );
        disp(err)
    end

    %SOLUTION STEPS END

end %while

% GENERATE OUTPUTS
CPIhat(1:DC,1)=prod(Phihat.^gammaL,2);
rGDPhat=incomehat./CPIhat;
rWhat=What./CPIhat;

trade_mask = ones(DC, DC, DJ);
for i = 1:DC
    trade_mask(i, i, :) = 0.0;
end
    
tmp1=Xhat.*d.X.*trade_mask;
tmp2=d.X.*trade_mask;
Exporthat=reshape(sum(tmp1,[1, 3])./sum(tmp2, [1, 3]),DC,1); 
Importhat=reshape(sum(tmp1,[2, 3])./sum(tmp2, [2, 3]),DC,1); 

%percent change    
cGDP = 100*(rGDPhat- 1.0);
cIncome = 100*(rWhat- 1.0);
cExport = 100*(Exporthat./CPIhat- 1.0);
cImport = 100*(Importhat./CPIhat- 1.0);

r.Xhat=Xhat; %Trade
r.Ymj_hat=Ymj_hat; %Output
r.wagehat=wagehat; %Wages
r.Lhat=Lhat; %Labor allocation
r.What=What; %Expected income (wage bill)
r.costMhat=costMhat; %Cost of materials
r.pihat=pihat; %Import share
r.scale=scale; 
r.Phihat=Phihat; %Price index for consumers
r.CPIhat=CPIhat; %CPI
r.err=err; %final error
r.iter=iter; %number of iters
r.rGDPhat=rGDPhat; %real GDP change in hats
r.cIncome=cIncome; %disp. income change in percent
r.cGDP=cGDP; %real GDP change in percent
r.cExport=cExport; %in percent
r.cImport=cImport; %in percent

end

%END FUNCTION FOR SOLUTION
%----------------------------------------------------------------------


