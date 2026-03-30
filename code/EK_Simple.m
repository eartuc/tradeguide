%MAIN SCRIPT
clear

% Citation: Artuc, Erhan and Ortega, Johan (2026). "International Trade Policy and Quantitative Models:
% A Practitioner’s Guide," World Bank Policy Research Working Paper Series.

%X: trade matrix --> row: n(importer) x column: m(exporter)

d=load("data_tiva25_simple.mat");

d.Niter=5000; %iterations
d.err_tol=1e-6; %tolerance
d.rho=0.2; %speed
d.theta=4.0; %trade elasticity

%shocks
s.deltahat=ones(d.DC,d.DC);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%prepare for 10% trade cost reduction
INDX(1,:)=1:d.DC;

for i=INDX

    i_=1:d.DC;
    i_(i)=[];

    s.deltahat(i,i_) = 1.1;
end


% RUN SIMULATION


tic
r=Solve_EK_simple(s,d);
toc
r.err

% SAVE RESULTS
resultTable = table(d.countries, r.cGDP, 'VariableNames', {'Country','GDP'}); 
writetable(resultTable,'results_ek.csv'); 

%----------------------------------------------------------------------
%DEFINE FUNCTION FOR SOLUTION

function r = Solve_EK_simple(s,d)

%number of countries
DC=d.DC;

%shock
deltahat=s.deltahat;

%parameters and shares
Niter=d.Niter;
err_tol=d.err_tol;
theta=d.theta;

Xbar=d.X;
pibar=d.pibar;

%speed variables
rho1=d.rho;

%intialize
Xhat=(ones(DC,DC)); %trade matrix --> row: n(importer) x column m(exporter)
wagehat=(ones(DC,1));


iter=0;
err_big=true;
err=100;

while (iter <  Niter && err_big )

    iter=iter+1;

    %SOLUTION STEPS BEGIN

    %Step 1: price including iceberg trade costs (Eq.1)
    costhat_nm=permute(wagehat,[2,1]);
    pnmhat= deltahat.*costhat_nm; 

    %Step 2: price index (Eq.2)
    Phihat = ( sum( pibar.*(pnmhat.^(-theta)) , 2 ) ).^(-(1/theta)); 

    %Step 3: trade flows (Eq.3)
    pihat = ( pnmhat./Phihat ).^(-theta); 

    %Step 4: market clearing and demand (Eq.4)
    incomehat = repmat(wagehat,1,DC); 
    Xhat = incomehat.*pihat; 

    %Step 5: new wages (Eq.5)
    wagehat_star=reshape(sum(Xbar.*Xhat,1)./sum(Xbar,1),[DC,1]);

    %Step 6: scale - need to scale guessed nominal variable
    scale=sum(Xhat(:).*Xbar(:))/sum(Xbar(:));
    wagehat_star=wagehat_star./scale;

    %Step 7: error computation then update variables
    err= sum( abs(wagehat_star(:) - wagehat(:)) );
    wagehat = rho1*wagehat_star + (1-rho1)*wagehat;

    %Step 8: check convergence
    if mod(iter,20)==0        
        err_big=(err > err_tol );
        disp(err)
    end

    %SOLUTION STEPS END

end %while

% GENERATE OUTPUTS
rGDPhat=wagehat./Phihat;

%percent changes  
cGDP = 100*(rGDPhat- 1.0);

r.Xhat=Xhat; %trade flows
r.wagehat=wagehat; %wages or production cost
r.pihat=pihat;  %import share
r.Phihat=Phihat; %price index
r.CPIhat=Phihat; %cpi same as above here
r.err=err; %final error
r.iter=iter; %iteration count
r.scale=scale; %normalization otherwise price can go to infinity or zero
r.rGDPhat=rGDPhat; %real gdp in hat
r.cGDP=cGDP; %in percent

end

%END FUNCTION FOR SOLUTION
%----------------------------------------------------------------------

 