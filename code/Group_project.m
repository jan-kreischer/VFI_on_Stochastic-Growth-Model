clear all
format long

% structural parameters of the economy

global alpha beta delta eta sigma rho zgrid kgrid hmat

alpha	= 0.27;                 %elasticity of production with respect to capital 
beta	= 0.994;                %discount factor
sigma	= 0.05;                 %standard deviation of the innovations of the productivity shock
delta	= 0.011;                %depretiation rat eof capital 
eta     = 2;                    %elasticity of marginal utility
rho     = 0.9;                  %autoregressive parameter

nz=9;                           %number of grid points for the productivity shock 
grid_size=5.5;                  %size of the grid for the productivity shock 
kmin_g=0.60;                    %lower bound of the grid for the capital stock
kmax_g=1.60;                    %upper bound of the grid for the capital stock
cumulative=1;                   %in successive computations, compute total run time, =0 do not

%Parameters to for the program

global VI_IP VI_xvec VI_ymat VI_zvec VI_pmat VI_beta VI_xex VI_zex VI_eps VI_nc VI_Max VI_BS VI_clear
VI_IP=1;                        % =1 linear interpolation, =2, cubic interpolation, =0 without interpolation
VI_xvec=0;                      % stores the x values used in interpolation
VI_ymat=0;                      % stores the y values uses in interpolation
VI_zvec=0;                      % stores the z values used in interpolation
VI_pmat=0;                      % stores the transition matrix used in interpolation
VI_beta=0;                      % stores information uses for linear interpolation
VI_xex=0;                       % stores information to write the value function as a function of x[i] alone
VI_zex=0;                       % stores information to write the value function as a function of x[i] alone

VI_eps=0.01;                    % stopping criterium
VI_nc=50;                       % stop if number of consecutive iterations with unchanged indices of policy function exceeds this number
VI_Max=5000;                    % maximal number of iterations
VI_BS=1;                        % binary search instead of sequential search
VI_clear=1;                     % clear screen for printing

%%Parameters for the computation of Euler equation residuals

kmin_e=0.8;                     % kmin_e*kstar is the lower bound
kmax_e=1.2;                     % kmax_e*kstar is the upper bound
nobs_e=200;                     % the number of residuals to be computed

%% Step 1: Markov chain 
%﻿The approximation by a Markov chain of the continuous valued AR(1).
% By using Markov chain, we force the realizations of zt to be taken from a finite grid of points.
% For the markov chain a "MarkovAR" fucntion was written
[zgrid,pmat]=MarkovAR(grid_size,nz,rho,sigma);

zgrid=exp(zgrid);

zmin=zgrid(1);
zmax=zgrid(nz);

kmin=((1-beta*(1-delta))/(alpha*beta*zmin))^(1/(alpha-1));
kmax=((1-beta*(1-delta))/(alpha*beta*zmax))^(1/(alpha-1));

%% Step 2: Compute stationary solution of deterministic model and intialize the value function

kstar	= ((1-beta*(1-delta))/(alpha*beta))^(1/(alpha-1)); %stationary capital stock
cstar   = kstar^alpha - delta*kstar;                       %stationary level of consumption


kmin_g	= kmin_g*kmin;
kmax_g	= kmax_g*kmax;

kmin_e=kmin_e*kstar;
kmax_e=kmax_e*kstar;

%initialize v0 (inital value function)
%Vector with different values of nk
nvec=[5];

nk=nvec(1);
v0=rf(1,kstar,kstar)/(1-beta);  % stationary solution
v0=ones(nk,nz).*v0;

%% Iterations of nvec start here
lmax=length(nvec);

policy=zeros(nobs_e,nobs_e,lmax);
emat=zeros(nobs_e,nobs_e,lmax);

%Iterations over different nk start here

s2=0;        % stores time needed to obtain initial v from previous v

tottime=zeros(lmax,1);

%% Step 3: Main Loop
% Compute the policy function for a stochastic DGE model with
% one endogenous state variable x and one exogenous shock z.
for l=1:lmax
    nk=nvec(l);
    kgrid=(kmin_g:(kmax_g-kmin_g)/(nk-1):kmax_g);

    % Solve for the policy function 
    s1=seconds(timeofday(datetime('now')));
    [v1,hmati]=SolveVIS(beta,kgrid,zgrid,pmat,v0);
    s1=seconds(timeofday(datetime('now')))-s1;
    
    if VI_IP==1
     hmat=hmati;
    else
        hmat=zeros(nk,nz);
        for i=1:nk
            for j=1:nz
                hmat(i,j)=kgrid(hmati(i,j));
            end
        end
    end

    % Computation of Euler equation residuals
    
    % For the Euler equation a "Euler" fucntion was written
    kvec=(kmin_e:(kmax_e-kmin_e)/(nobs_e-1):kmax_e);
    zvec=(0.95:0.1/(nobs_e-1):1.05);
    clear z0 k1;
    eer=Euler(kvec,zvec);
    emat(1:nobs_e,1:nobs_e)=eer;
    emax=max(max(abs(eer)));
    tottime(l)=s1+s2;
    
    file = fopen("Stochastic Ramsey Model.txt","a+");
    fprintf(file,"EER = %e\n",emax);
    fprintf(file,"Run time = %d minutes and %.2f seconds\n",floor((s1+s2)/60),rem((s1+s2),60));
    fclose(file);
    
    % computation of policy function 
    for i=1:nobs_e
        for j=1:nobs_e
            policy(i,j,l)=BLIP(kgrid,zgrid,hmat,kvec(i),zvec(j));
        end
    end
    % New initial v0    
    if l<lmax
        if nk==nvec(l+1)
            v0=v1;
        else
            s2=seconds(timeofday(datetime('now')));
            nk1=nvec(l+1);
            kgnew=kmin_g:(kmax_g-kmin_g)/(nk1-1):kmax_g;   
            kgnew(nk1)=kgrid(nk);
            kgnew(1)=kgrid(1); 
            v0=zeros(nk1,nz);    
            for j=1:nz
                v0(:,j)=LIP(kgrid,v1(:,j),kgnew);
            end
            s2=seconds(timeofday(datetime('now')))-s2;
        end
    end
    %optional plot for capital, tfp and cp
    cp=zeros(nk,nz);
    for i=1:nk
        for j=1:nz
            cp(i,j)=zgrid(j)*(kgrid(i)^alpha)+(1-delta)*kgrid(i)-hmat(i,j);
        end
    end
    surface(kgrid',zgrid,cp');
    %
end
save('emat.mat','emat');
save('policy.mat','policy');
save('kvec.mat','kvec');
save('zvec.mat','zvec');

%vector for tfp (tfp = Z_t, as described in the book)
tfp = (zmin:(zmax-zmin)/200:zmax);




