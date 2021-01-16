function [zt,p] = MarkovAR(size,nz,rho,sigma)
%   Purpose:    Markov function discribes the productivity shock
%
%   Input:      size := size of the grid for the productivity shock
%               nz := number of grid points for the productivity shock 
%               rho := autoregressive parameter
%               sigma := standard deviation of the innovations of the productivity shock
%
%   Output:     z:  m x 1 vector the grid approximating the process
%               p:  m x m matrix of transition probabilities

%STEP1 
% ﻿Compute the discrete approximation of the realizations
sigmaz=sqrt(sigma^2/(1-rho^2));
zbar=size*sigmaz;
zt=(-zbar:2*zbar/(nz-1):zbar);

p=zeros(nz);
%STEP2
% ﻿Compute the transition matrix
for i=1:nz
    p(i,1) = normcdf(((zt(1)-rho*zt(i))/sigma)+(zt(2)-zt(1))/(2*sigma));
    for j=2:nz-1
       p(i,j) = normcdf(((zt(j)-rho*zt(i))/sigma)+(zt(j)-zt(j-1))/(2*sigma)) - normcdf(((zt(j)-rho*zt(i))/sigma) - (zt(j)-zt(j-1))/(2*sigma));
    end
    p(i,nz)=1-sum(p(i,1:(nz-1)));
end
end

