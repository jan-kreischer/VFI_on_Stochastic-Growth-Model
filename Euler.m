function [eer] = Euler(kvec,zvec)
%   Purpose:  The stationary euler equation allows to the average capital stock  
%
%   Input:    kvec := n by 1 column vector with elements in ascending order that
%             represent valid elements of the variable k.
%             zvec := m by 1 column vector that represent the elements of the Markov chain
%             that represents the stochastic variable z
%
%   Output:   eer(i,j):  Euler equation residuals



    global z0 k1 alpha beta delta eta sigma rho
    n=length(kvec);
    m=length(zvec);
    eer=zeros(n,m);

    for i=1:n
        for j=1:m
            z0=rho*log(zvec(j));
            k1=PF(kvec(i),zvec(j));
            c0=zvec(j)*(kvec(i)^alpha)+(1-delta)*kvec(i)-k1;            
            rhs=GH_INT4(@GetRhs,sigma);
            rhs=rhs*beta;
            c1=(rhs^(-1/eta));
            eer(i,j)=(c1/c0)-1;
        end
    end
end

