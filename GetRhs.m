function f = GetRhs(x,z0,k1)
%   Purpose: rhs provides the compensation for an additional unit 
%            of savings: the increase in future output. 
%
%   Input:   z0 := scalar integer, the index of the current shock in the vector _VI_zvec;
%            x := scalar, the value of the endogenous state variable
%            k1 := scalar, the next-period vaulue of the capital stock
%   Output:  a pointer function that passes one number as an output each time, 
%            which is a part of the input for the gauss hermit integration.

%   Rhs of Euler equation
    global alpha delta eta z0 k1
    z1=z0+x;
    z1=exp(z1);

    k2=PF(k1,z1);
    c2=z1*(k1^alpha)+(1-delta)*k1-k2;
    
    f=(c2^(-eta))*(1.-delta+alpha*z1*(k1^(alpha-1)));
end

