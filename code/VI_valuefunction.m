function vi = VI_valuefunction(x1)
%   Purpose:    The value function as a function of x1 alone. z0 and x0 are
%               given in global variables. This is used to find x' via golden section
%               search. 
%   Usage:      vi=VI_valuefunction(x1)

    global VI_zex VI_xex;
    vi=rhs_bellman(VI_zex,VI_xex,x1);
end

