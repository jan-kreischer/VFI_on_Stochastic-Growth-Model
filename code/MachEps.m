function eps = MachEps()

%  Purpose: Computes the machine epsilon

%  Output:  the smallest number e so that (1+e)>1 is true

%------------------------------------------------------------------- 

    eps=1;
    while (1+eps)>1
        eps=eps/2;
    end
    eps=2*eps;
end

