function rhs = rhs_bellman(z0,x0,x1)
%  Purpose:     returns the rhs of the bellman equation rf(z0,x0,x1)+beta*pmat[z0,.]'v1[x1,.]'
%
%
%  Input:       z0 := scalar integer, the index of the current shock in the vector _VI_zvec;
%               x0 := scalar, the current value of the endogenous state variable
%               x1 := scalar, the next-period vaulue of the endogenous state variable
%
%  Output:      vij := ths value of the rhs of the Bellman equation at (z0,x0,x1)

    global VI_IP VI_beta VI_pmat VI_ymat VI_zvec VI_xvec
    m=size(VI_pmat,2);
    rhs=rf(VI_zvec(z0),x0,x1);
    if VI_IP==1
        for j=1:m
            rhs=rhs+VI_beta*VI_pmat(z0,j)*LIP(VI_xvec,VI_ymat(:,j),x1);
        end
    end
end

