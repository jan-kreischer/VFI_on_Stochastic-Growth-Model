function knext = PF(k,z)
%   Purpose:    Policy function via linear interpolation
%
%   Input:      kgrid:= n by 1, vector, the grid of variable k
%               zgrind:= m by 1, vector, the grid of variable z
%               hmat:=  n by m matrix, with tabulated function values
%
%   Output:     knext - the interpolated value of f(k,z)

    global kgrid zgrid hmat
    knext=BLIP(kgrid,zgrid,hmat,k,z);
    
end

