function z = BLIP(xvec,yvec,zmat,x,y)
%BLIP Bilinear Interpolation as given by Formulas (3.6.1) through (3.6.5) in Press et al. (1992), S. 116f.
% Input:   xvec = n by 1, vector, the grid of variable x
%          yvec = m by 1, vector, the grid of variable y
%          zmat = n by m matrix, with tabulated function values at
%                      z(i,j)=f(x(i),y(j), i=1, ..., n, j=1, ...,m
%             x = scalar,  the x coordinate
%             y = scalar,  the y coordinate
% Output:  z, the interpolated value of f(x,y)

n=length(xvec);
m=length(yvec);
%
z=zeros(4,1);

% first, locate the square that surrounds (x,y)
if (x<xvec(1) || x>xvec(n))
    disp("x outside of grid! Program stops. Press any key");
    pause;
    z=nan;
    return
end
if (y<yvec(1) || y>yvec(m))
    disp("y outside of grid! Program stops Press any key");
    pause;
    z=nan;
    return
end
% second, the optimal next-period asset is given by:
i=sum(xvec<=x);
j=sum(yvec<=y);

    if (i==n) && (j==m)
        z=zmat(n,m);
    elseif (i==n) && (j<m)
        u=(y-yvec(j))/(yvec(j+1)-yvec(j));
        z=(1-u)*zmat(n,j)+u*zmat(n,j+1);
    elseif (i<n) && (j==m)
        t=(x-xvec(i))/(xvec(i+1)-xvec(i));
        z=t*zmat(i+1,m)+(1.-t)*zmat(i,m);
    else
        t=(x-xvec(i))/(xvec(i+1)-xvec(i));
        u=(y-yvec(j))/(yvec(j+1)-yvec(j));
        z=(1.-t)*(1.-u)*zmat(i,j)+t*(1.-u)*zmat(i+1,j)+t*u*zmat(i+1,j+1)+(1.-t)*u*zmat(i,j+1);         
    end
end



