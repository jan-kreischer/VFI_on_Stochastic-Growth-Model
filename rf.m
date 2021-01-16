function c = rf(z,k1,k2)
%   Purpose:    rf defines the utility function
%
%   Input:      z  := scalar integer, the index of the shock;
%               k2 := scalar, the next-period value of the capital stock
%               k1 := scalar, the value of the capital stock
%
%   Output:     utility function


global alpha delta eta
%check for non-negative investment constraint
if(k2<(1-delta)*k1)
  disp("Error");
  pause;
end
c=z*(k1^alpha) + (1-delta)*k1 - k2;
  if c<0.0
     c=nan;
  end
     if eta==1
       c=ln(c);
     else
       c=(c^(1-eta))/(1-eta);
     end
end


