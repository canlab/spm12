function [ w , ierr ] = getsphwghts( g, g1, g2, g3 )
%
% compute weights for linear interpolation in g using g1,g2,g3
% pure matlab code
%
   ierr = 0;
   w = [0 0 0];
   w0 = sphtrarea1(g1,g2,g3);% in getspwghts.m
   w(1) = sphtrarea1(g,g2,g3);% in getspwghts.m
   w(2) = sphtrarea1(g,g1,g3);% in getspwghts.m
   w(3) = sphtrarea1(g,g1,g2);% in getspwghts.m
   if w0 < (sum(w) - 0.000001)
      ierr = 1;
   end   
   w = w/sum(w);
end

function area = sphtrarea1(g1,g2,g3)
%  Compute area of sherical triangle spanned by vectors 
%  g1,g2,g3 on unit sphere
%  use absolute values to identify opposite directions with each other
% pure matlab code
   c12 = abs(g1'*g2);% check ...
   c13 = abs(g1'*g3);
   c23 = abs(g2'*g3);
   s12 = sqrt(1-c12^2);
   s13 = sqrt(1-c13^2);
   s23 = sqrt(1-c23^2);
   a12 = acos(c12);
   a23 = acos(c23);
   a13 = acos(c13);
   b1 = min(1,max(-1,1+(c23-cos(a12-a13))/s12/s13));
   b2 = min(1,max(-1,1+(c13-cos(a12-a23))/s12/s23));
   b3 = min(1,max(-1,1+(c12-cos(a13-a23))/s13/s23));
   area = acos(b1)+acos(b2)+acos(b3)-pi;
end


