function [w,x] = neumann(in)

% NEUMANN returns the weights and nodes for Neumann's numerical integration
%
% HOW w     = neumann(x)   -- 1st Neumann method 
%     [w,x] = neumann(n)   -- 2nd Neumann method (Gauss quadrature)
%
% IN  x - base points (nodes) in the interval [-1;1]
%     n - number of weights and number of base points
% OUT w - quadrature weights
%     x - base points (nodes) in the interval [-1;1]
%
% NB1 In 1st N-method, length(x) should not become too high, 
%     since a linear system of this size is solved. Eg: 500.
% NB2 No use is made of potential symmetries of nodes.
%
% Nico Sneeuw                    Munich                      23/08/95 

% uses GRULE, PLM
% revision history:
%   25/11/98: - interchange names 1st and 2nd
%   30/06/99: - NB1 and NB2 added
%   14/02/00: - remove dependencies on isinteger, isvector, isscalar
%
% 1st N.-method: see Sneeuw (1994) GJI 118, pp 707-716, eq. 19.5
% 2nd N.-method: see GRULE

if length(in) == 1				% 2nd Neumann method
    
   if rem(in,1) ~= 0,      error('integer input argument required.'),     end
   if exist('grule') ~= 2, error('GRULE from the NIT-toolbox required.'), end
   [x,w] = grule(in);
   x = x(:); w = w(:);			% put vectors upright

elseif min(size(in)) == 1			% 1st Neumann method
   
   x  = in(:);
   th = acos(x)*180/pi;			% [deg]
   l  = 0:(length(x)-1);
   pp = plm(l,th)';				% normalized Legendre polynomials
   r  = [2;zeros(length(x)-1,1)];		% right-handside vector
   w  = pp \ r;				% solve system of equations
   if (size(x) ~= size(w)), w = w'; end
   
else
   
   error('input argument should be scalar or vector')

end
