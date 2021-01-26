function [p, dp] = plm(l,m,th)

% PLM Fully normalized associated Legendre functions for a selected order M
%
% HOW p      = plm(l,th)			- assumes M=0
%     p      = plm(l,m,th)
%     [p,dp] = plm(l,m,th)
%
% IN  l  - degree (vector). Integer, but not necessarily monotonic.
%          For l < m a vector of zeros will be returned.
%     m  - order (scalar). If absent, m=0 is assumed.
%     th - co-latitude [deg] (vector)
% OUT p  - Matrix with Legendre functions. The matrix has length(TH) rows
%          and length(L) columns, unless L or TH is scalar. Then the output
%          vector follows the shape of respectively L or TH. 
%    dp  - Matrix with first derivative of Legendre functions. The matrix 
%          has length(TH) rows and length(L) columns, unless L or TH is 
%          scalar. Then the output vector follows the shape of respectively 
%          L or TH. 
% 
% See also LEGPOL, YLM, IPLM

%-----------------------------------------------------------------------------
% Nico Sneeuw, IAPG, TU-Munich                                       08/08/94
%-----------------------------------------------------------------------------
% Uses none
%-----------------------------------------------------------------------------
% Revision history:
%  - NS09/06/97: help text brushed up
%  - NS13/07/98: Pmm non-recursive anymore 
%  - NS0299:     further help text brush-up
%  - MW13/08/04: extension for first derivative
%  - MW24/11/04: speed up calculation
%-----------------------------------------------------------------------------


% Some input checking.
if nargin == 2
   th = m;
   m  = 0;
end
if min(size(l)) ~= 1,  error('Degree l must be vector (or scalar)'), end
if any(rem(l,1) ~= 0), error('Vector l contains non-integers.'), end
if max(size(m)) ~= 1,  error('Order m must be scalar.'), end
if rem(m,1) ~=0,       error('Order m must be integer.'), end


% Preliminaries.
[lrow,lcol] = size(l);
[trow,tcol] = size(th);
lmax = max(l);
if lmax < m, error('Largest degree still smaller than order m.'), end
n    = length(th);				% number of latitudes
t    = th(:)*pi/180;
x    = cos(t);
y    = sin(t);
lvec = l(:)';					% l can be used now as running index.

% Recursive computation of the temporary matrix ptmp, containing the Legendre
% functions in its columns, with progressing degree l. The last column of
% ptmp will contain zeros, which is useful for assignments when l < m.
ptmp  = zeros(n,lmax-m+2);
if nargout == 2, dptmp = zeros(n,lmax-m+2); end

%--------------------------------------------------------------------
% sectorial recursion: PM (non-recursive, though)
%--------------------------------------------------------------------
if m == 0
   fac = 1;
else
   mm  = 2*(1:m);
   fac = sqrt(2*prod((mm+1)./mm));
end

ptmp(:,1) = fac*y.^m;                                      % The 1st column of ptmp.
if nargout == 2, dptmp(:,1) = m*fac*(y.^(m-1).*x); end     % The 1st column of dptmp.

%--------------------------------------------------------------------
% l-recursion: P
%--------------------------------------------------------------------
for l = m+1:lmax
   col   = l - m + 1;			% points to the next column of ptmp
   root1 = sqrt( (2*l+1)*(2*l-1)/((l-m)*(l+m)) ) ;
   root2 = sqrt( (2*l+1)*(l+m-1)*(l-m-1) / ( (2*l-3)*(l-m)*(l+m) ) );

   % recursion
   if l == m+1
       ptmp(:,col) = root1 *x.*ptmp(:,col-1);
   else
       ptmp(:,col) = root1 *x.*ptmp(:,col-1) - root2 *ptmp(:,col-2);
   end
       
   if nargout == 2, 
       if l == m+1
           dptmp(:,col) = root1 *(x.*dptmp(:,col-1)-y.*ptmp(:,col-1)); 
       else
           dptmp(:,col) = root1 *(x.*dptmp(:,col-1)-y.*ptmp(:,col-1)) - root2 *dptmp(:,col-2); 
       end
   end
end


% The Legendre functions have been computed. What remains to be done, is to
% extract the proper columns from ptmp, corresponding to the vector lvec. 
% If l or theta is scalar the output matrix p reduces to a vector. It should
% have the shape of respectively theta or l in that case.

p          = zeros(n,length(lvec));         % size declaration.
lind       = find(lvec < m);			    % index into l < m
pcol       = lvec - m + 1;			        % index into columns of ptmp
pcol(lind) = (lmax-m+2)*ones(size(lind));	% Now l < m points to last col.
p          = ptmp(:,pcol);			        % proper column extraction 
if nargout == 2, dp = dptmp(:,pcol); end    % proper column extraction

if max(size(lvec))==1  & min(size(th))==1 & (trow == 1), 
    p = p'; 
    if nargout == 2, dp = dp'; end
end
if max(size(th))==1 & min(size(lvec))==1  & (lcol == 1), 
    p = p'; 
    if nargout == 2, dp = dp'; end
end