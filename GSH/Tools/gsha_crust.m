function cs = gsha(f,method,grid,lmax)

% GSHA global spherical harmonic analysis
%
% HOW cs = gsha(f,method,grid)		
%
% IN  f      - global field of size (N+1)*2N or N*2N
%     method - string argument, defining the analysis method:
%              'ls'   - least squares
%              'wls'  - weighted least squares 
%              'aq'   - approximate quadrature 
%              'fnm'  - first neumann method
%              'snm'  - second neumann method
%              'mean' - block mean values (use of integrated Plm)
%     grid   - string argument, defining the grid:
%		       'pole' or 'mesh'     - equi-angular (N+1)*2N, including 
%		                              poles and Greenwich meridian.
%		       'block' or 'cell'    - equi-angular block midpoints N*2N
%		       'neumann' or 'gauss' - Gauss-Neumann grid (N+1)*2N
%     lmax   - maximum degree of development
% OUT cs     - Clm & Slm in |C\S| format
%
% NB  GRID argument is optional. Default: 'pole' if N+1, and 'block' if N.
%
% See also GSHS


%--------------------------------------------------------------------------
% Dimitris Tsoulis, IAPG, TU-Munich                                   11/98
%--------------------------------------------------------------------------
% uses PLM, IPLM, NEUMANN, SC2CS
%--------------------------------------------------------------------------
% revision history
%  - NS0299 - brush up (help text, layout, removal of unused commands, ...)
%           - restructuring of 'mean' method over 1st and 2nd analysis
%           - 'mean' method as quadrature instead of LS in 2nd step
%--------------------------------------------------------------------------
% remarks
%   TBD: - Zlm-functions option
%        - isotf, GRS80
%        - When 'pole' grid, m=1 yields singular Plm-matrix!
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% diagnostics and preliminaries
%--------------------------------------------------------------------------
error(nargchk(2,4,nargin))


%--------------------------------------------------------------------------
% Grid definition
%--------------------------------------------------------------------------
[rows,cols] = size(f);

if cols == 2*rows	 			% 'block' | 'cell' grid 

   if nargin < 4, lmax = rows; end      % default
   if nargin < 3, grid = 'block'; end	% default
   if ~strcmp(grid,'block') & ~strcmp(grid,'cell')
      error('Your GRID variable should be either block or cell')
   end
   
   n     = rows;
   dt    = 180 / n;
   theta = (dt/2:dt:180)';
   lam   = (dt/2:dt:360);           	% dt = dlam

elseif cols == 2*rows-2			% 'pole' | 'mesh' | 'neumann'
   
   if nargin < 4, lmax = rows-1; end    % default
   if nargin < 3, grid = 'pole'; end	% default
  
   n     = rows - 1;
   dt    = 180 / n;
   if strcmp(grid,'pole') | strcmp(grid,'mesh')
      theta = (0:dt:180)';
      lam   = (0:dt:360-dt);
   elseif strcmp(grid,'neumann') | strcmp(grid,'gauss')
      [gw,gx] = neumann(n+1);
      theta = acos(flipud(gx))*180/pi;
      lam   = (0:dt:360-dt);
   else
      error('The wrong type of GRID')
   end
   
else
   
   error('Invalid size of matrix F')
   
end	

%--------------------------------------------------------------------------
% further diagnostics
%--------------------------------------------------------------------------
if ~isstr(grid) | ~isstr(method)
   error('GRID and METHOD must be strings')
end
if strcmp(method,'snm') & ~(strcmp(grid,'neumann') | strcmp(grid,'gauss'))
   error('2nd Neumann method ONLY on a ''neumann''/''gauss'' GRID')
end
if strcmp(method,'mean') & ~(strcmp(grid,'block') | strcmp(grid,'cell'))
   error('Block mean method ONLY on a ''block''/''cell'' GRID')
end
if lmax > n
    error('Maximum degree of development is higher than number of rows of input.')
end

if strcmp(grid,'block') & lmax == n
    error('max. degree 180 is not possible for block grid, problem for m = 0')
end
if strcmp(grid,'block') & lmax == n
    error('max. degree 180 is not possible for polar grid, problem for m = 1')
end
%--------------------------------------------------------------------------
% Init.
%--------------------------------------------------------------------------
% L   = n;
L = lmax;
a   = zeros(rows,L+1);
b   = zeros(rows,L+1);
clm = zeros(L+1,L+1);
slm = zeros(L+1,L+1);

%--------------------------------------------------------------------------
% 1st step analysis: Am(theta) & Bm(theta)
%--------------------------------------------------------------------------
m   = 0:L;
c   = cos(lam'*m*pi/180);
s   = sin(lam'*m*pi/180);

% preserving the orthogonality (except for 'mean' case)
% we distinguish between 'block' and 'pole' type grids (in lambda)

if strcmp(grid,'block') | strcmp(grid,'cell')
   
   if strcmp(method,'mean')
      dl = dt;				% longitude block size
      c(:,1) = dl/360 * ones(2*n,1);	% ICm for m=0, including 1/(1+dm0)/pi
      m  = 1:L;
      ms = 2 ./ m .* sin(m*dl/2*pi/180) / pi;   
      c(:,2:L+1) = c(:,2:L+1) .* ms(ones(2*n,1),:);	% ICm
      s(:,2:L+1) = s(:,2:L+1) .* ms(ones(2*n,1),:);	% ISm
   else
      c = c/rows; s = s/rows;
      c(:,1)   = c(:,1)/2;			% / (1 + dm0 - dmL)
      s(:,L+1) = s(:,L+1)/2;		% / (1 - dm0 + dmL)
      c(:,L+1) = zeros(2*n,1); 		% CLL unestimable
      s(:,1)   = zeros(2*n,1);   		% Sl0    -"-
   end
   
else                                 	% 'pole' | 'mesh' | 'neumann'
   
   c = c/rows; s = s/rows;
   c(:,[1 L+1]) = c(:,[1 L+1])/2;		% / (1 + dm0 + dmL)
   s(:,[1 L+1]) = zeros(2*n,2);		% / (1 - dm0 - dmL), Sl0&SLL unestimable                 
   
end

a = f * c;
b = f * s;

%--------------------------------------------------------------------------
% 2nd step analysis: Clm & Slm
%--------------------------------------------------------------------------
% hwb = waitbar(0,'Percentage of orders m ready ...');
% set(hwb,'NumberTitle','off','Name','GSHA')

if strcmp(method,'ls')			% Least squares solution

   for m = 0:L
       p  = plm(m:L,m,theta);
%        if m == 0
%            ss = size(p)
%        end
       ai = a(:,m+1);
       bi = b(:,m+1);
       clm(m+1:L+1,m+1) = p \ ai;
       slm(m+1:L+1,m+1) = p \ bi;  
%        waitbar((m+1)/(L+1))			% Update the waitbar
   end

elseif strcmp(method,'wls')		% Weighted Least Squares
   
   si = sin(theta*pi/180);
   si = 2*si/sum(si);
   for m = 0:L
      p  = plm(m:L,m,theta);
      ai = a(:,m+1);
      bi = b(:,m+1);
      d   = 1:length(theta);
      pts = p' * sparse(d,d,si);
      clm(m+1:L+1,m+1) = (pts * p) \ pts * ai;
      slm(m+1:L+1,m+1) = (pts * p) \ pts * bi;  
%       waitbar((m+1)/(L+1))			% Update the waitbar
   end
   
elseif strcmp(method,'aq')			% Approximate Quadrature

   si = sin(theta*pi/180);
   si = 2*si/sum(si);
   for m = 0:L
       p  = plm(m:L,m,theta);
       ai = a(:,m+1);
       bi = b(:,m+1);
       clm(m+1:L+1,m+1) = (1 + (m==0))/4 * p' * (si.*ai);
       slm(m+1:L+1,m+1) = (1 + (m==0))/4 * p' * (si.*bi);   
%        waitbar((m+1)/(L+1))		% Update the waitbar
    end
    
elseif strcmp(method,'fnm')		% 1st Neumann method (exact up to L/2)
    
   w = neumann(cos(theta/180*pi));
   for m = 0:L
      p  = plm(m:L,m,theta);
      ai = a(:,m+1);
      bi = b(:,m+1);
      clm(m+1:L+1,m+1) = (1 + (m==0))/4 * p' * (w.*ai);
      slm(m+1:L+1,m+1) = (1 + (m==0))/4 * p' * (w.*bi);
%       waitbar((m+1)/(L+1))			% Update the waitbar
   end
   
elseif strcmp(method,'snm')		% 2nd Neumann method (exact)
   
   % weigths determined above already
   for m = 0:L
      p  = plm(m:L,m,theta);
      ai = a(:,m+1);
      bi = b(:,m+1);
      clm(m+1:L+1,m+1) = (1 + (m==0))/4 * p' * (gw.*ai);
      slm(m+1:L+1,m+1) = (1 + (m==0))/4 * p' * (gw.*bi);
%       waitbar((m+1)/(L+1))			% Update the waitbar
   end
   
elseif strcmp(method,'mean')		% block mean values
   
   for m = 0:L
      p  = iplm(m:L,m,theta);		% integrated Legendre
      ai = a(:,m+1);
      bi = b(:,m+1);
      clm(m+1:L+1,m+1) = (1 + (m==0))/4 * p' * ai;
      slm(m+1:L+1,m+1) = (1 + (m==0))/4 * p' * bi;   
%       waitbar((m+1)/(L+1))			% Update the waitbar
   end
      
else
   
   error('Choose a valid METHOD')
   
end

% close(hwb)					% close waitbar

%--------------------------------------------------------------------------
% Write the coefficients Clm & Slm in |C\S| format
%--------------------------------------------------------------------------
slm = fliplr(slm);
cs  = sc2cs([slm(:,1:L) clm]);
cs  = cs(1:lmax+1,1:lmax+1);