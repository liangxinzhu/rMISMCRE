% ---------------------------------------------------------- %
% [x,y,xy,bound] = domain(lx, ly)
%
% grid generator on a domain [0,1]^2
% inputs: lx = level of refinement in x direction
%         ly = level of refinement in y direction 
%
% outputs: x     = grid points in x
%          y     = grid points in y
%          xy    = grid points (x,y)
%          bound = index of boundary points
% ---------------------------------------------------------- %

function [x,y,xy,bound] = domain(lx,ly)
nx = 2^lx; ny = 2^ly; % time-step in x and y directions
x = [0:1/nx:1]; y = [0:1/ny:1]; % grid points in x and y directions

nxy = (nx+1)*(ny+1); % total number of grid points
[X,Y] = meshgrid(x,y);
xx = reshape(X',nxy,1); yy = reshape(Y',nxy,1);
xy = [xx(:),yy(:)]; % grid points (x,y)

% four boundary edges
k1 = find( xy(:,2)==0 ); % bottom 
k2 = find( xy(:,1)==1 & xy(:,2)<=1 & xy(:,2) >0 ); % right
k3 = find( xy(:,2)==1 & xy(:,1)<1  & xy(:,1) >0 ); % top
k4 = find( xy(:,1)==0 & xy(:,2)<=1 & xy(:,2) >0 ); % left

bound = sort([k1;k2;k3;k4]); % index of boundary points
end

