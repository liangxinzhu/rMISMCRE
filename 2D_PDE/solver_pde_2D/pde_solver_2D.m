% ------------------------------------------------------- %
% function [v, v_grid] = pde_solver_2D(lx,ly,XY,x)
%
% inputs: lx = level of refinement in x direction
%         ly = level of refinement in y direction
%         XY = interpolant points (x,y) vector
%         x  = parameter vector (R^2£©
%
% outputs: v      = solution at interpolant points
%          v_grid = solution at grid points
% ------------------------------------------------------- %
function [v, v_grid] = pde_solver_2D(lx,ly,XY,x)
if lx == 0 || ly == 0
    v = zeros(length(XY),1);
    nx = 2^lx; ny = 2^ly; % steps in x and y directions
    nxy = (nx+1)*(ny+1);  % total number of grid points
    v_grid = zeros(nxy,1);
else
% define geometry
[~,~,xy,bound] = domain(lx,ly);
% compute stiffness matrix and rhs without boundary
[A,f] = pde_solver_base(xy,x,lx,ly);
% compute stiffness matrix and rhs with boundary
[A,f] = boundary_conditions(A,f,xy,bound);
% compute solution at grid points
v_grid = A\f;
% interpolate solutions
[X,Y] = meshgrid(XY(:,1),XY(:,2));
v = griddata(xy(:,1),xy(:,2),v_grid,X,Y);
v = diag(v);
end
end

