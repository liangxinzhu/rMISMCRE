% ------------------------------------------------------- %
% [A, f] = boundary_conditions(A, f, xy, bound)
%
% inputs: A     = stiffness matrix without boundary
%         f     = rhs vector without boundary
%         xy    = grid coordiate vector
%         bound = boundary index vector
%
% output: A = stiffness matrix with boundary
%         f = rhs vector with boundary            
% ------------------------------------------------------- %
function [A, f] = boundary_conditions(A, f, xy, bound)
%% pre-settings
np  = length(f);     % number of solution points
nbd = length(bound); % number of points at the boundary

null_col = sparse(np,nbd);
null_row = sparse(nbd,np);
%% set boundary condition
% x-y coordinate of points at the boundary
xbd = xy(bound,1); ybd = xy(bound,2); 
% boundary conditions at different points
bc = zeros(size(xbd));
% reset stiffness matrix with boundary condition
dA = zeros(np,1); dA(bound) = ones(nbd,1);
A(:,bound) = null_col; A(bound,:) = null_row;
A = A + spdiags(dA,0,np,np);
f(bound) = bc;
end

