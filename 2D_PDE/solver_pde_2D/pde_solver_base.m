% --------------------------------------------------------------- %
% function [A,f] = pde_solver_base(xy,x,lx,ly)
%
% input: xy = grid points coordinate (x,y)
%        x  = value of parameter
%        lx = level of refinement in x direction
%        ly = level of refinement in y direction
%
% output: A = stiffness matrix without boundary condition
%         f = rhs without boundary condition
% --------------------------------------------------------------- %
function [A,f] = pde_solver_base(xy,x,lx,ly)
np = length(xy); % number of points
% initialise global matrix
A = sparse(np,np);
% compute stiffness matrix for each element
Ae = element_matrix(xy,x,lx,ly);
% corresponding the element index to global index
[elab] = element_label(lx,ly);
% perform assembly of global matrix
for row = 1:4
    nrow = sort(elab(:,row));
    for col = 1:4
        ncol = sort(elab(:,col));
        % global stiffness matrix
        A = A + sparse(nrow,ncol,Ae(:,row,col),np,np);
    end
end
% rhs
hlx = 2^-lx; hly = 2^-ly;
% non-liear PDE with rhs = 100
f = 100*hlx*hly*ones(np,1);
end

