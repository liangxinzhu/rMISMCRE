% ------------------------------------------------------ %
% function [v] = pde_solver(z,l,x)
%
% inputs: z = interpolation points
%         l = level of refinement
%         x = the value of parameter
%
% outputs: v = solution at corresponding points
% ------------------------------------------------------ %
function [v] = pde_solver(z,l,x)
if l == 0
    v = 0;
else
    h = 2^-l;
    u_vector = FEM_u(l,x);
    z_vector = (0:h:1);
    v = zeros(length(z),1);
    for i = 1:length(z)
        if z(i) == 1 || z(i) == 0
            v(i) = 0;
        else
            ind = find(z(i) <= z_vector, 1);
            z_u = z_vector(ind);
            z_d = z_vector(ind - 1);
            v(i) = ( z_u - z(i) )/h*u_vector(ind-1) + ...
                ( z(i) - z_d )/h*u_vector(ind);
        end
    end
end
end