% -------------------------------------------- %
% function [u_vector] = FEM_u(l,x)
%
% finite-element u_vector of Au = f
%
% inputs: l = level of refinement
%         x = the value of parameter
%
% outputs: u_vector = FEM grid solution
% -------------------------------------------- %
function [u_vector] = FEM_u(l,x)    
    
    h = 2^-l;
    
    K = 1/h;
     
    invh = h^(-1);
    
    A_i_j   = 2*ones(K-1,1); % diagonal
    A_i_jm1 = -ones(K-1,1); % lower
    A_i_jp1 = -ones(K-1,1); % upper
    
    A = spdiags([A_i_jm1 A_i_j A_i_jp1],[-1 0 1],K-1,K-1);
    
    A = A.*invh;
    
    f = x*h*ones(K-1,1);
    
    u_vector = [0; A\f; 0];

end
       
    
    