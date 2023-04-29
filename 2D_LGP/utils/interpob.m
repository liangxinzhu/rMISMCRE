% ----------------------------------------------------- %
% function [fun] = interpob(gridps,Nx,Ny,xpos,ypos)
% linear interpolant the points
% input:  gridps = points at grids 
%         Nx     = 1/step length in x
%         Ny     = 1/step length in y
%         xpos   = x coordinates of interpolant points
%         ypos   = y coordinates of interpolant points
% output: fun    = interpolant values
% ----------------------------------------------------- %
function [fun] = interpob(gridps,Nx,Ny,xpos,ypos)
[k1,k2]=meshgrid([0:1/Nx:1],[0:1/Ny:1]);
fun = interp2(k1,k2,...
    [[gridps,gridps(:,1)];gridps(1,:),gridps(1,1)],xpos,ypos);
end

