% ------------------------------------------------------ %
% function [x2] = analytic_x2(sigma,data,z_data)
%
% inputs: sigma  = std deviation of the error of y - G
%         date   = observations (value of y) vector
%         z_data = corresponding observation points
%
% outputs: x2 = analytic solution of x^2
% ------------------------------------------------------ %
function [x2] = analytic_x2(sigma,data,z_data)
x = z_data;
temp1 = sum(data.^2);
temp2 = sum(data.*(x.^2-x)/2);
temp3 = sum((x.^2-x).^2/4);
fun1 = @(z) exp(-(temp1 + 2*temp2.*z + temp3.*z.^2)/2/sigma/sigma);
cont = integral(fun1,-1, 1,'AbsTol',1e-16);
fun2 = @(z) z.^2.*exp(-(temp1 + 2*temp2.*z + temp3.*z.^2)/2/sigma/sigma)/cont;
x2 = integral(fun2,-1, 1,'AbsTol',1e-16);
end