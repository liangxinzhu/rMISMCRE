% ------------------------------------------------------------------- %
% function [sums, cost] = opre_rmlsmc_l(l,N,sigma,data,z_data)
%
% single level routine
%
% inputs: l      = level of refinement
%         N      = number of samples
%         sigma  = std deviation of the error of y - G
%         data   = observations (values of y) vector
%         z_data = corresponding observation points
%
% outputs: sums(1) = self-normalized increments estimator at level l
%          sums(2) = unnormalised increments \phi at level l
%          sums(3) = unnormalised increments 1 at level l
%          cost    = computational cost for one sample
% ------------------------------------------------------------------- %
function [sums, cost] = opre_rmlsmc_l(l,N,sigma,data,z_data)
y = data;
sums(1:3) = 0;

% set up the temporing step
lambda_diff = 0.5;
Lambda = (lambda_diff:lambda_diff:1);

% initialisation
x = 2*rand(N,1)-1;
Z = 1;

% multi-level increments
H   = zeros(N,1);
G_k = zeros(N,1);
f_f = zeros(N,1);
f_c = zeros(N,1);

parfor i = 1:N
    g_f = pde_solver(z_data, l, x(i));
    l_f_temp = exp(-0.5*(g_f - y)'*(g_f - y)/sigma/sigma);
    if l == 1
        g_c = 0;
        l_c_temp = 0;
    else
        g_c = pde_solver(z_data, l-1, x(i));
        l_c_temp = exp(-0.5*(g_c - y)'*(g_c - y)/sigma/sigma);
    end
    G_k(i) = max( l_f_temp, l_c_temp );
    f_f(i) = l_f_temp;
    f_c(i) = l_c_temp;
    H(i) = log( max( l_f_temp, l_c_temp )^lambda_diff );
end


for j = 1:length(Lambda)
    
    lambda = Lambda(j);
    % normalising constant
    Z = Z*sum(exp(H))/N;
    if j ~= 1
        W = exp( H - min(H) )./sum(exp( H - min(H)) );
        % resampling step
        A = Multinomial_Resampling(W);
        x = x(A');
    end
    x_rate = 0;
    for k = 1:N
        rate = 0;
        for i = 1:9
            x_star = x(k) + 0.5*randn(1);
            if x_star <= 1 && x_star >= -1 
                g_f = pde_solver(z_data, l, x_star);
                l_f_temp = exp(-0.5*(g_f - y)'*(g_f - y)/sigma/sigma);
                if l == 1
                    g_c =0;
                    l_c_temp = 0;
                else
                    g_c = pde_solver(z_data, l-1, x_star);
                    l_c_temp = exp(-0.5*(g_c - y)'*(g_c - y)/sigma/sigma);
                end
                G_star = max( l_f_temp, l_c_temp );
                alpha_temp = (G_star/G_k(k))^lambda;
                alpha = min( 1, alpha_temp );
                uni = rand(1);
                if uni < alpha
                    x(k) = x_star;
                    G_k(k) = G_star;
                    
                    f_f(k) = l_f_temp;
                    f_c(k) = l_c_temp;
                    H(k) = log( G_star^lambda_diff );
                    rate = rate + 1;
                end
            end
        end
        rate = rate/9;
        x_rate = x_rate + rate/N;
    end
    fprintf(1,'level %d, tempering step %.1f, acceptance rate %.5f \n',l,lambda,x_rate)
end

if l == 1
    sums(1) = sum(x.^2)/N;
    sums(2) = Z*sum(x.^2)/N;
    sums(3) = Z;
    cost = 2^l;
else
    sums(1) = sum(x.^2 .* f_f./G_k)/sum(f_f./G_k) - ...
              sum(x.^2 .* f_c./G_k)/sum(f_c./G_k);
    sums(2) = Z*( sum(x.^2 .* f_f./G_k)- sum(x.^2 .* f_c./G_k) )/N;
    sums(3) = Z*( sum(f_f./G_k) - sum(f_c./G_k) )/N;
    cost = 2*2^l;
end

end
