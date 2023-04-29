% ------------------------------------------------------------------- %
% function [sums, cost] = opre_sl(l,N,sigma,data,z_data)
%
% single level routine
%
% inputs: l      = level of refinement
%         N      = number of samples
%         sigma  = std deviation of the error of y - G
%         data   = observations (values of y) vector
%         z_data = corresponding observation points
%
% outputs: sums = estimator at level l
%          cost = computational cost for one sample
% ------------------------------------------------------------------- %
function [sums, cost] = opre_sl(l,N,sigma,data,z_data)
y = data;
% set up the temporing step
lambda_diff = 0.5;
Lambda = (0:lambda_diff:1);

likelihood = zeros(1,N);

for j = 1:length(Lambda)
    lambda = Lambda(j); 
    if j == 1
        x = 2*rand(N,1)-1;
        H = ones(1,N);
    else
        % resampling step
        A = Multinomial_Resampling(W);
        x = x(A');
        % mutation step
        x_rate = 0;
        for k = 1:N
            rate = 0;
            for i = 1:16
                x_star = x(k) + 0.5*randn(1);
                if x_star <= 1 && x_star >= -1
                    g_f = pde_solver(z_data, l, x_star);
                    l_temp = exp(-0.5/sigma/sigma*(g_f - y)'*(g_f - y));
                    alpha_temp = (l_temp/likelihood(k))^lambda;
                    alpha = min( 1, alpha_temp );
                    uni = rand(1);
                    if uni < alpha
                        x(k)   = x_star;
                        likelihood(k) = l_temp;
                        H(k) = log( l_temp^lambda_diff );
                        rate = rate + 1;
                    end
                end
            end
            rate = rate/16;
            x_rate = x_rate + rate/N;
        end
        fprintf(1,'tempering step %.1f, acceptance rate %.5f \n',lambda,x_rate)
    end
    W = exp( H - min(H) )./sum(exp( H - min(H)) );
end
sums = sum(x.^2)/N;
cost = 2^l;
end
