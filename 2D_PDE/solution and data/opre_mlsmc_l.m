% ------------------------------------------------------------------ %
% [sums, cost] = opre_mlsmc_l(lx,ly,N,sigma,data,x_data,K)   
%
% single level routine
%
% inputs:  lx     = level x
%          ly     = level y
%          N      = number of samples
%          sigma  = std deviation of the error of y - G
%          data   = observations
%          x_data = corresponding points of observations
%          K      = actual starting refinement level
%
% output:  sums(1) = normalised difference
%          sums(2) = difference unnormalised integral
%          sums(3) = difference normalising constant
%          cost    = computational cost
% ------------------------------------------------------------------ %
function [sums, cost] = opre_mlsmc_l(lx,ly,N,sigma,data,x_data,K)
y = data;
sums(1:3) = 0;

% set up the temporing step
lambda_diff = 0.5;
Lambda = (lambda_diff:lambda_diff:1);

% initialisation
u_single = 2*rand(N,2)-1; 
u = u_single;
Z_single = 1;
Z = Z_single;

% actual refinement level
rlx = lx + K;
rly = ly + K;

if rlx == K || rly == K
% single level approximation
H_single = zeros(1,N);
f_single = zeros(1,N);
parfor i = 1:N
    [g_single, ~] = pde_solver_2D(rlx,rly,x_data,u_single(i,:));
    l_f_temp = exp(-0.5/sigma/sigma*(g_single - y)'*(g_single - y));
    f_single(i) = l_f_temp;
    H_single(i) = log(l_f_temp^lambda_diff);
end

for j = 1:length(Lambda)
    
    lambda = Lambda(j);
    
    % normalising constant
    Z_single = Z_single*sum(exp(H_single))/N;
    
    W = exp(H_single - min(H_single))./ ...
        sum(exp(H_single - min(H_single)));

    % resampling step
    A = Multinomial_Resampling(W);
    u_single = u_single(A',:);
    
    % mutation step
    parfor K = 1:N
%          rate = 0;
        for i = 1:8
            u_star = u_single(K,:) + 1.5*randn(1,2);
            
            [g_single, ~] = pde_solver_2D(rlx,rly,x_data,u_star);
            
            l_f_temp = exp(-0.5/sigma/sigma*(g_single - y)'*(g_single - y));
            
            if u_star(1) >= 1 || u_star(1) <= -1 || u_star(2) >= 1 || u_star(2) <= -1
                alpha_temp = 0;
            else
                alpha_temp = (l_f_temp/f_single(K))^lambda;
            end
            
            alpha = min( 1, alpha_temp );
            
            uni = rand(1);
            
            if uni < alpha
                u_single(K,:) = u_star;
                f_single(K) = l_f_temp;
                H_single(K) = log(l_f_temp^lambda_diff);
%                  rate = rate + 1;
            end
        end
    end
end

L2norm2 = u_single(:,1).^2 + u_single(:,2).^2;
sums(1) = sum(L2norm2)/N;
sums(2) = Z_single*sum(L2norm2)/N;
sums(3) = Z_single;

else
% multi-level increments
H = ones(1,N);
G_k = ones(1,N);
f_1 = ones(1,N);
f_2 = ones(1,N);

parfor i = 1:N
    [g_1, ~] = pde_solver_2D(rlx,rly,x_data,u(i,:));
    [g_2, ~] = pde_solver_2D(rlx-1,rly-1,x_data,u(i,:));
    l_1_temp = exp(-0.5/sigma/sigma*(g_1 - y)'*(g_1 - y));
    l_2_temp = exp(-0.5/sigma/sigma*(g_2 - y)'*(g_2 - y));
    G_k(i) = max( [l_1_temp,l_2_temp] );
    f_1(i) = l_1_temp;
    f_2(i) = l_2_temp;
    H(i) = log( G_k(i)^lambda_diff );
end

for j = 1:length(Lambda)
    
    lambda = Lambda(j);
    
    % normalising constant
    Z = Z*sum(exp(H))/N;
    
    W = exp( H - min(H) )./sum(exp( H - min(H)) );

    % resampling step
    A = Multinomial_Resampling(W);
    u = u(A',:);
    
    parfor K = 1:N
%              rate = 0;
        for i = 1:8
            u_star = u(K,:) + 1.5*randn(1,2);
            
            [g_1, ~] = pde_solver_2D(rlx,rly,x_data,u_star);
            [g_2, ~] = pde_solver_2D(rlx-1,rly-1,x_data,u_star);
            l_1_temp = exp(-0.5/sigma/sigma*(g_1 - y)'*(g_1 - y));
            l_2_temp = exp(-0.5/sigma/sigma*(g_2 - y)'*(g_2 - y));
            G_star = max( [l_1_temp,l_2_temp] );
            
            if u_star(1) >= 1 || u_star(1) <= -1 || u_star(2) >= 1 || u_star(2) <= -1
                alpha_temp = 0;
            else
                alpha_temp = (G_star/G_k(K))^lambda;
            end
            
            alpha = min( 1, alpha_temp ); 
            uni = rand(1);
            
            if uni < alpha
                u(K,:) = u_star;
                G_k(K) = G_star;
                
                f_1(K) = l_1_temp;
                f_2(K) = l_2_temp;
                
                H(K) = log( G_star^lambda_diff );
%                 rate = rate + 1;
            end
        end
    end  
end

L2norm2 = u(:,1).^2 + u(:,2).^2;
      
sums(1) = sum(L2norm2' .* f_1./G_k)/sum(f_1./G_k) - ...
          sum(L2norm2' .* f_2./G_k)/sum(f_2./G_k);   
sums(2) = Z*( sum(L2norm2' .* f_1./G_k) - ...
              sum(L2norm2' .* f_2./G_k))/N;
sums(3) = Z*( sum(f_1./G_k) - ...
              sum(f_2./G_k))/N;

end
% cost defined as number of fine timesteps
cost = 2*2^(rlx+rly);
end

