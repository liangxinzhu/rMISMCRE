% ------------------------------------------------------------------ %
% [sums, cost] = opre_rmismc_l(lx,ly,Nl,sigma,data,z_data,K)  
%
% single level routine
%
% inputs: lx     = level x
%         ly     = level y
%         Nl      = number of samples
%         sigma  = std deviation of the error of y - G
%         data   = observations
%         z_data = corresponding points of observations
%         K      = starting refinement level
%
% output: sums(1) = normalised difference
%         sums(2) = difference unnormalised integral
%         sums(3) = difference normalising constant
%         cost    = computational cost
% ------------------------------------------------------------------ %
function [sums, cost] = opre_rmismc_l(lx,ly,Nl,sigma,data,z_data,K)
y = data;
sums(1:3) = 0;

% set up the temporing step
lambda_diff = 0.5;
Lambda = (lambda_diff:lambda_diff:1);

% initialisation
x_single = 2*rand(Nl,2)-1; 
x = x_single;
Z_single = 1;
Z = Z_single;

% actual refinement level
rlx = lx;
rly = ly;

% multi-index increments
H = ones(1,Nl);
G_k = ones(1,Nl);
f_1 = ones(1,Nl);
f_2 = ones(1,Nl);
f_3 = ones(1,Nl);
f_4 = ones(1,Nl);

parfor i = 1:Nl
    [g_1, ~] = pde_solver_2D(rlx,rly,z_data,x(i,:));
    l_1_temp = exp(-0.5/sigma/sigma*(g_1 - y)'*(g_1 - y));
    if rly <= K
        g_2 = 0;
        l_2_temp = 0;
    else
        [g_2, ~] = pde_solver_2D(rlx,rly-1,z_data,x(i,:));
        l_2_temp = exp(-0.5/sigma/sigma*(g_2 - y)'*(g_2 - y));
    end
    if rlx <= K
        g_3 = 0;
        l_3_temp = 0;
    else
        [g_3, ~] = pde_solver_2D(rlx-1,rly,z_data,x(i,:));
        l_3_temp = exp(-0.5/sigma/sigma*(g_3 - y)'*(g_3 - y));
    end
    if rlx <= K || rly <= K
        g_4 = 0;
        l_4_temp = 0;
    else
        [g_4, ~] = pde_solver_2D(rlx-1,rly-1,z_data,x(i,:));
        l_4_temp = exp(-0.5/sigma/sigma*(g_4 - y)'*(g_4 - y));
    end
    G_k(i) = max( [l_1_temp,l_2_temp,l_3_temp,l_4_temp] );
    f_1(i) = l_1_temp;
    f_2(i) = l_2_temp;
    f_3(i) = l_3_temp;
    f_4(i) = l_4_temp;
    H(i) = log( G_k(i)^lambda_diff );
end

for j = 1:length(Lambda)
    
    lambda = Lambda(j);
    
    % normalising constant
    Z = Z*sum(exp(H))/Nl;
    W = exp( H - min(H) )./sum(exp( H - min(H)) );
    
    % resampling step
    A = Multinomial_Resampling(W);
    x = x(A',:);
    
    x_rate = 0;
    parfor n = 1:Nl
        rate = 0;
        for i = 1:8
            u_star = x(n,:) + 1.5*randn(1,2);
            [g_1, ~] = pde_solver_2D(rlx,rly,z_data,u_star);
            l_1_temp = exp(-0.5/sigma/sigma*(g_1 - y)'*(g_1 - y));
            
            if rly <= K
                g_2 = 0;
                l_2_temp = 0;
            else
                [g_2, ~] = pde_solver_2D(rlx,rly-1,z_data,u_star);
                l_2_temp = exp(-0.5/sigma/sigma*(g_2 - y)'*(g_2 - y));
            end
            
            if rlx <= K
                g_3 = 0;
                l_3_temp = 0;
            else
                [g_3, ~] = pde_solver_2D(rlx-1,rly,z_data,u_star);
                l_3_temp = exp(-0.5/sigma/sigma*(g_3 - y)'*(g_3 - y));
            end
            
            if rlx <= K || rly <= K
                g_4 = 0;
                l_4_temp = 0;
            else
                [g_4, ~] = pde_solver_2D(rlx-1,rly-1,z_data,u_star);
                l_4_temp = exp(-0.5/sigma/sigma*(g_4 - y)'*(g_4 - y));
            end

            G_star = max( [l_1_temp,l_2_temp,l_3_temp,l_4_temp] );
            
            if u_star(1) >= 1 || u_star(1) <= -1 || u_star(2) >= 1 || u_star(2) <= -1
                alpha_temp = 0;
            else
                alpha_temp = (G_star/G_k(n))^lambda;
            end
            
            alpha = min( 1, alpha_temp );
            
            uni = rand(1);
            
            if uni < alpha
                x(n,:) = u_star;
                G_k(n) = G_star; 
                f_1(n) = l_1_temp;
                f_2(n) = l_2_temp;
                f_3(n) = l_3_temp;
                f_4(n) = l_4_temp;
                H(n) = log( G_star^lambda_diff );
                rate = rate + 1;
            end
        end
        rate = rate/8;
        x_rate = x_rate + rate/Nl;
    end 
    fprintf(1,'lx %d, ly %d, tempering step %.1f, acceptance rate %.5f \n',lx,ly,lambda,x_rate)
end

L2norm2 = x(:,1).^2 + x(:,2).^2;

sninc1  = sum(L2norm2' .* f_1./G_k)/sum(f_1./G_k);
reinc1  = sum(L2norm2' .* f_1./G_k);
reinc11 = sum(f_1./G_k);
if rly <= K
    sninc2  = 0;
    reinc2  = 0;
    reinc12 = 0;
else
    sninc2  = sum(L2norm2' .* f_2./G_k)/sum(f_2./G_k);
    reinc2  = sum(L2norm2' .* f_2./G_k);
    reinc12 = sum(f_2./G_k);
end

if rlx <= K
    sninc3  = 0;
    reinc3  = 0;
    reinc13 = 0;
else
    sninc3  = sum(L2norm2' .* f_3./G_k)/sum(f_3./G_k);
    reinc3  = sum(L2norm2' .* f_3./G_k);
    reinc13 = sum(f_3./G_k);
end

if rlx <= K || rly <= K 
    sninc4  = 0;
    reinc4  = 0;
    reinc14 = 0;
else
    sninc4  = sum(L2norm2' .* f_4./G_k)/sum(f_4./G_k);
    reinc4  = sum(L2norm2' .* f_4./G_k);
    reinc14 = sum(f_4./G_k);
end

sums(1) = (sninc1 - sninc2 - sninc3 + sninc4);  
sums(2) = Z*( reinc1 - reinc2 - reinc3 + reinc4 );
sums(3) = Z*( reinc11 - reinc12 - reinc13 + reinc14);

cost = 4*2^(rlx+rly);
end

