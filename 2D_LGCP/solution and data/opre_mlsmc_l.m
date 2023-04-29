% ----------------------------------------------------------------- %
% function [sums, cost] = opre_mlsmc_l(lx,ly,N,data,rate,a,b,mu,K)
%
% single level routine
%
% inputs:  lx   = level of refinement in x direction
%          ly   = level of refinement in y direction
%          N    = number of samples
%          data = observation points
%          rate = r+1
%          a    = \theta_2
%          b    = \theta_3
%          mu   = \theta_1
%          K    = starting level
%
% outputs: sums(1) = normalised increments of increments wpt QOI
%          sums(2) = unnormalised increments of increments wpt QOI
%          sums(3) = unnormalised increments of increments wpt 1
%          cost    = cost of one sample
% ----------------------------------------------------------------- %
function [sums, cost] = opre_mlsmc_l(lx,ly,N,data,rate,a,b,mu,K)
sums(1:3) = 0;
ESSmin = round(N/2);
xpos = data(:,1);
ypos = data(:,2);

Nx = 2^lx; Ny = 2^ly;

u1 = zeros(N,Nx*Ny);
u2 = zeros(N,Nx*Ny/4);

G  = zeros(1,N);
f1 = zeros(1,N);
f2 = zeros(1,N);

Z = 0;
lambda = 0;
for i = 1:N
    % sampling from prior
    if lx == K && ly == K
        [prior1] = LGCprior(lx,ly,rate,a,b,mu,1);
        fun1 = interpob(prior1,Nx,Ny,xpos,ypos);
        llik1  = sum(fun1(:)-mu) -  sum(exp(prior1(:)))/Nx/Ny;
        f1(i) = llik1;
        u1(i,:) = prior1(:);
        G(i) = llik1; 
    else
        [prior1,prior2] = LGCprior(lx,ly,rate,a,b,mu,5);
        fun1  = interpob(prior1,Nx,Ny,xpos,ypos);
        fun2  = interpob(prior2,Nx/2,Ny/2,xpos,ypos);
        llik1 = sum(fun1(:)-mu) -  sum(exp(prior1(:)))/Nx/Ny;
        llik2 = sum(fun2(:)-mu) -  sum(exp(prior2(:)))/Nx/Ny*4;
        f1(i) = llik1;
        f2(i) = llik2;
        u1(i,:) = prior1(:);
        u2(i,:) = prior2(:);
        G(i) = log( max( [exp(llik1),exp(llik2)] ) ); 
    end
        
end

count = 0;

while lambda < 1
    lambdaold = lambda;
    lambda = temstep(lambdaold,G,ESSmin);
    
    if count == 0 && lambda == 1
        lambda = 0.5;
    end
    
    count = count + 1;
    
    H = G.*(lambda-lambdaold); 
    Z = Z + log( mean(exp(H)) );
    W = exp( H - max(H) )./sum(exp( H - max(H)) );
    A = Multinomial_Resampling(W);
    u1 = u1(A',:);
    u2 = u2(A',:);
    
    x_rate = 0;
    for j = 1:N
        acc = 0;
        if lx == K && ly == K
            for i = 1:8
                [v1,llik1_temp] = ...
                    mutation(lx,ly,rate,a,b,mu,xpos,ypos,1,u1(j,:));
                G_star = llik1_temp;
                if log(rand) < lambda*(G_star-G(j))
                    u1(j,:) = v1(:);
                    G(j)   = G_star;
                    f1(j)  = llik1_temp;
                    acc = acc + 1;
                end  
            end
        else
            for i = 1:8
                [v1,llik1_temp,v2,llik2_temp] = ...
                    mutation(lx,ly,rate,a,b,mu,xpos,ypos,5,u1(j,:),u2(j,:));
                G_star = log( max( [exp(llik1_temp),exp(llik2_temp)] ) ); 
                if log(rand) < lambda*(G_star-G(j))
                    u1(j,:) = v1(:);
                    u2(j,:) = v2(:);
                    G(j)   = G_star;
                    f1(j)  = llik1_temp;
                    f2(j)  = llik2_temp;
                    acc = acc + 1;                    
                end
            end
        end
        acc = acc/8;
        x_rate = x_rate + acc/N;
    end
    fprintf(1,'lx: %d, ly: %d, tempering step: %.4f, acceptance rate: %.4f \n',lx,ly,lambda,x_rate)
end

QOI1 = sum(exp(u1), 2)/Nx/Ny;

sninc1  = sum(QOI1' .* exp(f1-G))/sum(exp(f1-G));
reinc1  = sum(QOI1' .* exp(f1-G));
reinc11 = sum(exp(f1-G));
cost = (lx+ly)*2^(lx+ly);

if lx <= K || ly <= K 
    sninc2  = 0;
    reinc2  = 0;
    reinc12 = 0;
else
    QOI2 = sum(exp(u2), 2)/Nx/Ny*4;
    sninc2  = sum(QOI2' .* exp(f2-G))/sum(exp(f2-G));
    reinc2  = sum(QOI2' .* exp(f2-G));
    reinc12 = sum(exp(f2-G));
    cost = cost + (lx+ly-2)*2^(lx+ly-2);
end

sums(1) = sninc1 - sninc2;  
sums(2) = exp(Z)*( reinc1 - reinc2 )/N;
sums(3) = exp(Z)*( reinc11 - reinc12)/N;

end


