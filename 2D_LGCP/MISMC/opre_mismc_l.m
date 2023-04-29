% ----------------------------------------------------------------- %
% function [sums, cost] = opre_mismc_l(lx,ly,N,data,rate,a,b,mu,K)
%
% single level routine
%
% inputs: lx   = level of refinement in x direction
%         ly   = level of refinement in y direction
%         N    = number of samples
%         data = observation points
%         rate = r+1
%         a    = \theta_2
%         b    = \theta_3
%         mu   = \theta_1
%         K    = actual starting refinement level
%
% outputs: sums(1) = normalised increments of increments wpt QOI
%          sums(2) = unnormalised increments of increments wpt QOI
%          sums(3) = unnormalised increments of increments wpt 1
%          cost    = cost of one sample
% ----------------------------------------------------------------- %
function [sums, cost] = opre_mismc_l(lx,ly,N,data,rate,a,b,mu,K)
sums(1:3) = 0;
ESSmin = round(N/2);
xpos = data(:,1);
ypos = data(:,2);

Nx = 2^lx; Ny = 2^ly;

u1 = zeros(N,Nx*Ny);
u2 = zeros(N,Nx*Ny/2);
u3 = zeros(N,Nx*Ny/2);
u4 = zeros(N,Nx*Ny/4);

G  = zeros(1,N);
f1 = zeros(1,N);
f2 = zeros(1,N);
f3 = zeros(1,N);
f4 = zeros(1,N);

Z = 0;
lambda = 0;
for i = 1:N
    % sampling from prior
    if lx == K && ly == K
        [prior1] = LGCprior(lx,ly,rate,a,b,mu,1);
        fun1 = interpob(prior1,Nx,Ny,xpos,ypos);
        llik1 = sum(fun1(:)-mu) - sum(exp(prior1(:)))/Nx/Ny;
        f1(i) = llik1;
        u1(i,:) = prior1(:);
        G(i) = llik1; 
    elseif ly == K && lx > K
        [prior1,prior2] = LGCprior(lx,ly,rate,a,b,mu,2);
        
        fun1 = interpob(prior1,Nx,Ny,xpos,ypos);
        fun2 = interpob(prior2,Nx/2,Ny,xpos,ypos);
        
        llik1  = sum(fun1(:)-mu) - sum(exp(prior1(:)))/Nx/Ny;
        llik2  = sum(fun2(:)-mu) - sum(exp(prior2(:)))/Nx/Ny*2;

        f1(i) = llik1;
        f2(i) = llik2;
        
        u1(i,:) = prior1(:);
        u2(i,:) = prior2(:);
        
        G(i) = log( max( [exp(llik1),exp(llik2)] ) ); 
    elseif lx == K && ly > K
        [prior1,prior3] = LGCprior(lx,ly,rate,a,b,mu,3);
        
        fun1 = interpob(prior1,Nx,Ny,xpos,ypos);
        fun3 = interpob(prior3,Nx,Ny/2,xpos,ypos);
        
        llik1 = sum(fun1(:)-mu) - sum(exp(prior1(:)))/Nx/Ny;
        llik3 = sum(fun3(:)-mu) - sum(exp(prior3(:)))/Nx/Ny*2;

        f1(i) = llik1;
        f3(i) = llik3;
        
        u1(i,:) = prior1(:);
        u3(i,:) = prior3(:);
        
        G(i) = log( max( [exp(llik1),exp(llik3)] ) ); 
    else
        [prior1,prior2,prior3,prior4] = LGCprior(lx,ly,rate,a,b,mu,4);
        
        fun1 = interpob(prior1,Nx,Ny,xpos,ypos);
        fun2 = interpob(prior2,Nx/2,Ny,xpos,ypos);
        fun3 = interpob(prior3,Nx,Ny/2,xpos,ypos);
        fun4 = interpob(prior4,Nx/2,Ny/2,xpos,ypos);
        
        llik1  = sum(fun1(:)-mu) - sum(exp(prior1(:)))/Nx/Ny;
        llik2  = sum(fun2(:)-mu) - sum(exp(prior2(:)))/Nx/Ny*2;
        llik3  = sum(fun3(:)-mu) - sum(exp(prior3(:)))/Nx/Ny*2;
        llik4  = sum(fun4(:)-mu) - sum(exp(prior4(:)))/Nx/Ny*4;

        f1(i) = llik1;
        f2(i) = llik2;
        f3(i) = llik3;
        f4(i) = llik4;
        
        u1(i,:) = prior1(:);
        u2(i,:) = prior2(:);
        u3(i,:) = prior3(:);
        u4(i,:) = prior4(:);
        
        G(i) = log( max( [exp(llik1),exp(llik2),exp(llik3),exp(llik4)] ) ); 
    end
end
count = 0;
% tempering
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
    u3 = u3(A',:);
    u4 = u4(A',:);
    % mutation
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
            acc = acc/8;
        elseif ly <= K && lx > K
            for i = 1:8
                [v1,llik1_temp,v2,llik2_temp] = ...
                    mutation(lx,ly,rate,a,b,mu,xpos,ypos,2,u1(j,:),u2(j,:));
                G_star = log( max( [exp(llik1_temp), exp(llik2_temp)] ) );
                if log(rand) < lambda*(G_star-G(j))
                    u1(j,:) = v1(:);
                    u2(j,:) = v2(:);
                    G(j)   = G_star;
                    f1(j)  = llik1_temp;
                    f2(j)  = llik2_temp;
                    acc = acc + 1;
                end
            end
        elseif lx <= K && ly > K
            for i = 1:8
                [v1,llik1_temp,v3,llik3_temp] = ...
                    mutation(lx,ly,rate,a,b,mu,xpos,ypos,3,u1(j,:),u3(j,:));
                G_star = log( max( [exp(llik1_temp), exp(llik3_temp)] ) );
                if log(rand) < lambda*(G_star-G(j))
                    u1(j,:) = v1(:);
                    u3(j,:) = v3(:);
                    G(j)   = G_star;
                    f1(j)  = llik1_temp;
                    f3(j)  = llik3_temp;
                    acc = acc + 1;
                end
            end
        else
            for i = 1:8
                [v1,llik1_temp,v2,llik2_temp,v3,llik3_temp,v4,llik4_temp] = ...
                    mutation(lx,ly,rate,a,b,mu,xpos,ypos,4,u1(j,:),u2(j,:),u3(j,:),u4(j,:));
                G_star = log( max( [exp(llik1_temp), exp(llik2_temp),...
                    exp(llik3_temp), exp(llik4_temp)] ) );
                if log(rand) < lambda*(G_star-G(j))
                    u1(j,:) = v1(:);
                    u2(j,:) = v2(:);
                    u3(j,:) = v3(:);
                    u4(j,:) = v4(:);
                    G(j)   = G_star;
                    f1(j)  = llik1_temp;
                    f2(j)  = llik2_temp;
                    f3(j)  = llik3_temp;
                    f4(j)  = llik4_temp;
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

if lx == K
    sninc2  = 0;
    reinc2  = 0;
    reinc12 = 0;
else
    QOI2 = sum(exp(u2), 2)/Nx/Ny*2;
    sninc2  = sum(QOI2' .* exp(f2-G))/sum(exp(f2-G));
    reinc2  = sum(QOI2' .* exp(f2-G));
    reinc12 = sum(exp(f2-G));
    cost = cost + (lx+ly-1)*2^(lx+ly-1);
end

if ly == K
    sninc3  = 0;
    reinc3  = 0;
    reinc13 = 0;
else
    QOI3 = sum(exp(u3), 2)/Nx/Ny*2;
    sninc3  = sum(QOI3' .* exp(f3-G))/sum(exp(f3-G));
    reinc3  = sum(QOI3' .* exp(f3-G));
    reinc13 = sum(exp(f3-G));
    cost = cost + (lx+ly-1)*2^(lx+ly-1);
end

if lx == K || ly == K 
    sninc4  = 0;
    reinc4  = 0;
    reinc14 = 0;
else
    QOI4 = sum(exp(u4), 2)/Nx/Ny*4;
    sninc4  = sum(QOI4' .* exp(f4-G))/sum(exp(f4-G));
    reinc4  = sum(QOI4' .* exp(f4-G));
    reinc14 = sum(exp(f4-G));
    cost = cost + (lx+ly-2)*2^(lx+ly-2);
end

sums(1) = sninc1 - sninc2 - sninc3 + sninc4;  
sums(2) = exp(Z)*( reinc1 - reinc2 - reinc3 + reinc4 )/N;
sums(3) = exp(Z)*( reinc11 - reinc12 - reinc13 + reinc14)/N;

end

