% --------------------------------------------------------------------------- %
% function ml_empirical(mlsmc_l,N,Lmax,sigma,data,z_data)
%
% multilevel Monte Carlo empirical test routine
%
% Input: mlmc_l = function for level l estimator
%        N      = number of samples for each level
%        Lmax   = maximum level for convergence test
%        sigma  = std deviation of the error of y - G
%        data   = observations (values of y) vector
%        z_data = corresponding observation points
%
% .mat: Lmax = maximum level we set
%       cost = cost at each level
%       we1  = difference with self-normalised increments estimator
%       we2  = difference unnormalised integral
%       we3  = difference normalising constant
%       se1  = variance of sn increments estimtor
%       se2  = unnormalised integral of squared target
%       se3  = unnormalised integral of 1
%       se4  = square of difference with self-normalised increments estimator
%       se5  = square of difference unnormalised integral
%       se6  = square of difference normalising constant
% --------------------------------------------------------------------------- %
function ml_empirical(mlsmc_l, N, Lmax, sigma, data, z_data)

we1  = zeros(Lmax,1);
we2  = zeros(Lmax,1);
we3  = zeros(Lmax,1);
se1  = zeros(Lmax,1);
se2  = zeros(Lmax,1);
se3  = zeros(Lmax,1);
se4  = zeros(Lmax,1);
se5  = zeros(Lmax,1);
se6  = zeros(Lmax,1);
cost = zeros(Lmax,1);

for l = 1:Lmax  

    sums = [];
    cst  = 0;

    M = 20;
    
    parfor j = 1:M  
        [sums_j, cst_j] = mlsmc_l(l,N,sigma,data,z_data);
        sums(j,:) = sums_j;
        cst  = cst  + cst_j/M;   
    end
    
    cost(l) = cst;
    we1(l) = abs( sum(sums(:,1))/M );
    we2(l) = abs( sum(sums(:,2))/M );
    we3(l) = abs( sum(sums(:,3))/M );
    se1(l) = N*sum( (sums(:,1) - sum(sums(:,1))/M).^2 )/M;
    se2(l) = N*sum( (sums(:,2) - sum(sums(:,2))/M).^2 )/M;
    se3(l) = N*sum( (sums(:,3) - sum(sums(:,3))/M).^2 )/M;
    se4(l)  = sum(sums(:,1).^2)/M;
    se5(l)  = sum(sums(:,2).^2)/M;
    se6(l)  = sum(sums(:,3).^2)/M;
        
end

we = {we1,we2,we3};
se = {se1,se2,se3,se4,se5,se6};

[filepart,~,~] = fileparts(pwd); 
savepath = fullfile(filepart,'Results','empirical results','mlsmc_empirical.mat');
save(savepath,'Lmax','cost','we','se')

end