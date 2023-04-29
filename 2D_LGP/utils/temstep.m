% ------------------------------------------------------------------ %
% [lambda] = temstep(lambda_pre,llik,ESSmin)    
% choosen tempering size such that ESS = ESSmin
%
% inputs:  lambda_pre = previous adaptive tempering size
%          llik       = log likelihood
%          ESSmin     = effective sample size threshold
%
% output:  lambda     = new adaptive tempering size
% ------------------------------------------------------------------ %

function [lambda] = temstep(lambda_pre,llik,ESSmin)

f = @(x) sum(exp(llik.*x - max(llik)*x),'all')^2/...
    sum((exp(llik.*x - max(llik)*x)).^2,'all')...
    - ESSmin;

if f(1 - lambda_pre) > 0
    lambda = 1;
else
    delta = fzero(f,[1e-12, 1-lambda_pre]);
    lambda = lambda_pre + delta;
end

end

