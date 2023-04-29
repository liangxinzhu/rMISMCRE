% -------------------------------------------------------------------------- %
% function [varargout] = mutation(lx,ly,rate,a,b,mu,xpos,ypos,type,varargin)
% inputs:  lx           = level of refinement in x direction
%          ly           = level of refinement in y direction
%          rate         = \beta^{\prime}
%          a            = \theta_2
%          b            = \theta_3
%          mu           = \theta_1
%          xpos         = x coordinate of observation points
%          ypos         = y coordinate of observation points
%          type         = 1: single level
%          type         = 2: increments of increments along x
%          type         = 3: increments of increments along y
%          type         = 4: increments of increments
%          type         = 5: increments
%          varargin{1}  = u1
%          varargin{2}  = u2
%          varargin{3}  = u3
%          varargin{4}  = u4
% ouputs:  varargout    = v1-4, loglikelihood 1-4
% -------------------------------------------------------------------------- %

function [varargout] = ...
    mutation(lx,ly,rate,a,b,mu,xpos,ypos,type,varargin)
% pCN rate
beta = 0.65;%0.65;%0.95;
n = 126;
% mutate samples with pCN
switch type
    case 1
        [samp1] = LGprior(lx,ly,rate,a,b,mu,1);
        v1 = sqrt(1-beta^2).*(varargin{1}(:) - mu) + mu + beta.*(samp1(:)-mu);
%         v1 = sqrt(1-beta^2).*varargin{1}(:) + beta.*samp1(:);
        fun1 = interpob(reshape(v1,size(samp1)),2^lx,2^ly,xpos,ypos);
        llik1 = sum(fun1(:)-mu) -  n*log(sum(exp(v1(:)))/2^lx/2^ly);
        
        varargout{1} = v1;
        varargout{2} = llik1;
    case 2
        [samp1,samp2] = LGprior(lx,ly,rate,a,b,mu,2);
        v1 = sqrt(1-beta^2).*(varargin{1}(:) - mu) + mu + beta.*(samp1(:)-mu);
        v2 = sqrt(1-beta^2).*(varargin{2}(:) - mu) + mu + beta.*(samp2(:)-mu);
%         v1 = sqrt(1-beta^2).*varargin{1}(:) + beta.*samp1(:);
%         v2 = sqrt(1-beta^2).*varargin{2}(:) + beta.*samp2(:);
        fun1 = interpob(reshape(v1,size(samp1)),2^lx,2^ly,xpos,ypos);
        llik1 = sum(fun1(:)-mu) -  n*log(sum(exp(v1(:)))/2^lx/2^ly);
        fun2 = interpob(reshape(v2,size(samp2)),2^lx/2,2^ly,xpos,ypos);
        llik2 = sum(fun2(:)-mu) -  n*log(sum(exp(v2(:)))/2^lx/2^ly*2);
        
        varargout{1} = v1;
        varargout{2} = llik1;
        varargout{3} = v2;
        varargout{4} = llik2;        
    case 3
        [samp1,samp2] = LGprior(lx,ly,rate,a,b,mu,3);
        v1 = sqrt(1-beta^2).*(varargin{1}(:) - mu) + mu + beta.*(samp1(:)-mu);
        v2 = sqrt(1-beta^2).*(varargin{2}(:) - mu) + mu + beta.*(samp2(:)-mu);
%         v1 = sqrt(1-beta^2).*varargin{1}(:) + beta.*samp1(:);
%         v2 = sqrt(1-beta^2).*varargin{2}(:) + beta.*samp2(:);
        fun1 = interpob(reshape(v1,size(samp1)),2^lx,2^ly,xpos,ypos);
        llik1 = sum(fun1(:)-mu) -  n*log(sum(exp(v1(:)))/2^lx/2^ly);
        fun2 = interpob(reshape(v2,size(samp2)),2^lx,2^ly/2,xpos,ypos);
        llik2 = sum(fun2(:)-mu) -  n*log(sum(exp(v2(:)))/2^lx/2^ly*2);
        
        varargout{1} = v1;
        varargout{2} = llik1;
        varargout{3} = v2;
        varargout{4} = llik2;  
    case 4
        [samp1,samp2,samp3,samp4] = LGprior(lx,ly,rate,a,b,mu,4);
        v1 = sqrt(1-beta^2).*(varargin{1}(:) - mu) + mu + beta.*(samp1(:)-mu);
        v2 = sqrt(1-beta^2).*(varargin{2}(:) - mu) + mu + beta.*(samp2(:)-mu);
        v3 = sqrt(1-beta^2).*(varargin{3}(:) - mu) + mu + beta.*(samp3(:)-mu);
        v4 = sqrt(1-beta^2).*(varargin{4}(:) - mu) + mu + beta.*(samp4(:)-mu);
%         v1 = sqrt(1-beta^2).*varargin{1}(:) + beta.*samp1(:);
%         v2 = sqrt(1-beta^2).*varargin{2}(:) + beta.*samp2(:);
%         v3 = sqrt(1-beta^2).*varargin{3}(:) + beta.*samp3(:);
%         v4 = sqrt(1-beta^2).*varargin{4}(:) + beta.*samp4(:);
        fun1 = interpob(reshape(v1,size(samp1)),2^lx,2^ly,xpos,ypos);
        llik1 = sum(fun1(:)-mu) -  n*log(sum(exp(v1(:)))/2^lx/2^ly);
        fun2 = interpob(reshape(v2,size(samp2)),2^lx/2,2^ly,xpos,ypos);
        llik2 = sum(fun2(:)-mu) -  n*log(sum(exp(v2(:)))/2^lx/2^ly*2);
        fun3 = interpob(reshape(v3,size(samp3)),2^lx,2^ly/2,xpos,ypos);
        llik3 = sum(fun3(:)-mu) -  n*log(sum(exp(v3(:)))/2^lx/2^ly*2);
        fun4 = interpob(reshape(v4,size(samp4)),2^lx/2,2^ly/2,xpos,ypos);
        llik4 = sum(fun4(:)-mu) -  n*log(sum(exp(v4(:)))/2^lx/2^ly*4);
        
        varargout{1} = v1;
        varargout{2} = llik1;
        varargout{3} = v2;
        varargout{4} = llik2;  
        varargout{5} = v3;
        varargout{6} = llik3;
        varargout{7} = v4;
        varargout{8} = llik4;          
    case 5
        [samp1,samp2] = LGprior(lx,ly,rate,a,b,mu,5);
        v1 = sqrt(1-beta^2).*(varargin{1}(:) - mu) + mu + beta.*(samp1(:)-mu);
        v2 = sqrt(1-beta^2).*(varargin{2}(:) - mu) + mu + beta.*(samp2(:)-mu);
%         v1 = sqrt(1-beta^2).*varargin{1}(:) + beta.*samp1(:);
%         v2 = sqrt(1-beta^2).*varargin{2}(:) + beta.*samp2(:);
        fun1 = interpob(reshape(v1,size(samp1)),2^lx,2^ly,xpos,ypos);
        llik1 = sum(fun1(:)-mu) -  n*log(sum(exp(v1(:)))/2^lx/2^ly);
        fun2 = interpob(reshape(v2,size(samp2)),2^lx/2,2^ly/2,xpos,ypos);
        llik2 = sum(fun2(:)-mu) -  n*log(sum(exp(v2(:)))/2^lx/2^ly*4);
        
        varargout{1} = v1;
        varargout{2} = llik1;
        varargout{3} = v2;
        varargout{4} = llik2;  
    otherwise
        fprintf(1,'Wrong Prior Type')
end
end

