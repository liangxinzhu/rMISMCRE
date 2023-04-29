% ---------------------------------------------------------%
% function A = Multinomial_Resampling(w)
%
% Multinomial_Resampling is used for resampling
%
% input: W = weights of samples
%
% output: A = sample indices with categorical distribution
% ---------------------------------------------------------%
function A = Multinomial_Resampling(W)
     N = length(W);
     
     s = W(1);
     m = 1;
     A = zeros(1,N);
     
     U = rand(1,N);
     U = sort(U);
     
     for  n = 1:N
         while s < U(n)
             m = m + 1;
             s = s + W(m);
         end
         A(n) = m;
     end
end