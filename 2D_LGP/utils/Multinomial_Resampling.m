% --------------------------------------- %
% function A = Multinomial_Resampling(W)
% multinomial resampling
% inputs:  W = weights of samples
% outputs: A = categorical distribution
% --------------------------------------- %
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
% --------------------------------------- %
% function A = Multinomial_Resampling(W)
% multinomial resampling
% inputs:  W = weights of samples
% outputs: A = categorical distribution
% --------------------------------------- %

% function A = Systematic_Resampling(W)
%     N = length(W);
%     A = zeros(N,1);
%     W = N * W;
%     j = 1;
%     csw = W(1);
%     U = rand(1);
%     for n = 1:N
%         while (csw < U)
%             j = j + 1;
%             csw = csw + W(j);
%         end
%         A(n) = j;
%         U = U + 1;
%     end
% 
% end