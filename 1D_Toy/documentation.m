% ----------------------------------------------------------- %
% Results are obtained with the following settings:
%   * PDE: - \Delta^2 u(x) ) = x,   on (0,1)
%                       u(x) = 0,   on \partial [0,1]
%   * U[-1,1] prior
%   * Noise N(0,0.2^2) of statistical model
%   * Quantity of interest: x^2 x ~ \pi(dx) (posterior)
%   * 10 observations on linspace(0.1,1,10)
% ----------------------------------------------------------- %

% --------------------------------------------------------------------------------------- %
% Before carrying on tests:
% Please create a folder named 'Results' in the 1D_Toy folder first, 
% and within the 'Results' folder, please create a folder named 'complexity results' 
% and a folder named 'empirical results'
%
% First time: before run complexity test
%   * generate data --> run generator.m for observation data
%   * empirical results --> run runMLregularity.m for empirical mean and var in MLSMC folder
%
% Complexity tests:
%   * MLSMC --> run runMLcomplexity.m for complexity test
%   * rMLSMC --> run runrMLcomplexity.m for complexity test
%   * SMC --> run runSLcomplexity.m for complexity test
% ----------------------------------------------------------------------------------------- %