% ---------------------------------------------------------------------- %
% Results are obtained with the following settings:
%   * PDE: - \Delta \cdot a(x) \Delta u(x) ) = 100, on (0,1)^2
%                                       u(x) = 0,   on \partial [0,1]^2
%   * a(X)(Z) = 3 + x_1cos(3z_1)sin(3z_2) + x_2cos(z_1)sin(z_2)
%   * U[-1,1]^2 prior
%   * Noise N(0,0.5^2) of statistical model
%   * Quantity of interest: x_1^2 + x_2^2, x_1, x_2 ~ \pi(dx) (posterior)
%   * 4 observations on [0.25, 0.25; 0.25, 0.75; 0.75, 0.25; 0.75, 0.75]
% ----------------------------------------------------------------------- %

% --------------------------------------------------------------------------------------- %
% Before carrying on tests:
% Please create a folder named 'Results' in folder named '2D_PDE' first, 
% and within the 'Results' folder, please create a folder named 'complexity results' 
% and a folder named 'empirical results'
%
% First time: before run complexity test
%   * generate data --> run generator.m for observation data in solution and data folder
%   * empirical results --> run runMLregularity.m for empirical mean and var in solution and data folder
%   * solution --> run solution.m for approximate solution in solution and data folder
%   * empirical results --> run runMIregularity.m for empirical mean and var in MISMC folder
%
% Complexity tests:
%   * MISMC --> run runMIcomplexity.m for complexity test
%   * rMISMC --> run runrMIcomplexity.m for complexity test
%   * SMC --> run runSLcomplexity.m for complexity test
% ----------------------------------------------------------------------------------------- %