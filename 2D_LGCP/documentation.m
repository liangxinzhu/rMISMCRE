% ----------------------------------------------------------------------- %
% For detailed explanation of the process, please refer to the paper.
% Results are obatained with the following parameter settings:
%   * the log-Gaussian Cox process:
%     Parameters: (\theta_1, \theta_2, \theta_3) = (0, 1, 110.339)
% ----------------------------------------------------------------------- %

% --------------------------------------------------------------------------------------- %
% Before carrying on tests:
% Please create a folder named 'Results' in folder named '2D_LGCP' first, 
% and within the results folder, please create a folder named 'complexity results' 
% and a folder named 'empirical results'
%
% First time: before run complexity test
%   * solution --> run solution.m for approximate solution
%   * empirical results --> run runMIregularity.m for empirical mean and var in MISMC folder
%
% Complexity tests:
%   * MISMC --> run runMIcomplexity.m for complexity test
%   * rMISMC --> run runrMIcomplexity.m for complexity test
% ----------------------------------------------------------------------------------------- %