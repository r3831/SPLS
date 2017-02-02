%% UPDATESGD_CXY updates the SVD of the cross-covariance matrix Cxy given a 
% new sample-pair (x,y) using the stochastic gradient descent(SGD) algorithm
% Usage:
% Inputs:  U,S,V - current SVD  
%          x,y   - new sample pair
%          RANK  - maximum SVD-rank desired
% Outputs: U,S,V - updated SVD
%
%%

function [U,V]=updateSGD_Cxy(U,V,x,y,etax,etay) 

% m=length(x);         % Number of new columns
% n=length(y);         % Number of new columns
% if(m~=size(U,1)||(n~=size(V,1))) % Check dimensionality
%     error('Dimensions of SVD are incompatible with data-dimensionality\n');
% end
U=U+etax*(x*(y'*V));
V=V+etay*(y*(x'*U));
end % End-function
