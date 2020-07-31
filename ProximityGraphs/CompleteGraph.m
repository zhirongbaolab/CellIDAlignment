function [CG]=CompleteGraph(N)

% CompleteGraph - Create a complete graph
% 
% CALL:
% [CG]=CompleteGraph(N)
%
% INPUT:
% N: number of vertices
%
% OUTPUT:
% CG: N cells containing all vertices as neighbors except itself
%
% Author   : Michael Aupetit
%            Qatar Computing Research Institute (QCRI)
%            Hamad Bin Khalifa University (HBKU)
%            Doha, Qatar
%            maupetit@qf.org.qa/michael.aupetit@gmail.com
% Last Rev : Feb 08 2016 by Michael Aupetit


CG{1}=[2:N];
CG{N}=[1:N-1];
for k=2:N-1
    CG{k}=[1:k-1 k+1:N];
end