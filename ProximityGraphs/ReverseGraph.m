function [RG]=ReverseGraph(G)

% ReverseGraph - Reverse a graph, every non neighbor becomes neighbor and vice-versa
% 
% CALL:
% [RG]=ReverseGraph(G)
%
% INPUT:
% G: a graph as N cells lists of neighbors
%
% OUTPUT:
% RG: N cells containing all vertices as neighbors except itself and
% vertices in G. RG{i}=S\{G{i} U i} where S is the set 1:N
%
% Author   : Michael Aupetit
%            Qatar Computing Research Institute (QCRI)
%            Hamad Bin Khalifa University (HBKU)
%            Doha, Qatar
%            maupetit@qf.org.qa/michael.aupetit@gmail.com
% Last Rev : Feb 08 2016 by Michael Aupetit

N=length(G);
RG=cell(N,1);
for k=1:N
    RG{k}=setdiff(1:N,[G{k} k]);
end