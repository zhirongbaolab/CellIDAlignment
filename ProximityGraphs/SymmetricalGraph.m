function [SG]=SymmetricalGraph(G)

% SymmetricalGraph - Replaces any directed edges by undirected edges
% 
% CALL:
% [SG]=SymmetricalGraph(G)
%
% INPUT:
% G: N cells containing neighbors need not to be a symetric graph
%
% OUTPUT:
% SG: N cells containing symetrized neighbors: if j in G{i} then 
% put j in SG{i} and i in SG{j}  (logical OR)
%
% NEEDS: addNeighbor; uniqueNeighbor.
%
% Author   : Michael Aupetit
%            Qatar Computing Research Institute (QCRI)
%            Hamad Bin Khalifa University (HBKU)
%            Doha, Qatar
%            maupetit@qf.org.qa/michael.aupetit@gmail.com
% Last Rev : Feb 08 2016 by Michael Aupetit

N=length(G);
SG=cell(N,1);
    
for i=1:N
    SG=addNeighbor(SG,i,G{i});
    for j=1:length(G{i})
       SG=addNeighbor(SG,G{i}(j),i); 
    end
end
SG=uniqueNeighbor(SG);