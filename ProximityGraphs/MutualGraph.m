function [MG]=MutualGraph(G)

% MutualGraph - Replaces two-ways directed edges by undirected edges
% 
% CALL:
% [MG]=MutualGraph(G)
%
% INPUT:
% G: N cells containing neighbors need not to be a symetric graph
%
% OUTPUT:
% MG: N cells containing mutual neighbors: if j in G{i} and i in G{j} then
% put j in MG{i} and i in MG{j}  (logical AND)
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
MG=cell(N,1);
    
for i=1:N
    
    for j=1:length(G{i})
        vois=G{i}(j);
        if ismember(i,G{vois});
            MG=addNeighbor(MG,i,vois);
            MG=addNeighbor(MG,vois,i);
        end
    end
end
MG=uniqueNeighbor(MG);