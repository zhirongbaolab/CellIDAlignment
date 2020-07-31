function [indNeighbor]=mat2cell(adjacencyMat)

% mat2cell - Creates a list of neighbors from an adjacency matrix
%
% CALL:
% [indNeighbor]=mat2cell(adjacencyMat)
%
% INPUT:
% adjacencyMat: (NxN) adjacency matrix containing 1 in (i,j) cell if i and j
%               are neighbors, and 0 otherwise
%
% OUTPUT:
% indNeighbor: N cells containing list of neighbors
%             indNeighbor{i} = [. . . j . ] if adjacencyMat(i,j)~=0 
%
% Author   : Michael Aupetit
%            Qatar Computing Research Institute (QCRI)
%            Hamad Bin Khalifa University (HBKU)
%            Doha, Qatar
%            maupetit@qf.org.qa/michael.aupetit@gmail.com
% Last Rev : Feb 08 2016 by Michael Aupetit

N=length(adjacencyMat(:,1));

for i=1:N
    indNeighbor{i}=find(adjacencyMat(i,:));
end

