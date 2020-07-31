function [adjacencyMat]=cell2mat(indNeighbor)

% cell2mat - Creates an adjacency matrix from a list of neighbors
%
% CALL:
% [adjacencyMat]=cell2mat(indNeighbor)
%
% INPUT:
% indNeighbor: N cells containing list of neighbors
%             indNeighbor{i} = [. . . j . ] if adjacencyMat(i,j)~=0 
%
% OUTPUT:
% adjacencyMat: (NxN) non-symmetrical adjacency matrix containing 1 in (i,j) cell
%                                       indNeighbor{i} = [. . . j . ]
%
% Author   : Michael Aupetit
%            Qatar Computing Research Institute (QCRI)
%            Hamad Bin Khalifa University (HBKU)
%            Doha, Qatar
%            maupetit@qf.org.qa/michael.aupetit@gmail.com
% Last Rev : Feb 08 2016 by Michael Aupetit


N=length(indNeighbor);

adjacencyMat=zeros(N,N);
for i=1:N
    adjacencyMat(i,indNeighbor{i})=1;
end

