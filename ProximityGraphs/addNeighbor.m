function [indNeighbor]=addNeighbor(indNeighbor,indRef,indToAdd)

% addNeighbor - Create/Append neighbors to a vertex 
%
% CALL:
% [indNeighbor]=addNeighbor(indNeighbor,indRef,indToAdd)
% 
% INPUT:
% indNeighbor: list of neighbors, indNeighbor{k} : list of neighbors of vertex k
% indRef : vertex to which neighbors indToAdd have to be added
% indToAdd : list of vertices to append to the list of neighbors of the indRef vertex
%
% OUTPUT:
% indNeighbor: updated list of neighbors
%
% Author   : Michael Aupetit
%            Qatar Computing Research Institute (QCRI)
%            Hamad Bin Khalifa University (HBKU)
%            Doha, Qatar
%            maupetit@qf.org.qa/michael.aupetit@gmail.com
% Last Rev : Feb 08 2016 by Michael Aupetit


if indRef<=length(indNeighbor)
    indNeighbor{indRef}=[indNeighbor{indRef} indToAdd];
else
    indNeighbor{indRef}=indToAdd;
end
