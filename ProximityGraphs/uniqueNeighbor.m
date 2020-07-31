function [UG]=uniqueNeighbor(G)

% uniqueNeighbor - Remove neighbors replica from lists of neighbors
% 
% CALL:
% [UG]=uniqueNeighbor(G)
%
% INPUT:
% G: N cells containing neighbors 
%
% OUTPUT:
% UG: N cells containing same content as G without replica
%     G{2}=[1 4 3 2 4 3 3 5] -> UG{2}=[1 2 3 4 5] 
%
% Author   : Michael Aupetit
%            Qatar Computing Research Institute (QCRI)
%            Hamad Bin Khalifa University (HBKU)
%            Doha, Qatar
%            maupetit@qf.org.qa/michael.aupetit@gmail.com
% Last Rev : Feb 08 2016 by Michael Aupetit


N=length(G);
UG=cell(N,1);
for i=1:N
    UG{i}=unique(G{i});
end
