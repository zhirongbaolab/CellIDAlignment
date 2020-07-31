function [NNG]=NearestNeighborGraph(dist)

% NearestNeighborGraph - Computes the Nearest Neighbor graph 
%
% Computes the Nearest Neighbor Graph  (NNG) of a set of points 
% NNG is a DIRECTED graph
%
% CALL:
% [NNG]=NearestNeighborGraph(dist)
% 
% INPUT:
% dist: NxN Euclidean distance matrix between each pair of points
%               dist(i,k)=norm(data(i,:)-data(k,:), (dist(k,k)=0)
%              
% OUPUT:
% NNG: N cells NNG{i}=[a b...f] set of index in data rows of the NNG
% neighbors of i
% 
% Author   : Michael Aupetit
%            Qatar Computing Research Institute (QCRI)
%            Hamad Bin Khalifa University (HBKU)
%            Doha, Qatar
%            maupetit@qf.org.qa/michael.aupetit@gmail.com
% Last Rev : Feb 08 2016 by Michael Aupetit

[N,N]=size(dist);
dist(eye(N)==1)=inf;
[val, indMin]=min(dist);
NNG=cell(N,1);

for i=1:N
    NNG{i}=indMin(i);
end