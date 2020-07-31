function [RNG]=RelativeNeighborhoodGraph(dist)

% RelativeNeighborhoodGraph - Computes the Relative Neighborhood Graph 
%
% Computes the Relative Neighborhood Graph of a set of points based on the Euclidean norm
% The Relative Neighborhood Graph (RNG) [1] is an empty region graphs in the sense defined in [2].
% RNG is a subgraph of the Gabriel Graph which connects two points p and q
% if no other point s has a distance to p and q smaller than the distance 
% between p and q. In other words, p and q are connected if the
% intersection of two balls, one centered at p intercepting q and the other 
% one centered at q intercepting p, is empty of any other point.
% RNG is UNDIRECTED
%
% CALL:
% [RNG]=RelativeNeighborhoodGraph(dist)
% 
% INPUT:
% dist: NxN Euclidean distance matrix between each pair of points
%               dist(i,k)=norm(data(i,:)-data(k,:), (dist(k,k)=0)
%
% OUPUT:
% RNG: N cells RNG{i}=[a b...f] set of index in data rows of the RNG neighbors of i
%
% NEEDS: addNeighbor.
%
% [1]  Toussaint, G. T. (1980), "The relative neighborhood graph of a
% finite planar set", Pattern Recognition 12 (4): 261–268, 
% doi:10.1016/0031-3203(80)90066-7.
% [2] Jean Cardinal, Sébastien Collette, Stefan Langerman, Empty region
%     graphs, Computational Geometry, Volume 42, Issue 3, April 2009, 
%     Pages 183-195, ISSN 0925-7721,
%     http://dx.doi.org/10.1016/j.comgeo.2008.09.003.
%
% Author   : Michael Aupetit
%            Qatar Computing Research Institute (QCRI)
%            Hamad Bin Khalifa University (HBKU)
%            Doha, Qatar
%            maupetit@qf.org.qa/michael.aupetit@gmail.com
% Last Rev : Feb 08 2016 by Michael Aupetit

[N,N]=size(dist);
RNG=cell(N,1);
for i=1:N-1
    for j=i+1:N
        neighbi=find(dist(i,:)<=dist(i,j));
        neighbj=find(dist(j,neighbi)<=dist(i,j));
        if length(neighbj)<=2
            RNG=addNeighbor(RNG,i,j);
            RNG=addNeighbor(RNG,j,i);
        end
    end
end

