function [SIG]=SphereOfInfluenceGraph(dist)

% SphereOfInfluenceGraph - Computes the Sphere of Influence Graph
%
% Computes the Sphere of Influence Graph of a set of points 
% The Sphere of Influence Graph (SIG) [1] connects
% two points p and q if the distance d(p,q) is lower than the sum
% of the distances to their own nearest neighbors. In other words, if
% the balls centered at points p and q and passing through their respective
% nearest neighbors intersect, then p and q are
% neighbors. Compared to the EBG, the radius of the ball
% is “automatically” adapted to the local density of the points.
% SIG is an UNDIRECTED graph
%
% CALL:
% [SIG]=SphereOfInfluenceGraph(dist)
%
% INPUT:
% dist: NxN  distance matrix between each pair of points
%              
% OUTPUT:
% SIG: N cells SIG{i}=[a b...f] set of index in data rows of the SIG
% neighbors of i
% 
% [1]  G. T. Toussaint, “Pattern recognition and geometrical complexity,” 
% Proceedings Fifth International Conference
% on Pattern Recognition, Miami Beach, December 1980, pp. 1324-1347.
%
% Author   : Michael Aupetit
%            Qatar Computing Research Institute (QCRI)
%            Hamad Bin Khalifa University (HBKU)
%            Doha, Qatar
%            maupetit@qf.org.qa/michael.aupetit@gmail.com
% Last Rev : Feb 08 2016 by Michael Aupetit

[N,N]=size(dist);

% recherche du ppv de chaque point
dist(eye(N)==1)=inf;
[radiusMin,indMin]=min(dist);

[X,Y] = meshgrid(radiusMin);
sumRadius=X+Y; 

diffRadius=sumRadius-dist;

SIG=cell(N,1);
for i=1:N
    SIG{i}=find(diffRadius(i,:)>=0);
end





