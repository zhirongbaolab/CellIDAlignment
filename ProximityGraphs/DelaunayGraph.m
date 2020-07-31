function [DG]=DelaunayGraph(data)

% DelaunayGraph - Computes the Delaunay graph
%
% Computes the Delaunay Graph (DG) of a set of points based on the
% Euclidean norm
% DG is one of the empty region graphs [1]
% DG connects two points p and q whose Voronoi cells are adjacent.
% Points are supposed to be in general position.
% DG is an UNDIRECTED graph
%
% CALL:
% [DG]=DelaunayGraph(data)
% 
% INPUT:
% data: NxD N point in D-dimension real space
%                
% OUTPUT:
% DG: N cells DG{i}=[a b...f] set of index in data rows of the DG neighbors of i
% 
% Does not work for data points in dimension greater than 6 (complexity issue)
% NEEDS: delaunayn (Matlab built-in) and uniqueNeighbor functions.
%
% [1] Jean Cardinal, Sébastien Collette, Stefan Langerman, Empty region
%     graphs, Computational Geometry, Volume 42, Issue 3, April 2009, 
%     Pages 183-195, ISSN 0925-7721,
%     http://dx.doi.org/10.1016/j.comgeo.2008.09.003.
% [2] R. Veltkamp. The gamma-neighborhood graph. Computational Geometry,1:227–246, 1991.
%
% Author   : Michael Aupetit
%            Qatar Computing Research Institute (QCRI)
%            Hamad Bin Khalifa University (HBKU)
%            Doha, Qatar
%            maupetit@qf.org.qa/michael.aupetit@gmail.com
% Last Rev : Feb 08 2016 by Michael Aupetit

[N,D]=size(data);

% compute the Delaunay triangulation
tri=delaunayn(data);

% extract the edges of the triangulation
DG=cell(length(unique(tri)),1);

[lt,ct]=size(tri);
for i=1:lt
    for j=1:ct
        indx=setdiff(1:ct,j);
        DG{tri(i,j)}=[DG{tri(i,j)} tri(i,indx)];
    end
end

[DG]=uniqueNeighbor(DG);


