function [ASG]=AlphaShapeGraph(data,dist,epsilon)

% AlphaShapeGraph - Computes the Alpha Shape graph
%
% Computes the Alpha Shape Graph (ASG) of a set of points based on the
% Euclidean norm
% The Alpha Shape Graph (ASG) is the 1-skeleton of the Alpha-Shape [] for
% alpha > 0. It connects two points p and q whose Voronoi cells are adjacent 
% and whose distance is lower or equal than 2/alpha (alpha>0).
% ASG is an UNDIRECTED graph
%
% CALL:
% [ASG]=AlphaShapeGraph(data)
% 
% INPUT:
% data: NxD N point in D-dimension real space
% dist: NxN Euclidean distance matrix between each pair of points
%               dist(i,k)=norm(data(i,:)-data(k,:), (dist(k,k)=0)
% epsilon: neighborhood radius (epsilon>0) (epsilon=alpha/2)
%
% OUTPUT:
% ASG: N cells ASG{i}=[a b...f] set of index in data rows of the ASG neighbors of i
% 
% Does not work for data points in dimension greater than 6 (complexity issue)
% NEEDS: delaunayn (Matlab built-in) and uniqueNeighbor functions.
%
% [1] Edelsbrunner, Herbert; Kirkpatrick, David G.; Seidel, Raimund (1983),
% "On the shape of a set of points in the plane", 
% IEEE Transactions on Information Theory 29 (4): 551–559, doi:10.1109/TIT.1983.1056714.
%
% Author   : Michael Aupetit
%            Qatar Computing Research Institute (QCRI)
%            Hamad Bin Khalifa University (HBKU)
%            Doha, Qatar
%            maupetit@qf.org.qa/michael.aupetit@gmail.com
% Last Rev : Feb 08 2016 by Michael Aupetit

[N,D]=size(data);
ASG=cell(N,1);

[DG]=DelaunayGraph(data);

for i=1:N
   neighbi=DG{i};
   for k=1:length(neighbi)
       if dist(i,neighbi(k))<=epsilon
           ASG=addNeighbor(ASG,i,neighbi(k));
           ASG=addNeighbor(ASG,neighbi(k),i);
       end
   end
end


