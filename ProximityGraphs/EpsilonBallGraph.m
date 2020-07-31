function [EBG]=EpsilonBallGraph(dist,epsilon)

% EpsilonBallGraph - Computes the Epsilon-Ball graph 
%
% Computes the Epsilon-Ball graphof a set of points 
% The Epsilon-Ball Graph (EBG) [1] connects points p and q
% if their distance d(p,q)<=epsilon (epsilon > 0). Any point in the
% epsilon-ball centered at p is neighbor of p.
% EBG is an UNDIRECTED graph
%
% CALL:
% [EBG]=EpsilonBallGraph(dist,epsilon)
%
% INPUT:
% dist: NxN  distance matrix between each pair of points
% epsilon: neighborhood radius (epsilon>0) epsilon=0: EBG is the empty graph
%              
% OUTPUT:
% EBG: N cells EBG{i}=[a b...f] set of index in data rows of the EBG
% neighbors of i
% 
% NEEDS: mat2cell.
% 
% [1] R. Veltkamp. The gamma-neighborhood graph. Computational Geometry, 1:227–246, 1991.
%
% Author   : Michael Aupetit
%            Qatar Computing Research Institute (QCRI)
%            Hamad Bin Khalifa University (HBKU)
%            Doha, Qatar
%            maupetit@qf.org.qa/michael.aupetit@gmail.com
% Last Rev : Feb 08 2016 by Michael Aupetit



[N,N]=size(dist);

epsmat=zeros(N,N);
if epsilon==0
    for k=1:N
        EBG{k}=[];
    end
else
    epsmat(dist<=epsilon)=1;
    EBG=mat2cell(epsmat-eye(N));  
end

