function [MSTG]=MinimumSpanningTree(dist)

% MinimumSpanningTree - Computes the Euclidean Minimum Spanning Tree
%
% Computes the Minimum Spanning Tree of a set of points 
% The Minimum Spanning Tree (MSTG) [1] is a subgraph of the RNG 
% and is the unique connected subgraph of the Delaunay graph forming 
% a tree and whose total edge length is minimum. It uses Prim algorithm.
% MSTG is an UNDIRECTED graph
%
% CALL:
% [MSTG,edges]=MinimumSpanningTree(dist)
%
% INPUT:
% dist: NxN  distance matrix between each pair of points
%              
% OUTPUT:
% MSTG: N cells MSTG{i}=[a b...f] set of index in data rows of the MSTG
% neighbors of i
% 
% NEEDS: addNeighbor.
% 
% [1]  Prim, R. C. (November 1957), "Shortest connection networks And some
% generalizations", Bell System Technical Journal 36 (6): 1389–1401, 
% doi:10.1002/j.1538-7305.1957.tb01515.x.
%
% Author   : Michael Aupetit
%            Qatar Computing Research Institute (QCRI)
%            Hamad Bin Khalifa University (HBKU)
%            Doha, Qatar
%            maupetit@qf.org.qa/michael.aupetit@gmail.com
% Last Rev : Feb 08 2016 by Michael Aupetit

[N,N]=size(dist);
MSTG=cell(N,1);

MST=1;  % explored list: index of explored points already in the MST
notMST=2:N; % candidate list: index of points still to explore

LnotMST=length(notMST);
LMST=length(MST);

while LMST<N
    
    indppvnotMST=zeros(1,LMST);  % index in notMST of the Nearest Neirghbor of each point in MST, among the notMST points
    indminMST=zeros(1,LMST);  % index in dist of Nearest Neighbor of each point in MST, among the notMST points
    valminMST=zeros(1,LMST);  % distance to the Nearest Neighbor of each point in MST, among the notMST points
    for i=1:LMST
        [valmin,indmin]=min(dist(MST(i),notMST));
        indppvnotMST(i)=indmin;  
        indminMST(i)=notMST(indmin); 
        valminMST(i)=valmin;
    end
    
    [val,ind]=min(valminMST);
    
    % update explored list
    MST=[MST indminMST(ind)];
    LMST=LMST+1;
    
    % add neighbors
    MSTG=addNeighbor(MSTG,MST(ind),indminMST(ind));
    MSTG=addNeighbor(MSTG,indminMST(ind),MST(ind));
    
    % update candidate list
    if LnotMST>1
        if indppvnotMST(ind)<LnotMST
            notMST(indppvnotMST(ind))=notMST(end);
        end
        notMST=notMST(1:end-1);
    else
        notMST=[];
    end
    LnotMST=LnotMST-1;
    
end

