function [CBSG]=CircleBetaSkeletonGraph(data,dist,beta)

% CircleBetaSkeletonGraph - Computes the circle-based beta-skeleton graph
%
% Computes the circle-based beta-skeleton graph of a set of points based on the
% Euclidean norm
% CBSG is one of the empty region graphs [1]
% The Circle-based beta-Skeleton Graph (CBSG) [2]
% connects two points p and q if for any other point s the angle formed
% by the lines joining s to p and q is lower than a threshold angle
% pi*(1+beta)/2. 
% For beta<=0: in the plane the empty region is the UNION of two discs of equal radius
% intercepting points p and q on each side of the line (pq): the larger
% this radius, the larger the forbidden empty region enclosing the edge [pq],
% the lower the number of edges in the graph.
% For beta>=0: in the plane the empty region is the INTERSECTION of two discs of equal radius
% intercepting points p and q on each side of the line (pq): the larger
% this radius, the smaller the forbidden empty region enclosing the edge [pq],
% the greater the number of edges in the graph.
% The parameter beta lies in the range of  [-1;1]: If -1,
% CBSG is the empty graph; if 0, CBSG is identical to the Gabriel
% graph; and if 1, then CBSG is the complete graph.
% 
% Computation complexity is O(N^3)
% Memory complexity is O(N^2)
% 
% CALL:
% [CBSG]=CircleBetaSkeletonGraph(data,dist,beta)
%
% INPUT:
% data: NxD N point in D-dimension real space
% dist: NxN Euclidean distance matrix between each pair of point dist(i,k)=norm(data(i,:)-data(k,:)||
% beta: [-1,1] parameter to tune the size of the neighborhood
%                empty graph    : beta=-1;
%                Gabriel graph  : beta=0; 
%                complete graph : beta=1;  
%                beta_i<beta_k => CBSG_i \subset CBSG_k
%                CBSG is an UNDIRECTED graph
%                
% OUTPUT:
% CBSG: N cells CBSG{i}=[a b...f] set of index in data rows of the CBSG neighbors of i
% 
%
% NEEDS: addNeighbor.
% 
% [1] Jean Cardinal, Sébastien Collette, Stefan Langerman, Empty region
% graphs, Computational Geometry, Volume 42, Issue 3, April 2009, 
% Pages 183-195, ISSN 0925-7721, http://dx.doi.org/10.1016/j.comgeo.2008.09.003.
% [2] R. Veltkamp. The gamma-neighborhood graph. Computational Geometry, 1:227–246, 1991.
%
% Author   : Michael Aupetit
%            Qatar Computing Research Institute (QCRI)
%            Hamad Bin Khalifa University (HBKU)
%            Doha, Qatar
%            maupetit@qf.org.qa/michael.aupetit@gmail.com
% Last Rev : Feb 08 2016 by Michael Aupetit

[N,D]=size(data);

angleMax=pi*0.5*(1+beta);

CBSG=cell(N,1);

if beta==-1 % empty graph
    
elseif beta==1 % complete graph
    [CBSG]=CompleteGraph(N);
else % other cases
    for i=1:N-1
        for j=i+1:N
            
            % test if forbidden region is empty
            novois=0;
            for k=1:N
                if (k~=i) && (k~=j)
                    Vki=data(i,:)-data(k,:);
                    Vkj=data(j,:)-data(k,:);
                    angleK=acos((Vki*Vkj')/(dist(k,i)*dist(k,j)));
                    if angleK>angleMax
                        novois=1;
                        break;
                    end
                end
            end
            if novois==0
                CBSG=addNeighbor(CBSG,i,j);
                CBSG=addNeighbor(CBSG,j,i);
            end
        end
    end
end


