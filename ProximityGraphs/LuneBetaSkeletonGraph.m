function [LBSG]=LuneBetaSkeletonGraph(data,dist,beta)

% LuneBetaSkeletonGraph - Computes the lune-based beta-skeleton graph
%
% Computes the lune-based beta-skeleton graph of a set of points based on the
% Euclidean norm
% LBSG is one of the empty region graphs [1]
% The Lune-based beta-Skeleton Graph (LBSG) [2]
% is identical to the Circle-Based Beta Skeleton (CBSG) for beta <=0
% For beta<0: in the plane the empty region is the INTERSECTION of two discs of equal
% radius depending on beta both centered on the line (pq), one intercepting p and containing q, 
% the other one intercepting q and containing p. The larger
% this radius, the larger the forbidden empty region enclosing the edge [pq],
% the lower the number of edges in the graph.
% For beta>=0: identically to CBSG, in the plane the empty region is the INTERSECTION of two discs of equal radius
% intercepting points p and q on each side of the line (pq): the larger
% this radius, the smaller the forbidden empty region enclosing the edge [pq],
% the greater the number of edges in the graph.
% 
% Computation complexity is O(N^3)
% Memory complexity is O(N^2)
%
% CALL:
% [LBSG]=LuneBetaSkeletonGraph(data,dist,beta)
%
% INPUT:
% data: NxD N point in D-dimension real space
% dist: NxN Euclidean distance matrix between each pair of point dist(i,k)=norm(data(i,:)-data(k,:)||
% beta: [-1,1] parameter to tune the size of the neighborhood
%                Infinite-strip neighborhood graph : beta=-1;
%                Relative Neighborhood graph : beta=-1/3;
%                Gabriel graph  : beta=0;
%                Complete graph : beta=1; 
%                beta_i<beta_k => LBSG_i \subset LBSG_k
%                LBSG is an UNDIRECTED graph
%                
% OUTPUT:
% LBSG: N cells CBSG{i}=[a b...f] set of index in data rows of the LBSG neighbors of i
% 
%
% NEEDS: addNeighbor; InfiniteStripBandGraph; CircleBetaSkeletonGraph.
% 
% [1] Jean Cardinal, Sébastien Collette, Stefan Langerman, Empty region
% graphs, Computational Geometry, Volume 42, Issue 3, April 2009, 
% Pages 183-195, ISSN 0925-7721, http://dx.doi.org/10.1016/j.comgeo.2008.09.003.
% [2] Kirkpatrick, David G.; Radke, John D. (1985), "A framework for
% computational morphology", Computational Geometry, 
% Machine Intelligence and Pattern Recognition 2, Amsterdam: North-Holland, pp. 217–248.
%
% Author   : Michael Aupetit
%            Qatar Computing Research Institute (QCRI)
%            Hamad Bin Khalifa University (HBKU)
%            Doha, Qatar
%            maupetit@qf.org.qa/michael.aupetit@gmail.com
% Last Rev : Feb 08 2016 by Michael Aupetit

[N,D]=size(data);

if beta>=0
    % run circle-based Beta-Skeleton
    [LBSG]=CircleBetaSkeletonGraph(data,dist,beta);
else
    % un lien i j existe dans LBSG si aucun point n'existe dans
    % l'intersection des disques Ci et Cj passant par i et j et centrés sur la
    % droite (ij) et de rayon ||ij||*beta/2
    LBSG=cell(N,1);
    
    if beta==-1
        % infinite-strip neighborhood
        [LBSG]=InfiniteStripBandGraph(data);
    else % other cases
        alpha=(1-beta)/(1+beta); 
        for i=1:N-1
            for j=i+1:N
                
                % on teste si la région interdite est vide
                novois=0;
                for k=1:N
                    if (k~=i) && (k~=j)
                        ci=data(i,:)+(data(j,:)-data(i,:))*(alpha/2);
                        cj=data(j,:)+(data(i,:)-data(j,:))*(alpha/2);
                        
                        distci=norm(data(k,:)-ci);
                        distcj=norm(data(k,:)-cj);
                        radiusci=norm(data(i,:)-ci);
                        radiuscj=norm(data(j,:)-cj);
                        if distci<=radiusci && distcj<=radiuscj
                            novois=1;
                            break;
                        end
                    end
                end
                if novois==0
                    LBSG=addNeighbor(LBSG,i,j);
                    LBSG=addNeighbor(LBSG,j,i);
                end
            end
        end
    end
end

