function [GG]=GabrielGraph(data,indall)

% GabrielGraph - Computes the Gabriel Graph 
%
% Computes the Gabriel Graph of a set of points based on the Euclidean norm
% The Gabriel Graph (GG) [1] is an empty region graphs in the sense defined in [2] and used
% in high-dimensional data analysis in [3].
% GG connects a point p to another point q if the ball with diameter [pq] 
% contains no other point than p and q. The Gabriel Graph is
% a connected Delaunay subgraph whose edges cross the
% shared Voronoi boundary of their endpoints at the center of these
% balls.
% GG is UNDIRECTED
%
% CALL:
% [GG]=GabrielGraph(data,indall)
% 
% INPUT:
% data: NxD N point in D-dimension real space
% indall: NxN indall(i,k) index of the ith nearest neighbor to k amond all
%                         the points (index are relative to the rows in
%                         data) (euclidean distance)
%                         indall can be obtained from 
%                               dataDist=distFast(data,data);
%                               [dall,indall]=sort(dataDist,'ascend');
%
%
% OUPUT:
% GG: N cells GG{i}=[a b...f] set of index in data rows of the Gabriel neighbors of i
%
% [1] Matula, D. W.; Sokal, R. R. (1980), "Properties of Gabriel graphs
% relevant to geographic variation research and clustering of points in the plane", 
% Geogr. Anal. 12 (3): 205–222, doi:10.1111/j.1538-4632.1980.tb00031.x
% [2] Jean Cardinal, Sébastien Collette, Stefan Langerman, Empty region
%     graphs, Computational Geometry, Volume 42, Issue 3, April 2009, 
%     Pages 183-195, ISSN 0925-7721,
%     http://dx.doi.org/10.1016/j.comgeo.2008.09.003.
% [3] M. Aupetit and T. Catz. High-dimensional labeled data analysis with
% topology representing graphs. Neurocomputing, 63:139–169, 2005.
%
% Author   : Michael Aupetit
%            Qatar Computing Research Institute (QCRI)
%            Hamad Bin Khalifa University (HBKU)
%            Doha, Qatar
%            maupetit@qf.org.qa/michael.aupetit@gmail.com
% Last Rev : Feb 08 2016 by Michael Aupetit

[N,D]=size(data);
GG=cell(N,1);

for i=1:N
    neighbi=ones(1,N); % Every points are neighbors of i
    neighbi(i)=0;
    
    indi=indall(:,i); % neighbors of i sorted by increasing Euclidean distance to i
    jj=0;
    while jj<N
        stop=0;
        jj=jj+1;
        j=indi(jj); % from nearest to farthest neighbor
        if j~=i 
            if ~ismember(i,GG{j}) % if i is not already neighbor of j
                % Compute the intermediary point Mi between i and j
                Mi=0.5*(data(i,:)+data(j,:));
                % compute Euclidean distance between Mi and any point j
                dMiNi=sum((Mi-data(j,:)).^2);
                kk=0;
                k=0;
                while stop==0 && k~=j
                    kk=kk+1;
                    k=indi(kk);  % consider neighbors in increasing distance order
                    if k~=j && k~=i
                        dMiNk=sum((Mi-data(k,:)).^2);
                        if dMiNk<=dMiNi
                            % k prevent j from being a neighbor of i
                            neighbi(j)=0;
                            stop=1;
                        end
                    end
                end
            end
        end
    end
    GG{i}=find(neighbi); % store neighbors of i in the cell
end
