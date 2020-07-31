function [ISBG]=InfiniteStripBandGraph(data)

% InfiniteStripBandGraph - Computes the Infinite Strip-Band Graph 
%
% Computes the Infinite Strip-Band Graph of a set of points
% The Infinite Strip-Band Graph (ISBG) [1] is an empty region graphs in the sense defined in [2].
% ISBG is a subgraph of the Gabriel Graph which connects two points p and q
% if no other point s projects orthogonally to the segment [pq].
% ISBG is UNDIRECTED
%
% CALL:
% [ISBG]=InfiniteStripBandGraph(data)
% 
% INPUT:
% data: NxD N point in D-dimension real space
% 
% OUPUT:
% ISBG: N cells ISBG{i}=[a b...f] set of index in data rows of the ISBG neighbors of i
%
% NEEDS: addNeighbor.
%
% [1] Kirkpatrick, David G.; Radke, John D. (1985), "A framework for
% computational morphology", Computational Geometry, 
% Machine Intelligence and Pattern Recognition 2, Amsterdam: North-Holland,
% pp. 217–248.
% [2] Jean Cardinal, Sébastien Collette, Stefan Langerman, Empty region
% graphs, Computational Geometry, Volume 42, Issue 3, April 2009, 
% Pages 183-195, ISSN 0925-7721, http://dx.doi.org/10.1016/j.comgeo.2008.09.003.
%
% Author   : Michael Aupetit
%            Qatar Computing Research Institute (QCRI)
%            Hamad Bin Khalifa University (HBKU)
%            Doha, Qatar
%            maupetit@qf.org.qa/michael.aupetit@gmail.com
% Last Rev : Feb 08 2016 by Michael Aupetit

[N,D]=size(data);
ISBG=cell(N,1);
% infinite-strip neighborhood
for i=1:N-1
    for j=i+1:N
        
        % on teste si la région interdite est vide
        novois=0;
        for k=1:N
            if (k~=i) && (k~=j)
                Vij=data(j,:)-data(i,:);
                Vik=data(k,:)-data(i,:);
                Vjk=data(k,:)-data(j,:);
 
                if Vik*Vij'>0 && Vjk*Vij'<0 
                    novois=1;
                    break;
                end
            end
        end
        if novois==0
            ISBG=addNeighbor(ISBG,i,j);
            ISBG=addNeighbor(ISBG,j,i);
        end
    end
end
