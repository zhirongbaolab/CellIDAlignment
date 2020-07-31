function [KNCG]=KNearestCenterOfGravity(data,K)

% KNearestCenterOfGravity - Computes the K Nearest Center of Gravity graph 
%
% Computes the K-Nearest Center of Gravity graph (KNCG) of a set of points 
% The KNCG [1] allows generate equal size neighborhoods as the KNNG but with
% the goal to put each point the closest possible to the center of gravity
% of its K neighbors. Initially, it connects each point p to its nearest
% neighbor for k=1. Then, iteratively at step k>1, it connects each
% point p to the candidate among the N-k remaining points, for
% which the center of gravity of this candidate together with the k-1 p’s 
% neighbors is the nearest to p. Stop when k==K. The KNCG with K=1 is identical
% to the NNG.
% KNCG is a DIRECTED graph
%
% CALL:
% [KNCG]=KNearestCenterOfGravity(data,K)
% 
% INPUT:
% data: (NxD) N point in D-dimension real space
% K: neighborhood size (K integer >=1) 
%              
% OUTPUT:
% KNCG: N cells KNCG{i}=[a b...f] set of index in data rows of the KNCG
% neighbors of i
% 
% NEEDS: mat2cell
% 
% [1] B. B. Chaudhuri. A new definition of neighborhood of a point in multidimensional
%     space. Pattern Recognition Letters, 17(1):11–17, 1996.
%
% Author   : Michael Aupetit
%            Qatar Computing Research Institute (QCRI)
%            Hamad Bin Khalifa University (HBKU)
%            Doha, Qatar
%            maupetit@qf.org.qa/michael.aupetit@gmail.com
% Last Rev : Feb 08 2016 by Michael Aupetit

if K==1
    dist=distFast(data,data);
    KNCG=NearestNeighborGraph(dist);
else

    [N,D]=size(data);
    KNCG=cell(N,1);
    
    for i=1:N
        candidates=ones(1,N);
        
        % Take v in data
        v=data(i,:);
        % Prevent v from being its own neighbor
        candidates(i)=0;
        
        % Seek the 1-NCG
        dist=inf(1,N);
        for j=find(candidates)
            dist(j)=norm(v-data(j,:));
        end
        [valmin,indmin]=min(dist);
        neighbi(1)=indmin;  % le 1-GPPV = 1-NN
        candidates(indmin)=0;
        dist(indmin)=inf;
        
        % Seek the 2nd neighbor (2-NCG)
        CGi=data(neighbi(1),:);
        for j=find(candidates)
            CG=(CGi+data(j,:))/2;
            dist(j)=norm(v-CG);
        end
        [valmin,indmin]=min(dist);
        neighbi(2)=indmin;
        candidates(indmin)=0;
        dist(indmin)=inf;
        
        % Seek the K-2 next neighbors
        for k=3:K
            CGi=(k-1)*mean(data(neighbi(1:k-1),:));
            for j=find(candidates)
                CG=(CGi+data(j,:))/k;
                dist(j)=norm(v-CG);
            end
            dist(i)=inf;
            [valmin,indmin]=min(dist);
            neighbi(k)=indmin;
            candidates(indmin)=0;
            dist(indmin)=inf;
        end
        
        KNCG{i}=neighbi;
    end
end
