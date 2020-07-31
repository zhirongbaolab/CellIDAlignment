function [GONG]=GammaObservableNeighborhoodGraph(data,indall,gamma,G)

% Compute the Gamma-Observable Graph of a set of points based on the
% Euclidean norm
% 
% Input:
% data: NxD N point in D-dimension real space
% indall: NxN indall(i,k) index of the ith nearest neighbor to k amond all
%                         the points (index are relative to the rows in
%                         data) (euclidean distance)
% gamma: [0...1] parameter to tune the size of the neighborhood
%                for gamma=0, GON == Nearest Neighbor Graph, 
%                for gamma=0.5, GON == Gabriel Graph
%                for gamma=1: GON == complete graph
%                gamma_i<gamma_k => GONG_i \subset GONG_k
%                0<=gamma<=0.5: GONG is subgraph of Delaunay graph, and is connected 
%                GONG is a DIRECTED graph
% G: (opt) N cells, G{i} set of GONG neighbors of i for gammaG>gamma
%          this option is to be used if GONG of data has already been computed
%          for some gammaG>gamma. It allows to reduce the search for
%          candidate neighbors of i to the set G{i}=GONG_gammaG{i}
%
% Output:
% GONG: N cells GONG{i}=[a b...f] set of index in data rows of the GON neighbors of i
% 
% Ref:Michaël Aupetit, Pierre Couturier, Pierre Massotte.
% gamma-Observable neighbours for vector quantization  
% Neural Networks, Volume 15, Issues 8-9, October-November 2002, Pages 1017-1027. Elsevier
%
% Requires distFast function.



[N,D]=size(data);
GONG=cell(N,1);

flagAll=0;
if nargin==3
    flagAll=1;
else
    if isempty(G)
        flagAll=1;
    end
end

if flagAll==1
    % consider all points as candidate neighbors in GONG
    for i=1:N
        voisi=ones(1,N); % Every points are neighbors of i
        voisi(i)=0;
        
        indi=indall(:,i); % neighbors of i sorted by increasing Euclidean distance to i
        jj=0;
        while jj<N
            stop=0;
            jj=jj+1;
            j=indi(jj); % from nearest to farthest neighbor
            if j~=i
                % Compute the intermediary point Mi between i and j
                Mi=(1-gamma)*data(i,:)+gamma*data(j,:);
                % compute Euclidean distance between Mi and any point j
                dMiNi=sum((Mi-data(j,:)).^2);
                kk=0;
                k=0;
                while stop==0 & k~=j
                    kk=kk+1;
                    k=indi(kk);  % consider neighbors in increasing distance order
                    if k~=j & k~=i
                        dMiNk=sum((Mi-data(k,:)).^2);
                        if dMiNk<dMiNi
                            % k prevent j from being a neighbor of i
                            voisi(j)=0;
                            stop=1;
                        end
                    end
                end
            end
        end
        GONG{i}=find(voisi); % store neighbors of i in the cell
    end
else    
    % consider only the neighbors in G as candidate neighbors in GONG
    for i=1:N
        if ~isempty(G{i})
            voisi=zeros(1,N);
            voisi(G{i})=1; % Every neighbors of i in G are GONG neighbors of i by default
            
            
            distvoisi=distFast(data(i,:),data(G{i},:));
            [valsortvois,indsortvoisi]=sort(distvoisi);
            nbvoisi=length(indsortvoisi);
            jj=0;
            while jj<nbvoisi
                stop=0;
                jj=jj+1;
                j=G{i}(indsortvoisi(jj));
                if j~=i
                    % Compute the intermediary point Mi between i and j
                    Mi=(1-gamma)*data(i,:)+gamma*data(j,:);
                    % compute Euclidean distance between Mi and any point j
                    dMiNi=sum((Mi-data(j,:)).^2);
                    k=0; % pour rentrer dans la boucle
                    kk=0;
                    
                    while stop==0 & k~=j
                        kk=kk+1;
                        
                        k=indall(kk,i);
                        if k~=j & k~=i
                            dMiNk=sum((Mi-data(k,:)).^2);
                            if dMiNk<dMiNi
                                % k prevent j from being a neighbor of i
                                voisi(j)=0;
                                stop=1;
                            end
                        end
                    end
                end
            end
            GONG{i}=find(voisi); % store neighbors of i in the cell
        end
    end
end
