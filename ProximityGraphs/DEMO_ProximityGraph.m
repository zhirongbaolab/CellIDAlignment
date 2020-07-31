function []=DEMO_ProximityGraph(graphName,ndata,parameter)


% DEMO_ProximityGraph - Runs and displays proximity graphs
%
% CALL:
% testProximityGraph(graphName,ndata,parameter)
%
% INPUT:
% graphName: name of the proximity graph to run: 
%      'GONG' : Gamma Observable Neighborhood Graph
%      'CBSG' : Circle-Based Beta Skeleton Graph
%      'LBSG' : Lune-Based Beta Skeleton Graph
%      'RNG'  : Relative Neighborhood Graph
%      'GG'   : Gabriel Graph
%      'DG'   : Delaunay Graph
%      'NNG'  : Nearest Neighbor Graph
%      'KNNG' : K-Nearest Neighbor Graph
%      'KNCG' : K-Nearest Cented of Gravity Graph
%      'EBG'  : Epsilon-Ball Graph
%      'MSTG' : Minimum Spanning Tree Graph
%      'SIG'  : Sphere of Influence Graph
%      'ASG'  : Alpha-Shape Graph
%      'ISBG' : Infinite Strip-Band Graph
%      (SYMMETRICAL and MUTUAL tested with KNNG)
%      (REVERSE tested with GG)
% ndata: number of data points (test with ndata=10 or 100)
% parameter (opt): parameter value for parametric proximity graphs:
%                  GONG, CBSG, LBSG, KNNG, KNCG, EBG, ASG
%
% OUTPUT:
% figures
%
% 
% NEEDS:
% AlphaShapeGraph                   
% CircleBetaSkeletonGraph           
% CompleteGraph                     
% DelaunayGraph                     
% EpsilonBallGraph                 
% GabrielGraph                    
% GammaObservableNeighborhoodGraph  
% InfiniteStripBandGraph         
% KNearestCenterOfGravity          
% KNearestNeighborGraph           
% LuneBetaSkeletonGraph             
% MinimumSpanningTree              
% MutualGraph                  
% NearestNeighborGraph             
% RelativeNeighborhoodGraph       
% ReverseGraph              
% SphereOfInfluenceGraph
% SymmetricalGraph
% distFast
% mat2cell
% cell2mat
% uniqueNeighbor
% addNeighbor
% 
% 
% EXAMPLES:
% testProximityGraph('GONG',100,0.8)
% testProximityGraph('GG',10)
% testProximityGraph('KNNG',10,5)
% testProximityGraph('KNCG',10,5)
%
% Author   : Michael Aupetit
%            Qatar Computing Research Institute (QCRI)
%            Hamad Bin Khalifa University (HBKU)
%            Doha, Qatar
%            maupetit@qf.org.qa/michael.aupetit@gmail.com
% Last Rev : Feb 08 2016 by Michael Aupetit


if nargin==2
    parameter=[];
end

G1=[];
parameter1=[];
G2=[];
parameter2=[];
    

% GENERATING DATA
rand('state',1);  % to get always the same set of random points
randn('state',1); % to get always the same set of random points
ndata1=ceil(0.5*ndata);
ndata2=ndata-ndata1;
data=[randn(ndata1,2); rand(ndata2,2)];

% COMPUTING DISTANCES
fprintf('COMPUTE DISTANCES...\n')
dist=distFast(data,data);
fprintf('SORT DISTANCES...\n')
[dall,indall]=sort(dist,'ascend');

switch graphName
    
% -------------------------------- GONG ------------------------------
    case 'GONG'
        
        % Computing the proximity graph
        fprintf('COMPUTE GRAPH...\n')
        tic
        [G]=GammaObservableNeighborhoodGraph(data,indall,parameter);
        fprintf('Time to compute GONG_0.5 from scratch: %f s...\n',toc);
        tic
        % get GONG_0.3 from scratch
        graphName1='From scratch: GONG';
        parameter1=parameter/2;
        [G1]=GammaObservableNeighborhoodGraph(data,indall,parameter/2);
        fprintf('Time to compute GONG_0.3 from scratch: %f s...\n',toc);
        tic
        % get GONG_0.3 knowing GONG_0.5
        graphName2='From GONG_0.5 : GONG';
        parameter2=parameter/2;
        [G2]=GammaObservableNeighborhoodGraph(data,indall,parameter/2,G);
        fprintf('Time to compute GONG_0.3 from GONG_0.5: %f s...\n',toc);
                
        
% -------------------------------- CBSG ------------------------------        
    case 'CBSG'
        [G]=CircleBetaSkeletonGraph(data,dist,parameter);
        
% -------------------------------- LBSG ------------------------------        
    case 'LBSG'
        [G]=LuneBetaSkeletonGraph(data,dist,parameter);
        
% -------------------------------- RNG ------------------------------        
    case 'RNG'
         [G]=RelativeNeighborhoodGraph(dist);

% -------------------------------- GG ------------------------------        
    case 'GG'
         [G]=GabrielGraph(data,indall);
         if ndata<20
             graphName1='REVERSE GG';
             parameter1=[];
             [G1]=ReverseGraph(G);
         end
% -------------------------------- DG ------------------------------
    case 'DG'
         [G]=DelaunayGraph(data);
         
% -------------------------------- NNG ------------------------------
    case 'NNG'
         [G]=NearestNeighborGraph(dist);
         
% -------------------------------- KNNG ------------------------------
    case 'KNNG'
         [G]=KNearestNeighborGraph(dist,parameter);
         graphName1='SYMMETRICAL KNNG';
         parameter1=[];
         [G1]=SymmetricalGraph(G);
         graphName2='MUTUAL KNNG';
         parameter2=[];
         [G2]=MutualGraph(G);
         
% -------------------------------- KNCG ------------------------------
    case 'KNCG'
         [G]=KNearestCenterOfGravity(data,parameter);
         
% -------------------------------- EBG ------------------------------
    case 'EBG'
         [G]=EpsilonBallGraph(dist,parameter);
         
% -------------------------------- MSTG ------------------------------
    case 'MSTG'
         [G]=MinimumSpanningTree(dist);

% -------------------------------- SIG ------------------------------
    case 'SIG'
         [G]=SphereOfInfluenceGraph(dist);
         
% -------------------------------- ASG ------------------------------
    case 'ASG'
         [G]=AlphaShapeGraph(data,dist,parameter);    
         
% -------------------------------- ISBG ------------------------------
    case 'ISBG'
         [G]=InfiniteStripBandGraph(data);    
end



% DISPLAY
figure(10)
plot(data(:,1),data(:,2),'k.');
hold on
for t=1:length(G)
    neighbors=G{t};
    numn=length(neighbors);
    for nt=1:numn
        midpoint=0.5*(data(t,:)+data(neighbors(nt),:));
        plot([data(t,1) midpoint(1)],[data(t,2) midpoint(2)],'r-');
        hold on
    end
end
plot(data(:,1),data(:,2),'k.');
axis equal
if isempty(parameter)
   title(sprintf('%s',graphName));
else
   title(sprintf('%s %f',graphName,parameter));
end
hold off

% OPTIONNAL DISPLAY 1
if ~isempty(G1)
    figure(11)
    plot(data(:,1),data(:,2),'k.');
    hold on
    for t=1:length(G1)
        neighbors=G1{t};
        numn=length(neighbors);
        for nt=1:numn
            midpoint=0.5*(data(t,:)+data(neighbors(nt),:));
            plot([data(t,1) midpoint(1)],[data(t,2) midpoint(2)],'b-');
            hold on
        end
    end
    plot(data(:,1),data(:,2),'k.');
    axis equal
    if isempty(parameter1)
        title(sprintf('%s',graphName1));
    else
        title(sprintf('%s %f',graphName1,parameter1));
    end
    hold off
end

% OPTIONNAL DISPLAY 2
if ~isempty(G2)
    figure(12)
    plot(data(:,1),data(:,2),'k.');
    hold on
    for t=1:length(G2)
        neighbors=G2{t};
        numn=length(neighbors);
        for nt=1:numn
            midpoint=0.5*(data(t,:)+data(neighbors(nt),:));
            plot([data(t,1) midpoint(1)],[data(t,2) midpoint(2)],'m-');
            hold on
        end
    end
    plot(data(:,1),data(:,2),'k.');
    axis equal
    if isempty(parameter2)
        title(sprintf('%s',graphName2));
    else
        title(sprintf('%s %f',graphName2,parameter2));
    end
    hold off
end
