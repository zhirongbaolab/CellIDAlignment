function [matches1to2,cost]=matchPointCloundsFast(pos_1,pos2,constraints)

%match two roughly aligned point clouds using distance

costmatrix=distance(pos_1',pos2');

if exist('constraints','var')
    costmatrix=max(costmatrix,constraints);
end

%costmatrix(costmatrix>40)=inf;
 
 
 %extend cost matrix with col for no link costs
%{
 nolink=ones(length(pos_1),length(pos_1)).*inf;
 for i=1:length(nolink)
     nolink(i,i)=40;
 end
 costmatrix=[costmatrix,nolink];
%costmatrix=[costmatrix,ones(length(pos_1),length(pos_1)).*35];
 %}
 
 %tic
 %[matches1to2,cost]=munkres(costmatrix);%matching
%toc
%slower even dammit
%tic
%[matches1to2j,costj]=linear_sum_assignment(costmatrix);
%toc

%note this implementation seems to be order of magnitude faster than

%munkres though puts output in weird format that I think this fixes


[matches1to2j,costj,u,v]=lapjv(costmatrix);%matching
if(size(costmatrix,1)>size(costmatrix,2))
    matches1to2jc=zeros(1,size(costmatrix,1));

    for i=1:length(matches1to2j)
        if(matches1to2j(i)~=0)
            matches1to2jc(matches1to2j(i))=i;
        end
    end
    matches1to2=matches1to2jc;
else
    matches1to2=matches1to2j;
end
%toc

cost=costj;
%}
 
% matches1to2(matches1to2>length(pos2))=0;% matches that are no link
% options are set to zero 
%}
end

