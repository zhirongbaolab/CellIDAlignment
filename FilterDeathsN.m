function [ bestconstraints,bestscore,bestpos2points,bestpos2names,bestmatches,bestmatchnames,changelog,changelogrejected] = FilterDeathsN(iterations,constraintnames,constraintN,bestpos2names,pos1_cur,bestpos2points,bestconstraints,deathlist,deaththresh,...
    confidentpos1,confidentnames2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%do deaths last

[bestmatches,cost3]=matchPointClounds(pos1_cur,bestpos2points,bestconstraints);


if (exist('confidentpos1','var'))
 dataDist=distFast( [pos1_cur;confidentpos1], [pos1_cur;confidentpos1]);
    [dall,indall]=sort(dataDist,'ascend');
    grapha=GabrielGraph( [pos1_cur;confidentpos1],indall);
    
   
else    
 dataDist=distFast(pos1_cur,pos1_cur);
    [dall,indall]=sort(dataDist,'ascend');
    grapha=GabrielGraph(pos1_cur,indall);
end
%baseline violations
if (exist('confidentpos1','var'))
    [v]=countConstraintViolationsNeighborhood(constraintnames,constraintN,...
        [bestmatches,linspace(length(bestpos2names)+1,length(bestpos2names)+length(confidentnames2),length(confidentnames2))],...
        {bestpos2names{:},confidentnames2{:}},[pos1_cur;confidentpos1],grapha);
else
    [v]=countConstraintViolationsNeighborhood(constraintnames,constraintN,bestmatches,bestpos2names,pos1_cur,grapha);
end

%[v]=countConstraintViolationsNeighborhood(constraintnames,constraintN,bestmatches,bestpos2names,pos1_cur);
bestscore=sum(v);


changelog={};
c=1;
changelogrejected={};
c2=1;
bestmatchnames=bestpos2names;%default if none are filtered
%cost is one proxy %violations might be a better one but for now use cost
%now iterate over all deaths and if they exist remove them and retest;
for i=1:length(deathlist)
    deathcomparison=strcmpi(deathlist{i},bestpos2names);
    deathindex=find(deathcomparison);
    if ~isempty(deathindex)
        %this death is present in our data set remove it rematch and
        %compare performance
        testpoints=bestpos2points(~deathcomparison,:);%all points that are not the death
        testnames=bestpos2names(~deathcomparison);
        testmatches=bestmatches;
        testmatches(testmatches==deathindex)=0; %this point no longer exists
        testmatches(testmatches>deathindex)=testmatches(testmatches>deathindex)-1; %update matches to match changed point list
        
        %havent checked if have dimensionality right here
        testconstraints=bestconstraints(:,find(~deathcomparison));
        
 %[v]=countConstraintViolationsNeighborhood(constraintnames,constraintN,bestmatches,bestpos2names,pos1_cur);
       
 
 if (exist('confidentpos1','var'))
    [v]=countConstraintViolationsNeighborhood(constraintnames,constraintN,...
        [bestmatches,linspace(length(bestpos2names)+1,length(bestpos2names)+length(confidentnames2),length(confidentnames2))],...
        {bestpos2names{:},confidentnames2{:}},[pos1_cur;confidentpos1],grapha);
else
    [v]=countConstraintViolationsNeighborhood(constraintnames,constraintN,bestmatches,bestpos2names,pos1_cur,grapha);
 end

        if(exist('confidentpos1','var'))
             [violations,matchestest,bestviolation,bestconstrainttest]=iterativeMatchWithConstraintsNoAnswerNPreGraph(constraintnames,constraintN,testmatches,testnames,pos1_cur,testpoints,iterations,testconstraints,grapha,...
                 confidentpos1,confidentnames2);

        else
        [violations,matchestest,bestviolation,bestconstrainttest]=iterativeMatchWithConstraintsNoAnswerNPreGraph(constraintnames,constraintN,testmatches,testnames,pos1_cur,testpoints,iterations,testconstraints,grapha);
        end
        %for some reason line below though essentially the same results in
       %completely different scores for permutations all looking much worse
       
        %[matchestest,cost3]=matchPointClounds(pos1_cur,testpoints,testconstraints);
       %bestconstrainttest=testconstraints;%not updating this anymore
 
        
        %note right now matchtestest is the new to be tested matches while
        %testmatches is the adapted to new points last best match
        % visualizeMatches( pos1_cur, testpoints,  unnamednames, testnames, matchestest );title('answer post 1 death permutation bad');
         
%        [v]=countConstraintViolationsNeighborhood(constraintnames,constraintN,matchestest,testnames,pos1_cur);
        if (exist('confidentpos1','var'))
    [v]=countConstraintViolationsNeighborhood(constraintnames,constraintN,...
        [matchestest,linspace(length(testnames)+1,length(testnames)+length(confidentnames2),length(confidentnames2))],...
        {testnames{:},confidentnames2{:}},[pos1_cur;confidentpos1],grapha);
else
     [v]=countConstraintViolationsNeighborhood(constraintnames,constraintN,matchestest,testnames,pos1_cur,grapha);
end
        costtest=sum(v);
        
        if(costtest+deaththresh<bestscore)
            bestmatches=matchestest;
            bestconstraints=bestconstrainttest;
            bestpos2names=testnames;
            bestpos2points=testpoints;
            
          %  visualizeMatches( pos1_cur, bestpos2points,  unnamednames, bestpos2names, bestmatches );title('answer post 1 death permutation');
            
            matchnames=[];
            for ij=1:length(matchestest)
                if(matchestest(ij)~=0)
                    matchnames{ij}=testnames{matchestest(ij)};
                else
                    matchnames{ij}='';
                end
            end
            bestmatchnames=matchnames;
            changelog{c}=['Removing death ',deathlist{i},' ',num2str(i),' ', num2str(bestscore), ' ',num2str(costtest),' ',num2str(bestscore-costtest)];
            c=c+1;
            bestscore=costtest;
        else
            changelogrejected{c2}=['Rejecting death ',deathlist{i},' ',num2str(i),' ', num2str(bestscore), ' ',num2str(costtest),' ',num2str(bestscore-costtest)];
            c2=c2+1;
            %implicitly the named template will retain this death
        end
    end
end


changelog=changelog';
changelogrejected=changelogrejected';
