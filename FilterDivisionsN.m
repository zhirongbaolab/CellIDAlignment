function [ bestconstraints,bestscore,bestpos2points,bestpos2names,bestmatches,bestmatchnames,changelog,changelogrejected ] = FilterDivisionsN(iterations,constraintnames,constraintN,bestpos2names,pos1_cur,bestpos2points,bestconstraints,divisionlist,divisionthresh,...
     confidentpos1,confidentnames2)
[bestmatches,cost3]=matchPointClounds(pos1_cur,bestpos2points,bestconstraints);

if (exist('confidentpos1','var'))
     dataDist=distFast([pos1_cur;confidentpos1],[pos1_cur;confidentpos1]);
    [dall,indall]=sort(dataDist,'ascend');
    grapha=GabrielGraph([pos1_cur;confidentpos1],indall);
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

bestscore=sum(v);

c=1;
changelog={};
c2=1;
changelogrejected={};
matchnames=[];
for ij=1:length(bestmatches)
    if(bestmatches(ij)~=0)
        matchnames{ij}=bestpos2names{bestmatches(ij)};
    else
        matchnames{ij}='';
    end
end
bestmatchnames=matchnames;
%
%if(length(bestpos2names)>8) %hack to skip early cases for now
for i=1:size(divisionlist,1)
    divcomparison=strcmpi(divisionlist{i,1},bestpos2names);%is division parent here?
    divisionindex=find(divcomparison);
    if ~isempty(divisionindex)
        
        % parent exists in this time point consider duplicating it to create adivision
        testpoints=bestpos2points(~divcomparison,:);
        testnames={bestpos2names{~divcomparison}};
        %Remove parent from comparison set and replace with two copies of its position with daughter names
        testnames={testnames{:},divisionlist{i,2:3}};
        testpoints=[testpoints;bestpos2points(divcomparison,:);bestpos2points(divcomparison,:)];
        
        %totally wrong
       % testmatches=[bestmatches(~divcomparison),bestmatches(divcomparison),0];%old matches and additional cell is unmatched

       testmatches=[bestmatches];
       testmatches(testmatches==divisionindex)=size(testpoints,1);%set the one matching to updated position the other is un reverse matched by definition
        testmatches(testmatches>divisionindex)=testmatches(testmatches>divisionindex)-1; %update matches to match changed point list
        
        
        %note not sure about row column thing
        testconstraints=bestconstraints(:,find(~divcomparison));
        testconstraints=[testconstraints,bestconstraints(:,divcomparison),zeros(size(bestconstraints(:,divcomparison)))];
        
        
        %[matchestest,costtestm]=matchPointCloundsFast(pos1_cur,testpoints(:,1:3));
        % iterations=20;
        if(exist('confidentpos1','var'))
        [violations,matchestest,bestviolation,bestconstrainttest]=iterativeMatchWithConstraintsNoAnswerNPreGraph(constraintnames,constraintN,testmatches,testnames,pos1_cur,testpoints,iterations,testconstraints,...
             grapha,confidentpos1,confidentnames2);
        else
             [violations,matchestest,bestviolation,bestconstrainttest]=iterativeMatchWithConstraintsNoAnswerNPreGraph(constraintnames,constraintN,testmatches,testnames,pos1_cur,testpoints,iterations,testconstraints,grapha);
       
        end
         

       %  [matchestest,cost3]=matchPointClounds(pos1_cur,testpoints,testconstraints);
       %bestconstrainttest=testconstraints;%not updating this anymore
 
        %[violationsAP,violationsDV,violationsLR]=countConstraintViolationsNamed(constraintnames,constraintAP,constraintLR,constraintDV,matchnames,testnames,pos1_cur);

        if (exist('confidentpos1','var'))
            [v]=countConstraintViolationsNeighborhood(constraintnames,constraintN,...
                [matchestest,linspace(length(testnames)+1,length(testnames)+length(confidentnames2),length(confidentnames2))],...
                {testnames{:},confidentnames2{:}},[pos1_cur;confidentpos1],grapha);
        else
            [v]=countConstraintViolationsNeighborhood(constraintnames,constraintN,matchestest,testnames,pos1_cur,grapha);
        end
        
       % [v]=countConstraintViolationsNeighborhood(constraintnames,constraintN,matchestest,testnames,pos1_cur);
        costtest=sum(v);
        
        if(costtest+divisionthresh<bestscore)
            bestmatches=matchestest;
            bestconstraints=bestconstrainttest;
            bestpos2names=testnames;
            bestpos2points=testpoints;
            matchnames=[];
            for ij=1:length(matchestest)
                if(matchestest(ij)~=0)
                    matchnames{ij}=testnames{matchestest(ij)};
                else
                    matchnames{ij}='';
                end
            end
            bestmatchnames=matchnames;
            
            changelog{c}=['splitting parent ', divisionlist{i,1},' ', num2str(bestscore), ' ',num2str(costtest),' ',num2str(bestscore-costtest)];
           c=c+1;
            bestscore=costtest;
        else
            changelogrejected{c2}=['rejecting splitting parent ', divisionlist{i,1},' ', num2str(bestscore), ' ',num2str(costtest),' ',num2str(bestscore-costtest)];
           c2=c2+1;
            %implicitly the named template will retain this parent rather
            %than divided children
        end
        %}
    else
        
        %here the cells have divided already consider merging them
        % division happened consider undoing it
        divcomparison1=strcmpi(divisionlist{i,2},bestpos2names);%is division
        divcomparison2=strcmpi(divisionlist{i,3},bestpos2names);%is division
        if~isempty(find(divcomparison1))&&~isempty(find(divcomparison2)) %both daughers should be present happens more with death filtering
            %daughter 2 becomes parent
            testpoints=bestpos2points;
            testnames=bestpos2names;
            testnames{divcomparison2}=divisionlist{i,1};
            testpoints(divcomparison2,:)=(testpoints(divcomparison2,:)+testpoints(divcomparison1,:))/2;%replace with mean
            
            %removed daughter 1
            testpoints=testpoints(~divcomparison1,:);
            testnames=testnames(~divcomparison1);
            
            %            [matchestest,costtestm]=matchPointCloundsFast(pos1_cur,testpoints(:,1:3));
            %          iterations=20;
            
            %totally wrong way to update
           % testmatches=bestmatches(~divcomparison1); %arbitrarily keep match and constraints for 2
            testmatches=bestmatches;
            div1index=find(divcomparison1);
            testmatches(testmatches==div1index)=0;
        testmatches(testmatches>div1index)=testmatches(testmatches>div1index)-1; %update matches to match changed point list
            
            
            %havent checked if have dimensionality right here
            testconstraints=bestconstraints(:,find(~divcomparison1));
            if(exist ('confidentpos1','var'))
            
            [violations,matchestest,bestviolation,bestconstrainttest]=iterativeMatchWithConstraintsNoAnswerNPreGraph(constraintnames,constraintN,testmatches,testnames,pos1_cur,testpoints,iterations,testconstraints,...
                 grapha,confidentpos1,confidentnames2);
            else
                  [violations,matchestest,bestviolation,bestconstrainttest]=iterativeMatchWithConstraintsNoAnswerNPreGraph(constraintnames,constraintN,testmatches,testnames,pos1_cur,testpoints,iterations,testconstraints,grapha);
      
            end
            
           if (exist('confidentpos1','var'))
            [v]=countConstraintViolationsNeighborhood(constraintnames,constraintN,...
                [matchestest,linspace(length(testnames)+1,length(testnames)+length(confidentnames2),length(confidentnames2))],...
                {testnames{:},confidentnames2{:}},[pos1_cur;confidentpos1],grapha);
        else
            [v]=countConstraintViolationsNeighborhood(constraintnames,constraintN,matchestest,testnames,pos1_cur,grapha);
        end
        
           % [v]=countConstraintViolationsNeighborhood(constraintnames,constraintN,matchestest,testnames,pos1_cur);
            costtest=sum(v);
            
            if(costtest+divisionthresh<bestscore)
                bestmatches=matchestest;
                bestconstraints=bestconstrainttest;
                
                bestpos2names=testnames;
                bestpos2points=testpoints;
                matchnames=[];
                for ij=1:length(matchestest)
                    if(matchestest(ij)~=0)
                        matchnames{ij}=testnames{matchestest(ij)};
                    else
                        matchnames{ij}='';
                    end
                end
                bestmatchnames=matchnames;
                changelog{c}=['merging daughters to parent ', divisionlist{i,1},' ', num2str(bestscore), ' ',num2str(costtest),' ',num2str(bestscore-costtest)];
          c=c+1;
                bestscore=costtest;
            else
                 changelogrejected{c2}=['rejectingmerging daughters to parent ', divisionlist{i,1},' ', num2str(bestscore), ' ',num2str(costtest),' ',num2str(bestscore-costtest)];
          c2=c2+1;
                %implicitly the named template will retain this parent rather
                %than divided children
            end
            %}
        end
    end
end
%}
changelog=changelog';
changelogrejected=changelogrejected';
end

