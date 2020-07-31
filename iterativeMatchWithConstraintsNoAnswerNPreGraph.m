function [violations,bestanswer,bestviolation,bestconstraint]=iterativeMatchWithConstraintsNoAnswerNPreGraph(constraintnames,constraintN,matches1to2,names_2,pos_1,pos_2,iterations,constraints,grapha,...
    confidentpos1,confidentnames2)



if (exist('confidentpos1','var'))
    [v]=countConstraintViolationsNeighborhood(constraintnames,constraintN,...
        [matches1to2,linspace(length(names_2)+1,length(names_2)+length(confidentnames2),length(confidentnames2))],...
        {names_2{:},confidentnames2{:}},[pos_1;confidentpos1],grapha);
     v(length(pos_1):end)=0; %these are fixed and cant be picked to be changed
else
    [v]=countConstraintViolationsNeighborhood(constraintnames,constraintN,matches1to2,names_2,pos_1,grapha);
end

%[v]=countConstraintViolationsNeighborhood(constraintnames,constraintN,matches1to2,names_2,pos_1);

%correctness1to2=correctmatches(matches1to2, names_1,names_2);
%figure;scatter(correctness1to2,violationsDV+violationsLR+violationsAP);

%pick worst and make it illegal rematch
matches1to2c=matches1to2;

violations=[];
%constraints is optional variable if is provided you continue starting with
%existing constraint matrix logically breaks unless it matches the
%constraints used to create the current answer
if(~exist('constraints','var'))
    constraints12=zeros(length(pos_1),length(pos_2));
    
else
    constraints12=constraints;
end
bestanswer=matches1to2;
bestviolation=sum(v);
bestconstraint=constraints12;
iterationswoimprovement=0;
for i=1:iterations
    
    % length(find(correctmatches(matches1to2c, names_1,names_2)))/length(matches1to2)
    % accuracy(i)=length(find(correctmatches(matches1to2c, names_1,names_2)))/length(matches1to2);
    
    %remove worst violation
    %{
    [mviolationsnum,worstindex]=max(violationsDV+violationsLR+violationsAP);
    constraints12(worstindex,matches1to2c(worstindex))=inf;%make the most violating constraint illegal
    %}
    % if(violations(i)>0&&iterationswoimprovement<400)
    if(sum(v)>0)
        %remove all bad violations (bad is  >= (worst -10)
        maxviol=max(v);
        minviol=min(v);
       
        %sanity insurance on bad range for small number of examples
        %at least removes the worst otherwise the min of 10 or 1/10th the
        %total range of violation scores, ensures you never make everything
        %illegal in one go which breaks the method which assumes you just
        %make the worst thing bad (using 10 instead of 1 is a speedup
        %specific to us)
        %range=min(10,max(1,(maxviol-minviol)/10));%10;%max(violationsDV+violationsLR+violationsAP)/3; %alternate definition of range of bad cases
        range=1;
        %worstones=find((violationsDV+violationsLR+violationsAP)>=max(violationsDV+violationsLR+violationsAP)-10);
        worstones=find((v)>=(max(1,max(v)-range)));
        
        for j=1:length (worstones)
            constraints12(worstones(j),matches1to2c(worstones(j)))=inf;%make the most violating constraint illegal
        end
        %    constraints(i)=length(find(constraints12==inf));
        
        [matches1to2c,cost12c]=matchPointClounds(pos_1,pos_2(:,1:3),constraints12);
        
       if (exist('confidentpos1','var'))
            [v]=countConstraintViolationsNeighborhood(constraintnames,constraintN,...
                [matches1to2c,linspace(length(names_2)+1,length(names_2)+length(confidentnames2),length(confidentnames2))],...
                {names_2{:},confidentnames2{:}},[pos_1;confidentpos1],grapha);
             v(length(pos_1):end)=0; %these are fixed and cant be picked to be changed
        else
           [v]=countConstraintViolationsNeighborhood(constraintnames,constraintN,matches1to2c,names_2,pos_1,grapha);
        end
        
       % [v]=countConstraintViolationsNeighborhood(constraintnames,constraintN,matches1to2c,names_2,pos_1);
    else
        'warning exiting iterations early bc of convergence'
        break;
    end
    %note I reversed these to avoid the kind of unclear initial score boost
    %from just permuting point set, shouldnt change anythign else
    violations(i)=sum(v);
    if (violations(i)<bestviolation)
        bestviolation=violations(i);
        bestanswer=matches1to2c;
        bestconstraint=constraints12;%store constraints used in best answer
        iterationswoimprovement=0;
    else
        iterationswoimprovement=iterationswoimprovement+1;
       %'short circut to always return next iteration for debugging  '
       %     bestanswer=matches1to2c;
        
    end
    
end
end

