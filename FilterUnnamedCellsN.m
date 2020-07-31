function  [bestconstraintsperm,bestmatchnames,bestmatchesremoval,changelog,changelogrejected]=FilterUnnamedCellsN(iterations,constraintnames,constraintN,bestpos2names,pos1_cur,bestpos2points,bestconstraints,filterthreshold,unnamednames,...
    confidentpos1,confidentnames2)
[bestmatches,cost3]=matchPointClounds(pos1_cur,bestpos2points,bestconstraints);

%baseline violations
%[v]=countConstraintViolationsNeighborhood(constraintnames,constraintN,bestmatches,bestpos2names,pos1_cur);

if (exist('confidentpos1','var'))
    [v]=countConstraintViolationsNeighborhood(constraintnames,constraintN,...
        [bestmatches,linspace(length(bestpos2names)+1,length(bestpos2names)+length(confidentnames2),length(confidentnames2))],...
        {bestpos2names{:},confidentnames2{:}},[pos1_cur;confidentpos1]);
else
    [v]=countConstraintViolationsNeighborhood(constraintnames,constraintN,bestmatches,bestpos2names,pos1_cur);
end


bestscore=sum(v);
       changelog={};
       c=1;
       changelogrejected={};
       c2=1;
%for each point in the unknown set consider the possiblity it is
%unmachable, remove it from the constraint matrix and the point set and
%compare if this has greatly reduced violations
bestconstraintsperm=bestconstraints; %initialize if never gets updated
bestgoodpos1=ones(length(pos1_cur),1);
bestmatchesremoval=bestmatches;
bestviolationsbrokenout=v;
for i=1:length(pos1_cur)
        goodpos1test=bestgoodpos1;
        goodpos1test(i)=0;
       %{
        %by removal
        pos1_curtest=pos1_cur(logical(goodpos1test),:);%all points that are not the death
        testmatches=bestmatches(logical(goodpos1test));
        %mod constraint
        testconstraints=bestconstraints(logical(goodpos1test),:);
        %} 
        %by inf weights
        pos1_curtest=pos1_cur;
        testmatches=bestmatches;
        testmatches(~goodpos1test)=0;
        testconstraints=bestconstraints;
        testconstraints(~goodpos1test,:)=inf;
        
        %because we're not doing this in the context of constraint
        %optimized matrix can use fast?
       % [violations,matchestest,bestviolation,bestconstrainttest]=iterativeMatchWithConstraintsNoAnswer(constraintnames,constraintAP,constraintLR,constraintDV,testmatches,bestpos2names,pos1_curtest,bestpos2points,iterations,testconstraints);
       [matchestest,cost3]=matchPointClounds(pos1_curtest,bestpos2points,testconstraints);
       bestconstrainttest=testconstraints;%not updating this anymore
       %[matchestest,cost3]=matchPointClounds(pos1_curtest,bestpos2points,bestconstraints);
 %
 %
% visualizeMatches( pos1_curtest, bestpos2points, unnamednames, bestpos2names, matchestest );title('match changes after perm');
%scatter3(pos1_curtest(i,1),pos1_curtest(i,2),pos1_curtest(i,3),'g*');
  
%if(i==435||i==468)
%    'stop'
%end
        %note right now matchtestest is the new to be tested matches while
        %testmatches is the adapted to new points last best match
        %[v]=countConstraintViolationsNeighborhood(constraintnames,constraintN,matchestest,bestpos2names,pos1_curtest);
        
        if (exist('confidentpos1','var'))
    [v]=countConstraintViolationsNeighborhood(constraintnames,constraintN,...
        [matchestest,linspace(length(bestpos2names)+1,length(bestpos2names)+length(confidentnames2),length(confidentnames2))],...
        {bestpos2names{:},confidentnames2{:}},[pos1_curtest;confidentpos1]);
else
    [v]=countConstraintViolationsNeighborhood(constraintnames,constraintN,matchestest,bestpos2names,pos1_curtest);
end

        costtest=sum(v);
        currentremoval=zeros(length(pos1_cur),1);
        currentremoval(i)=1;
        %gainalternative=sum(violationsAP(logical(currentremoval))+violationsDV(logical(currentremoval))+violationsLR(logical(currentremoval))); %alternative cost of prior matching minus directly removed violations
        
       % if(exist('bestviolationsbrokenout','var'))
            %its violations of the original match if not but if have a last
            %best match use that as baseline conflicts of the removed
            %element
            gainalternative=bestviolationsbrokenout(logical(currentremoval));
        %end
        
        %gainalternative=sum(violationsAP(logical(~goodpos1test))+violationsDV(logical(~goodpos1test))+violationsLR(logical(~goodpos1test))); %alternative cost of prior matching minus directly removed violations
        gain=bestscore-costtest;
        percellold=bestscore/(sum(goodpos1test)+1);
        percellnew=costtest/sum(goodpos1test);
        %logic is improvement of >15 and fractional improvement > fraction
        %of points lost
      %  if(gain>50&&gain>gainalternative&&percellnew<percellold)%/bestscore>1/sum(goodpos1test))
             if(gain>filterthreshold)%/bestscore>1/sum(goodpos1test))
      
           bestgoodpos1=goodpos1test; %elimination vector
           bestviolationsbrokenout=v;
           bestmatchesremoval=matchestest;
             %this info is only here in test environment
             if(exist('unnamednames','var'))
            changelog{c}=['Removing cell ',num2str(i),' ',unnamednames{i},' ', num2str(bestscore), ' ',num2str(costtest),' ',num2str(gain),' ',num2str(gainalternative),' ',num2str(percellold),' ', num2str(percellnew)];
             else
                 changelog{c}=['Removing cell ',num2str(i),' ', num2str(bestscore), ' ',num2str(costtest),' ',num2str(gain),' ',num2str(gainalternative),' ',num2str(percellold),' ', num2str(percellnew)];
          
             end
            c=c+1;
            bestscore=costtest;
            bestconstraintsperm=bestconstrainttest;%includes infs to make some illegal
        else
            %implicitly the unnamed data set will retain point
              if(exist('unnamednames','var'))
         
             changelogrejected{c2}=['Reject removing cell ',num2str(i),' ',unnamednames{i},' ', num2str(bestscore), ' ',num2str(costtest),' ',num2str(gain),' ',num2str(gainalternative),' ',num2str(percellold),' ', num2str(percellnew)];
              else
                   changelogrejected{c2}=['Reject removing cell ',num2str(i),' ', num2str(bestscore), ' ',num2str(costtest),' ',num2str(gain),' ',num2str(gainalternative),' ',num2str(percellold),' ', num2str(percellnew)];
           
              end
             c2=c2+1;
        end 
end
%at this point have best matches and used points for them reconstruct full
%size matches


for ij=1:length(bestmatchesremoval)
    if(bestmatchesremoval(ij)~=0)
        matchnames{ij}=bestpos2names{bestmatchesremoval(ij)};
    else
        matchnames{ij}='';
    end
end
bestmatchnames=matchnames;
changelog=changelog';
changelogrejected=changelogrejected';
