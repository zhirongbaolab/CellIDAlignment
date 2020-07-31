function [matches1to2names,changelog,bestviolationd,accuracyinfo]=matchPointCloundsWithPerturbationUpdatingConstraintsN(pos1_cur,pos2_cur,names2_cur,divisionlist,deathlist,ICP,PERMUTE,constraintnames,constraintN,iterations,initialiterations,finaliterations,divisionthresh,deaththresh,filterthresh,...
    unnamednames)
%confident pos1,2 and names 2 are confident cases that should not be
%matched anew, but which need to be inserted before scoring of matches is
%done against neighbor model.

accuracyinfo=[];
%given a set of points and an embryo match against temporal match
%then permute each division in local vicinity and each dead cell present

%really means do prealignment
if ICP
    
       %first do ICP to refine initial affine alignment %always
    [TR,TT,pos2transto1ICP] = icp(pos1_cur,pos2_cur(:,1:3));
    pos2transto1ICP=pos2transto1ICP';
    pos2_cur=pos2transto1ICP;
 
   if (length(pos1_cur)>25) %if more than minimal points
        %being close to optimal then do non-rigid warp
        unknown=pos1_cur;
        known=pos2_cur;
        config = gmmreg_load_config('D:\gmmreg_cpd\gmmreg_cpd/face.ini');
        config.motion = 'grbf';
       % config.percellsigma=percellvariance;
     %  config.init_sigma=.375;%was .5 in original param file
        %config.motion='tps';
        config.init_param = zeros(size(unknown,1),3);
        config.ctrl_pts=unknown;
        config.model=config.ctrl_pts;
        config.scene=known;
        [fp,fm] = gmmreg_cpd(config);
       % [fp,fm] = gmmreg_cpd_withmultisigmas(config);
   
        pos1_cur=fm;
    end
end

%{
[TR,TT,pos2transto1ICP] = icp(pos1_cur,pos2_cur(:,1:3));
pos2transto1ICP=pos2transto1ICP';
if ICP
    pos2_cur=pos2transto1ICP;
end
  %}
% pretransform match set qualityto that matches (constraints,
% cost?
%repeat
%best match is answer
bestscore=inf;
bestpos2points=[];
bestpos2names=names2_cur;

%initial match
[matches1to2,cost12]=matchPointClounds(pos1_cur,pos2_cur(:,1:3));
 
initialconstraints=zeros(length(pos1_cur),length(pos2_cur));
 
   % visualizeMatches( pos1_cur, pos2_cur,  unnamednames, names2_cur, matches1to2 );
   %  visualizeMatches( pos1_cur, fm,  unnamednames,unnamednames, zeros(1,length(fm)) );
    %[ correct,correctreverse,realmatch,realmatchreverse ] = computeCorrectnessAndRealMatches( unnamednames,names2_cur,matches1to2);
% saved P from alignment
%test code to see if ihas useful indications on unmachable
%{
%matchless
%p values
matchlessrows=P(realmatch==0,:);
matchedrowscorrect=P(realmatch~=0&realmatch==matches1to2,:);
matchedrowswrong=P(realmatch~=0&realmatch~=matches1to2,:);
bins=linspace(0,1,100);
figure
plot(histc(max(matchlessrows'),bins));
hold on
plot(histc(max(matchedrowscorrect'),bins));
plot(histc(max(matchedrowswrong'),bins));
title ('cases final warping probability matrix');

Psecond=P';
[V,ind]=max(Psecond);
for i=1:length(ind)
Psecond(ind(i),i)=-inf;
end
[V2,ind2]=max(Psecond);
bins=linspace(0,.000125,100);
figure
plot(histc(V2(realmatch==0)./max(matchlessrows'),bins));
hold on
plot(histc(V2(realmatch~=0&realmatch==matches1to2)./max(matchedrowscorrect'),bins));
plot(histc(V2(realmatch~=0&realmatch~=matches1to2)./max(matchedrowswrong'),bins));
title ('cases ratio of best and second final warping probability matrix');




costmatrix=distance(pos1_cur',pos2_cur');
matchlessrows=costmatrix(realmatch==0,:);
matchedrowscorrect=costmatrix(realmatch~=0&realmatch==matches1to2,:);
matchedrowswrong=costmatrix(realmatch~=0&realmatch~=matches1to2,:);
bins=linspace(0,50,50);
figure
plot(histc(min(matchlessrows'),bins));
hold on
plot(histc(min(matchedrowscorrect'),bins));
plot(histc(min(matchedrowswrong'),bins));
title ('cases warped distance matrix');

Psecond=costmatrix';
[V,ind]=std(Psecond);
for i=1:length(ind)
Psecond(ind(i),i)=inf;
end
[V2,ind2]=min(Psecond);
bins=linspace(0,1,100);
figure
plot(histc(V2(realmatch==0)./max(matchlessrows'),bins));
hold on
plot(histc(V2(realmatch~=0&realmatch==matches1to2)./max(matchedrowscorrect'),bins));
plot(histc(V2(realmatch~=0&realmatch~=matches1to2)./max(matchedrowswrong'),bins));
title ('cases ratio of best and second distance final warping');

Psecond=costmatrix';
V=std(Psecond);
bins=linspace(0,100,100);
figure
plot(histc(V2(realmatch==0),bins));
hold on
plot(histc(V2(realmatch~=0&realmatch==matches1to2),bins));
plot(histc(V2(realmatch~=0&realmatch~=matches1to2),bins));
title ('std dev distance');




%}
if(exist('unnamednames','var'))
accuracyinfo(1)=length(find(correctmatches(matches1to2, unnamednames,names2_cur)))/length(matches1to2);
 end
 
%initialiterations=1400;%initial different from number of iterations with each test
%[violations,matches1to2,bestviolation,bestconstraints]=iterativeMatchWithConstraintsNoAnswer(constraintnames,constraintAP,constraintLR,constraintDV,matches1to2,bestpos2names,pos1_cur,pos2_cur,initialiterations,initialconstraints);
% if(exist('unnamednames','var'))
%accuracyinfo(2)=length(find(correctmatches(matches1to2, unnamednames,names2_cur)))/length(matches1to2);
% end
 bestconstraints=initialconstraints;
 
matchnames=[];
for i=1:length(matches1to2)
    if(matches1to2(i)~=0)
        matchnames{i}=names2_cur{matches1to2(i)};
    else
        matchnames{i}='';
    end
end
%[violationsAP,violationsDV,violationsLR]=countConstraintViolationsNamed(constraintnames,constraintAP,constraintLR,constraintDV,matchnames,names2_cur,pos1_cur);


%[violationsAP,violationsDV,violationsLR]=countConstraintViolations(constraintnames,constraintAP,constraintLR,constraintDV,matches1to2,names2_cur,pos1_cur);
 %graph for this timepoint
 %precompute gabriel graph once and use it  for whole run 
    dataDist=distFast(pos1_cur,pos1_cur);
    [dall,indall]=sort(dataDist,'ascend');
    grapha=GabrielGraph(pos1_cur,indall);

[violations]=countConstraintViolationsNeighborhood(constraintnames,constraintN,matches1to2,names2_cur,pos1_cur,grapha);

bestmatches=matches1to2;
bestmatchnames=matchnames;
bestpos2points=pos2_cur;
%bestscore=sum(violationsAP+violationsDV+violationsLR)
bestscore=sum(violations);

%baseline performance
%length(find(correctmatchesnames(bestmatchnames, unnamednames)))/length(bestmatchnames)
%length(find(correctmatches(bestmatches, unnamednames,names2_cur)))/length(bestmatches)
%visualizeMatches( pos1_cur, pos2_cur,  unnamednames, names2_cur, bestmatches );title('answer post optimization');

%{
%filter divisions
[ bestconstraintsdiv,bestscorediv,bestpos2pointsdiv,bestpos2namesdiv,bestmatchesdiv,bestmatchnamesdiv,changelogdiv,changelogdivrejected ] = FilterDivisions(iterations,bestscore,constraintnames,constraintAP,constraintLR,constraintDV,bestmatches,bestpos2names,pos1_cur,bestpos2points,bestconstraints,divisionlist);
 bestmatches=bestmatchesdiv;
 bestpos2names=bestpos2namesdiv;
 bestpos2points=bestpos2pointsdiv;
 bestconstraints=bestconstraintsdiv;
 bestscore=bestscorediv;
 

%filter deaths
 [ bestconstraintsdeath,bestscoredeath,bestpos2pointsdeath,bestpos2namesdeath,bestmatchesdeath,bestmatchnamesdeath,changelogdeath,changelogdeathrejected ] = FilterDeaths(iterations,bestscore,constraintnames,constraintAP,constraintLR,constraintDV,bestmatches,bestpos2names,pos1_cur,bestpos2points,bestconstraints,deathlist,unnamednames);
%with death
 %length(find(correctmatchesnames(bestmatchnamesdeath, unnamednames)))/length(bestmatchnamesdeath)
% length(find(correctmatches(bestmatchesdeath, unnamednames,bestpos2namesdeath)))/length(bestmatchesdeath)
%visualizeMatches( pos1_cur, bestpos2pointsdeath,  unnamednames, bestpos2namesdeath, bestmatchesdeath );title('answer post death permutation');
%
%finaliterations=800;
[violationsd,matches1to2d,bestviolationd,bestconstraintsd]=iterativeMatchWithConstraintsNoAnswer(constraintnames,constraintAP,constraintLR,constraintDV,bestmatchesdeath,bestpos2namesdeath,pos1_cur,bestpos2pointsdeath,finaliterations,bestconstraintsdeath);
%visualizeMatches( pos1_cur, bestpos2pointsdeath,  unnamednames, bestpos2namesdeath, matches1to2d );title('answer post death permutation and reoptimization 800 rounds');
% length(find(correctmatches(matches1to2d, unnamednames,bestpos2namesdeath)))/length(matches1to2d)

 %make so death results used in next step
 %bestmatches=bestmatchesdeath;
 bestmatches=matches1to2d;%use reoptimized 
 bestpos2names=bestpos2namesdeath;
 bestpos2points=bestpos2pointsdeath;
 bestconstraints=bestconstraintsd;%use reoptimized
 bestscore=bestviolationd;%use reoptimized
 
%do reverse filtering
%unnamednames=names2_cur;
[bestconstraintsperm,bestmatchnamesperm,bestmatchesperm,changelogperm,changelogpermrejected]=FilterUnnamedCells(iterations,bestscore,constraintnames,constraintAP,constraintLR,constraintDV,bestmatches,bestpos2names,pos1_cur,bestpos2points,bestconstraints, unnamednames);
%reoptimize post reverse filtering
[violationsd,matches1to2d,bestviolationd,bestconstraintsd]=iterativeMatchWithConstraintsNoAnswer(constraintnames,constraintAP,constraintLR,constraintDV,bestmatchesperm,bestpos2namesdeath,pos1_cur,bestpos2pointsdeath,finaliterations,bestconstraintsperm);
length(find(correctmatches(matches1to2d, unnamednames,bestpos2namesdeath)))/length(matches1to2d)
visualizeMatches( pos1_cur, bestpos2pointsdeath,  unnamednames, bestpos2namesdeath, matches1to2d );title('answer post death and drop permutation and re-reoptimization 800 rounds');

%length(find(correctmatchesnames(bestmatchnamesperm, unnamednames)))/length(bestmatchnamesperm)


%test post optimization %dosent seem to do anything
%finaliterations=800;%initial different from number of iterations with each test
%[violationsp,matches1to2p,bestviolationp,bestconstraintsp]=iterativeMatchWithConstraintsNoAnswer(constraintnames,constraintAP,constraintLR,constraintDV,bestmatches,bestpos2names,pos1_cur,bestpos2points,finaliterations,bestconstraints);

%length(find(correctmatchesnames(bestmatchnames, names2_cur)))/length(bestmatchnames)
%length(find(correctmatches(matches1to2p, names2_cur,bestpos2names)))/length(matches1to2p)

matches1to2names=matches1to2d;

%}
if (PERMUTE)
%'initial match done'
%alternate path where do it div,death,perm,reoptimize
[ bestconstraintsdiv,bestscorediv,bestpos2pointsdiv,bestpos2namesdiv,bestmatchesdiv,bestmatchnamesdiv,changelogdiv,changelogdivrejected ] =...
    FilterDivisionsN(iterations,constraintnames,constraintN,bestpos2names,pos1_cur,bestpos2points,bestconstraints,divisionlist,divisionthresh);
    
bestmatches=bestmatchesdiv;
 bestpos2names=bestpos2namesdiv;
 bestpos2points=bestpos2pointsdiv;
 bestconstraints=bestconstraintsdiv;
 bestscore=bestscorediv;
  if(exist('unnamednames','var'))
 accuracyinfo(2)=length(find(correctmatches(bestmatchesdiv, unnamednames,bestpos2namesdiv)))/length(bestmatchesdiv);
  end
 %'divs done'
%filter deaths
 %[ bestconstraintsdeath,bestscoredeath,bestpos2pointsdeath,bestpos2namesdeath,bestmatchesdeath,bestmatchnamesdeath,changelogdeath,changelogdeathrejected ] = FilterDeaths(iterations,bestscore,constraintnames,constraintAP,constraintLR,constraintDV,bestmatches,bestpos2names,pos1_cur,bestpos2points,bestconstraints,deathlist,unnamednames);
[ bestconstraintsdeath,bestscoredeath,bestpos2pointsdeath,bestpos2namesdeath,bestmatchesdeath,bestmatchnamesdeath,changelogdeath,changelogdeathrejected ] = FilterDeathsN...
    (iterations,constraintnames,constraintN,bestpos2names,pos1_cur,bestpos2points,bestconstraints,deathlist,deaththresh);
%'deaths done'
  %make so death results used in next step
 %bestmatches=bestmatchesdeath;
 bestmatches=bestmatchesdeath;%use reoptimized 
 bestpos2names=bestpos2namesdeath;
 bestpos2points=bestpos2pointsdeath;
 bestconstraints=bestconstraintsdeath;%use death
 bestscore=bestscoredeath;%use reoptimized
  if(exist('unnamednames','var'))
  accuracyinfo(3)=length(find(correctmatches(bestmatches, unnamednames,bestpos2namesdeath)))/length(bestmatchesdiv);
  end
%dont do reverse filtering since doesnt effect score and is time consuming
PERMUNKNOWN=false;
if (PERMUNKNOWN)
    [bestconstraintsperm,bestmatchnamesperm,bestmatchesperm,changelogperm,changelogpermrejected]=FilterUnnamedCellsN(iterations,constraintnames,constraintN,bestpos2names,pos1_cur,bestpos2points,bestconstraints,filterthresh, unnamednames);
    if(exist('unnamednames','var'))
        accuracyinfo(4)=length(find(correctmatches(bestmatchesperm, unnamednames,bestpos2namesdeath)))/length(bestmatchesperm);
    end
else
    bestconstraintsperm=bestconstraintsdeath;
    bestmatchesperm=bestmatchesdeath;
    changelogperm=[];
    changelogpermrejected=[];
end
 
%reoptimize post reverse filtering
%finaliterations=800;
finaliterations=min(finaliterations,prod(size(bestconstraints))/3);

%note its fiddly to pass grapha as optional variable because its alerady
%using confident as optional variable to signal could do it with checking
%for empty etc but not now
[violationsd,matches1to2d,bestviolationd,bestconstraintsd]=iterativeMatchWithConstraintsNoAnswerNPreGraph(constraintnames,constraintN,bestmatchesperm,bestpos2namesdeath,pos1_cur,bestpos2pointsdeath,finaliterations,bestconstraintsperm,grapha);
%[violationsd,matches1to2d,bestviolationd,bestconstraintsd]=iterativeMatchWithConstraintsNoAnswer(constraintnames,constraintAP,constraintLR,constraintDV,bestmatches,bestpos2namesdeath,pos1_cur,bestpos2pointsdeath,finaliterations,bestconstraints);
%'rematch done'
%length(find(correctmatches(matches1to2d, unnamednames,bestpos2namesdeath)))/length(matches1to2d)
%length(find(correctmatches(bestmatchesperm, unnamednames,bestpos2namesdeath)))/length(bestmatchesperm)
 if(exist('unnamednames','var'))
 accuracyinfo(5)=length(find(correctmatches(matches1to2d, unnamednames,bestpos2namesdeath)))/length(matches1to2d);
 end
[v]=countConstraintViolationsNeighborhood(constraintnames,constraintN,matches1to2d,bestpos2namesdeath,pos1_cur,grapha);
bestviolationd=v; %
changelog=[];
changelog.div=changelogdiv;
changelog.divr=changelogdivrejected;
changelog.death=changelogdeath;
changelog.deathr=changelogdeathrejected;
changelog.perm=changelogperm;
changelog.permr=changelogpermrejected;
else
    matches1to2d=matches1to2;
    bestviolations=violations;
    changelog=[];
    bestpos2namesdeath=names2_cur;
end

matches1to2names=[];
for i=1:length(matches1to2d)
    if(matches1to2d(i)~=0)
        matches1to2names{i}=bestpos2namesdeath{matches1to2d(i)};
    else
       matches1to2names{i}='';
    end
end