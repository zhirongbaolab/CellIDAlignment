function [matches1to2names,changelog,bestviolationd,accuracyinfo]=matchPointCloundsWithPerturbationUpdatingConstraintsNIt(pos1_cur,pos2_cur,names2_cur,divisionlist,deathlist,ICP,PERMUTE,constraintnames,constraintN,iterations,initialiterations,finaliterations,divisionthresh,deaththresh,filterthresh,...
    confidentpos1,confidentnames2,unnamednames)
%confident pos1,2 and names 2 are confident cases that should not be
%matched anew, but which need to be inserted before scoring of matches is
%done against neighbor model.
%I made this a separate function rather than making these 3 optional
%entries because there are many existing calls to the basic function
%without them but with the optional answer key variable. 

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
        config = gmmreg_load_config('G:\My Documents\MATLAB\gmmreg_cpd\gmmreg_cpd/face.ini');
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

bestscore=inf;
bestpos2points=[];
bestpos2names=names2_cur;

%initial match
[matches1to2,cost12]=matchPointClounds(pos1_cur,pos2_cur(:,1:3));
 initialconstraints=zeros(length(pos1_cur),length(pos2_cur));
 
 
if(exist('unnamednames','var'))
accuracyinfo(1)=length(find(correctmatches(matches1to2, unnamednames,names2_cur)))/length(matches1to2);
 end
 

 bestconstraints=initialconstraints;
 
matchnames=[];
for i=1:length(matches1to2)
    if(matches1to2(i)~=0)
        matchnames{i}=names2_cur{matches1to2(i)};
    else
        matchnames{i}='';
    end
end

%revised for iterative (all else is same)
%concatenate certain cases which have been rearranged to be in same order
%onto end of points and names and matches before evaluating 
%***
[violations]=countConstraintViolationsNeighborhood(constraintnames,constraintN,...
    [matches1to2,linspace(length(names2_cur)+1,length(names2_cur)+length(confidentnames2),length(confidentnames2))],...
    {names2_cur{:},confidentnames2{:}},[pos1_cur;confidentpos1]);


bestmatches=matches1to2;
bestmatchnames=matchnames;
bestpos2points=pos2_cur;
bestscore=sum(violations);

if (PERMUTE)
%'initial match done'
%alternate path where do it div,death,perm,reoptimize
%revise
%****
[ bestconstraintsdiv,bestscorediv,bestpos2pointsdiv,bestpos2namesdiv,bestmatchesdiv,bestmatchnamesdiv,changelogdiv,changelogdivrejected ] =...
    FilterDivisionsN(iterations,constraintnames,constraintN,bestpos2names,pos1_cur,bestpos2points,bestconstraints,divisionlist,divisionthresh,...
      confidentpos1,confidentnames2);
    
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
%revise
[ bestconstraintsdeath,bestscoredeath,bestpos2pointsdeath,bestpos2namesdeath,bestmatchesdeath,bestmatchnamesdeath,changelogdeath,changelogdeathrejected ] = FilterDeathsN...
    (iterations,constraintnames,constraintN,bestpos2names,pos1_cur,bestpos2points,bestconstraints,deathlist,deaththresh,...
      confidentpos1,confidentnames2);
  
%'deaths done'
  %make so death results used in next step
 bestmatches=bestmatchesdeath;%use reoptimized 
 bestpos2names=bestpos2namesdeath;
 bestpos2points=bestpos2pointsdeath;
 bestconstraints=bestconstraintsdeath;%use death
 bestscore=bestscoredeath;%use reoptimized
  if(exist('unnamednames','var'))
  accuracyinfo(3)=length(find(correctmatches(bestmatches, unnamednames,bestpos2namesdeath)))/length(bestmatchesdiv);
  end
%dont do reverse filtering since doesnt effect score and is time consuming
%revise

%this was turned on which is why initial run was longer than expected (and
%had more graph calls) but doesnt explain difference in runtime of 0th
%iteration on other computer, reprofile on mine?
PERMUNKNOWN=false;
if (PERMUNKNOWN)
    [bestconstraintsperm,bestmatchnamesperm,bestmatchesperm,changelogperm,changelogpermrejected]=FilterUnnamedCellsN(iterations,constraintnames,constraintN,bestpos2names,pos1_cur,bestpos2points,bestconstraints,filterthresh, unnamednames,...
        confidentpos1,confidentnames2);
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

%revise
[violationsd,matches1to2d,bestviolationd,bestconstraintsd]=iterativeMatchWithConstraintsNoAnswerN(constraintnames,constraintN,bestmatchesperm,bestpos2namesdeath,pos1_cur,bestpos2pointsdeath,finaliterations,bestconstraintsperm,...
      confidentpos1,confidentnames2);

if(exist('unnamednames','var'))
 accuracyinfo(5)=length(find(correctmatches(matches1to2d, unnamednames,bestpos2namesdeath)))/length(matches1to2d);
end

 %revise
%[v]=countConstraintViolationsNeighborhood(constraintnames,constraintN,matches1to2d,bestpos2namesdeath,pos1_cur);

[v]=countConstraintViolationsNeighborhood(constraintnames,constraintN,...
    [matches1to2d,linspace(length(bestpos2namesdeath)+1,length(bestpos2namesdeath)+length(confidentnames2),length(confidentnames2))],...
    {bestpos2namesdeath{:},confidentnames2{:}},[pos1_cur;confidentpos1]);


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
    bestviolationd=violations;
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