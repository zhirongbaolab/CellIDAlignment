function [matches1to2consensus,support1to2consensus,allnames]=matchPointCloundsConsensusMultiembOptimizedIteratedN(pos1_cur,emb2s,temb2matchs,window,divisionlist,deathlist,ICP,PERMUTE,constraintnames,neighborconstraints,iterations,initialiterations,finaliterations,divisionthresh,deaththresh,filterthresh,consensusmatch,consensussupport,supportthreshold,unnamednames)

%given a set of points and an embryo match against temporal match and
%window to either side extract consensus name for each pos1 point
%match is now a name allowing say a name that does not exist in temporal
%match time but does in majority of neihboring ones to be an answer
allnames={};
c=1;
for embs=1:length(emb2s)
    for win=-window:window
        currindex=temb2matchs(embs)+win;
        emb2=emb2s{embs};
        if(currindex>=1&currindex<size(emb2,2))
            pos2_cur=emb2(currindex).finalpoints;
            names2_cur=emb2(currindex).names;
          
            
            reliable3con=consensussupport>=supportthreshold;
             LM3=[];
            LM1=[];
            confidentp3=zeros(1,length(pos1_cur));
 confidentp1=zeros(1,length(pos2_cur));
confidentp1index=[];
for i=1:length(pos1_cur)
    if(reliable3con(i))
        inde1=find(strcmp(names2_cur,consensusmatch{i}));
        if~isempty(inde1)
            confidentp3(i)=1;
            confidentp1(inde1)=1;
            confidentp1index(end+1)=inde1; %record inex to reorder confident templates
            LM3=[LM3;pos1_cur(i,:)];
            LM1=[LM1;pos2_cur(inde1,:)];
        end
    end
end
confidentp3=logical(confidentp3);
confidentp1=logical(confidentp1);

pos1tps=TPS3D(LM3, LM1,pos1_cur);
%now rematch 
%alternative approach having warped based on confident
%filter p1 and p2 based on confident and rematch
%this seems to work in terms of improving the single case from .819 to
%.842 but unclear if iterating the whole consensus will get us ove .87
%divisionlist=divisionlist{embs};
ICP=false;%test for this iteration
%for neighbor version I'm passing in the confident ones and names as well,
%because they need to be reinserted before neighbor scoring is done
[matches3to1ptpsitv2,changelog31tpsitv2,violations31tpsitv2,accuracyinfo31tpsitv2]=matchPointCloundsWithPerturbationUpdatingConstraintsNIt...
    (pos1tps(~confidentp3,:),pos2_cur(~confidentp1,:),names2_cur(~confidentp1),divisionlist{embs},deathlist,ICP,PERMUTE,constraintnames,neighborconstraints,iterations,initialiterations,finaliterations,divisionthresh,deaththresh,filterthresh,...
    pos1tps(confidentp3,:),names2_cur(confidentp1index),unnamednames(~confidentp3));

%length(find(correctmatchesnames(matches3to1ptpsitv2, names3_cur(~confident3))))/length(matches3to1ptpsitv2)
matchnames=consensusmatch;
matchnames(~confidentp3)=matches3to1ptpsitv2; %assemble full result into result for this match


%[matchnames,changelog]=matchPointCloundsWithPerturbationUpdatingConstraints(pos1_cur,pos2_cur,names2_cur,divisionlist{embs},deathlist,ICP,constraintnames,constraintAP,constraintLR,constraintDV,iterations,initialiterations,finaliterations,divisionthresh,deaththresh,filterthresh,unnamednames);

allnames{c}=matchnames;
c=c+1;
        end
    end
end

%at this point we have n lists of name matches
%for each position we find consensus
matches1to2consensus={};
support1to2consensus=[];

for i=1:size(pos1_cur,1)
        currstrings={};
    for j=1:length(allnames)
        currstrings{j}=allnames{j}{i};
        
    end
    uniquelist=unique(currstrings);
    counts=[];
    for j=1:length(uniquelist)
        counts(j)=length(find(strcmpi(uniquelist{j},currstrings)));
    end
    [v,ind]=max(counts);
    matches1to2consensus{i}=uniquelist{ind};
    support1to2consensus(i)=counts(ind);
end

end

