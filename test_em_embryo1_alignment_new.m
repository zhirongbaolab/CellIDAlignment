%assume have generated spatial and temporal alignments and constraint
%matrix

load('datafiles\Global_Model_Fully_Loaded.mat');
load('datafiles\embryodata\e1names.mat');
load('datafiles\embryodata\e1pos.mat');

alldispim={emb_e1,transformedembs{:}};
allstart=[t4emb1,t4emb2,t4emb3];
allend=[temb1,temb2,temb3];

i=temb1;
names1_cur=emb_e1(i).names;
pos1_cur=emb_e1(i).finalpoints;

alpha1=(i-t4emb1)/(temb1-t4emb1);
temb2match=round(t4emb2+(temb2-t4emb2)*alpha1); 
temb3match=round(t4emb3+(temb3-t4emb3)*alpha1);

pos2_cur=transformedembs{1}(temb2match).finalpoints;
pos3_cur=transformedembs{2}(temb3match).finalpoints;
names2_cur=transformedembs{1}(temb2match).names;
names3_cur=transformedembs{2}(temb3match).names;

landmarks={'ABplaaappa';'ABplaaappp';'ABarppaaap';'ABarppapaa';'ABarppapap';'ABplappapa';'ABarppappa';'ABplapapaa';'ABarppappp';'ABarpapppp';'ABarpppaap';'ABarppppaa';'ABarppppap';'ABprappapa';'ABarpppppa';'ABprapapaa';'ABarpppppp'};
%gut
landmarks2={'Ealaad';'Earaad';'Ealaav';'Earaav';'Earpa';'Ealpa';'Earap';'Ealap';'Epraa';'Eplaa';'Earpp';'Ealpp';'Eprap';'Eplap';'Eprpa';'Eplpa';'Eprppa';'Eplppa';'Eplppp';'Eprppp'}
%headlm exc hyp6 7,4 6x3 4x2 hyp7x2
landmarks3={'ABplpappaap','ABplaaaapp','ABarpaapap','ABarpapapa','ABplaaaapa','ABarpaapaa','ABarpapapp','ABplaappaa','ABpraappaa','ABplaapppp','Abpraapppp'};
%tail lm %pvqr pvql p11,p12
landmarks4={'ABprapppaaa','ABplapppaaa','ABplapappa','ABprapappa'};

%landmarks={landmarks{:},landmarks2{:},landmarks3{:},landmarks4{:}};

landmarks={landmarks{:},landmarks3{:},landmarks4{:}};


e1names=emb_e1(temb1).names;
e1pos=emb_e1(temb1).finalpoints;
LMWG=[];
LMEM=[];
for i=1:length(landmarks)
    indwg=find(strcmpi(landmarks{i},e1names));
    indem=find(strcmpi(landmarks{i},EMNAMES));
    if~isempty(indwg)&&~isempty(indem)
        LMWG=[LMWG;e1pos(indwg,:)];
        LMEM=[LMEM;EMPOS(indem,:)];
    end
end

%affine
xreal=[LMWG';LMEM'];%sort into matched pairs
transform=computeAffineLeastSquares( xreal);
EMtrans=(transform*[EMPOS,ones(length(EMPOS),1)]')';
EMtrans=EMtrans(:,1:3);



numberEMNamed=0;

for i=1:length(EMNAMES)
    if~isempty(EMNAMES{i})
        numberEMNamed=numberEMNamed+1;
    end
end
%at this EM and two other refs have all been matched to ref emb 1


PERMUTE=true;
iterations=200;%10;
initialiterations=0;%1400;
finaliterations=1400;%800;
divisionthresh=3;%55;
deaththresh=3;%55;
filterthresh=3;

%iterations=100;
%initialiterations=0;
%finaliterations=700;




%test run iterations
%{
iterations=2;
initialiterations=4;
finaliterations=4;
%}
%do optimized match for each of 3 targets
%{
divisionlist=divisionliste3all{temb3match};%use known embryo
[matchesEMto3p,changelogEM3,violationsEMto3,accuracyinfoEM3]=matchPointCloundsWithPerturbationUpdatingConstraints(EMtrans,pos3_cur,names3_cur,divisionlist,deathlist,ICP,constraintnames,constraintAP,constraintLR,constraintDV,iterations,initialiterations,finaliterations,divisionthresh,deaththresh,filterthresh,EMNAMES);
length(find(correctmatchesnames(matchesEMto3p, EMNAMES)))/length(matchesEMto3p)
length(find(correctmatchesnames(matchesEMto3p, EMNAMES)))/numberEMNamed

confident=[]
for i=1:length(EMCONFIDENCE)
    if(isempty(EMCONFIDENCE{i}))
        confident(i)=0;
    else
        confident(i)=1;
    end
end
length(find(correctmatchesnames(matchesEMto3p, EMNAMES)&confident))/sum(confident)


 %[matches3,cost3]=matchPointCloundsFast(EMtrans,pos3_cur(:,1:3));
% iterations=1000;
% [violations3t,bestanswer3t,bestviolation3t,bestconstraint3t]=iterativeMatchWithConstraintsNoAnswer(constraintnames,constraintAP,constraintLR,constraintDV,matches3,names_3,EMtrans,pos3_cur(:,1:3),iterations);
%length(find(correctmatches(bestanswer3, EMNAMES,names_3)))/length(bestanswer3)
%visualizeMatches( EMtrans, pos3_cur,  EMNAMES, names3_cur,matchesEMto3p );title('long simple iteration 3');

matchsimple=[];
for j=1:length(matchesEMto3p)
    v=find(strcmp(matchesEMto3p{j},names3_cur));
    if(isempty(v))
        matchsimple(j)=0;
    else
        matchsimple(j)=v;
    end
end
visualizeMatches( EMtrans, pos3_cur,  EMNAMES, names3_cur, matchsimple );title(' 1 2 2 emb3');

divisionlist=divisionliste2all{temb2match};
[matchesEMto2p,changelogEM2,violationsEMto2,accuracyinfoEM2]=matchPointCloundsWithPerturbationUpdatingConstraints(EMtrans,pos2_cur,names2_cur,divisionlist,deathlist,ICP,constraintnames,constraintAP,constraintLR,constraintDV,iterations,initialiterations,finaliterations,divisionthresh,deaththresh,filterthresh,EMNAMES);
length(find(correctmatchesnames(matchesEMto2p, EMNAMES)))/length(matchesEMto2p)
length(find(correctmatchesnames(matchesEMto2p, EMNAMES)))/numberEMNamed
length(find(correctmatchesnames(matchesEMto2p, EMNAMES)&confident))/sum(confident)



matchsimple=[];
for j=1:length(matchesEMto2p)
    v=find(strcmp(matchesEMto2p{j},names2_cur));
    if(isempty(v))
        matchsimple(j)=0;
    else
        matchsimple(j)=v;
    end
end
visualizeMatches( EMtrans, pos2_cur,  EMNAMES, names2_cur, matchsimple );title(' 1 2 2 emb2');


divisionlist=divisionliste1all{temb1};
[matchesEMto1p,changelogEM1,violationsEMto1,accuracyinfoEM1]=matchPointCloundsWithPerturbationUpdatingConstraints(EMtrans,pos1_cur,names1_cur,divisionlist,deathlist,ICP,constraintnames,constraintAP,constraintLR,constraintDV,iterations,initialiterations,finaliterations,divisionthresh,deaththresh,filterthresh,EMNAMES);
length(find(correctmatchesnames(matchesEMto1p, EMNAMES)))/length(matchesEMto1p)
length(find(correctmatchesnames(matchesEMto1p, EMNAMES)))/numberEMNamed
length(find(correctmatchesnames(matchesEMto1p, EMNAMES)&confident))/sum(confident)


matchsimple=[];
for j=1:length(matchesEMto1p)
    v=find(strcmp(matchesEMto1p{j},names1_cur));
    if(isempty(v))
        matchsimple(j)=0;
    else
        matchsimple(j)=v;
    end
end
visualizeMatches( EMtrans, pos1_cur,  EMNAMES, names1_cur, matchsimple );title(' 1 2 2 emb1');


%}

%one possible match with 1 as unknown etc
window=6;
[matches2consensus,support2consensus,allnamesinit]=matchPointCloundsConsensusMultiembOptimizedN(EMtrans,alldispim,allend,window,alldivisionlist,deathlist,ICP,PERMUTE,constraintnames,neighborconstraints,iterations,initialiterations,finaliterations,divisionthresh,deaththresh,filterthresh,EMNAMES);
%length(find(correctmatchesnames(matches2consensus, EMNAMES)))/length(matches2consensus)
length(find(correctmatchesnames(matches2consensus, EMNAMES)))/(numberEMNamed)
%length(find(correctmatchesnames(matches2consensus, EMNAMES)&EMCONFIDENCE'))/sum(EMCONFIDENCE)
[accuracyinfo,tissuetally]=tissueLevelAccuracy(matches2consensus, EMNAMES,partlist)

%{
%note imposing exclusive matches on consensus names does not seem to help
matches2consensus_ex=computeExclusiveConsensusNames(allnamesinit);
length(find(correctmatchesnames(matches2consensus_ex, EMNAMES)))/(numberEMNamed)
length(find(correctmatchesnames(matches2consensus_ex, EMNAMES)&EMCONFIDENCE'))/sum(EMCONFIDENCE)
[accuracyinfo,tissuetally]=tissueLevelAccuracy(matches2consensus_ex, EMNAMES,partlist)
visualizeMatchesNamed( EMtrans, pos1_cur,  EMNAMES, names1_cur, matchesEMconsensus );title('EM e1 to WG consensus shown on emb 1)');
%}

'done with round 1'

ICPIT=false
thresh=17;%8
[matches2consensusIT,support2consensusIT,allnamesit1]=matchPointCloundsConsensusMultiembOptimizedIteratedN(EMtrans,alldispim,allend,window,alldivisionlist,deathlist,ICPIT,PERMUTE,constraintnames,neighborconstraints,iterations,initialiterations,finaliterations,divisionthresh,deaththresh,filterthresh,matches2consensus,support2consensus,thresh,EMNAMES);
%length(find(correctmatchesnames(matches2consensusIT, EMNAMES)))/length(matches2consensusIT)
%compute exculsive names and see if any better
length(find(correctmatchesnames(matches2consensusIT, EMNAMES)))/(numberEMNamed)
%length(find(correctmatchesnames(matches2consensusIT, EMNAMES)&EMCONFIDENCE'))/sum(EMCONFIDENCE)
[accuracyinfo,tissuetally]=tissueLevelAccuracy(matches2consensusIT, EMNAMES,partlist)

%{
matches2consensusIT_ex=computeExclusiveConsensusNames(allnamesit1);
length(find(correctmatchesnames(matches2consensusIT_ex, EMNAMES)))/(numberEMNamed)
length(find(correctmatchesnames(matches2consensusIT_ex, EMNAMES)&EMCONFIDENCE'))/sum(EMCONFIDENCE)
[accuracyinfo,tissuetally]=tissueLevelAccuracy(matches2consensusIT_ex, EMNAMES,partlist)
%}

%it 2
%thresh=8;
[matches2consensusIT2,support2consensusIT2 ,allnamesit2]=matchPointCloundsConsensusMultiembOptimizedIteratedN(EMtrans,alldispim,allend,window,alldivisionlist,deathlist,ICPIT,PERMUTE,constraintnames,neighborconstraints,iterations,initialiterations,finaliterations,divisionthresh,deaththresh,filterthresh,matches2consensusIT,support2consensusIT,thresh,EMNAMES);
%length(find(correctmatchesnames(matches2consensusIT2, EMNAMES)))/length(matches2consensusIT2)
length(find(correctmatchesnames(matches2consensusIT2, EMNAMES)))/(numberEMNamed)
%length(find(correctmatchesnames(matches2consensusIT2, EMNAMES)&EMCONFIDENCE'))/sum(EMCONFIDENCE)


%iteration 2 seems to be best on test data set so we are calling that our
%final result
[ terminalfinalmatches ] = SulstontoTerminal(matches2consensusIT2,partlist,deathlist );

%it 3
%thresh=8
[matches2consensusIT3,support2consensusIT3,allnamesit3]=matchPointCloundsConsensusMultiembOptimizedIteratedN(EMtrans,alldispim,allend,window,alldivisionlist,deathlist,ICPIT,PERMUTE,constraintnames,neighborconstraints,iterations,initialiterations,finaliterations,divisionthresh,deaththresh,filterthresh,matches2consensusIT2,support2consensusIT2,thresh,EMNAMES);
%length(find(correctmatchesnames(matches2consensusIT3, EMNAMES)))/length(matches2consensusIT3)
length(find(correctmatchesnames(matches2consensusIT3, EMNAMES)))/(numberEMNamed)
%length(find(correctmatchesnames(matches2consensusIT3, EMNAMES)&EMCONFIDENCE'))/sum(EMCONFIDENCE)

return
%{
numberEMNamed=0;
numberEMNamedConfident=0;

for i=1:length(EMNAMES)
    if~isempty(EMNAMES{i})
        numberEMNamed=numberEMNamed+1;
        if (EMCONFIDENCE(i)==1)
            numberEMNamedConfident=numberEMNamedConfident+1;
        end
    end
end
%}
[test1,EMNAMES]=xlsread('D:\EM_embryo\cellnames_6_5_2020.xlsx');
EMCONFIDENCE=xlsread('D:\EM_embryo\cellconfidences_6_5_2020.xlsx');
%EMCONFIDENCE=cell2mat(EMCONFIDENCE);
EMCONFIDENCE(isnan(EMCONFIDENCE))=0;

return
%do consensus match
window=4;
alldivisionlist={divisionliste1all,divisionliste2all,divisionliste2all};
%one possible match with 1 as unknown etc
[matchesEMconsensus,supportEMconsensus]=matchPointCloundsConsensusMultiembOptimized(EMtrans,alldispim,allend,window,alldivisionlist,deathlist,ICP,constraintnames,constraintAP,constraintLR,constraintDV,iterations,initialiterations,finaliterations,divisionthresh, deaththresh,filterthresh,EMNAMES);
length(find(correctmatchesnames(matchesEMconsensus, EMNAMES)))/length(matchesEMconsensus)
length(find(correctmatchesnames(matchesEMconsensus, EMNAMES)))/numberEMNamed
length(find(correctmatchesnames(matchesEMconsensus, EMNAMES)&confident))/sum(confident)

%support
correctsupport=histc(supportEMconsensus(correctmatchesnames(matchesEMconsensus, EMNAMES)&confident),linspace(0,15,15));
wrongsuport=histc(supportEMconsensus(~correctmatchesnames(matchesEMconsensus, EMNAMES)&confident),linspace(0,15,15));
figure; plot(correctsupport(1:end-1)/(correctsupport(1:end-1)+wrongsuport(1:end-1)))
return
%below is test code to see what the result of base matching for each emb
%was (not so good)
    allmatches={}
    for j=1:length(alldispim)
        temb2match=allend(j);
        pos2_cur=alldispim{j}(temb2match).finalpoints;
        names2_cur=alldispim{j}(temb2match).names;
        
        constraintnames=allconstraints(temb2match).constraintnames;
        constraintAP=allconstraints(temb2match).constraintAP;
        constraintLR= allconstraints(temb2match).constraintLR;
        constraintDV=allconstraints(temb2match).constraintDV;
        %divisionlist=divisionlistall{temb2match};
        
        [TR,TT,pos2transto1ICP] = icp(pos1_cur,pos2_cur(:,1:3));
        pos2transto1ICP=pos2transto1ICP';
        [matches1to2ICP,cost12ICP]=matchPointCloundsFast(pos1_cur,pos2transto1ICP(:,1:3));
        matchnames=[];
            for i=1:length(matches1to2ICP)
                if(matches1to2ICP(i)~=0)
                    matchnames{i}=names2_cur{matches1to2ICP(i)};
                else
                    matchnames{i}='';
                end
            end
        allmatches{j}=matchnames';        
    end
    ICP=true;
    window=2;
    %[matches1to2consensus,support1to2consensus]=matchPointCloundsConsensus(pos1_cur,transformedembs{1},temb2match,window,ICP);
    [matches1to2consensus,support1to2consensus]=matchPointCloundsConsensusMultiemb(pos1_cur,alldispim,allend,window,ICP);
    matches1to2consensus=matches1to2consensus';
    support1to2consensus=support1to2consensus';
    
    %[matches1to2consensus]=matchPointCloundsWithPerturbation(pos1_cur,pos2_cur,names2_cur,divisionlist,deathlist,ICP,constraintnames,constraintAP,constraintLR,constraintDV);

