%assume have generated spatial and temporal alignments and constraint
%matrix

load('datafiles\Global_Model_Fully_Loaded.mat');
load('datafiles\embryodata\e1names.mat');
load('datafiles\embryodata\e1pos.mat');


alldispim={emb_e1,transformedembs{:}};

allend=[300,290,270];
i=allend(1);
names1_cur=emb_e1(i).names;
pos1_cur=emb_e1(i).finalpoints;

pos2_cur=transformedembs{1}(allend(2)).finalpoints;
pos3_cur=transformedembs{2}(allend(3)).finalpoints;
names2_cur=transformedembs{1}(allend(2)).names;
names3_cur=transformedembs{2}(allend(3)).names;
%seam cells w.o tr/l and wo h0r (I think) for no opbvious reason
landmarks1={'ABplaaappa';'ABplaaappp';'ABarppaaap';'ABarppapaa';'ABarppapap';'ABplappapa';'ABarppappa';'ABplapapaa';'ABarppappp';'ABarpapppp';'ABarpppaap';'ABarppppaa';'ABarppppap';'ABprappapa';'ABarpppppa';'ABprapapaa';'ABarpppppp'};

%seam cells (still w/0 tr/l but also without h0 h1 bc unsure about them.
landmarks1={'ABarppaaap';'ABarppapaa';'ABarppapap';'ABplappapa';'ABarppappa';'ABplapapaa';'ABarppappp';'ABarpppaap';'ABarppppaa';'ABarppppap';'ABprappapa';'ABarpppppa';'ABprapapaa';'ABarpppppp'};

%gut% early generation
landmarks2={'Ealaa';'Earaa';'Earpa';'Ealpa';'Earap';'Ealap';'Epraa';'Eplaa';'Earpp';'Ealpp';'Eprap';'Eplap';'Eprpa';'Eplpa';'Eprpp';'Eplpp'}
%headlm hyp3(early) exc hyp6 7,4 6x3 4x2 hyp7x2
landmarks3={'ABplaapaaaa','ABplpappaap','ABplaaaapp','ABarpaapap','ABarpapapa','ABplaaaapa','ABarpaapaa','ABarpapapp','ABplaappaa','ABpraappaa','ABplaapppp','Abpraapppp'};
%tail lm %pvqr pvql p11,p12
landmarks4={'ABprapppaaa','ABplapppaaa','ABplapappa','ABprapappa'};

%landmarks={landmarks{:},landmarks2{:},landmarks3{:},landmarks4{:}};
%  hyp4,   hyp7(2)  hyp8/9x2   spikex2, hyp6x2,AVL,RMDDR/L,RIH, hyp1, hyp3, exc cell tl tr p12 p11
landmarks5={'ABarpapapa','caappd','cpappd','Abplpppapap','Abprpppapap','Abprppppppa','Abplppppppa','ABplaaaapa','Abarpaapaa','ABprpappaap','ABarappapaa','Abalpapapaa','ABprpappaaa','ABarappaapa','ABplaapaaaa','ABplpappaap','ABplappppp','ABprappppp','ABprapappa','ABplapappa'};
landmarks={landmarks1{:},landmarks2{:},landmarks5{:}};


e1names=emb_e1(allend(1)).names;
e1pos=emb_e1(allend(1)).finalpoints;
LMWG=[];
LMEM=[];
for i=1:length(landmarks)
    indwg=find(strcmpi(landmarks{i},e1names));
    indem=find(strcmpi(landmarks{i},EMNAMES));
    if(length(indem)>1)
    'typo'
    end
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

numberEMNamed=0;

for i=1:length(EMNAMES)
    if~isempty(EMNAMES{i})
        numberEMNamed=numberEMNamed+1;
    end
end

iterations=200;%10;
initialiterations=0;%1400;
finaliterations=1400;%800;
divisionthresh=3;%25;%55;
deaththresh=3;%25;%55;
filterthresh=3;%25;

%test run iterations
iterations=2;
initialiterations=0;
finaliterations=2;


ICP=true;
PERMUTE=true;


   constraintnames=allconstraints(allend(1)).constraintnames;
   neighborconstraints=allconstraints(allend(1)).neighborconstraints;
   %    constraintAP=allconstraints(allend(1)).constraintAP;
%    constraintLR= allconstraints(allend(1)).constraintLR;
%    constraintDV=allconstraints(allend(1)).constraintDV;
 

[transformedembs]=coalignNamedEmbryosPerTime(emb_e1,temb1,t4emb1,anisotropy_1,{emb_e2,emb_e3},{temb2,temb3},{t4emb2,t4emb3},anisotropy_2);

    %do consensus match
window=6;
alldivisionlist={divisionliste1all,divisionliste2all,divisionliste2all};
[matchesEMconsensus,supportEMconsensus,allnamesinit]=matchPointCloundsConsensusMultiembOptimizedN(EMtrans,alldispim,allend,window,alldivisionlist,deathlist,ICP,PERMUTE,constraintnames,neighborconstraints,iterations,initialiterations,finaliterations,divisionthresh,deaththresh,filterthresh,EMNAMES);
length(find(correctmatchesnames(matchesEMconsensus, EMNAMES)))/numberEMNamed
%length(find(correctmatchesnames(matchesEMconsensus, EMNAMES)&EMCONFIDENCE'))/sum(EMCONFIDENCE)


visualizeMatchesNamed( EMtrans, pos1_cur,  EMNAMES, names1_cur, matchesEMconsensus );title('EM e2 to WG consensus shown on emb 1)');
terminalEMmatch={};
terminalEMmatch=SulstontoTerminal(matchesEMconsensus,partlist,deathlist);
    
%iteration


ICPIT=false
thresh=33;%8
[matchesEMconsensusIT,support2consensusIT,allnamesit]=matchPointCloundsConsensusMultiembOptimizedIteratedN(EMtrans,alldispim,allend,window,alldivisionlist,deathlist,ICPIT,PERMUTE,constraintnames,neighborconstraints,iterations,initialiterations,finaliterations,divisionthresh,deaththresh,filterthresh,matchesEMconsensus,supportEMconsensus,thresh,EMNAMES);
%length(find(correctmatchesnames(matchesEMconsensusIT, EMNAMES)))/length(matchesEMconsensusIT)
length(find(correctmatchesnames(matchesEMconsensusIT, EMNAMES)))/(numberEMNamed)

%length(find(correctmatchesnames(matchesEMconsensusIT, EMNAMES)&EMCONFIDENCE'))/sum(EMCONFIDENCE)


%tissue level accuracy
[accuracyinfo,tissuetally]=tissueLevelAccuracy(matchesEMconsensusIT, EMNAMES,partlist);

%length(find(EMCONFIDENCE'&correctmatchesnames(matchesEMconsensusIT, EMNAMES)))/(numberEMNamedConfident)

terminalEMmatch=SulstontoTerminal(matchesEMconsensusIT,partlist,deathlist);

%it 2
%thresh=8;
[matches2consensusIT2,support2consensusIT2,allnamesit2]=matchPointCloundsConsensusMultiembOptimizedIteratedN(EMtrans,alldispim,allend,window,alldivisionlist,deathlist,ICPIT,PERMUTE,constraintnames,neighborconstraints,iterations,initialiterations,finaliterations,divisionthresh,deaththresh,filterthresh,matchesEMconsensusIT,support2consensusIT,thresh,EMNAMES);
%length(find(correctmatchesnames(matches2consensusIT2, EMNAMES)))/length(matches2consensusIT2)
length(find(correctmatchesnames(matches2consensusIT2, EMNAMES)))/(numberEMNamed)
%length(find(correctmatchesnames(matches2consensusIT2, EMNAMES)&EMCONFIDENCE'))/sum(EMCONFIDENCE)


%it 3
%thresh=8
[matchesEMconsensusIT3,support2consensusIT3,allnamesit3]=matchPointCloundsConsensusMultiembOptimizedIteratedN(EMtrans,alldispim,allend,window,alldivisionlist,deathlist,ICPIT,PERMUTE,constraintnames,neighborconstraints,iterations,initialiterations,finaliterations,divisionthresh,deaththresh,filterthresh,matches2consensusIT2,support2consensusIT2,thresh,EMNAMES);
%length(find(correctmatchesnames(matchesEMconsensusIT3, EMNAMES)))/length(matchesEMconsensusIT3)
length(find(correctmatchesnames(matchesEMconsensusIT3, EMNAMES)))/(numberEMNamed)
%length(find(correctmatchesnames(matchesEMconsensusIT3, EMNAMES)&EMCONFIDENCE'))/sum(EMCONFIDENCE)


%{
%note these lose leading blanks
EMCONFIDENCE=xlsread('D:\EM_embryo\EM_emb2\textbasedconfidence642020.xlsx')
EMCONFIDENCE(isnan(EMCONFIDENCE))=0;
test=xlsread('D:\EM_embryo\EM_emb2\textbasedanswerkey642020.xlsx')
EMNAMES=test;

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

return


window=1;
alldivisionlist={divisionliste1all,divisionliste2all,divisionliste2all};
%one possible match with 1 as unknown etc
[matchesEMconsensus,supportEMconsensus]=matchPointCloundsConsensusMultiembOptimized(EMtrans,alldispim,allend,window,alldivisionlist,deathlist,ICP,PERMUTE,constraintnames,constraintAP,constraintLR,constraintDV,iterations,initialiterations,finaliterations,divisionthresh, deaththresh,filterthresh,EMNAMES);
length(find(correctmatchesnames(matchesEMconsensus, EMNAMES)))/numberEMNamed
return

%
    
%}
%do optimized match for each of 3 targets

divisionlist=divisionliste1all{allend(1)};
[matchesEMto1p,changelogEM1,violationsEMto1,accuracyinfoEM1]=matchPointCloundsWithPerturbationUpdatingConstraints(EMtrans,pos1_cur,names1_cur,divisionlist,deathlist,ICP,PERMUTE,constraintnames,constraintAP,constraintLR,constraintDV,iterations,initialiterations,finaliterations,divisionthresh,deaththresh,filterthresh,EMNAMES);
length(find(correctmatchesnames(matchesEMto1p, EMNAMES)))/numberEMNamed
%length(find(correctmatchesnames(matchesEMto1p, EMNAMES)&confident))/sum(confident)
visualizeMatchesNamed( EMtrans, pos1_cur,  EMNAMES, names1_cur, matchesEMto1p );title('EM e2 to WG emb 1)');
















return

divisionlist=divisionliste3all{allend(3)};%use known embryo
[matchesEMto3p,changelogEM3,violationsEMto3,accuracyinfoEM3]=matchPointCloundsWithPerturbationUpdatingConstraints(EMtrans,pos3_cur,names3_cur,divisionlist,deathlist,ICP,PERMUTE,constraintnames,constraintAP,constraintLR,constraintDV,iterations,initialiterations,finaliterations,divisionthresh,deaththresh,filterthresh,EMNAMES);
%length(find(correctmatchesnames(matchesEMto3p, EMNAMES)))/length(matchesEMto3p)
length(find(correctmatchesnames(matchesEMto3p, EMNAMES)))/numberEMNamed

visualizeMatchesNamed( EMtrans, pos3_cur,  EMNAMES, names3_cur, matchesEMto3p );title('EM e2 to WG emb 3)');


%{
confident=[]
for i=1:length(EMCONFIDENCE)
    if(isempty(EMCONFIDENCE{i}))
        confident(i)=0;
    else
        confident(i)=1;
    end
end
length(find(correctmatchesnames(matchesEMto3p, EMNAMES)&confident))/sum(confident)
%}
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
%visualizeMatches( EMtrans, pos3_cur,  EMNAMES, names3_cur, matchsimple );title(' 1 2 2 emb3');

divisionlist=divisionliste2all{allend(2)};
[matchesEMto2p,changelogEM2,violationsEMto2,accuracyinfoEM2]=matchPointCloundsWithPerturbationUpdatingConstraints(EMtrans,pos2_cur,names2_cur,divisionlist,deathlist,ICP,PERMUTE,constraintnames,constraintAP,constraintLR,constraintDV,iterations,initialiterations,finaliterations,divisionthresh,deaththresh,filterthresh,EMNAMES);
length(find(correctmatchesnames(matchesEMto2p, EMNAMES)))/numberEMNamed
%length(find(correctmatchesnames(matchesEMto2p, EMNAMES)&confident))/sum(confident)
visualizeMatchesNamed( EMtrans, pos2_cur,  EMNAMES, names2_cur, matchesEMto2p );title('EM e2 to WG emb 2)');




matchsimple=[];
for j=1:length(matchesEMto2p)
    v=find(strcmp(matchesEMto2p{j},names2_cur));
    if(isempty(v))
        matchsimple(j)=0;
    else
        matchsimple(j)=v;
    end
end
%visualizeMatches( EMtrans, pos2_cur,  EMNAMES, names2_cur, matchsimple );title(' 1 2 2 emb2');




return
matchsimple=[];
for j=1:length(matchesEMto1p)
    v=find(strcmp(matchesEMto1p{j},names1_cur));
    if(isempty(v))
        matchsimple(j)=0;
    else
        matchsimple(j)=v;
    end
end
%visualizeMatches( EMtrans, pos1_cur,  EMNAMES, names1_cur, matchsimple );title(' 1 2 2 emb1');

%length(find(correctmatchesnames(matchesEMconsensus, EMNAMES)&confident))/sum(confident)

return

%support
correctsupport=histc(supportEMconsensus(correctmatchesnames(matchesEMconsensus, EMNAMES)&confident),linspace(0,15,15));
wrongsuport=histc(supportEMconsensus(~correctmatchesnames(matchesEMconsensus, EMNAMES)&confident),linspace(0,15,15));
figure; plot(correctsupport(1:end-1)/(correctsupport(1:end-1)+wrongsuport(1:end-1)))




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

