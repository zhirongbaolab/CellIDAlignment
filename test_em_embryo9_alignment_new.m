%assume have generated spatial and temporal alignments and constraint
%matrix
clear
load('datafiles\Global_Model_Fully_Loaded.mat');
load('datafiles\partlist_withtissuelevelannotation.mat');
[~,EMNAMES]=xlsread('datafiles\embryodata\namesformatlab_e9.xlsx');
[EMPOS,~]=xlsread('datafiles\embryodata\positionsformatlab_e9.xlsx');
EMPOS(:,3)=EMPOS(:,3).*6; %compensate for anisotropy didnt do this in previous ones and relied on alignment to isotropic florescent data but think with few landmarks this is problematic
if (length(EMNAMES)<length(EMPOS))
    EMNAMES{length(EMPOS)}='';
end


%e2 was 360 350 looks better test after debug basic alignment
allend=[375,360,322];%[300,290,270];
allend=[360,350,315];%[300,290,270];
allend=[360,340,315];%[300,290,270];
i=allend(1);
names1_cur=emb_e1(i).names;
pos1_cur=emb_e1(i).finalpoints;

pos2_cur=transformedembs{1}(allend(2)).finalpoints;
pos3_cur=transformedembs{2}(allend(3)).finalpoints;
names2_cur=transformedembs{1}(allend(2)).names;
names3_cur=transformedembs{2}(allend(3)).names;

%seam cells h0r,         h0l,           h1l,        h1r         h2l             h2r         v4r v4l
landmarks1={'ABarpapppa','ABplaaappa','ABplaaappp','ABarpapppp','ABarppaaap','ABarpppaap','ABarpppppa','ABarppappa',...
    'ABarppapaa','ABarppppaa','ABarppapap','ABarppppap','ABplappapa','Abprappapa',}; 
%v1l v1r                    v2l v2r                     v3l, v3r

%tail seam cells
landmarks1b={'ABplapapaa','ABprapapaa','ABprappppp','Abplappppp','ABarpppppp','ABarppappp'};
%seam cells (still w/0 tr/l but also without h0 h1 bc unsure about them.
%landmarks1={'ABarppaaap';'ABarppapaa';'ABarppapap';'ABplappapa';'ABarppappa';'ABplapapaa';'ABarppappp';'ABarpppaap';'ABarppppaa';'ABarppppap';'ABprappapa';'ABarpppppa';'ABprapapaa';'ABarpppppp'};

%gut% early generation
%landmarks2={'Ealaa';'Earaa';'Earpa';'Ealpa';'Earap';'Ealap';'Epraa';'Eplaa';'Earpp';'Ealpp';'Eprap';'Eplap';'Eprpa';'Eplpa';'Eprpp';'Eplpp'}
%gut terminal generation
landmarks2={'Ealaad';'Earaad';'Ealaav';'Earaav';'Earpa';'Ealpa';'Earap';'Ealap';'Epraa';'Eplaa';'Earpp';'Ealpp';'Eprap';'Eplap';'Eprpa';'Eplpa';'Eprppa';'Eplppa';'Eplppp';'Eprppp'}

%headlm hyp3(early) exc hyp6 7,4 6x3 4x2 hyp7x2
%landmarks3={'ABplaapaaaa','ABplpappaap','ABplaaaapp','ABarpaapap','ABarpapapa','ABplaaaapa','ABarpaapaa','ABarpapapp','ABplaappaa','ABpraappaa','ABplaapppp','Abpraapppp'};
%tail lm %pvqr pvql p11,p12
%landmarks4={'ABprapppaaa','ABplapppaaa','ABplapappa','ABprapappa'};

%landmarks={landmarks{:},landmarks2{:},landmarks3{:},landmarks4{:}};
%hyp4,   hyp7(2)  hyp8/9x2   spikex2, hyp6x2,AVL,RMDDR/L,RIH, hyp1, hyp3, exc cell tl tr p12 p11

% ala exc cell, avl,                                    
landmarks5={'ABalapppaaa','ABplpappaap','ABprpappaap',...
    'ABarpapapa','ABarpaappa','ABarpapapp','Abpraappaa','Abplaappaa','Abpraapppa','Abplaapppa',...
    'ABplappapp','ABprappapp','ABarappapapp','ABalpapapapp','Abarpapapa','MSappaap','MSpppaap'...
    'ABprpaapapa','ABplpaapapa','ABprappppa', 'ABplappppa','ABarpappap','ABplaaapap',...
    'ABplaapppp','ABpraapppp','ABplppaaaaa','ABprppaaaaa','Abprpppapap','Abplpppapap'};
%hyp4mindline, hyp7%% midline  hyp 6  hyp4 x2 hyp7 ventralx2, 
%p7/8l r smbdr smbdl,midline hyp 4 z4, z1
%amsor amsol, hyp7tail,l/r hyp5 l/r
%hyp7 above p's sibdl,r hyp8/9 in tail

%landmarks={landmarks1{:},landmarks5{:}};

landmarks={landmarks1{:},landmarks1b{:},landmarks2{:},landmarks5{:}};
%landmarks={landmarks1{:},landmarks1b{:},landmarks5{:}};


e1names=emb_e1(allend(1)).names;
e1pos=emb_e1(allend(1)).finalpoints;
LMWG=[];
LMEM=[];
for i=1:length(landmarks)
    indwg=find(strcmpi(landmarks{i},e1names));
    indem=find(strcmpi(landmarks{i},EMNAMES));
    if(length(indem)>1)
    ['typo ',landmarks{i},' ',num2str(indem')]
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
%try visualizing to find error in tps landmarks that I'm sure is still
%there
dummymatches=zeros(length(EMtrans),1);

%visualizeMatches( EMtrans, pos1_cur,  EMNAMES, names1_cur, dummymatches );


%tps prealign
EMtrans=TPS3D(LMEM,LMWG, EMPOS);



%at this EM and two other refs have all been matched to ref emb 1
numberEMNamed=0;
numberEMNamedConfident=0;

for i=1:length(EMNAMES)
    if~isempty(EMNAMES{i})
        numberEMNamed=numberEMNamed+1;
%        if (EMCONFIDENCE(i)==1)
%            numberEMNamedConfident=numberEMNamedConfident+1;
%        end
    end
end

iterations=200;%10;
initialiterations=0;%1400;
finaliterations=1400;%800;
divisionthresh=3;%25;%55;
deaththresh=3;%25;%55;
filterthresh=3;%25;

%test run iterations
%iterations=5;
%initialiterations=0;
%finaliterations=5;

ICP=true;
PERMUTE=false;

%'test without permute'

   constraintnames=allconstraints(allend(1)).constraintnames;
   neighborconstraints=allconstraints(allend(1)).neighborconstraints;%
 

    %do consensus match
window=6; %was 6 for other embs maybe smaller better bc fast elongation?

alldispim={emb_e1,transformedembs{:}};

%replace alldispim linearl coalinged to e1 with tps coaligned to EM
alldispim=coalignNamedEmbryosToEM(EMtrans,EMNAMES,landmarks,{emb_e1,emb_e2,emb_e3},allend-window,allend+window);


pos1_cur=alldispim{1}(allend(1)+3).finalpoints;
pos2_cur=alldispim{2}(allend(2)+3).finalpoints;
pos3_cur=alldispim{3}(allend(3)+3).finalpoints;
names1_cur=alldispim{1}(allend(1)+3).names;
names2_cur=alldispim{2}(allend(2)+3).names;
names3_cur=alldispim{3}(allend(3)+3).names;

%debug initial alignment which seems to be broken partially

visualizeMatches( EMtrans, pos1_cur,  EMNAMES, names1_cur, dummymatches );
axis vis3d
visualizeMatches( EMtrans, pos2_cur,  EMNAMES, names2_cur, dummymatches );
axis vis3d
visualizeMatches( EMtrans, pos3_cur,  EMNAMES, names3_cur, dummymatches );
axis vis3d

%this was not how divisions were handled in prior scripts, not sure why
%possibly a bug but wouldn't make much difference since should be same
alldivisionlist={divisionliste1all{allend(1)},divisionliste2all{allend(2)},divisionliste3all{allend(3)}};
[matchesEMconsensus,supportEMconsensus,allnamesinit]=matchPointCloundsConsensusMultiembOptimizedN(EMtrans,alldispim,allend,window,alldivisionlist,deathlist,ICP,PERMUTE,constraintnames,neighborconstraints,iterations,initialiterations,finaliterations,divisionthresh,deaththresh,filterthresh,EMNAMES);
length(find(correctmatchesnames(matchesEMconsensus, EMNAMES)))/numberEMNamed
%length(find(correctmatchesnames(matchesEMconsensus, EMNAMES)&EMCONFIDENCE'))/sum(EMCONFIDENCE)

visualizeMatchesNamed( EMtrans, pos1_cur,  EMNAMES, names1_cur, matchesEMconsensus );title('EM e9 to WG consensus shown on emb 1)');
visualizeMatchesNamed( EMtrans, pos2_cur,  EMNAMES, names2_cur, matchesEMconsensus );title('EM e9 to WG consensus shown on emb 2)');
visualizeMatchesNamed( EMtrans, pos3_cur,  EMNAMES, names3_cur, matchesEMconsensus );title('EM e9 to WG consensus shown on emb 3)');

terminalEMmatch={};
terminalEMmatch=SulstontoTerminal(matchesEMconsensus,partlist,deathlist);
    
%iteration

ICPIT=false

thresh=28; %same 85 percent threshold %33 of 40 or 17 of 21?
[matchesEMconsensusIT,support2consensusIT,allnamesit]=matchPointCloundsConsensusMultiembOptimizedIteratedN(EMtrans,alldispim,allend,window,alldivisionlist,deathlist,ICPIT,PERMUTE,constraintnames,neighborconstraints,iterations,initialiterations,finaliterations,divisionthresh,deaththresh,filterthresh,matchesEMconsensus,supportEMconsensus,thresh,EMNAMES);
%length(find(correctmatchesnames(matchesEMconsensusIT, EMNAMES)))/length(matchesEMconsensusIT)
length(find(correctmatchesnames(matchesEMconsensusIT, EMNAMES)))/(numberEMNamed)

%length(find(correctmatchesnames(matchesEMconsensusIT, EMNAMES)&EMCONFIDENCE'))/sum(EMCONFIDENCE)

%tissue level accuracy
%[accuracyinfo,tissuetally]=tissueLevelAccuracy(matchesEMconsensusIT, EMNAMES,partlist);

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
