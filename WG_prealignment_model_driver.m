%a prealignment driver preprocessing the WG embryos to make model 
%timings are set up for full data set, but results of this should now be
%sucessor to everything in old files (and equivalent other than loading all
%timepoints)  the emnames and emconfidences variables are what are missing
%to actually run a script
%though note the details of interior temporal alignment are actually
%different because of different endpoint


%note uses endpoint alignment like prior single timepoint scripts, rather
%than per timepoint alingment like the time experiments within WG embs
anisotropy_2=1;
%emb_2='L:\santella\nih_emb_qc\emb2\Decon_emb1_beforeMGedits.zip'; %e2
emb_2='L:\santella\nih_emb_qc\emb2\Decon_emb1_MGedits.zip'; %e2

temb2=360; %last edited frame and best overall shape match %329;%e2
t4emb2=4;


%emb_1='L:\santella\nih_emb_qc\new_new_version\Decon_emb1_edited_v13_mirrored.zip';%e1
emb_1='L:\santella\nih_emb_qc\new_new_version\Decon_emb1_edited_v13_11202018namingfix.zip';
temb1=380; %380 last frame 375 340;%e1 %synced timepoint matching EM
t4emb1=16;
anisotropy_1=1;
%emb_3='L:\santella\nih_emb_qc\emb3\Decon_emb1_updated_v2.zip';%e3
emb_3='L:\santella\nih_emb_qc\emb3\Decon_emb1_updated.zip';%e3

temb3=329;%322 seems to be best tail length match 296;%way worse310;%;296;%299 worse%%; %e3 %entirely possible that 310 is better
t4emb3=3;
anisotropy_3=1;

templocation='temp_unzip\';
%unzip zipfile to temp file
if ~exist(emb_1,'file')
    errors.zipfilemissing=1;
    return
end
try
    unzip(emb_1,templocation);
catch exception
    errors.zipfilecorrupted=1;
    return
end

[ cells_1,emb_e1] = loadcells_unnamed(templocation,temb1,anisotropy_1,false );
rmdir(templocation,'s');


if ~exist(emb_2,'file')
    errors.zipfilemissing=1;
    return
end
try
    unzip(emb_2,templocation);
catch exception
    errors.zipfilecorrupted=1;
    return
end

[ cells_2,emb_e2] = loadcells_unnamed(templocation,temb2,anisotropy_2,false );
rmdir(templocation,'s');



if ~exist(emb_3,'file')
    errors.zipfilemissing=1;
    return
end
try
    unzip(emb_3,templocation);
catch exception
    errors.zipfilecorrupted=1;
    return
end

[ cells_3,emb_e3] = loadcells_unnamed(templocation,temb3,anisotropy_3,false );
rmdir(templocation,'s');



names_1=emb_e1(temb1).names;
names_2=emb_e2(temb2).names;
names_3=emb_e3(temb3).names;
pos_1=emb_e1(temb1).finalpoints;

%aligns each timepoint on emb list to first template embryo using fixed
%landmarks 
%for testing matching aside from initial alignment/temporal,spatial
[transformedembs]=coalignNamedEmbryos(emb_e1,temb1,t4emb1,anisotropy_1,{emb_e2,emb_e3},{temb2,temb3},{t4emb2,t4emb3},anisotropy_2);

'coalignment done'

window=6; %2;%5;
allconstraints=[];
'note only computing constraints for temoral neighborhood needed for naming e9'
for i=temb1-40:temb1
    metanames={};
    metapos={};
    for h=-window:window
        if(i+h>=1&i+h<length(emb_e1))
            names1_cur=emb_e1(i+h).names;
            pos1_cur=emb_e1(i+h).finalpoints;
            
            alpha1=((i+h)-t4emb1)/(temb1-t4emb1);
            temb2match=round(t4emb2+(temb2-t4emb2)*alpha1);
            temb3match=round(t4emb3+(temb3-t4emb3)*alpha1);
            temb2match=max(1,temb2match);
            temb3match=max(1,temb3match);
            pos2_cur=transformedembs{1}(temb2match).finalpoints;
            pos3_cur=transformedembs{2}(temb3match).finalpoints;
            names2_cur=transformedembs{1}(temb2match).names;
            names3_cur=transformedembs{2}(temb3match).names;
            metanames={metanames{:},names1_cur,names2_cur,names3_cur};
            metapos={metapos{:},pos1_cur,pos2_cur(:,1:3),pos3_cur(:,1:3)};
        end
    end
    
    %now building constraint marix each timepoint
    [constraintnames,neighborconstraints,concount]=buildConstraintMatrixNeighborhood(metanames,metapos);
    
    allconstraints(i).constraintnames=constraintnames;
    allconstraints(i).concount=concount;
    allconstraints(i).neighborconstraints=neighborconstraints;
    i
end

%this is not how divisions were done for e1,e2 but this makes more sense
%basically e1,e2 seemed to be using one over and over
 for i=t4emb1:temb1
  divisionwindow=10;
    divisionlist=findNearbyDivisions(emb_e1,i,divisionwindow); %extract list of nearby divisions in template that might need to try splitting or merging %returns list of parents
    divisionliste1all{i}=divisionlist;
 end
    
 for i=t4emb2:temb2
   divisionwindow=10;
    divisionlist=findNearbyDivisions(emb_e2,i,divisionwindow); %extract list of nearby divisions in template that might need to try splitting or merging %returns list of parents
    divisionliste2all{i}=divisionlist;
 end
 
 for i=t4emb3:temb3
   divisionwindow=10;
    divisionlist=findNearbyDivisions(emb_e3,i,divisionwindow); %extract list of nearby divisions in template that might need to try splitting or merging %returns list of parents
    divisionliste3all{i}=divisionlist;
 end
 
 
 [~,deathlist]=xlsread('D:\EM_embryo\CellDeaths.csv');
