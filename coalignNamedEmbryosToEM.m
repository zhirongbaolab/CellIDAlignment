function [transformedembs]=coalignNamedEmbryosToEM(EMPOS,EMNAMES,lmtargetnames,embryostoalign,starttimes,endtimes)

%UNTITLED given a template embryo and list of n embryos to be aligned to it
%along with start and end time for computing linear temporal alignment 
%for end frame compute transform for each emb 
%do this anew at each frame
transformedembs={};
        names_1=EMNAMES;
        pos_1=EMPOS;
        
for e=1:length(embryostoalign)
    embt=embryostoalign{e};
    for frames=starttimes(e):endtimes(e)
        lmpositions1=[];
        lmpositions2=[];
        lmnames={};
        names_2=embryostoalign{e}(frames).names;
        pos_2=embryostoalign{e}(frames).finalpoints;
        c=1;
        for i=1:length(lmtargetnames)
            matchpoint1=[];
            matchpoint2=[];
            
            for j=1:length(names_1)
                %if same or if current cell is ancestor of target
                if strcmpi(lmtargetnames{i},names_1{j})||~isempty(strfind(lmtargetnames{i},names_1{j}))
                    indmatch1=j;
                    matchpoint1=pos_1(j,:);
                end
            end
            for j=1:length(names_2)
                %if same or if current cell is ancestor of target
                if strcmpi(lmtargetnames{i},names_2{j})||~isempty(strfind(lmtargetnames{i},names_2{j}))
                    indmatch2=j;
                    matchpoint2=pos_2(j,:);
                end
            end
            if (~isempty(matchpoint1)&~isempty(matchpoint2))
                lmnames{c}=lmtargetnames{i};
                c=c+1;
                lmpositions1=[lmpositions1;matchpoint1];
                lmpositions2=[lmpositions2;matchpoint2];
            end
        end
        
        postrans=TPS3D(lmpositions2,lmpositions1,embt(frames).finalpoints);
        embt(frames).finalpoints=postrans;
    end
    transformedembs{e}=embt;
end