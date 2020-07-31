function [transformedembs]=coalignNamedEmbryos(emb_e1,temb1,t4emb1,anisotropy_1,embryostoalign,tsembstoalign,t4embstoalign,anisotropy_2)

%UNTITLED given a template embryo and list of n embryos to be aligned to it
%along with start and end time for computing linear temporal alignment 
%for end frame compute transform for each emb 
%apply to each frame for now
%optionally run ICP to realign each roughly aligned frame


%seam cell landmarks
landmarks={'ABplaaappa';'ABplaaappp';'ABarppaaap';'ABarppapaa';'ABarppapap';'ABplappapa';'ABarppappa';'ABplapapaa';'ABarppappp';'ABarpapppp';'ABarpppaap';'ABarppppaa';'ABarppppap';'ABprappapa';'ABarpppppa';'ABprapapaa';'ABarpppppp'};
%gut
landmarks2={'Ealaad';'Earaad';'Ealaav';'Earaav';'Earpa';'Ealpa';'Earap';'Ealap';'Epraa';'Eplaa';'Earpp';'Ealpp';'Eprap';'Eplap';'Eprpa';'Eplpa';'Eprppa';'Eplppa';'Eplppp';'Eprppp'};
%headlm exc hyp6 7,4 6x3 4x2 hyp7x2
landmarks3={'ABplpappaap','ABplaaaapp','ABarpaapap','ABarpapapa','ABplaaaapa','ABarpaapaa','ABarpapapp','ABplaappaa','ABpraappaa','ABplaapppp','ABpraapppp'};
%tail lm %pvqr pvql p11,p12
landmarks4={'ABprapppaaa','ABplapppaaa','ABplapappa','ABprapappa','Cappppv','Cpppppv'};

lmtargetnames={landmarks{:},landmarks3{:},landmarks4{:}};

transformedembs={};

for e=1:length(embryostoalign)

lmpositions1=[];
lmpositions2=[];
lmnames={};
names_1=emb_e1(temb1).names;
names_2=embryostoalign{e}(tsembstoalign{e}).names;
pos_1=emb_e1(temb1).finalpoints;
pos_2=embryostoalign{e}(tsembstoalign{e}).finalpoints;
c=1;
for i=1:length(lmtargetnames)
    matchpoint1=[];
    matchpoint2=[];
      
    for j=1:length(names_1)
        if strcmp(lmtargetnames{i},names_1{j})
            indmatch1=j;
            matchpoint1=pos_1(j,:);
        end
    end
    for j=1:length(names_2)
        if strcmp(lmtargetnames{i},names_2{j})
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

%compute transformation

lmpositions1(:,3)=lmpositions1(:,3).*anisotropy_1;
lmpositions2(:,3)=lmpositions2(:,3).*anisotropy_2;

transform2to1=[lmpositions1,ones(length(lmpositions1),1)]'/[lmpositions2,ones(length(lmpositions2),1)]';
alltransforms{e}=transform2to1;
end
%applytransforms to all timepoints
for i=1:length(embryostoalign)
    e=embryostoalign{i};
    for t=1:length(e)-1
        e(t).finalpoints=(alltransforms{i}*[e(t).finalpoints,ones(size(e(t).finalpoints,1),1)]')';
        e(t).finalpoints=e(t).finalpoints(:,1:3);
    end
    transformedembs{i}=e;
    
end


    %optionally apply icp to do