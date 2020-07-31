function  divisionlist=findNearbyDivisions(emb_e2,temb2match,divisionwindow)

%extract list of nearby divisions in template that might need to try splitting or merging %returns list of parents
divisionlist={};
c=1;
for i=-divisionwindow:divisionwindow
    tpoint=temb2match+i;
    if(tpoint>=1&tpoint<length(emb_e2)-1)
        for j=1:size(emb_e2(tpoint).suc,1)
            if(emb_e2(tpoint).suc(j,2)~=-1)
                divisionlist{c,1}=emb_e2(tpoint).names{j}; %store parent name if there is a division
                divisionlist{c,2}=emb_e2(tpoint+1).names{emb_e2(tpoint).suc(j,1)}; %store daughter 1 name if there is a division
                divisionlist{c,3}=emb_e2(tpoint+1).names{emb_e2(tpoint).suc(j,2)}; %store daughter 2 name if there is a division
                 
                c=c+1;
            end
        end
    
    end

end