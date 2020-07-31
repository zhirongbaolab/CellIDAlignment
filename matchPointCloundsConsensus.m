function [matches1to2consensus,support1to2consensus]=matchPointCloundsConsensus(pos1_cur,emb2,temb2match,window,ICP)
   
%given a set of points and an embryo match against temporal match and
%window to either side extract consensus name for each pos1 point
%match is now a name allowing say a name that does not exist in temporal
%match time but does in majority of neihboring ones to be an answer
allnames={};
c=1;
for win=-window:window
    currindex=temb2match+win;
    if(currindex>=1&currindex<size(emb2,2))
        pos2_cur=emb2(currindex).finalpoints;
        names2_cur=emb2(currindex).names;
         [TR,TT,pos2transto1ICP] = icp(pos1_cur,pos2_cur(:,1:3));
         pos2transto1ICP=pos2transto1ICP';
        if ICP
            pos2_cur=pos2transto1ICP;
        end
        
        [matches1to2,cost12]=matchPointCloundsFast(pos1_cur,pos2_cur(:,1:3));
        matchnames=[];
        for i=1:length(matches1to2)
            if(matches1to2(i)~=0)
                matchnames{i}=names2_cur{matches1to2(i)};
            else
                matchnames{i}='';
            end
        end
        allnames{c}=matchnames;
        c=c+1;
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

