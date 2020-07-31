function [constraintnames,c,ccount]=buildConstraintMatrixNeighborhood(allnamesa,allposa)
%takes a cell array of name cell arrays one per embryo and a positions
%arraay of matricies
%all embryos are aligned already  and returned is a constraint matrix where
%entry i,j asks is i anterior of j in all data sets etc

%extract list of unique cell names for all embryos
constraintnames=allnamesa{1};
for i=2:length(allnamesa)
    constraintnames=unique({allnamesa{i}{:},constraintnames{:}});
end

%note my initial implementation initialized this to ones
%and relied on updating it to 0 if the pair of names is found but not
%contacting
%this has the bug of setting all non-co-occuring pairs to 1
%likely just created a tare of 'missing' links on all occurances since
%these are connectins between parents and daughters primarily 
%but fixing it just in case
%c=ones(length(constraintnames),length(constraintnames));

c=zeros(length(constraintnames),length(constraintnames));

for j=1:length(allposa)
     dataDist=distFast(allposa{j},allposa{j});
    [dall,indall]=sort(dataDist,'ascend');
    grapha=GabrielGraph(allposa{j},indall);
    
    for h=1:length(constraintnames) %iterate over list of all names
        for i=1:length(constraintnames) %iterate over list of all names
            
            
            %find cell h if in this emb
            cellindh=find(strcmp(constraintnames(h),allnamesa{j}));
            
            %find cell i if in this emb
            cellindi=find(strcmp(constraintnames(i),allnamesa{j}));
            
            if ~isempty(cellindh)&&~isempty(cellindi)
                neighbors=grapha{cellindh};
                
                if(~isempty(find(neighbors==cellindi, 1)))
                 %previously this checked if it was empty and set to 0 if
                 %was
                    c(h,i)=c(h,i)+1;
                end
            end
        end
    end
    
end
ccount=c;
%these are newconstraint exists if present in all otherwise not
c(c<length(allposa))=0;
%think this was a bug in my bug fix that would make a lot (all?)
%constraints be lost, not sure what is going on here, yes I am, it leaves
%the equal cases as counts, and so fail ==1 test so always zero constraint
%violations
%c(c>length(allposa))=1;
c(c>=length(allposa))=1;



