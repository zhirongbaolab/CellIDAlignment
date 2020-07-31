function[v]=countConstraintViolationsNeighborhood(constraintnames,c,matches,names,cellpos,grapha)
%given constraint matrix on sorting of points, matches, names for the
%destination set and points for the source set
%count for each matched name, if all matches are correct how many
%violations of constrains are present. 
%note matches,cellpos are the unnamed embryo, names are from the named
%embryo

v=zeros(size(matches));
if(~exist('grapha','var'))
    %graph for this timepoint
    dataDist=distFast(cellpos,cellpos);
    [dall,indall]=sort(dataDist,'ascend');
    grapha=GabrielGraph(cellpos,indall);
end

%compute match names
matchnames=cell(length(matches),1);
for i=1:length(matches)
    if(matches(i)~=0)
        matchnames{i}=names{matches(i)};
    end
end

%precompute matched name locations in constraint
namelookupindexC=zeros(1,length(constraintnames));%index of each name in names in the constraint matrix compute this outside loop bc it is a hot spot 50%of calculation
namelookupindex=zeros(1,length(matchnames));%index of each name in names in the constraint matrix compute this outside loop bc it is a hot spot 50%of calculation
for i=1:length(matches)
    if ~isempty(matchnames{i})
        constraintindi=find(strcmp(matchnames{i},constraintnames));
        if isempty(constraintindi)
           % namelookupindex(i)=0;
        else
            namelookupindexC(constraintindi)=i;     
            namelookupindex(i)=constraintindi;
        end
    end
end

%}
%{
%this should be faster by 2x but theres some type bug so I give up
[~,namelookupindex]=(ismember(matchnames,constraintnames'));
namelookupindexC=zeros(size(constraintnames));
for i=1:length(namelookupindex)
    if(namelookupindex(i)~=0)
        namelookupindexC(namelookupindex(i))=i;
    end
end
%}

%now we iterate over each cell, iterate over its row in constraint matrix
%and if find a 1 and find that other cell in our data set and their not
%neighbors in grapha add a violation
%{
for i=1:length(matchnames)
    indi=namelookupindex(i); %i is cell in question in current data set indi is in constraints
    if(indi~=0)
           %iterate over constraints
        for j=1:length(constraintnames)
            %our cell and something is a neighbor
            if(c(indi,j)==1)
                indj=namelookupindexC(j); %location in current data set of match of this constraint cell
                if(indj~=0)
                    %exists in current data set
                      neighbors=grapha{i};
                
                    if(max(neighbors==indj)==0)
                        v(i)=v(i)+1;
                    end
                end
            end
        end
    
    end
end
%}
%alternative loop, is find more efficient than iteration above where
%if(c(indi,j)==1) was crazy hot spot
for i=1:length(matchnames)
    indi=namelookupindex(i); %i is cell in question in current data set indi is in constraints
    if(indi~=0)
        %iterate over constraints
        targets=find(c(indi,:)==1);
        for j=1:length(targets)
            indj=namelookupindexC(targets(j));
            if(indj~=0)
                %exists in current data set
                neighbors=grapha{i}; 
                if(max(neighbors==indj)==0)
                    v(i)=v(i)+1;
                end
            end
        end
    end
end
