function [correct ] = correctmatches(matches1to2, names_1,names_2)
%array of whether match names match real names

rearrangednames2={};
%{names_1{matches1to2}}';
for i=1:length(names_1)
    if(matches1to2(i)~=0)
        rearrangednames2{i}=names_2{matches1to2(i)};
    else
        rearrangednames2{i}=[];
    end
end


correct=(strcmpi(rearrangednames2',names_1));


end

