function [correct ] = correctmatchesnames(matches1to2names, names_1)
%array of whether match names match real names

for i=1:length(names_1)
correct(i)=(strcmpi(matches1to2names{i},names_1{i})&~isempty(names_1{i}));
end

end

