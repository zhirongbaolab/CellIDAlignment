function [accuracyinfo,tissuetally]=tissueLevelAccuracySimplified(matches,answers,partlist)

tissuetally={'neuron';'muscle',;'amphid';'hyp';'seam';'gut';'pharynx';'other'};%,unique({partlist{notempt,4}});
for i=1:size(tissuetally,1)
    tissuetally{i,2}=0;%cases col
    tissuetally{i,3}=0;%cases correct at cell level
    tissuetally{i,4}=0;%non correct cases correct at tissue level
end

accuracyinfo=-1*ones(size(matches));

correctness=correctmatchesnames(matches, answers);
for i=1:length(matches)
    
    %look up tissue index
    indmatch=[];
    indanswer=[];
    %partlist match of automated
    if ~isempty(matches{i})
        indmatch=find(strcmpi(matches{i},partlist(:,2)));
    end
    %parlist match of answer
    if ~isempty(answers{i})
        indanswer=find(strcmpi(answers{i},partlist(:,2)));
    end
    tissueindex=[];
    if ~isempty(indanswer)
        tissueindex=find(strcmpi(partlist{indanswer,4},tissuetally(:,1)));
    end
    if(isempty(tissueindex))
        tissueindex=size(tissuetally,1);
    end
    if~isempty(answers{i})
       %increment found in category
        tissuetally{tissueindex,2}=tissuetally{tissueindex,2}+1;
        if(correctness(i))
          accuracyinfo(i)=1;
        tissuetally{tissueindex,3}=tissuetally{tissueindex,3}+1;
    else
        if ~isempty (indanswer)&~isempty(indmatch)
            if strcmpi(partlist{indanswer,4},partlist{indmatch,4})
                accuracyinfo(i)=2;
                tissuetally{tissueindex,4}=tissuetally{tissueindex,4}+1;
            end
        end
    end
    end
    
    
end

%{
for i=1:length(matches)
    indmatch=[];
    indanswer=[];
    %partlist match of automated
    if ~isempty(matches{i})
        indmatch=find(strcmpi(matches{i},partlist(:,2)));
    end
    if(strcmpi(answers{i},'Abarpaapaa'))
        'test case'
    end
    %parlist match of answer
    if ~isempty(answers{i})
        indanswer=find(strcmpi(answers{i},partlist(:,2)));
    end
    %if both real and predicted are in part list
    if ~isempty(indanswer)&~isempty(indmatch)
        %find row in tally corresponding to actual tissue
        tissueindex=find(strcmpi(partlist{indanswer,4},tissuetally(:,1)));
        %        tissueindexnswer=find(strcmpi(tissuetally(:,1),partlist{indmatch,4}));
        if isempty(tissueindex)
            tissueindex=size(tissuetally,1);
        end
        %increment found in category
        tissuetally{tissueindex,2}=tissuetally{tissueindex,2}+1;
        
        %match by identical
        if strcmpi(matches{i},answers{i}) %indanswer==indmatch
            accuracyinfo(i)=1;
            tissuetally{tissueindex,3}=tissuetally{tissueindex,3}+1;
        else
            if strcmpi(partlist{indanswer,4},partlist{indmatch,4})
                accuracyinfo(i)=2;
                tissuetally{tissueindex,4}=tissuetally{tissueindex,4}+1;
            end
        end
        
    else %is nonterminal/death and counts as other\
        if(~isempty(answers{i}))
           tissueindex=size(tissuetally,1);  
            %increment found in category
            tissuetally{tissueindex,2}=tissuetally{tissueindex,2}+1;        
        end
        if ~isempty(answers{i})&~isempty(matches{i})      
            %uncovered case of non terminal cells
            if strcmpi(matches{i},answers{i})
                accuracyinfo(i)=1;
                tissuetally{tissueindex,3}=tissuetally{tissueindex,3}+1;    
            end
        end
    end
end
%}