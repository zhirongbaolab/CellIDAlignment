function [ correct,correctreverse,realmatch,realmatchreverse ] = computeCorrectnessAndRealMatches( names1,names2,matches)
%given two sets of names and a set of computed matches between them
%computes correctness and  real matches for each in both directions

correct=[];
realmatch=[];
for i=1:length(names1)
    realmatchi=find(strcmpi(names1{i},names2));
       if(matches(i)~=0)
        correct(i)=strcmpi(names1{i},names2{matches(i)});
    else %no match given
        if(~isempty(realmatchi))
            correct(i)=0;
        else
            correct(i)=1;%no match exists so is correct to not match it
        end
    end
    
    if ~isempty(realmatchi)
        realmatch(i)=realmatchi;
    else
        realmatch(i)=0;
    end
    
end
%matches==linspace(1,length(matches),length(matches)); % 70% correct

realmatchreverse=[];
correctreverse=[]
for i=1:length(names2)
      
    realmatchi=find(strcmpi(names2{i},names1));
    computedmatchi=find(matches==i);
    if length(computedmatchi)>1
        computedmatchi=computedmatchi(1);
    end
    if isempty(computedmatchi)
        if isempty(realmatchi)
            correctreverse(i)=1;
        else
            correctreverse(i)=0;
        end
    else
        
        if(~isempty(realmatchi)&realmatchi==computedmatchi)
            correctreverse(i)=1;
        else
            correctreverse(i)=0;
        end
    end
  
    if ~isempty(realmatchi)
        realmatchreverse(i)=realmatchi;
    else
        realmatchreverse(i)=0;
    end

end

end

