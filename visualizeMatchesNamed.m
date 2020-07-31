function visualizeMatchesNamed( p1, p2,  names1, names2, matches )
%visualize the incorrect and missling links from match over point cloud
matchsimple=[];
for j=1:length(matches)
    v=find(strcmp(matches{j},names2));
    if(isempty(v))
        matchsimple(j)=0;
    else
        matchsimple(j)=v;
    end
end
matches=matchsimple;
[ correct,correctreverse,realmatch,realmatchreverse ] = computeCorrectnessAndRealMatches( names1,names2,matches);
%{
length(find(correct==1))/length(correct)
length(find(correctreverse==1))/length(correctreverse)

length(find(correct==0))
length(find(correctreverse==0))

length(find(correct==0&realmatch==0))
length(find(correctreverse==0&realmatchreverse==0))

realmatchdistances=[];
for i=1:length(realmatch)
if(realmatch(i)~=0)
    realmatchdistances=[realmatchdistances;distance(p1(i,1:3)',p2(realmatch(i),1:3)')];
end
end
%}

%2d plot


figure
hold on
scatter3(p2(:,1),p2(:,2),p2(:,3),'.b')
scatter3(p1(:,1),p1(:,2),p1(:,3),'.g')
for i=1:length(p1)
    
    if~correct(i) %wrong cases
        if(matches(i)~=0) %match exists
            if (realmatch(i)~=0)
                plot3([p2(matches(i),1),p1(i,1)],[p2(matches(i),2),p1(i,2)],[p2(matches(i),3),p1(i,3)],'r'); %wrong match
           end
        else
                 scatter3(p1(i,1),p1(i,2),p1(i,3),'*m');%wrong no match made 
           
        end
        
        if(realmatch(i)~=0)
            
            plot3([p1((i),1),p2(realmatch(i),1)],[p1(i,2),p2(realmatch(i),2)],[p1((i),3),p2(realmatch(i),3)],'g'); %missing correct
            if(distance(p1',p2')>30)
                ['weird long distance match for ' ,names2{realmatch(i)}]
            end
        else
           % scatter3(p1(i,1),p1(i,2),p1(i,3),'*k'); %no real match exists for set 1
        end
    else %correct cases
        %     plot3([p1(matchi),1),p2(i,1)],[p1(matches(i),2),p2(i,2)],[p1(matches(i),3),p2(i,3)],'b'); %correct and chosen
        if (matches(i)~=0)
            plot3([p2(matches(i),1),p1(i,1)],[p2(matches(i),2),p1(i,2)],[p2(matches(i),3),p1(i,3)],'b'); %wrong match
        else
          %  'correct match to nothing'
            scatter3(p1(i,1),p1(i,2),p1(i,3),'*b'); %no real match exists for set 1 and none is made
   
        end
    end
    if (realmatch(i)~=0)
      if(distance(p1(i,:)',p2(realmatch(i),:)')>35)
                ['weird long distance match ',num2str(distance(p1(i,:)',p2(realmatch(i),:)')), ' for ' ,names1{i}, names2{realmatch(i)}]
            end
    end
end

for i=1:length(p2)
    if~correctreverse(i)
       
        if(realmatchreverse(i)==0)
          %  scatter3(p2(i,1),p2(i,2),p2(i,3),'*c'); %no real match exists for set2
        else
            
            % plot3([p2((i),1),p1(realmatchreverse(i),1)],[p2(i,2),p1(realmatchreverse(i),2)],[p2((i),3),p1(realmatchreverse(i),3)],'g');
            % %missing correct alerady covered

        end 
    else
        %correct cases a
       % plot3([p1(matches(i),1),p2(i,1)],[p1(matches(i),2),p2(i,2)],[p1(matches(i),3),p2(i,3)],'b'); %correct and chosen
    end
     
end



axis equal
end

