function [ m ] = computeAffineLeastSquares( x )
%given x which consists of [x...;y,...;z...;x2..;y2...;z2..]
%   compute best fit transform of xyz1 to xyz2


x1=[x(1:3,:);ones(1,size(x,2))];
x2=[x(4:6,:);ones(1,size(x,2))];
m=x1/x2;

end

%{
figure; scatter3(x1(1,:),x1(2,:),x1(3,:))
hold on;
 scatter3(x2(1,:),x2(2,:),x2(3,:),'.k')
x2transform=m*x2;
%}