function [ndiscp,samplev]=Discretize(ndp,xmax,ymax)

% pic2=imread('1.png');
% BW2= im2bw(pic2,0.9);
% m=size(BW2,1);
% n=size(BW2,2);
discp=zeros(ndp,2);
idp=0;

while idp<ndp
idp=idp+1;
x=xmax*rand(1,1);
y=ymax*rand(1,1);
discp(idp,:)=[x,y];
end
ndiscp=discp(1:ndp,:);

samplev=zeros(34*34,2);
A=0.5:3:100;
samplev(1:34,:)=[A',repmat(0.5,34,1)];
for i=1:33
    j=i*34+1;
samplev(j:j+33,:)=[A',repmat(i*3+0.5,34,1)];   
end

