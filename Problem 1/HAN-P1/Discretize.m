function ndiscp=Discretize(ndp,xmin,ymin,xmax,ymax)

discp=zeros(ndp,2);
idp=0;

while idp<ndp
idp=idp+1;

x=xmin+(xmax-xmin)*rand(1,1);
y=ymin+(ymax-ymin)*rand(1,1);
discp(idp,:)=[x,y];

end
ndiscp=discp(1:ndp,:);

