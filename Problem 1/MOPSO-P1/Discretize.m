function ndiscp=Discretize(ndp,xmax,ymax)

discp=zeros(1000,2);
idp=0;

while idp<ndp
idp=idp+1;

x=xmax*rand(1,1);
y=ymax*rand(1,1);
discp(idp,:)=[x,y];

end
ndiscp=discp(1:ndp,:);

