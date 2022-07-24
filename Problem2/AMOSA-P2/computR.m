function  R=computR(xn,xc,Arc)

for i=1:length(Arc)
Ctem(i,:)=Arc(i).cost';
end
Ctem(length(Arc)+1,:)=xc';
Ctem(length(Arc)+2,:)=xn';


R(1)=max(Ctem(:,1))-min(Ctem(:,1));
R(2)=max(Ctem(:,2))-min(Ctem(:,2));
