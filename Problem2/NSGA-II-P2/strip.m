

function indices=strip(X,k)

T2 = clusterdata(X,'linkage','single','MaxClust',k);
Indix=[0];
j=1;

for i=1:length(X)
if ~ismember(T2(i),Indix)
indices(i)=i;
j=j+1;
Indix(j+1)=T2(i);
end
end

indices=nonzeros(indices)';