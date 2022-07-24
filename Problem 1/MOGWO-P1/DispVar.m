function  gamRR=DispVar(ndiscp,a,c0,nu)

% find mean variogram between discretization points
F=squareform(pdist(ndiscp));
disp=F(:)';

% find corresponding variogram values

[~,n]=size(disp);
varp=zeros(1,n);
for i=1:n
varp(i)=sv(a,disp(i),c0,nu);
end
gamRR=mean(varp);

end