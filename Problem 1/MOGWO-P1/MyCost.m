 
function Z=MyCost(x,ndiscp,a,c0,nu,gamRR)
yy=[x;[0,0]];
yy=unique(yy,'rows');
[sizr,~]=size(yy);
[~,Z(2)]=tspsearch(yy,sizr);

% find mean variogram between sample points
F=squareform(pdist(x));
ss=F(:)';

% find corresponding variogram values

[~,n]=size(ss);
varss=zeros(1,n);
for i=1:n
varss(i)=sv(a,ss(i),c0,nu);
end
gamSS=mean(varss);

% find mean variogram between sample points and
% discretization points


dm=pdist2(x,ndiscp);
[NN,~]=size(ndiscp);
gamxR=0;
[N,~]=size(x);
variR=zeros(N,NN);
for i=1:N
    
    for j=1:NN
     variR(i,j)=sv(a,dm(i,j),c0,nu);
    end
gamxR=gamxR+mean(variR(i,:));
end
gamxR=gamxR/N;


Z(1)=(2*gamxR)-(gamSS)-gamRR;
end



