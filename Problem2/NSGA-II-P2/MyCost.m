 
function Z=MyCost(x,samplev,a,c0,nu)
  yy=[x;[0,0]];
  [sizr,~]=size(yy);
  [~,Z(2)]=tspsearch(yy,sizr);


%% Relation (4.3)
F=squareform(pdist(x));
ss=F(:)';
[~,n]=size(ss);
varss=zeros(1,n);
for i=1:n
varss(i)=cv(a,ss(i),c0,nu);
end
Cmat=vec2mat(varss,sqrt(n)); %matrix C without one's
A=[Cmat ones(sqrt(n),1); [ones(1,sqrt(n)),0]]; %matrix C and one's
invA=pinv(A);
invCm=pinv(Cmat);
% d, lambda and sigma_ok w.r.t s_0 and equation (4.4)
for i=1:size(samplev,1)
d(i,:)=pdist2(samplev(i,:),x);
  [~,nd]=size(d(i,:));
  varss0=zeros(1,nd);
   for j=1:nd
   varss0(i,j)=cv(a,d(i,j),c0,nu);
   end
%Lambda(i,:)=(inv(A)* [varss0(i,:),1]')';
Lambda(i,:)=(invA*[varss0(i,:),1]')';
sigmaok(i)=1-(Lambda(i,:)*[varss0(i,:),1]');
end

%% covariance between the varigram parameters equation (4.7)

rhoCnu=zeros(1,n); 
rhoCa=zeros(1,n); 
for i=1:n
[rhoCnu(i),rhoCa(i)]=parti(a,ss(i),nu);
end
RhoCnu=vec2mat(rhoCnu,sqrt(n));%\partial C / \partal nu
RhoCa=vec2mat(rhoCa,sqrt(n));%\partial C / \partal a



%covariance 4.7
covcnuca=1/(0.5*trace((invCm*RhoCnu)*(invCm*RhoCa)));

DAnu=[RhoCnu zeros(sqrt(n),1); [zeros(1,sqrt(n)),0]];%\partial A / \partial nu
DAa=[RhoCa zeros(sqrt(n),1); [zeros(1,sqrt(n)),0]]; %\partial A / \partial a
%\partial d
for i=1:size(samplev,1)
  [~,nd]=size(d(i,:));
   dnuvarss0=zeros(1,nd);
   davarss0=zeros(1,nd);
   for j=1:nd
   [dnuvarss0(i,j),davarss0(i,j)]=parti(a,d(i,j),nu);
   end
  
   Dnuvarss0(i,:)=[dnuvarss0(i,:),0];%\partial d / \partial nu
   Davarss0(i,:)=[davarss0(i,:),0]; %\partial d / \partial a
   DLamnu(i,:)=(invA*(Dnuvarss0(i,:)-(DAnu*Lambda(i,:)')')')';
   DLama(i,:)=(invA*(Davarss0(i,:)-(DAa*Lambda(i,:)')')')';
end

Cmat2=[Cmat zeros(sqrt(n),1); [zeros(1,sqrt(n)),0]];
for i=1:size(samplev,1)
  ET(i)=covcnuca*(DLamnu(i,:)*Cmat2*DLama(i,:)');
  sigma2p(i)=sigmaok(i)+ET(i);
end

Sigma2p=mean(sigma2p);



 Z(1)=Sigma2p;

  end
 


