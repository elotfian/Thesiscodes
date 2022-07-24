function  crosspop=crossover(crosspop,pop,ncross,nvar,samplev,a,c0,nu,lbc2)

npop=length(pop);


for n=1:2:ncross
    
%     id3=1;
%     while(id3==1)
    jj=0;
    id=1;
    while (id==1)
    i1=randi([1 npop]);
    i2=randi([1 npop]);
    
    
    p1=pop(i1).pos;
    p2=pop(i2).pos;
    
    
    r=rand(1,nvar);
    o1=([r.*p1(:,1)';r.*p1(:,2)']+[(1-r).*p2(:,1)';(1-r).*p2(:,2)'])';
    o2=([r.*p2(:,1)';r.*p2(:,2)']+[(1-r).*p1(:,1)';(1-r).*p1(:,2)'])';
    
     Ftem1=min(pdist(o1));
     Ftem2=min(pdist(o2));
     if Ftem1>=lbc2 && Ftem2>=lbc2
         id=0;
crosspop(n).pos=o1;
crosspop(n).cost=MyCost(o1,samplev,a,c0,nu)';

crosspop(n+1).pos=o2;
crosspop(n+1).cost=MyCost(o2,samplev,a,c0,nu)';
     
     else
          jj=jj+1;
          if jj==3
             if Ftem1>lbc2
              crosspop(n).pos=o1;
              crosspop(n).cost=MyCost(o1,samplev,a,c0,nu)'; 
              crosspop(n+1).pos=o2;
              crosspop(n+1).cost=[10^6;10^6];
             end
             if Ftem2>lbc2
              crosspop(n).pos=o1;
              crosspop(n).cost=[10^6;10^6];
              crosspop(n+1).pos=o2;
              crosspop(n+1).cost=MyCost(o2,samplev,a,c0,nu)';
             end
             if Ftem1<lbc2 && Ftem2<lbc2 
              crosspop(n).pos=o1;
              crosspop(n).cost=[10^6;10^6];
              crosspop(n+1).pos=o2;
              crosspop(n+1).cost=[10^6;10^6];  
             end
         id=0;
         end
    end
    end


end





