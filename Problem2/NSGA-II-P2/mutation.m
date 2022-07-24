function  mutpop=mutation(mutpop,pop,nmut,nvar,samplev,a,c0,nu,lbc2,xmax,ymax)

npop=length(pop);


for n=1:nmut
    

   jj=0;    
   id=1;
   while (id==1)
   ii=randi([1 npop]);
   p=pop(ii).pos;
   
   j=randi([1 nvar]);
 
   if rand<0.5
       
    p(j,:)=p(j,:)-rand*0.1*(p(j,:));  
   else
    p(j,:)=p(j,:)+rand*0.1*(p(j,:));       
   end
   p(j,1)=min(max(0,p(j,1)),xmax);
   p(j,2)=min(max(0,p(j,2)),ymax);
   
   Ftem=min(pdist(p));
    
     if Ftem>lbc2
         id=0;
         mutpop(n).pos=p;
         mutpop(n).cost=MyCost(p,samplev,a,c0,nu)';
                        
     else
         jj=jj+1;
         if jj==3
         mutpop(n).pos=p;
         mutpop(n).cost=[10^6;10^6];
         id=0;
         end
     end
   end
   
  
   
     
end






