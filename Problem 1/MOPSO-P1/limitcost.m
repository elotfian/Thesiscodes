function pop=limitcost(pop,cost1max,cost2max,ndiscp,nvar,a,C0,C1,gamRR,pcen,lb,ub,prad)

n=size(pop,1);
emp.pos=[];
emp.cost=[];
emp.rank=[];      
emp.cdis=[];      % crowding distance


pop2=repmat(emp,n,1);

for i=1:n
    if pop(i).cost(1)<=cost1max && pop(i).cost(2)<=cost2max
    pop2(i)=pop(i);
    else
     id1=0;
 while id1==0
    pop2(i).pos=datasample(ndiscp,nvar,'Replace',false);
    pop2(i).cost=MyCost(pop2(i).pos,ndiscp,a,C0,C1,gamRR)';
   if pop2(i).cost(1) <=cost1max && pop2(i).cost(2)<=cost2max
       id1=1;
   end
 end
    end
end
pop=pop2;
end

 
   