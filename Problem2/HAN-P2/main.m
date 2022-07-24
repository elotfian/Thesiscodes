clc;
clear;
close all;

%% Feasible set
% Specify variogram model
a=30;
C0=0;
nu=0.2;


% Specify length of domain sides
xmax=100;
ymax=100;
xmin=0;
ymin=0;
VarMin=0;
VarMax=100;
% Provide coordinates of pond centre, and radius
% discretize domain: form a set of ndp
% random locations across the domain for numerical integration
% discretize domain: form a set of ndp
% random locations across the domain for numerical integration
ndp=1000;
[ndiscp,samplev]=Discretize(ndp,xmax,ymax);
lbc2=1; %a minimum distance of between any two points
%plot(ndiscp(:,1),ndiscp(:,2),'ko');



% Specify the number of points
nvar=50;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%lb=zeros(1,nvar);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ub=100*ones(1,nvar);

%% Archive size limits and initialization at length CL

SL=70;  % Soft limit on archive size
HL=40;  % Hard limit on archive size	

CL=70;   %Current archive size 
M=2;    %Number of objective functions
maxiter=2000;
pc=0.8;       % percent of crossover
ncross=2*round(HL*pc/2);  % number of crossover offspring

pm=0.2;        %  percent of mutation
nmut=round(HL*pm);  % number of mutation offspring

% Arc is a matrix which holds the initial values of the objective functions
% for the current solutions. 

emp.pos=[];
emp.cost=[];
emp2.pos=[];
emp2.cost=[];
emp2.rank=[];      
emp2.cdis=[];
Arc=repmat(emp,CL,1);
pop=repmat(emp2,CL,1);

for i=1:CL
  i
    Arc(i).pos=Position(xmax,ymax,lbc2,nvar);
    Arc(i).cost=MyCost(Arc(i).pos,samplev,a,C0,nu)';  
end


for i=1:CL
    pop(i).pos=Arc(i).pos;
    pop(i).cost=Arc(i).cost;
end
[pop,F]=non_dominated_sorting(pop);
pop=crowding_distance(pop,F);
pop=sorting(pop);
      
      



%use the first layout to initial Arc
[AArc,F]=non_dominated_sorting(Arc);
sizeF1=size(F{1},2);
AArc=crowding_distance(AArc,F);
AArc=sorting(AArc);
AArc=AArc(1:sizeF1);
[CL,~]=size(AArc);
Arc=repmat(emp,CL,1);
for i=1:CL
    Arc(i).pos=AArc(i).pos;
    Arc(i).cost=AArc(i).cost;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(CL>HL)
Ctem=[Arc.cost]';
indices=strip(Ctem,HL-2);
%add best solution with respect to each OF to the set to retain
[imin,~]=find(Ctem==min(Ctem(:,1)));
[jmin,~]=find(Ctem==min(Ctem(:,2)));
indices=[indices,imin(1),jmin(1)];
indices=unique(indices);
Arc=Arc(indices);
end
[CL,~]=size(Arc);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5



propacc=zeros(1,maxiter); %Vector to store proportion of dominated 
   					   %changes accepted in each iteration

% Cooling schedule

tau=1.0;  	%Initial temperature
alp=0.95;    %Change in temperature at the end of each iteration
maxjump=5;	%Maximum distance to move a sample point
     
%R=[cost1max,cost2max];
xn=zeros(1,M)';

%pick a point from the archive at random 

nxc=randsample(CL,1);
xc=Arc(nxc).cost;  
xyc=Arc(nxc).pos;  % current solution (coordinates of sample points)

%% Start of the SA iterations.  To run, select from the
% following line to "This is the end of the SA iterations"
% below.

for iter=1:maxiter

nrandom=0;
nacc=0;
disp(['Iteration ' num2str(iter) ': Number of Rep Members = ' num2str(length(Arc))]);
for j=1:nvar

       jj=0;
       
       while true
%             
       jump = [unifrnd(-maxjump,maxjump),unifrnd(-maxjump,maxjump)];
       prop=xyc(j,:)+jump;
       prop(1)=min(max(xmin,prop(1)),xmax);
       prop(2)=min(max(ymin,prop(2)),ymax);
       
       xycheck=xyc;
       xycheck(j,:)=xycheck(j,:)+jump;
       if min(pdist(xycheck))>lbc2
       FF=MyCost(xycheck,samplev,a,C0,nu)';     
       break;
       else 
              jj=jj+1;
              if jj==3
                 FF=[10^6,10^6]';
                break
              end
       end        
       end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Proposed change is permissible.  Now evaluate 
        % dominance (function Dom)with respect to current point and archive
        % points, and accept or reject and modify the archive 
        % accordingly according to the rules of Bandyopadhyay et al
        xyn=xycheck;  %vector for proposed new solution
        xn=FF;
        
        R=computR(xn,xc,Arc);
        [op,klist]=Dom(xn,xc,Arc,tau,R,CL);
        
         if(op(1)==1)
          %Bandyopadhyay et al. Case 1
          nrandom=nrandom+1;
          pacc=op(2);
             if(unifrnd(0,1)<pacc)
             xc=xn;
             xyc=xyn;
             nacc=nacc+1;
             end
         end
 
if(op(1)==2)
%Bandyopadhyay et al. Case 2
        if(op(2)==1)
               nrandom=nrandom+1;
               pacc=op(3);
                 if(unifrnd(0,1)<pacc)
                    xc=xn;
                    xyc=xyn;
                    nacc=nacc+1;
                  end
        else
            if(op(4)>0)
               newsol.cost=xn;
               newsol.pos=xyn;
               %Arc(op(5:(op(4)+4)))=[];
               Arc(klist)=[];
               Arc=[Arc
                          newsol];
               [CL,~]=size(Arc);
             else
                newsol.cost=xn;
                newsol.pos=xyn;
                Arc=[Arc
                          newsol];
                [CL,~]=size(Arc);
            end
        end

end

if(op(1)==3)
%Bandyopadhyay et al. Case 3
if(op(2)==1)
%Case 3.1 
imax=length(op)-2;
pacc=op(imax-1);
kmindo=op(imax);
if(unifrnd(0,1)<pacc)
xc=Arc(kmindo).cost;
xyc=Arc(kmindo).pos;
end
end

if(op(2)==2)
%Case 3.2 
xc=xn;
xyc=xyn;

%remove the (new) current point from the archive if already in there

lop=length(op);
if(op(lop-1)==1)
xcrow=op(lop);
Arc(xcrow)=[];
[CL,~]=size(Arc);
end

% add current point to archive (previous step avoids duplication)
newsol.cost=xn;
newsol.pos=xyn;
Arc=[Arc
    newsol];
[CL,~]=size(Arc);

end
if(op(2)==3)
%Case 3.3 
newsol.cost=xn;
newsol.pos=xyn;
Arc(klist)=[];
Arc=[Arc
    newsol];
[CL,~]=size(Arc);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% thin the archive if its length exceeds the soft limit
%
if(CL>SL)
    disp('Strip');
    
% indices=sort(datasample(1:CL,HL-2,'Replace',false));
clear indices Ctem

Ctem=[Arc.cost]';
indices=strip(Ctem,HL-2);
%add best solution with respect to each OF to the set to retain
[imin,~]=find(Ctem==min(Ctem(:,1)));
[jmin,~]=find(Ctem==min(Ctem(:,2)));
indices=[indices,imin(1),jmin(1)];
indices=unique(indices);
Arc=Arc(indices);
end
[CL,~]=size(Arc);

end


 Arctem=Arc;
 for i=1: size(Arc,1)
    Arctem(i).rank=[];
    Arctem(i).cdis=[];
 end

% crossover
crosspop=repmat(emp2,ncross,1);
crosspop=crossover(crosspop,pop,ncross,nvar,samplev,a,C0,nu,lbc2);

mutpop=repmat(emp2,nmut,1);
mutpop=mutation(mutpop,pop,nmut,nvar,samplev,a,C0,nu,lbc2,xmax,ymax);
[pop]=[pop;crosspop;mutpop;Arctem];
    
     [pop,F]=non_dominated_sorting(pop);
      pop=crowding_distance(pop,F);
      pop=sorting(pop);
      pop=pop(1:HL);
      [pop,F]=non_dominated_sorting(pop);
      pop=crowding_distance(pop,F);
      pop=sorting(pop);


%  Camosa=[Arctem.cost]';
%  CNsagaii=[pop.cost]';
% 
% figure(1)
% plot(Camosa(:,1),Camosa(:,2),'r*')
% hold on
% plot(CNsagaii(:,1),CNsagaii(:,2),'bo')
% hold on
% plot(Cnsg100(:,1),Cnsg100(:,2),'ko','markerfacecolor','k','MarkerSize',6)
% hold off
% pause(0.01);


if(iter==1) 
    mtem=nacc/nrandom;
display(['' num2str(mtem)]);
end

propacc(iter)=nacc/nrandom;

tau=tau*alp;

%One iteration completed
end
% This is the end of the SA iterations
% Plot the proportion of dominated changes for successive iterations
% figure(1);
% plot(Ctem(:,1),Ctem(:,2),'o')
%  xlim([0 350])
%  ylim([0 500])
% figure(2)
% plot([1:iter],propacc)

% if necessary, strip the final set of solutions to a set length HL

[CL,~]=size(Arc);
ktem=0;
indtem=0;
for i=1:CL
if min(pdist(Arc(i).pos)<lbc2)
ktem=ktem+1;
indtem(ktem)=i;
end
end
if ktem>0
Arc(indtem)=[];
[CL,~]=size(Arc);
end

if(CL>HL)
    disp('Strip');
% indices=sort(datasample(1:CL,HL-2,'Replace',false));
clear indices Ctem

Ctem=[Arc.cost]';
indices=strip(Ctem,HL-2);
%add best solution with respect to each OF to the set to retain
[imin,~]=find(Ctem==min(Ctem(:,1)));
[jmin,~]=find(Ctem==min(Ctem(:,2)));
indices=[indices,imin(1),jmin(1)];
indices=unique(indices);
Arc=Arc(indices);
end
clear Ctem
Ctem=zeros(length(Arc),2);
for i=1:length(Arc)
Ctem(i,:)=Arc(i).cost';
end

C=Ctem;

 for i=1: size(Arc,1)
    Arc(i).rank=[];
    Arc(i).cdis=[];
end

Arc3=[Arc;pop];

[CL3,~]=size(Arc3);
ktem=0;
indtem=0;
for i=1:CL3
if min(pdist(Arc3(i).pos)<lbc2)
ktem=ktem+1;
indtem(ktem)=i;
end
end
if ktem>0
Arc3(indtem)=[];
[CL3,~]=size(Arc3);
end

if(CL3>HL)
    disp('Strip');
Ctem3=[Arc3.cost]';    
indices=0;
indices=strip(Ctem3,HL-2);
[imin,~]=find(Ctem3==min(Ctem3(:,1)));
[jmin,~]=find(Ctem3==min(Ctem3(:,2)));
indices=[indices,imin(1),jmin(1)];
indices=unique(indices);
Arc3=Arc3(indices);
end


[AArc3,F1]=non_dominated_sorting(Arc3);
sizeF1=size(F1{1},2);
AArc3=crowding_distance(AArc3,F1);
AArc3=sorting(AArc3);
AArc3=AArc3(1:min(sizeF1,HL));

[CL3,~]=size(AArc3);
Arc3=repmat(emp,CL3,1);
for i=1:CL3
    Arc3(i).pos=AArc3(i).pos;
    Arc3(i).cost=AArc3(i).cost;
end

Ctem3=[Arc3.cost]';

figure(2)
plot(Ctem(:,1),Ctem(:,2),'r*')
hold on
plot(Ctem3(:,1),Ctem3(:,2),'bo')

C=0;
C=Ctem3;

[n,~]=size(C(:,1));
f1best=min(C(:,1));
f2best=min(C(:,2));
f1min=min(C(:,1));
f2min=min(C(:,2));
f1max=max(C(:,1));
f2max=max(C(:,2));

Mid=0;
for i=1:n
   Mid=Mid+sqrt(((C(i,1)-f1best)/(f1max-f1min))^2+((C(i,2)-f2best)/(f2max-f2min))^2);    
end
Mid=Mid/n;

row=size(C,1);
for i=1:size(C,1)
    DD1=repmat(C(i,:),row,1);
    DD2=DD1-C;
    for j=1:row
        if (j~=i)
        DD3(j)=norm(DD2(j,:),1);
        end
    end
    DDistance(i)=min(DD3);
      
end

d_average=(1/size(C,1))*sum(DDistance);

total_distance2=0;
for i=1:size(C,1)
    total_distance2=total_distance2+abs((DDistance(i)-d_average));
end

SM=total_distance2/((size(C,1)-1)*d_average)

%%%%%%%%%%%%%%%%%%%%%
DM=sqrt((f1min-f1max)^2+(f2min-f2max)^2)
%%%%%%%%%%%%
CC=zeros(size(C,1),1);
SN1=0;
for i=1:size(C,1)
    CC(i)=sqrt(C(i,1)^2+C(i,2)^2);
    SN1=SN1+(Mid-CC(i))^2;
end
SNS=sqrt(SN1/(size(C,1)-1))


