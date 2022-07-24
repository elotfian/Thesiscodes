clear
clc


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
Varmima=[VarMin,VarMax];


% discretize domain: form a set of ndp
% random locations across the domain for numerical integration
ndp=1000;
[ndiscp,samplev]=Discretize(ndp,xmax,ymax);
lbc2=1; %a minimum distance of between any two points
%plot(ndiscp(:,1),ndiscp(:,2),'ko');



% Specify the number of points
N=50;
nvar=N;
lb=zeros(1,nvar);
ub=100*ones(1,nvar);
%% NSGA-II Parameters
npop=40;     % number of population


pc=0.7;       % percent of crossover
ncross=2*round(npop*pc/2);  % number of crossover offspring

pm=0.15;        %  percent of mutation
nmut=round(npop*pm);  % number of mutation offspring


maxiter=2000;

%% initialization
tic
emp.pos=[];
emp.cost=[];
emp.rank=[];      
emp.cdis=[];      % crowding distance


pop=repmat(emp,npop,1);

for i=1:npop
    i
    pop(i).pos=Position(xmax,ymax,lbc2,nvar);
    pop(i).cost=MyCost(pop(i).pos,samplev,a,C0,nu)';

end

[pop,F]=non_dominated_sorting(pop);
pop=crowding_distance(pop,F);
pop=sorting(pop);
%% main loop

for iter=1:maxiter

    % crossover
     crosspop=repmat(emp,ncross,1);
     crosspop=crossover(crosspop,pop,ncross,nvar,samplev,a,C0,nu,lbc2);
                     
     mutpop=repmat(emp,nmut,1);
     mutpop=mutation(mutpop,pop,nmut,nvar,samplev,a,C0,nu,lbc2,xmax,ymax);
                           
     
     [pop]=[pop;crosspop;mutpop];
    
     [pop,F]=non_dominated_sorting(pop);
      pop=crowding_distance(pop,F);
      pop=sorting(pop);
      pop=pop(1:npop);
      
      [pop,F]=non_dominated_sorting(pop);
      pop=crowding_distance(pop,F);
      pop=sorting(pop);
      
      
      C=[pop.cost]';
      
      
      
      
      disp([ ' iter =   '  num2str(iter) ' N Pareto = '  num2str(length(F{1})) ])

end
toc
t=toc
%% results

pareto=pop(F{1},:);
[CL,~]=size(pareto);
ktem=0;
indtem=0;
for i=1:CL
if min(pdist(pareto(i).pos)<lbc2)
ktem=ktem+1;
indtem(ktem)=i;
end
end
if ktem>0
pareto(indtem)=[];
[CL,~]=size(pareto);
end

if(CL>npop)
    disp('Strip');
indices=strip([pareto.cost]',npop-2);
%add best solution with respect to each OF to the set to retain
Ctem=zeros(length(pareto),2);
for i=1:length(pareto)
Ctem(i,:)=pareto(i).cost';
end
[imin,~]=find(Ctem==min(Ctem(:,1)));
[jmin,~]=find(Ctem==min(Ctem(:,2)));
indices=[indices,imin(1),jmin(1)];
indices=unique(indices);
CL=npop;
pareto=pareto(indices);
end
Ctem=zeros(length(pareto),2);
for i=1:length(pareto)
Ctem(i,:)=pareto(i).cost';
end

C=Ctem;


[n,~]=size(C(:,1));
f1best=min(C(:,1));
f2best=min(C(:,2));
f1min=min(C(:,1));
f2min=min(C(:,2));
f1max=max(C(:,1));
f2max=max(C(:,2));

%Mean Ideal Distance (MID)
Mid=0;
for i=1:n
   Mid=Mid+sqrt(((C(i,1)-f1best)/(f1max-f1min))^2+((C(i,2)-f2best)/(f2max-f2min))^2);    
end
Mid=Mid/n

%Spacing Metric (SM)
row=size(C,1);
for i=1:size(C,1)
    DD1=repmat(C(i,:),row,1);
    DD2=DD1-C;
    for j=1:row
        if (j==i)
        DD3(j)=10^6;
        else
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
%Diversification Metric (DM)
DM=sqrt((f1min-f1max)^2+(f2min-f2max)^2)
%%%%%%%%%%%%
%Spread of Non-dominance Solution (SNS)
CC=zeros(size(C,1),1);
SN1=0;
for i=1:size(C,1)
    CC(i)=sqrt(C(i,1)^2+C(i,2)^2);
    SN1=SN1+(Mid-CC(i))^2;
end
SNS=sqrt(SN1/(size(C,1)-1))


figure(1)
plotpareto(F,C)
%xlim([0 350])
%ylim([150 500])
xlabel({'Mean total error'})
ylabel({'Total distance /m'})
%set(gca, 'XTick', [0:50:300],'YTick', [150 :50:500])



Mid
SM
DM
SNS

for i=1:length(pareto)
    AA(i)=pareto(i).cost(1);
   end
find(AA==min(AA))
find(AA==max(AA))

for i=1:length(pareto)
    figure(i)
    plot(pareto(i).pos(:,1), pareto(i).pos(:,2),'.k')
    xlim([0 100])
    ylim([0 100])
end
 





