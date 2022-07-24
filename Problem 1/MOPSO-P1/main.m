clc
clear
close

%% Feasible set
% Specify variogram model
a=30;
C0=0;
nu=0.2;


% Specify length of domain sides
xmax=100;
ymax=100;
VarMin=0;
VarMax=100;
% discretize domain: form a set of ndp
% random locations across the domain for numerical integration

ndp=1000;
ndiscp=Discretize(ndp,xmax,ymax);
lbc2=1; %a minimum distance of between any two points

gamRR=DispVar(ndiscp,a,C0,nu);

% Specify the number of points
N=50;
nVar=N;
VarSize=[2 nVar];
lb=zeros(1,nVar);
ub=100*ones(1,nVar);

%% MOPSO Parameters

MaxIt=500;           % Maximum Number of Iterations

nPop=40;            % Population Size

nRep=40;            % Repository Size

w=0.9;              % Inertia Weight
wdamp=0.99;         % Intertia Weight Damping Rate
c1=2;               % Personal Learning Coefficient
c2=1.5;               % Global Learning Coefficient
% c1=1;
% c2=1;

nGrid=7;            % Number of Grids per Dimension
alpha=0.15;          % Inflation Rate

beta=3;             % Leader Selection Pressure
gamma=2;            % Deletion Selection Pressure

mu=0.1;             % Mutation Rate
%% Initialization
tic
empty_particle.Position=[];
empty_particle.Velocity=[];
empty_particle.Cost=[];
empty_particle.Best.Position=[];
empty_particle.Best.Cost=[];
empty_particle.IsDominated=[];
empty_particle.GridIndex=[];
empty_particle.GridSubIndex=[];

pop=repmat(empty_particle,nPop,1);

for i=1:nPop
    pop(i).Position=Position(xmax,ymax,lbc2,nVar);
    pop(i).Cost=MyCost(pop(i).Position,ndiscp,a,C0,nu,gamRR)'; 
    pop(i).Velocity=zeros(VarSize);
    pop(i).Best.Position=pop(i).Position;
    pop(i).Best.Cost=pop(i).Cost;
end

pop=DetermineDomination(pop);

rep=pop(~[pop.IsDominated]);

Grid=CreateGrid(rep,nGrid,alpha);

for i=1:numel(rep)
    rep(i)=FindGridIndex(rep(i),Grid);
end

%% main loop
for it=1:MaxIt
    for i=1:nPop
       
       leader=SelectLeader(rep,beta);
       
         v1=pop(i).Velocity;
         v2=pop(i).Position;
         id3=0;
         jj=0;
        while id3==0
        
       pop(i).Velocity = w*pop(i).Velocity ...
            +c1*rand(VarSize).*(pop(i).Best.Position-pop(i).Position)' ...
            +c2*rand(VarSize).*(leader.Position-pop(i).Position)';
        
        pop(i).Velocity=min(max(pop(i).Velocity,[lb;lb]),[ub;ub]);
        
        pop(i).Position = pop(i).Position + pop(i).Velocity';
       
        pop(i).Position=min(max(pop(i).Position,[lb;lb]'),[ub;ub]');
        
        Ftem=min(pdist(pop(i).Position));
     
     if Ftem>lbc2 
          pop(i).Cost=MyCost(pop(i).Position,ndiscp,a,C0,nu,gamRR)';
          id3=1;
     else
         jj=jj+1;
         pop(i).Velocity=v1;
         pop(i).Position=v2;
           if jj==3
           pop(i).Cost=[10^6;10^6]; 
           id3=1;
           end
           
     end 
       end
            
             
          
       
     % Apply Mutation
        pm=(1-(it-1)/(MaxIt-1))^(1/mu);
        if rand<pm
            
           v3=pop(i).Position;
             id2=0;
             jj=0;
         while id2==0
        NewSol.Position=Mutate(pop(i).Position,pm,VarMin,VarMax);
        
        Ftem2=min(pdist(NewSol.Position));
     
     if Ftem2>lbc2 
         NewSol.Cost=MyCost(NewSol.Position,ndiscp,a,C0,nu,gamRR)';
         id2=1;
     else
         jj=jj+1;
         pop(i).Position=v3;
         if jj==3
         NewSol.Cost=[10^6;10^6];
         id2=1;
         end
     end 
        end
          
          
        
            
            if Dominates(NewSol,pop(i))
                pop(i).Position=NewSol.Position;
                pop(i).Cost=NewSol.Cost;

            elseif Dominates(pop(i),NewSol)
                % Do Nothing

            else
                if rand<0.5
                    pop(i).Position=NewSol.Position;
                    pop(i).Cost=NewSol.Cost;
                end
            end
        end
        
        if Dominates(pop(i),pop(i).Best)
            pop(i).Best.Position=pop(i).Position;
            pop(i).Best.Cost=pop(i).Cost;
            
        elseif Dominates(pop(i).Best,pop(i))
            % Do Nothing
            
        else
            if rand<0.5
                pop(i).Best.Position=pop(i).Position;
                pop(i).Best.Cost=pop(i).Cost;
            end
        end
     
     
    end

     % Add Non-Dominated Particles to REPOSITORY
    rep=[rep
         pop(~[pop.IsDominated])]; %#ok
    
    % Determine Domination of New Resository Members
    rep=DetermineDomination(rep);
    
    % Keep only Non-Dminated Memebrs in the Repository
    rep=rep(~[rep.IsDominated]);
     
    %Update Grid
    Grid=CreateGrid(rep,nGrid,alpha);

    % Update Grid Indices
    for i=1:numel(rep)
        rep(i)=FindGridIndex(rep(i),Grid);
    end
    
    % Check if Repository is Full
    if numel(rep)>nRep
        
        Extra=numel(rep)-nRep;
        for e=1:Extra
            rep=DeleteOneRepMemebr(rep,gamma);
        end
        
    end
    
    % Plot Costs
    %figure(1);
    %PlotCosts(pop,rep);
    %pause(0.01);
    
    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Number of Rep Members = ' num2str(numel(rep))]);
    
    % Damping Inertia Weight
    w=w*wdamp;
    
    
end

toc
t=toc
%% results
% 

[CL,~]=size(rep);
ktem=0;
indtem=0;
for i=1:CL
if min(pdist(rep(i).Position)<lbc2)
ktem=ktem+1;
indtem(ktem)=i;
end
end
if ktem>0
rep(indtem)=[];
[CL,~]=size(rep);
end

if(CL>nRep)
    disp('Strip');
% indices=sort(datasample(1:CL,HL-2,'Replace',false));
indices=strip([rep.Cost]',nRep-2);
%add best solution with respect to each OF to the set to retain
Ctem=zeros(length(rep),2);
for i=1:length(rep)
Ctem(i,:)=rep(i).Cost';
end
% [imax,~]=find(Ctem==max(Ctem(:,1)));
% [jmax,~]=find(Ctem==max(Ctem(:,2)));
[imin,~]=find(Ctem==min(Ctem(:,1)));
[jmin,~]=find(Ctem==min(Ctem(:,2)));
%indices=[indices,imax(1),jmax(1)];
indices=[indices,imin(1),jmin(1)];
CL=nRep;
rep=rep(indices);
end
Ctem=zeros(length(rep),2);
for i=1:length(rep)
Ctem(i,:)=rep(i).Cost';
end


C=Ctem;
%C=[rep.Cost]';
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
Mid=Mid/n

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

Mid
SM
DM
SNS

for i=1:length(rep)
    figure(i)
    plot(rep(i).Position(:,1), rep(i).Position(:,2),'.k')
    xlim([0 100])
    ylim([0 100])
end
 
for i=1:40
    AA(i)=rep(i).Cost(1);
end
    
find(AA==min(AA))
find(AA==max(AA))


figure(1)
PlotCosts(pop,rep)
%PlotCosts(pop,Ctem')
%xlim([0 350])
%ylim([150 500])
% grid off
 xlabel({'Variance of mean'})
 ylabel({'Total distance /m'})
% set(gca, 'XTick', [0 50 100 150 200 250 300],'YTick', [150 200 250 300 350 400 450 500])


