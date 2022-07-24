clc;
clear;
close all;

%% Feasible set
% Specify variogram model
a=100;
C0=500;
C1=4500;


% Specify length of domain sides
xmax=100;
ymax=100;
VarMin=0;
VarMax=100;
cost1max=300;
cost2max=500;
% Provide coordinates of pond centre, and radius

pcen=[100,100];
prad=40;

% discretize domain: form a set of ndp
% random locations across the domain for numerical integration

ndp=1000;
ndiscp=Discretize(ndp,pcen,prad,xmax,ymax);
%plot(ndiscp(:,1),ndiscp(:,2),'ko');

gamRR=DispVar(ndiscp,a,C0,C1);

% Specify the number of points
N=20;
nVar=N;
VarSize=[2 nVar];
lb=zeros(1,nVar);
ub=100*ones(1,nVar);

%% MOPSO Parameters

MaxIt=2000;           % Maximum Number of Iterations

nPop=20;            % Population Size

nRep=40;            % Repository Size

w=0.3;              % Inertia Weight
wdamp=0.99;         % Intertia Weight Damping Rate
c1=1;               % Personal Learning Coefficient
c2=2;               % Global Learning Coefficient

nGrid=7;            % Number of Grids per Dimension
alpha=0.1;          % Inflation Rate

beta=2;             % Leader Selection Pressure
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
    id1=0;
    while id1==0
    pop(i).Position=datasample(ndiscp,nVar,'Replace',false);
    pop(i).Cost=MyCost(pop(i).Position,ndiscp,a,C0,C1,gamRR)';
    if pop(i).Cost(1) <=cost1max && pop(i).Cost(2)<=cost2max
    id1=1;
    end
    end
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

propor=zeros(1,MaxIt);
%% main loop
for it=1:MaxIt
    for i=1:nPop
       leader=SelectLeader(rep,beta);
       
       id3=0;
       while id3==0
        id=0;
        while id==0
        pop(i).Velocity = w*pop(i).Velocity ...
            +c1*rand(VarSize).*(pop(i).Best.Position-pop(i).Position)' ...
            +c2*rand(VarSize).*(leader.Position-pop(i).Position)';
        pop(i).Position = pop(i).Position + pop(i).Velocity';
       
       [nn,~]=size(pop(i).Position); 
       id2=0;
       for j=1:nn
       dpond=0;
       dpond=norm(pop(i).Position(j,:)-pcen);
       if pop(i).Position(j,1)<VarMin
           id2=1;
       elseif pop(i).Position(j,2)<VarMin
            id2=1;
       elseif pop(i).Position(j,1)>VarMax
           id2=1;
       elseif pop(i).Position(j,2)>VarMax
           id2=1;
       elseif dpond<=prad
           id2=1;
       end
       end
       if id2==0
           id=1;
       end    
       end
     pop(i).Cost=MyCost(pop(i).Position,ndiscp,a,C0,C1,gamRR)'; 
     if pop(i).Cost(1)<=cost1max && pop(i).Cost(2)<=cost2max
         id3=1;
     end
       end
     
     % Apply Mutation
        pm=(1-(it-1)/(MaxIt-1))^(1/mu);
        if rand<pm
           
            id3=0;
       while id3==0
        id=0;
        while id==0
            NewSol.Position=Mutate(pop(i).Position,pm,VarMin,VarMax);
        [nn,~]=size(NewSol.Position); 
        id2=0;    
       for j=1:nn
       dpond=0;
       dpond=norm(NewSol.Position(j,:)-pcen);
       if NewSol.Position(j,1)<VarMin
           id2=1;
       elseif NewSol.Position(j,2)<VarMin
            id2=1;
       elseif NewSol.Position(j,1)>VarMax
           id2=1;
       elseif NewSol.Position(j,2)>VarMax
           id2=1;
       elseif dpond<=prad
           id2=1;
       end
       end
       
       
       
       
       %%%%%
       if id2==0
           id=1;
      
       end    
       end
          NewSol.Cost=MyCost(NewSol.Position,ndiscp,a,C0,C1,gamRR)';
           if NewSol.Cost(1)<=cost1max && NewSol.Cost(2)<=cost2max
           id3=1;
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

figure(1)
PlotCosts(pop,rep)
xlim([0 250])
ylim([0 500])
%xlabel({'Variance of mean'})
%ylabel({'Total distance /m'})
%set(gca, 'XTick', [0 50 100 150 200 250],'YTick', [0 100 200 300 400 500])


C=[rep.Cost]';
[n,~]=size(C(:,1));
f1best=max(C(:,1));
f2best=max(C(:,2));
f1min=min(C(:,1));
f2min=min(C(:,2));

Mid=0;
for i=1:n
   Mid=Mid+sqrt(((C(i,1)-f1best)/(f1best-f1min))^2+((C(i,2)-f2best)/(f2best-f2min))^2);    
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
DM=sqrt((f1min-f1best)^2+(f2min-f2best)^2)
%%%%%%%%%%%%
CC=zeros(size(C,1),1);
SN1=0;
for i=1:size(C,1)
    CC(i)=sqrt(C(i,1)^2+C(i,2)^2);
    SN1=SN1+(Mid-CC(i))^2;
end
SNS=sqrt(SN1/(size(C,1)-1))