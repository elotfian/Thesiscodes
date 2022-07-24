clear
clc

% Specify variogram model
aa=30;
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
ndiscp=Discretize(ndp,xmin,ymin,xmax,ymax);
lbc2=1; %a minimum distance of between any two points
%plot(ndiscp(:,1),ndiscp(:,2),'ko');

gamRR=DispVar(ndiscp,aa,C0,nu);

%%
% Specify the number of points
N=50;
nvar=N;

% Lower bound and upper bound
lb=zeros(1,nvar);
ub=100*ones(1,nvar);

GreyWolves_num=40;  %Population size (N-pop)
MaxIt=500;  % Maximum Number of Iterations
Archive_size=40;   % Repository Size

alpha=0.15;  % Grid Inflation Parameter
nGrid=7;   % Number of Grids per each Dimension
beta=1; %=4;    % Leader Selection Pressure Parameter
gamma=2;    % Extra (to be deleted) Repository Member Selection Pressure

% Initialization
tic
GreyWolves=CreateEmptyParticle(GreyWolves_num);


for i=1:GreyWolves_num
    GreyWolves(i).Velocity=0;
    GreyWolves(i).Position=Position(xmin,ymin,xmax,ymax,lbc2,nvar);
    GreyWolves(i).Cost=MyCost(GreyWolves(i).Position,ndiscp,aa,C0,nu,gamRR);
    GreyWolves(i).Best.Position=GreyWolves(i).Position;
    GreyWolves(i).Best.Cost=GreyWolves(i).Cost;
end

GreyWolves=DetermineDomination(GreyWolves);

Archive=GetNonDominatedParticles(GreyWolves);

Archive_costs=GetCosts(Archive);
G=CreateHypercubes(Archive_costs,nGrid,alpha);

for i=1:numel(Archive)
    [Archive(i).GridIndex, Archive(i).GridSubIndex]=GetGridIndex(Archive(i),G);
end

%% MOGWO main loop

for it=1:MaxIt
    a=2-it*((2)/MaxIt);
    for i=1:GreyWolves_num
       
       
        jj=0;
        id3=0;
        while id3==0
        clear rep2
        clear rep3
        
        % Choose the alpha, beta, and delta grey wolves
        Delta=SelectLeader(Archive,beta);
        Beta=SelectLeader(Archive,beta);
        Alpha=SelectLeader(Archive,beta);
        
        % If there are less than three solutions in the least crowded
        % hypercube, the second least crowded hypercube is also found
        % to choose other leaders from.
        if size(Archive,1)>1
            counter=0;
            for newi=1:size(Archive,1)
                if sum(Delta.Position~=Archive(newi).Position)~=0
                    counter=counter+1;
                    rep2(counter,1)=Archive(newi);
                end
            end
            Beta=SelectLeader(rep2,beta);
        end
        
        % This scenario is the same if the second least crowded hypercube
        % has one solution, so the delta leader should be chosen from the
        % third least crowded hypercube.
        if size(Archive,1)>2
            counter=0;
            for newi=1:size(rep2,1)
                if sum(Beta.Position~=rep2(newi).Position)~=0
                    counter=counter+1;
                    rep3(counter,1)=rep2(newi);
                end
            end
            Alpha=SelectLeader(rep3,beta);
        end
        
        % Eq.(3.4) in the paper
        c=2.*rand(1, nvar);
        % Eq.(3.1) in the paper
        D=abs([(c.*Delta.Position(:,1)')' ,(c.*Delta.Position(:,2)')']-GreyWolves(i).Position);
        % Eq.(3.3) in the paper
        A=2.*a.*rand(1, nvar)-a;
        % Eq.(3.8) in the paper
        X1=Delta.Position-[(A.*abs(D(:,1)'))',(A.*abs(D(:,2)'))'];
        
        
        % Eq.(3.4) in the paper
        c=2.*rand(1, nvar);
        % Eq.(3.1) in the paper
        D=abs([(c.*Beta.Position(:,1)')' ,(c.*Beta.Position(:,2)')']-GreyWolves(i).Position);
        % Eq.(3.3) in the paper
        A=2.*a.*rand(1, nvar)-a;
        % Eq.(3.9) in the paper
        X2=Beta.Position-[(A.*abs(D(:,1)'))',(A.*abs(D(:,2)'))'];
        
        
        % Eq.(3.4) in the paper
        c=2.*rand(1, nvar);
        % Eq.(3.1) in the paper
        D=abs([(c.*Alpha.Position(:,1)')' ,(c.*Alpha.Position(:,2)')']-GreyWolves(i).Position);
        % Eq.(3.3) in the paper
        A=2.*a.*rand(1, nvar)-a;
        % Eq.(3.10) in the paper
        X3=Alpha.Position-[(A.*abs(D(:,1)'))',(A.*abs(D(:,2)'))'];
        
        % Eq.(3.11) in the paper
         GreyWolves(i).Position=(X1+X2+X3)./3;
  %        GreyWolves(i).Position=X3;
        

        % Boundary checking
        GreyWolves(i).Position=min(max(GreyWolves(i).Position,[lb;lb]'),[ub;ub]'); 
       
       Ftem=min(pdist(GreyWolves(i).Position));
       if Ftem>lbc2
       GreyWolves(i).Cost=MyCost(GreyWolves(i).Position,ndiscp,aa,C0,nu,gamRR);
       id3=1;
       else 
              jj=jj+1;
              if jj==3
               GreyWolves(i).Cost=[10^6,10^6];
               id3=1;
              end
       end        
       end  

    end
    
    GreyWolves=DetermineDomination(GreyWolves);
    non_dominated_wolves=GetNonDominatedParticles(GreyWolves);
    
    Archive=[Archive
        non_dominated_wolves];
    
    Archive=DetermineDomination(Archive);
    Archive=GetNonDominatedParticles(Archive);
    
    for i=1:numel(Archive)
        [Archive(i).GridIndex, Archive(i).GridSubIndex]=GetGridIndex(Archive(i),G);
    end
    
    if numel(Archive)>Archive_size
        EXTRA=numel(Archive)-Archive_size;
        Archive=DeleteFromRep(Archive,EXTRA,gamma);
        
        Archive_costs=GetCosts(Archive);
        G=CreateHypercubes(Archive_costs,nGrid,alpha);
        
    end
    
    disp(['In iteration ' num2str(it) ': Number of solutions in the archive = ' num2str(numel(Archive))]);
    save results
    
 
    
end
toc
t=toc

 % Results
    

[CL,~]=size(Archive);
ktem=0;
indtem=0;
for i=1:CL
if min(pdist(Archive(i).Position)<lbc2)
ktem=ktem+1;
indtem(ktem)=i;
end
end
if ktem>0
Archive(indtem)=[];
[CL,~]=size(Archive);
end

if(CL>Archive_size)
    disp('Strip');
% indices=sort(datasample(1:CL,HL-2,'Replace',false));
indices=strip([Archive.Cost]',Archive_size-2);
%add best solution with respect to each OF to the set to retain
Ctem=zeros(length(Archive),2);
for i=1:length(Archive)
Ctem(i,:)=Archive(i).Cost';
end
% [imax,~]=find(Ctem==max(Ctem(:,1)));
% [jmax,~]=find(Ctem==max(Ctem(:,2)));
[imin,~]=find(Ctem==min(Ctem(:,1)));
[jmin,~]=find(Ctem==min(Ctem(:,2)));
%indices=[indices,imax(1),jmax(1)];
indices=[indices,imin(1),jmin(1)];
indices=unique(indices);
CL=Archive_size;
Archive=Archive(indices);
end
Ctem=zeros(length(Archive),2);
for i=1:length(Archive)
Ctem(i,:)=Archive(i).Cost';
end
 


%     costs=GetCosts(GreyWolves);
%     Archive_costs=GetCosts(Archive);
    
%C=[Archive_costs]';
C=Ctem;
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
    
        figure(1)
        plot(Ctem(:,1),Ctem(:,2),'ko','markerfacecolor','k','MarkerSize',6)
%         xlim([0 1])
%         ylim([200 1000])
        xlabel({'Variance of mean'})
        ylabel({'Total distance /m'})
       % set(gca, 'XTick', [0: 50: 300],'YTick',[150: 50: 500]);
       
Mid 
SM 
DM 
SNS
        
  for i=1:length(Archive)
    AA(i)=Archive(i).Cost(1);
end
    
find(AA==min(AA))
find(AA==max(AA))


 for i=1:length(Archive)
     figure(i)
     plot(Archive(i).Position(:,1), Archive(i).Position(:,2),'.k')
     xlim([0 100])
     ylim([0 100])
 end




% figure(2)
% subplot(2,2,1)
% plot(Archive(14).Position(:,1),Archive(14).Position(:,2),'ko','markerfacecolor','k','MarkerSize',4)
% xlim([0 100])
% ylim([0 100])
% %set(gca, 'XTick', [0 20 40 60 80 100],'YTick', [0 20 40 60 80 100])
% set(gcf,'Position',[10 10 900 900])
% subplot(2,2,2)
% plot(Archive(17).Position(:,1),Archive(17).Position(:,2),'ko','markerfacecolor','k','MarkerSize',4)
% %39
% xlim([0 100])
% ylim([0 100])
% %set(gca, 'XTick', [0 20 40 60 80 100],'YTick', [0 20 40 60 80 100])
% set(gcf,'Position',[10 10 900 900])
% subplot(2,2,3)
% plot(Archive(3).Position(:,1),Archive(3).Position(:,2),'ko','markerfacecolor','k','MarkerSize',4)
% %31
% xlim([0 100])
% ylim([0 100])
% %set(gca, 'XTick', [0 20 40 60 80 100],'YTick', [0 20 40 60 80 100])
% set(gcf,'Position',[10 10 900 900])
% subplot(2,2,4)
% plot(Archive(5).Position(:,1),Archive(5).Position(:,2),'ko','markerfacecolor','k','MarkerSize',4)
% xlim([0 100])
% ylim([0 100])
% %set(gca, 'XTick', [0 20 40 60 80 100],'YTick', [0 20 40 60 80 100])
% set(gcf,'Position',[10 10 900 900])
% save results

    


