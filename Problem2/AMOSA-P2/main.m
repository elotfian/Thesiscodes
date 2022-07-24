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

ndp=1000;
[ndiscp,samplev]=Discretize(ndp,xmax,ymax);
lbc2=1; %a minimum distance of between any two points
%plot(ndiscp(:,1),ndiscp(:,2),'ko');



% Specify the number of points
nvar=50;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%lb=zeros(1,nvar);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ub=100*ones(1,nvar);

%% Archive size limits and initialization at length CL

SL=60;  % Soft limit on archive size
HL=40;  % Hard limit on archive size	

CL=5;   %Current archive size 
M=2;    %Number of objective functions
maxiter=2000;

% Arc is a matrix which holds the initial values of the objective functions
% for the current solutions. 

emp.pos=[];
emp.cost=[];
Arc=repmat(emp,CL,1);

for i=1:CL
    i
    Arc(i).pos=Position(xmax,ymax,lbc2,nvar);
    Arc(i).cost=MyCost(Arc(i).pos,samplev,a,C0,nu)';                
end

[AArc,F]=non_dominated_sorting(Arc);
AArc=crowding_distance(AArc,F);
AArc=sorting(AArc);
AArc=AArc(F{1},:);
[CL,~]=size(AArc);
Arc=repmat(emp,CL,1);
for i=1:CL
    Arc(i).pos=AArc(i).pos;
    Arc(i).cost=AArc(i).cost;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(CL>HL)
indices=strip([Arc.cost]',HL-2);
%add best solution with respect to each OF to the set to retain
Ctem=zeros(length(Arc),2);
for i=1:length(Arc)
Ctem(i,:)=Arc(i).cost';
end
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
alp=0.99;    %Change in temperature at the end of each iteration
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
       id3=0;
       while id3==0
%             
       jump = [unifrnd(-maxjump,maxjump),unifrnd(-maxjump,maxjump)];
       prop=xyc(j,:)+jump;
       prop(1)=min(max(xmin,prop(1)),xmax);
       prop(2)=min(max(ymin,prop(2)),ymax);
       
       if min(pdist2(prop, xyc))>lbc2
       xycheck=xyc;
       xycheck(j,:)=xycheck(j,:)+jump;
       FF=MyCost(xycheck,samplev,a,C0,nu);           
       id3=1;
       else 
              jj=jj+1;
              if jj==3
                xycheck=xyc;
                xycheck(j,:)=xycheck(j,:)+jump;  
                FF=[10^6,10^6];
                id3=1;
              end
       end        
       end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Proposed change is permissible.  Now evaluate 
        % dominance (function Dom)with respect to current point and archive
        % points, and accept or reject and modify the archive 
        % accordingly according to the rules of Bandyopadhyay et al
        xyn=xyc;  %vector for proposed new solution
        xyn(j,:)=xycheck(j,:); % apply change 
        xn=FF';
        
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

% strip out any solutions that fall outside specified limits on each
% objective
% ktem=0;
% indtem=0;
% for i=1:CL
% if Arc(i).cost(1)>cost1max || Arc(i).cost(2)>cost2max
% ktem=ktem+1;
% indtem(ktem)=i;
% end
% end
% if ktem>0
% Arc(indtem)=[];
% [CL,~]=size(Arc);
% end

% check the archive against the initial points, remove	
% any archive points dominated by initial ones and put the initial point
% back in the archive


% checklist=ones(1,CL); %#list set to 1
% arclist=zeros(1,3);
% kn=0;
% for jj=1:3  % set up here with three initial points
% xin=urArc(jj).cost;
% for i=1:CL
% dv=ddom(xin,Arc(i).cost,R);
% if(dv(1)==1)
% if(sum(abs(xin-Arc(i).cost))~=0)
% kn=kn+1;
% arclist(i)=1;
% checklist(i)=0; %set checklist to 0 for dominated point in archive
% end
% end
% end 
% end
% arclist=[0,1,0];
% if (kn>0)
%     for i=1:length(checklist)
%         if checklist(i)~=1
%             Arc(i)=[];
%         end
%     end
% %     for i=1:length(arclist)
% %        if arclist(i)~=1
% %         urrArc(i)=urArc(i);
% %        end
% %     end
%     Arc=[Arc
%             urArc(find(arclist==1))];
%    [CL,~]=size(Arc);
% end
    
    



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% thin the archive if its length exceeds the soft limit
%
if(CL>SL)
    disp('Strip');
% indices=sort(datasample(1:CL,HL-2,'Replace',false));
indices=0;
indices=strip([Arc.cost]',HL-2);
%add best solution with respect to each OF to the set to retain
Ctem=zeros(length(Arc),2);
for i=1:length(Arc)
Ctem(i,:)=Arc(i).cost';
end
% [imax,~]=find(Ctem==max(Ctem(:,1)));
% [jmax,~]=find(Ctem==max(Ctem(:,2)));
[imin,~]=find(Ctem==min(Ctem(:,1)));
[jmin,~]=find(Ctem==min(Ctem(:,2)));
%indices=[indices,imax(1),jmax(1)];
indices=[indices,imin(1),jmin(1)];
indices=unique(indices);
CL=HL;
Arc=Arc(indices);
[CL,~]=size(AArc);
end


end




if(iter==1) 
    mtem=nacc/nrandom;
display(['' num2str(mtem)]);
end

propacc(iter)=nacc/nrandom;
% Ctem=zeros(length(Arc),2);
% for i=1:length(Arc)
% Ctem(i,:)=Arc(i).cost';
% end
% figure(1);
% plot(Ctem(:,1),Ctem(:,2),'o')
%  xlim([0 350])
%  ylim([0 500])

%plot(-Arc,xlim=c(0,350),ylim=c(0,500),pch=16,col="red",
%xlab="Variance of mean", ylab="Total distance /m")
%points(-urArc,pch="+",col="blue")
%text(340,390,round(propacc[iter],3))
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
indices=0;
indices=strip([Arc.cost]',HL-2);
%add best solution with respect to each OF to the set to retain
Ctem=zeros(length(Arc),2);
for i=1:length(Arc)
Ctem(i,:)=Arc(i).cost';
end
% [imax,~]=find(Ctem==max(Ctem(:,1)));
% [jmax,~]=find(Ctem==max(Ctem(:,2)));
[imin,~]=find(Ctem==min(Ctem(:,1)));
[jmin,~]=find(Ctem==min(Ctem(:,2)));
%indices=[indices,imax(1),jmax(1)];
indices=[indices,imin(1),jmin(1)];
indices=unique(indices);
CL=HL;
Arc=Arc(indices);
end
Ctem=zeros(length(Arc),2);
for i=1:length(Arc)
Ctem(i,:)=Arc(i).cost';
end

%  figure(1);
%  plot(Ctem(:,1),Ctem(:,2),'ko','markerfacecolor','k','MarkerSize',6)
%  xlim([0 350])
%  ylim([150 500])
%  xlabel({'Variance of mean'})
%  ylabel({'Total distance /m'})
% set(gca, 'XTick', [0 :50 :300],'YTick', [150 :50 :500])

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

Mid
SM
DM
SNS

figure(1)
plot(Ctem(:,1),Ctem(:,2),'ko','markerfacecolor','k','MarkerSize',6)
%PlotCosts(pop,Ctem')
%xlim([0 350])
%ylim([150 500])
%grid off
xlabel({'Mean total error'})
ylabel({'Total distance /m'})
% set(gca, 'XTick', [0 50 100 150 200 250 300],'YTick', [150 200 250 300 350 400 450 500])


for i=1:length(Arc)
    AA(i)=Arc(i).cost(1);
end
    
find(AA==min(AA))
find(AA==max(AA))


for i=1:length(Arc)
    figure(i)
    plot(Arc(i).pos(:,1), Arc(i).pos(:,2),'.k')
    xlim([0 100])
    ylim([0 100])
end

% 
%  figure(2)
%  plot([1:iter],propacc)
% indices=sort(datasample(1:CL,HL-2,'Replace',false));
% add best solution with respect to each OF to the set to retain
% Ctem=zeros(length(Arc),2);
% for i=1:length(Arc)
% Ctem(i,:)=Arc(i).cost';
% end
[imax,~]=find(Ctem==max(Ctem(:,1)));
[jmax,~]=find(Ctem==max(Ctem(:,2)));
% indices=[indices,imax,jmax];
% CL=HL;
% Arc=Arc(indices);


% figure(2)
% subplot(2,2,1)
% plot(Arc(1).pos(:,1), Arc(1).pos(:,2),'ko','markerfacecolor','k','MarkerSize',4)
% xlim([0 100])
% ylim([0 100])
% %set(gca, 'XTick', [0 20 40 60 80 100],'YTick', [0 20 40 60 80 100])
% set(gcf,'Position',[10 10 900 900])
% subplot(2,2,2)
% plot(Arc(3).pos(:,1), Arc(3).pos(:,2),'ko','markerfacecolor','k','MarkerSize',4)
% xlim([0 100])
% ylim([0 100])
% %set(gca, 'XTick', [0 20 40 60 80 100],'YTick', [0 20 40 60 80 100])
% set(gcf,'Position',[10 10 900 900])
% subplot(2,2,3)
% plot(Arc(11).pos(:,1), Arc(11).pos(:,2),'ko','markerfacecolor','k','MarkerSize',4)
% xlim([0 100])
% ylim([0 100])
% %set(gca, 'XTick', [0 20 40 60 80 100],'YTick', [0 20 40 60 80 100])
% set(gcf,'Position',[10 10 900 900])
% subplot(2,2,4)
% plot(Arc(40).pos(:,1), Arc(40).pos(:,2),'ko','markerfacecolor','k','MarkerSize',4)
% xlim([0 100])
% ylim([0 100])
% %set(gca, 'XTick', [0 20 40 60 80 100],'YTick', [0 20 40 60 80 100])
% set(gcf,'Position',[10 10 900 900])





% replot the Pareto front

%plot(-Arc,xlim=c(0,350),ylim=c(0,450),pch=16,col="red",
%xlab="Variance of mean", ylab="Total distance /m")
%points(-urArc,pch="+",col="blue")


%write.table(Arc,"ParetoSet.dat",quote=F,col.names=c("Variance","Distance"),
%row.names=F)
%write.table(ArcP,"ParetoSetDesigns.dat",quote=F)
%write.table(propacc,"Cooling.dat",quote=F)
%write.table(urArc,"InitialSet.dat",quote=F)
% for i=1:length(Arc)
%     figure(i)
%     plot(Arc(i).pos(:,1), Arc(i).pos(:,2),'o')
% end





