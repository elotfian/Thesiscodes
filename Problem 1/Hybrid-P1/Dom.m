function [op,klist]=Dom(xn,xc,Arc,tau,R,CL)

% Function to assess domination of xn,xc and Ar which is a set of objects 
% which are M-vectors like xn and xc.  Output is a vector.  
%
%The structure of the vector depends on the case (sensu Bandyopadhyay et al).


%case code 1, 2, or 3
%casecode=0;
% probdec  1 if a random decision on acceptance is to be made, with
% probacc the probability of acceptance.

%probdec=1;

% First question, dominance between xn and xc

dxnxc=ddom(xn,xc,R);
dcode=dxnxc(1);
domxnxc=dxnxc(2);

if dcode==-1

%Bandyo's case 1
%find mean dominance, Deldomav.  

casecode=1;
probdec=1;
klist=0;

md=domxnxc;
kd=0;
for i=1:CL
dv=ddom(xn,Arc(i).cost,R);
if(dv(1)==-1)
kd=kd+1;
md=md+dv(2);
end
end
Deldomav=md/(kd+1);
probacc=1/(1+exp(Deldomav/tau));

% set op for Bandyopadhyay et al. Case 1
%
op=[1,probacc];
elseif(dcode==0)

% Bandyo's case 2

% Compute mean dominance of xn by kd points in archive
% or identify kn points in archive dominated by xn (in klist) or
% note that xn and archive are non-dominating

%casecode=2;
klist=zeros(1,CL);
md=0;
kd=0;
kn=0;

for i=1:CL
dv=ddom(xn,Arc(i).cost,R);
if(dv(1)==-1)
kd=kd+1;
md=md+dv(2);
elseif (dv(1)==1)
kn=kn+1;
klist(kn)=i;
end 
end

if(kd>0)
probdec=1;
Deldomav=md/(kd);
probacc=1/(1+exp(Deldomav*tau));
elseif (kn>0)
probdec=0;
probacc=0;
list=find(klist~=0);
klist=klist(list);
else 
probdec=0;
probacc=0;
end

%set op for Bandyopadhyay et al. Case 2

op=[2,probdec,probacc,kn,klist];
else
  %Bandyo's case 3, xn dominates xc
% Compute dominance of xn by least dominating of kd points in archive
% or identify kn points in archive dominated by xn (in klist) or
% note that xn and archive are non-dominating.
casecode=3;
klist=zeros(1,CL);
md=0;
kd=0;
kn=0;
probdec=0;
probacc=0;
mindom=100000;
kmindom=0;
xcinarch=0; % says whether xc is in archive (case 3.2)
xcarch=0;	%index of xc in archive (case 3.2)

for i=1:CL
    tem=Arc(i).cost;
% if xc(1)==tem(1)
% if xc(2)==tem(2)
if abs(xc(1)-tem(1))<10e-4
if abs(xc(2)-tem(2))<10e-4
xcinarch=1;
xcarch=i;
end
end
dv=ddom(xn,Arc(i).cost,R);
if(dv(1)==-1)
kd=kd+1;
if(mindom>dv(2))
mindom=dv(2);
kmindom=i;
end
elseif (dv(1)==1) 
kn=kn+1;
klist(kn)=i;
end
end
Ban3case=2;

if(kd>0)
Ban3case=1;
probacc=1/(1+exp(-mindom));
elseif (kn>0)
Ban3case=3;   
list=find(klist~=0);
klist=klist(list);
end

%op=[3,Ban3case,kn,klist(1:kn),probacc,kmindom,xcinarch,xcarch];
op=[3,Ban3case,kn,klist,probacc,kmindom,xcinarch,xcarch];
end

end

