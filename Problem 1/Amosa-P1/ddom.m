function dxnxc=ddom(x,y,R)

% function to indicate domination of x and y (minimum problem)
% which are M-vectors.  Output is a vector.  First element is
% 0 (non-dominating pair), 1 (x dominates y) or -1 (y dominates x)
% Second element is total domination, given ranges in vector R


diff=y-x;

mind=min(diff);
maxd=max(diff);

if mind*maxd<0
do=0;
elseif mind>=0
do=1;
else
do=-1;
end

Del(1)=abs(diff(1))/R(1);
Del(2)=abs(diff(2))/R(2);
Del(find(Del==0))=1;
Tdom=prod(Del);

dxnxc=[do,Tdom];

