function pos=Position(xmax,ymax,lbc2,nvar)


x=xmax*rand(1,1);
y=ymax*rand(1,1);
pos(1,:)=[x,y];
i=1;

 while i<nvar
 x=xmax*rand(1,1);
 y=ymax*rand(1,1);
 a = bsxfun(@minus,[x y],pos);
 out = cellfun(@norm,num2cell(a,2));
 if sum(out>lbc2)>=i
 i=i+1;
 pos(i,:)=[x,y];
 end
 end
 
 
