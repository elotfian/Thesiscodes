function  semiv=sv(a,h,c0,nu)



% if h>a
%     gh=1;
% else
% gh=(1.5*h/a)-(0.5*(h/a)^3);
% end
% 
% if h==0 
%    semiv=0;
% else
% semiv=c0+c1*gh;
% end
% end

 if h==0 
    rho=1;
 else
  rho=((1-c0)/((2^(nu-1))*gamma(nu)))*(((2*(nu^(1/2))*h)/a)^nu)*besselk(nu,(2*(nu^(1/2))*h)/a);
 end
semiv=1-rho;
 







