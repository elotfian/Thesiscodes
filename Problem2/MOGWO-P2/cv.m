function  Cv=cv(a,h,c0,nu)


if h==0 
    rho=1;
 else
rho=((1-c0)/((2^(nu-1))*gamma(nu)))*(((2*(nu^(1/2))*h)/a)^nu)*besselk(nu,(2*(nu^(1/2))*h)/a);
 end
Cv=rho;