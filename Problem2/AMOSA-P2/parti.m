function  [rhoCnu,rhoCa]=parti(a,h,nu)


%\partial C / \partal nu
if h==0
    rhoCnu=0;
else
    rhoCnu=(2^(1 - nu)*(-1/2)*(besselk(nu-1, (2*h*nu^(1/2))/a)-besselk(nu+1, (2*h*nu^(1/2))/a))*((2*h*nu^(1/2))/a)^nu)/gamma(nu)...
        + (2^(1 - nu)*besselk(nu, (2*h*nu^(1/2))/a)*(log((2*h*nu^(1/2))/a)*((2*h*nu^(1/2))/a)^nu...
        + (h*nu^(1/2)*((2*h*nu^(1/2))/a)^(nu - 1))/a))/gamma(nu)...
        - (2^(1 - nu)*psi(nu)*besselk(nu, (2*h*nu^(1/2))/a)*((2*h*nu^(1/2))/a)^nu)/gamma(nu)...
        - (2^(1 - nu)*log(2)*besselk(nu, (2*h*nu^(1/2))/a)*((2*h*nu^(1/2))/a)^nu)/gamma(nu);

end

%\partial C / \partal a
if h==0 
    rhoCa=0;
else
    rhoCa=- (2^(1 - nu)*((nu*besselk(nu, (2*h*nu^(1/2))/a))/a - ...
        (2*h*nu^(1/2)*besselk(nu + 1, (2*h*nu^(1/2))/a))/a^2)*((2*h*nu^(1/2))/a)^nu)/gamma(nu)...
        - (2*2^(1 - nu)*h*nu^(3/2)*besselk(nu, (2*h*nu^(1/2))/a)*((2*h*nu^(1/2))/a)^(nu - 1))/(a^2*gamma(nu));
end


