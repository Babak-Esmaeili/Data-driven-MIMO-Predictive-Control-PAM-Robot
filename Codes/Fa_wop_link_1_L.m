function Fa_wop_1_L = Fa_wop_link_1_L(q1,r1)

    L = 185e-3;
%     r1 = 0.015;
    lambda_1_L = 1 + ((pi-q1)/L)*r1;
    
    A = 2.7e-3;
    B = 4.7e-3;
    
    alpha = 67;
    n = (L/((2*pi*B)*tan(alpha*(pi/180))));
    l = L.*lambda_1_L;

    C1 = 243;
    C2 = 0;

    E = C1.*(lambda_1_L.^-1) + C2.*lambda_1_L;
    b = (sqrt((((L^2).*(1-(lambda_1_L.^2))) + 4*(pi^2)*(n^2)*(B^2))/(4*(n^2)*(pi^2))));
    
    D = ((((b.^2)-(B^2).*(lambda_1_L.^-1)))+eps);
    a = sqrt(((b.^2)-(lambda_1_L.^-1).*(B^2-A^2)));

    X2 = E.*b.*log((D-(b.^2))./(D-(a.^2)));
    X3 = -2.*E.*b.*log(b./a);
    X4 = E.*D.*((-1./b)+(b./(a.^2)));
    
    Fz = ((log(b./a).*(4*pi.*D.*C2 - 2.*pi.*D.*E.*lambda_1_L))./lambda_1_L - log((- b.^2 + D)./(- a.^2 + D)).*(pi.*D.*E - 2.*pi.*D.*C1.*lambda_1_L) - 2.*E.*b.^2.*pi.*log(b./a) -(1/2).*(-a.^2+b.^2).*(-4.*pi.*a.^2.*C2.*lambda_1_L.^4-4.*pi.*a.^2.*C1.*lambda_1_L.^3-2.*D.*E.*pi.*lambda_1_L.^2+4.*pi.*a.^2.*C2.*lambda_1_L+4.*pi.*a.^2*C1)./(a.^2.*lambda_1_L.^2));
    
    Fteta = (l.*(X2 + X3 + X4));
    
    Fa_wop_1_L = Fz - ((-(Fteta.*(l)))./(2*pi.*b.*(n.^2)));
    
end