function deltaP_2_L_coeff = coeff_2_L(q2,r2)

    L = 185e-3;
%     r2 = 0.015;
    lambda_2_L = 1 + ((pi-q2)/L)*r2;
    
    A = 2.7e-3;    %#ok
    B = 4.7e-3;
    
    alpha = 67;
    n = (L/((2*pi*B)*tan(alpha*(pi/180))));
    l = L.*lambda_2_L;

    C1 = 243;
    C2 = 0;

    E = C1.*(lambda_2_L.^-1) + C2.*lambda_2_L;  %#ok
    b = (sqrt((((L^2).*(1-(lambda_2_L.^2))) + 4*(pi^2)*(n^2)*(B^2))/(4*(n^2)*(pi^2))));
    
    deltaP_2_L_coeff = pi*b.^2 - (l.^2*b)./(2*pi*b.*n.^2);

end