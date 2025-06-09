function qddot = PAM_twoLink_robot_dynamics(t,q,qdot,u,u_dist)
    
    x1 = q;
    x2 = qdot;

    params = getParams();

    g = params.g;

    L1 = params.L1;
    L2 = params.L2; %#ok

    Lc_1 = params.Lc1;
    Lc_2 = params.Lc2;

    m1 = params.m1;
    m2 = params.m2;

    r1 = params.r1;
    r2 = params.r2;

    n1 = params.n1; %#ok
    n2 = params.n2; %#ok

    J1 = params.J1;
    J2 = params.J2;
    
    a11 = m2*(2*L1*Lc_2*cos(x1(2)) + Lc_2^2 + L1^2) + m1*Lc_1^2 + J1 + J2;
    a12 = m2*(Lc_2^2 + L1*Lc_2*cos(x1(2))) + J2;
    a22 = m2*Lc_2^2 + J2;

    A = [a11 a12
         a12 a22];

    b11 = -m2*L1*Lc_2*x2(2)*sin(x1(2));
    b12 = -m2*L1*Lc_2*(x2(1) + x2(2))*sin(x1(2));
    b21 = m2*L1*Lc_2*x2(1)*sin(x1(2));
    b22 = 0;

    B = [b11 b12
         b21 b22];

    c11 = -m1*g*Lc_1*sin(x1(1)) - m2*g*(L1*sin(x1(1)) + Lc_2*sin(x1(1)+x1(2)));
    c21 = -m2*g*Lc_2*sin(x1(1)+x1(2));

    C = [c11
         c21];

    d11 = (-coeff_1_L(x1(1),r1) - coeff_1_R(x1(1),r1))*r1;
    d22 = (-coeff_2_L(x1(2),r2) - coeff_2_R(x1(2),r2))*r2;

    D = [d11   0
         0   d22];

    h11 = (Fa_wop_link_1_R(x1(1),r1) - Fa_wop_link_1_L(x1(1),r1))*r1;
    h22 = (Fa_wop_link_2_R(x1(2),r2) - Fa_wop_link_2_L(x1(2),r2))*r2;

    H = [h11
         h22];

    Bp = A\B;
    Cp = A\C;
    Hp = A\H;

    F = -Bp*[x2(1);x2(2)] - Cp + Hp;
    G = (10^5)*(A\D);
    G_dist = -(10^5)*A\eye(2);
        
%    u_dist = [2*sign(x2(1))
%              2*sign(x2(2))];
        
%    u_dist = [2*sin(t) + 0.5*sin(200*pi*t)
%              cos(2*t) + 0.5*sin(200*pi*t)];
        
%     u_dist = [5*cos(t)
%               5*cos(2*t)];

    qddot = F + G*u + 1*G_dist*u_dist;

end