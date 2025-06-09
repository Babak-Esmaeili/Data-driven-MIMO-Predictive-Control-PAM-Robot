function out = getParams()

    g = 9.81;

    L1 = 0.20;
    L2 = 0.20;

    Lc1 = L1/2;
    Lc2 = L2/2;

    m1 = 0.3;
    m2 = 0.3;

    r1 = 0.015;
    r2 = 0.015;

    n1 = 6;
    n2 = 3;

    J1 = m1*Lc1^2;
    J2 = m2*Lc2^2;
    
    out.g = g;
    out.L1 = L1;
    out.L2 = L2;
    out.Lc1 = Lc1;
    out.Lc2 = Lc2;
    out.m1 = m1;
    out.m2 = m2;
    out.r1 = r1;
    out.r2 = r2;
    out.n1 = n1;
    out.n2 = n2;
    out.J1 = J1;
    out.J2 = J2;

end