function [sys, x0, str, ts] = vssModel(t,x,u,flag)
    switch flag
    case 0
        [sys, x0, str,ts] = mdlInitializeSizes;
    case 1
        sys = mdlDerivatives(t,x,u);
    case 3
        sys = mdlOutputs(t,x,u);
    case {2,4,9}
        sys = [];
    otherwise
        error(['Unhandled flag = ', num2str(flag)]);
    end
end
function [sys,x0,str,ts]=mdlInitializeSizes
    sizes = simsizes;
    sizes.NumContStates = 8;
    sizes.NumDiscStates = 0;
    sizes.NumOutputs = 14;
    sizes.NumInputs = 4 ;
    sizes.DirFeedthrough = 1;
    sizes.NumSampleTimes = 0;
    sys = simsizes(sizes);
    x0 = [zeros(8,1)];
    str = [];
    ts = [];
end
function sys = mdlDerivatives(t,x,u)
    global vss;
    Ksf = vss.Ksf; Ktf = vss.Ktf; Csf = vss.Csf; Ksr = vss.Ksr;g = vss.g;
    Ktr = vss.Ktr; Csr = vss.Csr; Io = vss.Io; l1 = vss.l1; l2 = vss.l2;
    ms = vss.ms; muf = vss.muf; mur = vss.mur;
    
    x1 = x(1);x2=x(2);gd = [0.5*x1;0.5*x2];
    uf = u(1); ur = u(2);zrf = u(3);zrr = u(4);
    U = [uf;ur];w = [zrf;zrr];
    
    a1 = 1/ms + l1^2/Io; a2=1/ms-l1*l2/Io; a3=1/ms+l2^2/Io;
    a00 = zeros(4);
    a12 = [1 0 -1 0;
           0 1  0 -1;
           0 0 1 0;
           0 0 0 1];
    a21 = [-a1*Ksf -a2*Ksr 0 0;
           -a2*Ksf -a3*Ksr 0 0;
           Ksf/muf 0 -Ktf/muf 0;
           0 Ksr/mur 0 -Ktr/mur];
    a22 = [-a1*Csf -a2*Csr a1*Csf a2*Csr;
           -a2*Csf -a3*Csr a2*Csf a3*Csr;
           Csf/muf 0 -Csf/muf 0;
           0 Csr/mur 0 -Csr/mur];
    
    A = [a00 a12;
          a21 a22];
    B = [0 0 0 0 a1 a2 -1/muf 0;
          0 0 0 0 a2 a3 0 -1/mur]';
    B1 = [0 0 -1 0 0 0 0 0;
           0 0 0 -1 0 0 0 0]';
    
    dx = A*x + B*(U+gd) + B1*w;
    sys= dx;
end

function sys = mdlOutputs(t,x,u) 
    global vss;
    Ksf = vss.Ksf; Ktf = vss.Ktf; Csf = vss.Csf; Ksr = vss.Ksr;g = vss.g;
    Ktr = vss.Ktr; Csr = vss.Csr; Io = vss.Io; l1 = vss.l1; l2 = vss.l2;
    ms = vss.ms; muf = vss.muf; mur = vss.mur;
    zfMax = vss.zfMax; zrMax = vss.zrMax;
    
    x1 = x(1);x2=x(2);x5=x(5);x6=x(6);x7=x(7);x8=x(8);
    gd = [0.5*x1;0.5*x2];
    uf = u(1); ur = u(2);zrf = u(3);zrr = u(4);
    U = [uf;ur];w = [zrf;zrr];
    
    Ff = (ms*g*l1 + mur*g*(l1+l2))/(l1+l2);
    Fr = (ms+muf+mur)*g - Ff;
    
    C1 = [-Ksf/ms -Ksr/ms 0 0 -Csf/ms -Csr/ms Csf/ms Csr/ms;
           l1*Ksf/Io -l2*Ksr/Io 0 0 l1*Csf/Io -l2*Csr/Io -l1*Csf/Io l2*Csr/Io];
    C2 = [1/zfMax 0 0 0 0 0 0 0;
           0 1/zrMax 0 0 0 0 0 0;
           0 0 Ksf/Ff 0 0 0 0 0;
           0 0 0 Ksr/Fr 0 0 0 0];
    D1 = [1/ms 1/ms;
           -l1/Io l2/Io];
    z1 = C1*x + D1*(U+gd);
    z2 = C2*x;

    sys = [x;z1;z2];

end

