function [sys, x0, str, ts] = vssController(t,x,u,flag)
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
    sizes.NumContStates = 3;
    sizes.NumDiscStates = 0;
    sizes.NumOutputs = 3;
    sizes.NumInputs = 8 ;
    sizes.DirFeedthrough = 1;
    sizes.NumSampleTimes = 0;
    sys = simsizes(sizes);
    x0 = [zeros(3,1)];
    str = [];
    ts = [];
end
function sys = mdlDerivatives(t,x,u)
    global vss controlParams
    
    Ksf = vss.Ksf; Ktf = vss.Ktf; Csf = vss.Csf; Ksr = vss.Ksr;g = vss.g;
    Ktr = vss.Ktr; Csr = vss.Csr; Io = vss.Io; l1 = vss.l1; l2 = vss.l2;
    zfMax = vss.zfMax; zrMax = vss.zrMax;
    
    mss = [vss.msMin;vss.msMax];
    mufs = [vss.mufMin;vss.mufMax];
    murs = [vss.murMin;vss.murMax];
    
    states = reshape(u(1:8),[],1);
    
    for i = 1:1:2
        ms = mss(i);
        for j=1:1:2
            muf = mufs(j);
            for k=1:1:2
                mur = murs(k);
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
    
                Ai = [a00 a12;
                      a21 a22];
                Bi = [0 0 0 0 a1 a2 -1/muf 0;
                      0 0 0 0 a2 a3 0 -1/mur]';
             
                A{k+2*(j-1)+2*2*(i-1)} = Ai;
                B{k+2*(j-1)+2*2*(i-1)} = Bi;
            end   
        end
    end
    
    hs = membershipFunc();
    At = zeros(8,8);
    for i=1:1:8
        for j=1:1:8
            At = At + hs(i)*hs(j)*(A{i} + B{i}*controlParams.K{j});
        end
    end
    
    is = x(1:2); delta = x(3);
    s = controlParams.G*states - is;
    
    ds = controlParams.G*At*states;
    ddelta = controlParams.eta*norm(s,1)*norm(states,1);
    
    sys = [ds;ddelta];
end

function sys = mdlOutputs(t,x,u) 
    global vss controlParams;
    hs = membershipFunc();
    states = reshape(u(1:8),[],1);
    Kall = zeros(2,8);
    for i=1:1:8
        Kall = Kall + hs(i)*controlParams.K{i};
    end
    
    is = x(1:2); delta = x(3);
    s = controlParams.G*states - is;
    psi = controlParams.lambda + delta*norm(states,1); 
    control = Kall*states - psi*sign(s);
    sys = [control;delta];
end

function hs = membershipFunc()
    global vss controlParams;

    mss = [vss.msMin;vss.msMax];
    mufs = [vss.mufMin;vss.mufMax];
    murs = [vss.murMin;vss.murMax];
    hs = zeros(1,8);
    for i = 1:1:2
        mem1 = ((-1)^(i+1))*(1/mss(i)-1/controlParams.ms)/(1/mss(1)-1/mss(2));
        for j=1:1:2
            mem2 = ((-1)^(i+1))*(1/mufs(i)-1/controlParams.muf)/(1/mufs(1)-1/mufs(2));
            for k=1:1:2
                mem3 = ((-1)^(i+1))*(1/murs(i)-1/controlParams.mur)/(1/murs(1)-1/murs(2));
                hs(k+2*(j-1)+2*2*(i-1)) = mem1*mem2*mem3;
            end
        end
    end
end
