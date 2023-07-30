function [sys, x0, str, ts] = roadProfile(t,x,u,flag)
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
    sizes.NumContStates = 0;
    sizes.NumDiscStates = 0;
    sizes.NumOutputs = 2;
    sizes.NumInputs = 0 ;
    sizes.DirFeedthrough = 1;
    sizes.NumSampleTimes = 0;
    sys = simsizes(sizes);
    x0 = [];
    str = [];
    ts = [];
end

function sys = mdlOutputs(t,x,u) 
    global vss;
    A = 0.1;L=2.5;V=20/3.6;

    if (t <= L/V)
        zrf = A/2*(1-cos(2*pi*V*t/L));
        zrr = A/2*(1-cos(2*pi*V*(t-(vss.l1+vss.l2)/V)/L));
    else
        zrf = 0;
        zrr = 0;
    end

    sys = [zrf;zrr];
end

