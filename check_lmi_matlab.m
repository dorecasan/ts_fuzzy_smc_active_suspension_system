clear all;
%% Parameters
run('globalParams.m');
Ksf = vss.Ksf; Ktf = vss.Ktf; Csf = vss.Csf; Ksr = vss.Ksr;g = vss.g;
Ktr = vss.Ktr; Csr = vss.Csr; Io = vss.Io; l1 = vss.l1; l2 = vss.l2;
zfMax = vss.zfMax; zrMax = vss.zrMax;

%%  TS configuration
mss = [vss.msMin;vss.msMax];
mufs = [vss.mufMin;vss.mufMax];
murs = [vss.murMin;vss.murMax];

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

            Fr = (ms*g*l1 + mur*g*(l1+l2))/(l1+l2);
            Ff = (ms+muf+mur)*g - Fr;
            Ai = [a00 a12;
                  a21 a22];
            Bi = [0 0 0 0 a1 a2 -1/muf 0;
                  0 0 0 0 a2 a3 0 -1/mur]';
            B1i = [0 0 -1 0 0 0 0 0;
                   0 0 0 -1 0 0 0 0]';
            C1i = [-Ksf/ms -Ksr/ms 0 0 -Csf/ms -Csr/ms Csf/ms Csr/ms;
                   l1*Ksf/Io -l2*Ksr/Io 0 0 l1*Csf/Io -l2*Csr/Io -l1*Csf/Io l2*Csr/Io];
            C2i = [1/zfMax 0 0 0 0 0 0 0;
                   0 1/zrMax 0 0 0 0 0 0;
                   0 0 Ksf/Ff 0 0 0 0 0;
                   0 0 0 Ksr/Fr 0 0 0 0];
            D1i = [1/ms 1/ms;
                   -l1/Io l2/Io];

            A{k+2*(j-1)+2*2*(i-1)} = Ai;
            B{k+2*(j-1)+2*2*(i-1)} = Bi;
            B1{k+2*(j-1)+2*2*(i-1)} = B1i;
            C1{k+2*(j-1)+2*2*(i-1)} = C1i;
            C2{k+2*(j-1)+2*2*(i-1)} = C2i;
            D1{k+2*(j-1)+2*2*(i-1)} = D1i;
        end   
    end
end

nin = 8;
nout = 2;

setlmis([]);
global gamma;
gamma = 5.8483;
rho = 1;
X = lmivar(1,[nin,1]); 
for i =1:1:8
    M{i} = lmivar(2,[nout,nin]);
end
numIeq = 1;
lmiterm([-numIeq 1 1 X],1,1);

for i = 1:1:8
    numIeq = numIeq + 1;
    lmiterm([numIeq,1,1,X],A{i},1,'s');
    lmiterm([numIeq,1,1,M{i}],B{i},1,'s');
    lmiterm([numIeq,1,2,0],B1{i});
    lmiterm([numIeq,1,3,X],1,C1{i}');
    lmiterm([numIeq,1,3,-M{i}],1,D1{i}')
    lmiterm([numIeq,2,2,0],-gamma^2)
    lmiterm([numIeq,3,3,0],-1);

    for j=1:1:4
        numIeq = numIeq + 1;
        lmiterm([numIeq,1,1,X],-1,1);
        lmiterm([numIeq,1,2,X],sqrt(rho),C2{i}(j,1:end)');
        lmiterm([numIeq,2,2,0],-1);
    end
    
end

for i=1:1:7
    for j=i+1:1:8
        numIeq = numIeq + 1;
        lmiterm([numIeq,1,1,X],A{i},1,'s');
        lmiterm([numIeq,1,1,M{j}],B{i},1,'s');
        lmiterm([numIeq,1,1,X],A{j},1,'s');
        lmiterm([numIeq,1,1,M{i}],B{j},1,'s');

        lmiterm([numIeq,1,2,0],B1{i});
        lmiterm([numIeq,1,2,0],B1{j});

        lmiterm([numIeq,1,3,X],1,C1{i}');
        lmiterm([numIeq,1,3,-M{j}],1,D1{i}')
        lmiterm([numIeq,1,3,X],1,C1{j}');
        lmiterm([numIeq,1,3,-M{i}],1,D1{j}')

        lmiterm([numIeq,2,2,0],-gamma^2)
        lmiterm([numIeq,2,2,0],-gamma^2)

        lmiterm([numIeq,3,3,0],-1);
        lmiterm([numIeq,3,3,0],-1);
    end
end

lmisys = getlmis;
[tmin,xfeas] = feasp(lmisys);

X_sol = dec2mat(lmisys,xfeas,X);
P = inv(X_sol);
controlParams.K{1} = dec2mat(lmisys,xfeas,M{1})*inv(X_sol);
controlParams.K{2} = dec2mat(lmisys,xfeas,M{2})*inv(X_sol); 
controlParams.K{3} = dec2mat(lmisys,xfeas,M{3})*inv(X_sol);
controlParams.K{4} = dec2mat(lmisys,xfeas,M{4})*inv(X_sol);
controlParams.K{5} = dec2mat(lmisys,xfeas,M{5})*inv(X_sol);
controlParams.K{6} = dec2mat(lmisys,xfeas,M{6})*inv(X_sol); 
controlParams.K{7} = dec2mat(lmisys,xfeas,M{7})*inv(X_sol);
controlParams.K{8} = dec2mat(lmisys,xfeas,M{8})*inv(X_sol);
