global vss controlParams;

vss.Ksf = 18000;vss.Ktf = 200000;vss.Csf = 1000;vss.Ksr = 22000;
vss.Ktr = 200000;vss.Csr = 1000;vss.Io = 1222;vss.l1 = 1.3;vss.l2 = 1.5;vss.g=9.81;
vss.ms = 700; vss.muf = 40; vss.mur = 45;
vss.msMin = 621; vss.msMax = 759; vss.mufMin = 39.6; vss.mufMax = 40.4;
vss.murMin = 44.55; vss.murMax = 45.45;
vss.zfMax = 0.1; vss.zrMax = 0.1;

controlParams.G = [0 0 0 0 0 1 1 1;
                   0 0 0 0 1 0 1 1];
controlParams.lambda = 0.5;
controlParams.eta = 0.1;
controlParams.ms = 650;
controlParams.muf = 40.3;
controlParams.mur = 44.8;