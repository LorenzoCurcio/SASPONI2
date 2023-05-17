clear all
close all
clc

%% Import Input Audio Signal
[Vin,~] = audioread('ExpSweep.wav');

%% LTSpice Files for Ground-Truth
[OutLowSpice,~]=audioread('outlowsweep.wav');
[OutMidSpice,~]=audioread('outmidsweep.wav');
[OutHighSpice,FsLTSpice]=audioread('outhighsweep.wav');
TsLTSpice=1/FsLTSpice;

%% Sampling frequency (to be varied: FsLTSpice/downSampFact, with downSampFact={4,3,2})
downSampFact=2;
Fs =FsLTSpice/downSampFact; 

%% Downsample Input Signal
Vin=Vin([1:downSampFact:end]);
%% Sampling Period
Ts=1/Fs;
%% Number of Samples
Nsamp=length(Vin);
%% Simulated time
tstop=Nsamp*Ts;
%% Parameters of Dynamic Element
L1=0.35*10^(-3);
L2=0.35*10^(-3);
L3=3.5*10^(-3);
L4=3.5*10^(-3);
C1=2.8*10^(-6);
C2=2.8*10^(-6);
C3=28*10^(-6);
C4=4.7*10^(-6);
C5=28*10^(-6);
C6=47*10^(-6);
%% Resistive Parameters
R1=10;
RspkLow=8;
R2=10;
RspkMid=8;
RspkHigh=8;
%% WDF setting of free parameters (adaptation conditions)
Z12 = 2*L1/Ts;
Z11 = RspkHigh;
Z9 = Z11*Z12/(Z11+Z12);
Z10 = Z9;
Z8 = Ts/(2*C1);
Z5 = Z9+Z8;
Z7 = Z5;
Z24 = R2;
Z23 = Ts/(2*C6);
Z21 = Z24+Z23;
Z22 = Z21;
Z20 = RspkLow;
Z18 = Z20*Z21/(Z20+Z21);
Z19 = Z18;
Z17 = Ts/(2*C5);
Z15 = Z17*Z18/(Z17+Z18);
Z16 = Z15;
Z14 = 2*L4/Ts;
Z6 = Z14+15;
Z13 = Z6;
Z2 = Z6*Z5/(Z6+Z5);
Z4 = Z2;
Z42 = R1;
Z41 = Ts/(2*C4);
Z39 = Z42+Z41;
Z40 = Z39;
Z38 = RspkMid;
Z36 = Z39*Z38/(Z39+Z38);
Z37 = Z36;
Z35 = 2*L3/Ts;
Z33 = Z36*Z35/(Z36+Z35);
Z34 = Z33;
Z32 = Ts/(2*C3);
Z30 = Z32+Z33;
Z31 = Z30;
Z29 = Ts/(2*C2);
Z27 = Z29*Z30/(Z29+Z30);
Z28 = Z27;
Z26 = 2*L2/Ts;
Z3 = Z26+Z27;
Z25 = Z3;
Z1 = Z2*Z3/(Z2+Z3);

%% Computation of Scattering Matrices 
% parallels
S1 = 2/(1/Z1+1/Z2+1/Z3)*[1/Z1; 1/Z2; 1/Z3]*ones(1,3)-eye(3);
S4 = 2/(1/Z4+1/Z5+1/Z6)*[1/Z4; 1/Z5; 1/Z6]*ones(1,3)-eye(3);
S10 = 2/(1/Z10+1/Z11+1/Z12)*[1/Z10; 1/Z11; 1/Z12]*ones(1,3)-eye(3);
S16 = 2/(1/Z16+1/Z17+1/Z18)*[1/Z16; 1/Z17; 1/Z18]*ones(1,3)-eye(3);
S19 = 2/(1/Z19+1/Z20+1/Z21)*[1/Z19; 1/Z20; 1/Z21]*ones(1,3)-eye(3);
S28 = 2/(1/Z28+1/Z29+1/Z30)*[1/Z28; 1/Z29; 1/Z30]*ones(1,3)-eye(3);
S34 = 2/(1/Z34+1/Z35+1/Z36)*[1/Z34; 1/Z35; 1/Z36]*ones(1,3)-eye(3);
S37 = 2/(1/Z37+1/Z38+1/Z39)*[1/Z37; 1/Z38; 1/Z39]*ones(1,3)-eye(3);

%series
S7 = eye(3)-2/(Z7+Z8+Z9)*[Z7;Z8;Z9]*ones(1,3);
S13 = eye(3)-2/(Z13+Z14+Z15)*[Z13;Z14;Z15]*ones(1,3);
S22 = eye(3)-2/(Z22+Z23+Z24)*[Z22;Z23;Z24]*ones(1,3);
S25 = eye(3)-2/(Z25+Z26+Z27)*[Z25;Z26;Z27]*ones(1,3);
S31 = eye(3)-2/(Z31+Z32+Z33)*[Z31;Z32;Z33]*ones(1,3);
S40 = eye(3)-2/(Z40+Z41+Z42)*[Z40;Z41;Z42]*ones(1,3);

%% Initialization of Waves
a = zeros(Nsamp+1,42);
b = zeros(Nsamp+1,42);

%% Initialize Output Signals
% Low
VoutLow=zeros(size(Vin));
% Mid
VoutMid=zeros(size(Vin));
% High
VoutHigh=zeros(size(Vin));

ii=0;
while (ii<Nsamp)
    ii=ii+1;

    %% Manage Dynamic Elements

    a(ii+1,12) = -b(ii,12);
    a(ii+1,26) = -b(ii,26);
    a(ii+1,35) = -b(ii,35);
    a(ii+1,14) = -b(ii,14);

    a(ii+1,8) = b(ii,8);
    a(ii+1,29) = b(ii,29);
    a(ii+1,32) = b(ii,32);
    a(ii+1,41) = b(ii,41);
    a(ii+1,17) = b(ii,17);
    a(ii+1,22) = b(ii,22);

    
    
    %% Forward Scan
    %leaves
   
    a(ii+1,11) = 0;
    a(ii+1,20) = 0;
    a(ii+1,24) = 0;
    a(ii+1,42) = 0;
    a(ii+1,38) = 0;

    %first layer of adaptors
    b(ii+1,22) = S22(1,1)*a(ii+1,22)+S22(1,2)*a(ii+1,23)+S22(1,3)*a(ii+1,24);
    b(ii+1,10) = S10(1,1)*a(ii+1,10)+S10(1,2)*a(ii+1,11)+S22(1,3)*a(ii+1,12);
    b(ii+1,40) = S40(1,1)*a(ii+1,40)+S40(1,2)*a(ii+1,41)+S40(1,3)*a(ii+1,42);

    %second layer of adaptors
    a(ii+1,9) = b(ii+1,10);
    a(ii+1,21) = b(ii+1,22);
    a(ii+1,39) = b(ii+1,40);
    b(ii+1,7) = S7(1,1)*a(ii+1,7)+S7(1,2)*a(ii+1,8)+S7(1,3)*a(ii+1,9);
    b(ii+1,19) = S19(1,1)*a(ii+1,19)+S19(1,2)*a(ii+1,20)+S19(1,3)*a(ii+1,21);
    b(ii+1,37) = S37(1,1)*a(ii+1,37)+S37(1,2)*a(ii+1,38)+S37(1,3)*a(ii+1,39);

    %third layer
    a(ii+1,5) = b(ii+1,7);
    a(ii+1,18) = b(ii+1,19);
    a(ii+1,36) = b(ii+1,37);
    b(ii+1,16) = S16(1,1)*a(ii+1,16)+S16(1,2)*a(ii+1,17)+S16(1,3)*a(ii+1,18);
    b(ii+1,34) = S34(1,1)*a(ii+1,34)+S34(1,2)*a(ii+1,35)+S34(1,3)*a(ii+1,36);

    %4th layer
    a(ii+1,15)=b(ii+1,16);
    a(ii+1,33)=b(ii+1,34);
    b(ii+1,13) = S13(1,1)*a(ii+1,13)+S13(1,2)*a(ii+1,14)+S13(1,3)*a(ii+1,15);
    b(ii+1,31) = S31(1,1)*a(ii+1,31)+S31(1,2)*a(ii+1,32)+S31(1,3)*a(ii+1,33);

    %5th layer
    a(ii+1,6)=b(ii+1,13);
    a(ii+1,30)=b(ii+1,31);
    b(ii+1,4) = S4(1,1)*a(ii+1,4)+S4(1,2)*a(ii+1,5)+S4(1,3)*a(ii+1,6);
    b(ii+1,28) = S28(1,1)*a(ii+1,28)+S28(1,2)*a(ii+1,29)+S28(1,3)*a(ii+1,30);

    %6th layer
    a(ii+1,2) = b(ii+1,4);
    a(ii+1,27) = b(ii+1,28);
    b(ii+1,4) = S25(1,1)*a(ii+1,25)+S25(1,2)*a(ii+1,26)+S25(1,3)*a(ii+1,27);

    %7th layer
    a(ii+1,3) = b(ii+1,25);
    b(ii+1,1) = S1(1,1)*a(ii+1,1)+S1(1,2)*a(ii+1,2)+S1(1,3)*a(ii+1,3);


    %% Local Root Scattering
    a(ii+1,1) = 2*Vin(ii)-b(ii+1,1);

    %% Backward Scan

    %1st layer
    b(ii+1,2) = S1(2,1)*a(ii+1,1)+S1(2,2)*a(ii+1,2)+S1(2,3)*a(ii+1,3);
    b(ii+1,3) = S1(3,1)*a(ii+1,1)+S1(3,2)*a(ii+1,2)+S1(3,3)*a(ii+1,3);

    %2nd layer
    a(ii+1,4) = b(ii+1,2);
    a(ii+1,25) = b(ii+1,3);
    b(ii+1,5) = S4(2,1)*a(ii+1,4)+S4(2,2)*a(ii+1,5)+S4(2,3)*a(ii+1,6);
    b(ii+1,6) = S4(3,1)*a(ii+1,4)+S4(3,2)*a(ii+1,5)+S4(3,3)*a(ii+1,6);
    b(ii+1,26) = S25(2,1)*a(ii+1,25)+S25(2,2)*a(ii+1,26)+S25(2,3)*a(ii+1,27);
    b(ii+1,27) = S25(3,1)*a(ii+1,4)+S25(3,2)*a(ii+1,5)+S25(3,3)*a(ii+1,6);

    %3rd layer
    a(ii+1,13) = b(ii+1,6);
    a(ii+1,7) = b(ii+1,5);
    a(ii+1,28) = b(ii+1,27);
    b(ii+1,14) = S13(2,1)*a(ii+1,13)+S13(2,2)*a(ii+1,14)+S13(2,3)*a(ii+1,15);
    b(ii+1,15) = S13(3,1)*a(ii+1,13)+S13(3,2)*a(ii+1,14)+S13(3,3)*a(ii+1,15);
    b(ii+1,8) = S7(2,1)*a(ii+1,7)+S7(2,2)*a(ii+1,8)+S7(2,3)*a(ii+1,9);
    b(ii+1,9) = S7(3,1)*a(ii+1,7)+S7(3,2)*a(ii+1,8)+S7(3,3)*a(ii+1,9);
    b(ii+1,29) = S28(2,1)*a(ii+1,28)+S28(2,2)*a(ii+1,29)+S28(2,3)*a(ii+1,30);
    b(ii+1,30) = S28(3,1)*a(ii+1,28)+S28(3,2)*a(ii+1,29)+S28(3,3)*a(ii+1,30);
    
    %4th layer
    a(ii+1,16) = b(ii+1,15);
    a(ii+1,10) = b(ii+1,9);
    a(ii+1,31) = b(ii+1,30);
    b(ii+1,17) =  S16(2,1)*a(ii+1,16)+S16(2,2)*a(ii+1,17)+S16(2,3)*a(ii+1,18);
    b(ii+1,18) =  S16(3,1)*a(ii+1,16)+S16(3,2)*a(ii+1,17)+S16(3,3)*a(ii+1,18);
    b(ii+1,11) =  S10(2,1)*a(ii+1,10)+S10(2,2)*a(ii+1,11)+S10(2,3)*a(ii+1,12);
    b(ii+1,12) =  S10(3,1)*a(ii+1,10)+S10(3,2)*a(ii+1,11)+S10(3,3)*a(ii+1,12);
    b(ii+1,32) =  S31(2,1)*a(ii+1,31)+S31(2,2)*a(ii+1,32)+S31(2,3)*a(ii+1,33);
    b(ii+1,33) =  S31(3,1)*a(ii+1,31)+S31(3,2)*a(ii+1,32)+S31(3,3)*a(ii+1,33);

    %5th layer
    a(ii+1,19) = b(ii+1,18);
    a(ii+1,34) = b(ii+1,33);
    b(ii+1,20) = S19(2,1)*a(ii+1,19)+S19(2,2)*a(ii+1,20)+S19(2,3)*a(ii+1,21);
    b(ii+1,21) = S19(3,1)*a(ii+1,19)+S19(3,2)*a(ii+1,20)+S19(3,3)*a(ii+1,21);
    b(ii+1,35) = S34(2,1)*a(ii+1,34)+S34(2,2)*a(ii+1,35)+S34(2,3)*a(ii+1,36);
    b(ii+1,36) = S34(3,1)*a(ii+1,34)+S34(3,2)*a(ii+1,35)+S34(3,3)*a(ii+1,36);
    
    %6th layer
    a(ii+1,22) = b(ii+1,21);
    a(ii+1,37) = b(ii+1,36);
    b(ii+1,23) = S22(2,1)*a(ii+1,22)+S22(2,2)*a(ii+1,23)+S22(2,3)*a(ii+1,24);
    b(ii+1,24) = S22(3,1)*a(ii+1,22)+S22(3,2)*a(ii+1,23)+S22(3,3)*a(ii+1,24);
    b(ii+1,38) = S37(2,1)*a(ii+1,37)+S37(2,2)*a(ii+1,38)+S37(2,3)*a(ii+1,39);
    b(ii+1,39) = S37(3,1)*a(ii+1,37)+S37(3,2)*a(ii+1,38)+S37(3,3)*a(ii+1,39);

    %7th layer
    a(ii+1,40) = b(ii+1,39);
    b(ii+1,41) = S40(2,1)*a(ii+1,40)+S40(2,2)*a(ii+1,41)+S40(2,3)*a(ii+1,42);
    b(ii+1,42) = S40(3,1)*a(ii+1,40)+S40(3,2)*a(ii+1,41)+S40(3,3)*a(ii+1,42);


    %% Read Output
    VoutLow(ii) = (a(ii+1,20) + b(ii+1,20))/2;
    VoutMid(ii) = (a(ii+1,38) + b(ii+1,38))/2;
    VoutHigh(ii) = (a(ii+1,11) + b(ii+1,11))/2;
end

%% Output Plots
figure
set(gcf, 'Color', 'w');
subplot(311)
plot(TsLTSpice*[1:length(OutLowSpice)],OutLowSpice,'r','Linewidth',2); hold on;
plot(Ts*[1:Nsamp],VoutLow,'b--','Linewidth',1); grid on; xlim([0,tstop]); 
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$V_{\mathrm{outLow}}$ [V]','Fontsize',16,'interpreter','latex');
legend('LTspice','WDF','Fontsize',16,'interpreter','latex');
title('Output Signals','Fontsize',18,'interpreter','latex');
subplot(312)
plot(TsLTSpice*[1:length(OutMidSpice)],OutMidSpice,'r','Linewidth',2); hold on;
plot(Ts*[1:Nsamp],VoutMid,'b--','Linewidth',1); grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$V_{\mathrm{outMid}}$ [V]','Fontsize',16,'interpreter','latex');
subplot(313)
plot(TsLTSpice*[1:length(OutHighSpice)],OutHighSpice,'r','Linewidth',2); hold on;
plot(Ts*[1:Nsamp],VoutHigh,'b--','Linewidth',1); grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$V_{\mathrm{outHigh}}$ [V]','Fontsize',16,'interpreter','latex');

%% Error Plots
figure
set(gcf, 'Color', 'w');
subplot(311)
plot(Ts*[1:Nsamp],OutLowSpice([1:downSampFact:end])-VoutLow,'k','Linewidth',1);grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$E_{\mathrm{outLow}}$ [V]','Fontsize',16,'interpreter','latex');
title(['Error Signals. $F_{\mathrm{s}}=$ ',num2str(Fs),' Hz'],'Fontsize',18,'interpreter','latex');
subplot(312)
plot(Ts*[1:Nsamp],OutMidSpice([1:downSampFact:end])-VoutMid,'k','Linewidth',1);grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$E_{\mathrm{outMid}}$ [V]','Fontsize',16,'interpreter','latex');
subplot(313)
plot(Ts*[1:Nsamp],OutHighSpice([1:downSampFact:end])-VoutHigh,'k','Linewidth',1);grid on; xlim([0,tstop]);
xlabel('time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$E_{\mathrm{outHigh}}$ [V]','Fontsize',16,'interpreter','latex');

