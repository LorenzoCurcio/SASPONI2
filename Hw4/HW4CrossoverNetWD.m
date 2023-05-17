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
C1= 2.8*10^(-6);
C2= 2.8*10^(-6);
C3= 28*10^(-6);
C4= 4.7*10^(-6);
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
Z38 = RspkHigh;
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
Z26 = 2*L2/(Ts);
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
a = zeros(Nsamp,42);
b = zeros(Nsamp,42);

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
    if ii==1
        a(ii,12)=0;
        a(ii,26)=0;
        a(ii,35)=0;
        a(ii,8)=0;
        a(ii,29)=0;
        a(ii,32)=0;
    else
        a(ii,12) = -b(ii-1,12);
        a(ii,26) = -b(ii-1,26);
        a(ii,35) = -b(ii-1,35);
        a(ii,8) = b(ii-1,8);
        a(ii,29) = b(ii-1,29);
        a(ii,32) = b(ii-1,32);

    end
    
    
    %% Forward Scan
    a(ii,1) = Vin(ii);

    %% Local Root Scattering


    %% Backward Scan


    %% Read Output
  
    
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

