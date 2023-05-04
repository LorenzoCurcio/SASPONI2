function [y,y_lpf,y_hpf,y_hp_sdf] = leslie_o(x, Fs, freq)
%Leslie Speaker Emulation
%
% J. Pekonen et al. Computationally Efficient Hammond Organ Synthesis
% in Proc. of the 14th International Conference on Digital Audio
% Effects(DAFx-11), Paris, France, Sept. 19-23, 2011

% length of the input signal
N = length(x);

% global modulator parameters
alpha=0.9;
% tremble spectral delay filter parameter 
Ms_t=0.2;
Mb_t=-0.75;
N_sdf_t=4;
% bass spectral delay filter parameter 
Ms_b=0.04;
Mb_b=-0.92;
N_sdf_b=3;

% cross-over network design
fc=800;                 % cutoff frequency

%TODO: compute the coefficients for the two 4th order butterworth filters
%with cutoff frequency fc
[b_lp, a_lp]=butter(4,fc/Fs,'low');... %LPF design
[b_hp, a_hp]=butter(4,fc/Fs,'high');...  %HPF design

% allocate input and output buffers for IIR filters
% hp filter buffers
hpf.state=zeros(N+4,1);
hpf.in=zeros(N+4,1);
% lp filter buffers
lpf.state=zeros(N+4,1);
lpf.in=zeros(N+4,1);
% treble sdf filter buffers
sdf_h.state=zeros(N+N_sdf_t,1);
sdf_h.in=zeros(N+N_sdf_t,1);
% bass sdf filter buffers
sdf_b.state=zeros(N+N_sdf_t,1);
sdf_b.in=zeros(N+N_sdf_t,1);

% modulators
sinusoid = sin(2*pi*freq*linspace(0,N/Fs,N));
m_b=Ms_b*sinusoid+Mb_b; % bass modulator
m_t=Ms_t*sinusoid+Mb_t; % tremble modulator

%sample processing
x = [0;0;0;0;x];
hpf.in = x;
lpf.in = x;
y = hpf.in;
for n=5:N+4

    % compute crossover network filters outputs
    y_lpf= (a_lp(5)/b_lp(1))*lpf.in(n) + (a_lp(4)/b_lp(1))*lpf.in(n-1) + (a_lp(3)/b_lp(1))*lpf.in(n-2) + (a_lp(2)/b_lp(1))*lpf.in(n-3) + (a_lp(1)/b_lp(1))*lpf.in(n-4) - (b_lp(2)/b_lp(1))*lpf.state(n-1) - (b_lp(3)/b_lp(1))*lpf.state(n-2) - (b_lp(4)/b_lp(1))*lpf.state(n-3) - (b_lp(5)/b_lp(1))*lpf.state(n-4);
    lpf.state(n) = y_lpf;
    %y_lpf = [0;0;0;0;y_lpf];
    y_hpf= (a_hp(5)/b_hp(1))*hpf.in(n) + (a_hp(4)/b_hp(1))*hpf.in(n-1) + (a_hp(3)/b_hp(1))*hpf.in(n-2) + (a_hp(2)/b_hp(1))*hpf.in(n-3) + (a_hp(1)/b_hp(1))*hpf.in(n-4) - (b_hp(2)/b_hp(1))*hpf.state(n-1) - (b_hp(3)/b_hp(1))*hpf.state(n-2) - (b_hp(4)/b_hp(1))*hpf.state(n-3) - (b_hp(5)/b_hp(1))*hpf.state(n-4);
    hpf.state(n) = y_hpf;
    %y_hpf = [0;0;0;0;y_hpf];

    % compute bass SDF output
    sum = 0;
    sdf_b.in(n) = lpf.state(n);
    for i=0:1:N_sdf_b
        sum = sum + nchoosek(N_sdf_b,i)*m_b(n-4)^i*(sdf_b.in(n-(N_sdf_b-i))-sdf_b.state(n-i));
    end
    y_lp_sdf= sum;
    sdf_b.state(n) = y_lp_sdf;

    % compute treble SDF output
    sum = 0;
    sdf_h.in(n) = hpf.state(n);
    for i=0:1:N_sdf_t
        sum = sum + nchoosek(N_sdf_t,i)*m_t(n-4)^i*(sdf_h.in(n-(N_sdf_t-i))-sdf_h.state(n-i));
    end
    y_hp_sdf= sum;
    sdf_b.state(n) = y_hp_sdf;

    % implement AM modulation blocks
    y_lp_am= (1 + alpha*m_b(n-4))*y_lp_sdf;
    y_hp_am= (1 + alpha*m_b(n-4))*y_hp_sdf;

    y(n-4)= y_lp_am +y_hp_am;

end
y = y(1:N);
end

