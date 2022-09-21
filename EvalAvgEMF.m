function emf = EvalAvgEMF(tool_params, pipe_params, cfg_params)
%EvalAvgEMF  Calculate averaged EMF curve every ms.
% This function is based on the MTDe+ Job Planner forward modeling, 
%   fourth internal release (R4) on April 17, 2015, by Jin YingXin.
% Kc calculated from transendental equation by optimization function.
% ePDT V1 finished on May 11th 2017, by Yanxiang Yu.
%   Improved speed on June 30th 2018, by Kuang Qin.
% ePDT-B finished on May 18th 2020, by Kuang Qin.
% Generalized with input parameters on March 3rd 2021, by Kuang Qin.
%   Except those in cfg_params, all units should be SI units.
%   The unit of EMF is in volt.
%
% Input:
%  tool_params  tool parameters: geometry, core em properties, coil turns
%    h            half height of the simulation region, unit: m
%    No           number of points in simulation region
%    hc           half height of the magnetic core, unit: m
%    rc           radius of the magnetic core, unit: m
%    mur_cm       relative permeability of the core
%    hpx          half height of the pickup coil, unit: m
%    rp           outer radius of pickup coil, unit: m
%    Np           turns of pickup coil
%    he           half height of the excitation coil, unit: m
%    re           outer radius of excitation coil, unit: m
%    Ne           turns of excitation coil
%    rih          inner radius of titanium alloy housing, unit: m
%    roh          outer radius of titanium alloy housing, unit: m
%
%  pipe_params  pipe parameters: struct array of od, thickness, em properties
%    od           outer diameter of pipe, unit: m
%    tk           thickness of pipe, unit: m
%    sig          conductivity of pipe, unit: S/m
%    mur          relative permeability of pipe
%
%  cfg_paras    config parameters: input from config file
%    tx_cur       tx pulse amplitude, unit: mA
%    tx_dur       tx pulse duration, unit: ms
%    rx_delay     acq delay time, unit: ms
%    rx_time      acq total length, unit: ms
%
% Output:
%  emf          electromotive force (voltage in volt), 1 x rx_time


%% multi layer structure
%   |
% --|--|-|--|-|--|----|---|------|---|----|---|------|
%   |  | |  | |  |    | / |      | / |    | / |      |
%   |  | |  | |  |    | / |      | / |    | / |      |
% --|--|-|--|-|--|----|---|------|---|----|---|------|
%   |
%   0  1 2  3 4  5    6   7      8   9    10  11

% simulation range
h = tool_params.h;
No = tool_params.No;

% 0: center line
% 1: magnetic core
hc = tool_params.hc;           % half height of the magnetic core
rc = tool_params.rc;           % radius of the magnetic core
mur_cm = tool_params.mur_cm;   % permeability of the core

% 2: pickup coil
hp = tool_params.hp;           % half height of the pickup coil
rp = tool_params.rp;           % outer radius of pickup coil
Np = tool_params.Np;           % turns of pickup coil
rp1 = rc;
rp2 = rp;
zp1 = -hp;                     % lower boundary of pickup coil
zp2 = hp;                      % upper boundary of pickup coil

% 3: excitation coil
he = tool_params.he;           % half height of the pickup coil
re = tool_params.re;           % outer radius of pickup coil
Ne = tool_params.Ne;           % turns of pickup coil
re1 = rp;
re2 = re;
ze1 = -he;                     % lower boundary of excitation coil
ze2 = he;                      % upper boundary of excitation coil

% 4: air gap
% sig_a = 0;
% mur_a = 1;

% 5: titanium alloy housing
rih = tool_params.rih;
roh = tool_params.roh;
sig_h = 5.848e5; % S/m
mur_h = 1;

% 6: oil, gas or water
% the outer boundary of this layer is the inner boundary of the first pipe
sig_o = 0;
mur_o = 1;

% 7: tubing - 1st pipe
% 8: material between 1st and 2nd pipe
% 9: casing - 2nd pipe
% 10: material between 2nd and 3rd pipe
% 11: casing - 3rd pipe
% more pipes can be added from well schematic...
% total layer is twice (ID and OD) of each pipe plus housing
n_p = length(pipe_params);
n_layers = 2 * (n_p + 1);
r = zeros(1, n_layers);      % radius of each layer
r(1) = rih; r(2) = roh;
sig = zeros(1, n_layers);    % conductivity of each layer
sig(1) = sig_h; sig(2) = sig_o;
mur = zeros(1, n_layers);    % relative permeability of each layer
mur(1) = mur_h; mur(2) = mur_o;

% add more layers
for ip = 1:n_p
    r(2 * ip + 1) = pipe_params(ip).od / 2 - pipe_params(ip).tk;
    r(2 * ip + 2) = pipe_params(ip).od / 2;
    
    sig(2 * ip + 1) = pipe_params(ip).sig;
    sig(2 * ip + 2) = 0;
    
    mur(2 * ip + 1) = pipe_params(ip).mur;
    mur(2 * ip + 2) = 1;
end

%% Fourier transform
% units and constants
mA = 1.0e-3;
kHz = 1.0e3;
ms = 1.0e-3;
mu0 = 4 * pi * 1.0e-7;

% configurations
tx_cur = cfg_params.tx_cur * mA;
tx_dur = cfg_params.tx_dur;
rx_delay = cfg_params.rx_delay;
rx_time = cfg_params.rx_time;

% coefficient of the magnetic flux linkage
CT = Ne ./ ( (re - rp) * (ze2 - ze1) ) * Np / ( (rp - rc) * (zp2 - zp1) );

% calculate eigenvalues for the z-direction basis functions for the ferrite core region
% the general solution inside the magnetic core
Nc = No;
kc = Calculatekc(mur_cm, h, hc);
% load('kc_ePDT.mat'); % kc

% calculate the sampling points
k = (2 * ( 0:(No - 1) )' + 1) * pi / 2 / h;
fs = 22 * kHz;                         % Sampling frequency
Tm = 1536 * ms;                        % Total time

% Make sure Ns is an odd number,
% number of samples in time domain, and number of
% frequency components in frequency domain.
% One component for direct current,
% (Ns - 1) / 2 components for positive frequencies,
% (Ns - 1) / 2 components for negative frequencies.
Ns = int32(Tm * fs / 2 - 0.5) * 2 + 1;                                       
Ns = double(Ns);

% frequency vector for the fourier transform
Deltaf = 1 / Tm;    % Sampling frequency on the frequency domain
f = [0:(Ns - 1) / 2, -(Ns - 1) / 2:-1]' * Deltaf;
omega = 2.0 * pi * f;
Deltat = Tm / (Ns - 1);
t = (0:Ns - 1) * Deltat;

Te = (tx_dur - 1) * ms;               % input pulse duration
Nse = Te * fs;
I = [ones(1, Nse) * tx_cur, zeros(1, Ns - Nse) ];   % Input pulse current
FI = fft(I);

%% matrix calculation
% The following matrices UC, UD, ITZp, KTZp, ZKe, ZIe can be found in Dr.
% Jin's powerpoint slides on the foward modeling. Please check the original
% formula in the document
UC = zeros(No, Nc);
UD = zeros(No, Nc);

ITZp = zeros(1, No);
KTZp = zeros(1, No);

ZKe = zeros(No, 1);
ZIe = zeros(No, 1);

% This double loop is trying to match the expansion of z-direction basis
% function of pipe (and formation) layers and the ferrite core region.
% In the pipe and formation layer, the z-direction basis function is 
% exp(i*k1*z), exp(i*k2*z), exp(i*k3*z), ...
% But in the ferrite core region the basis function is a piecewise
% continous function that is described in Dr. Jin's power point slides.
% Please check Dr. Jin's power point slide for detailed description. 
% The theory of basis function matching is derived in the paper
% J. R. Bowler and T. P. Theodoulidis, "Eddy currents induced in a
% conducting rod of finite length by a coaxial encircling coil," J. Phys.
% D: Appl. Phys. 38(2005)2861-2868. 
% The pdf is on the share point. Please refer to this paper for
% mathematical details.
for i_ = 1:No
    ki = k(i_);
    
    kirc = ki * rc;
    
    for j_ = 1:Nc
        kcj = kc(j_);
         
        if k(i_) == kc(j_)
            kihc = ki * hc;
            kih = ki * h;
            sin2kihc = sin(2 * kihc);
            sinkih = sin(kih);
            sinkihc = sin(ki * hc);
                  
            Zcijm = (kihc + 0.5 * sin2kihc) * sin(kc(j_) * (h - hc) ) / kih;
            Zcija = (sinkih * (kih - kihc + 0.5 * (sin(2 * kih) - sin2kihc) ) - ...
                     cos(kih) * (sinkih * sinkih - sinkihc * sinkihc) ) * cos(kc(j_) * hc) / kih;
        else
            kimkcj = ki - kcj;
            kipkcj = ki + kcj;
            
            kimkcjhc = kimkcj * hc;
            kipkcjhc = kipkcj * hc;
            
            kimkcjh = kimkcj * h;
            kipkcjh = kipkcj * h;
            
            kcjh = kcj * h;
            kcjhc = kcj * hc;
            
            sinkimkcjhc = sin(kimkcjhc);
            sinkipkcjhc = sin(kipkcjhc);
            
            Zcijm = (sinkimkcjhc / kimkcj + sinkipkcjhc / kipkcj) * sin(kcjh - kcjhc) / h;
            Zcija = ( ( (sin(kimkcjh) - sinkimkcjhc) / kimkcj + ...
                        (sin(kipkcjh) - sinkipkcjhc) / kipkcj) * sin(kcjh) - ...
                      ( (cos(kimkcjh) - cos(kimkcjhc) ) / kimkcj + ...
                        (cos(kipkcjhc) - cos(kipkcjh) ) / kipkcj) * cos(kcjh) ) * cos(kcjhc) / h;
        end
        
        Zcij = Zcijm + Zcija;
        Zcmij = Zcijm / mur_cm + Zcija;
        
%        kirc = ki * rc;
        kcjrc = kcj * rc;
        
        kiZcijI1kcjrc = ki * Zcij * besseli(1, kcjrc);
        ZcmijkcjI0kcjrc = Zcmij * kcj * besseli(0, kcjrc);
        
        UC(i_, j_) = rc * (besselk(0, kirc) * kiZcijI1kcjrc + besselk(1, kirc) * ZcmijkcjI0kcjrc);  % Za(z) = 1/rc * I_n, I_n is an unit matrix
        UD(i_, j_) = rc * (besseli(0, kirc) * kiZcijI1kcjrc - besseli(1, kirc) * ZcmijkcjI0kcjrc);
        
    end

    invkiki = 1 / ki / ki;
    
    kirp2 = ki * rp2;
    kirp1 = ki * rp1;
    
    Zepi = 2.0 * sin(ki * hp) / ki; 
    
    ITZp(i_) = IntxI1(kirp2, kirp1) * invkiki * Zepi;
    KTZp(i_) = IntxK1(kirp2, kirp1) * invkiki * Zepi;
    
    kire2 = ki * re2;
    kire1 = ki * re1;
    
    ZKe(i_) = Zepi * IntxK1(kire2, kire1) * invkiki;
    ZIe(i_) = Zepi * IntxI1(kire2, kire1) * invkiki;
    
end

%% magnetic flux calculation
FPsi = zeros(1, Ns);

Nf = length(f);
gcp;
% Match boundary conditions of pipe and formation layers
% The extact organization of the code is still not sure
parfor i_f = 1:ceil(Nf/2)
    kappa = sqrt(repmat(k .* k, 1, n_layers) + ...
        repmat(1i * omega(i_f) * mur * mu0 .* sig, No, 1) );
    L = n_layers;
                    mur0 = 1;
                    
                    kappalm1rl = kappa(:, L - 1) * r(L);
                    kappalrl = kappa(:, L) * r(L);
                    invkappamur = mur(L) ./ kappa(:, L) .* kappa(:, L - 1) / mur(L - 1);
                    
                    
                    Rof = -invkappamur .* (besselk(1, kappalrl, 1) ./ besselk(0, kappalrl, 1));
                    Rofl12 = Rof .* besselk(0, kappalm1rl, 1) + besselk(1, kappalm1rl, 1);
                    Rofl22 = Rof .* besseli(0, kappalm1rl, 1) - besseli(1, kappalm1rl, 1);
                    Ror1 = exp(log(Rofl12 ./ Rofl22) - abs(real(kappalm1rl) ) - kappalm1rl);
                    
                    for el = (L - 1):-1:1
                        
                        if el > 1
                            kappalm1rl = kappa(:, el - 1) * r(el);
                            kappalrl = kappa(:, el) * r(el);
                            invkappamur = mur(el) ./ kappa(:, el) .* kappa(:, el - 1) / mur(el - 1);
                        else
                            kappalm1rl = k * r(el);
                            kappalrl = kappa(:, el) * r(el);
                            invkappamur = mur(el) ./ kappa(:, el) .* k / mur0;
                        end
                        
                        if kappalrl <= 700
                            expmabsrealkappalrl = exp(-abs(real(kappalrl) ) );
                            expkappalrl = exp(kappalrl);
                            
                            Rof = invkappamur .* (Ror1 .* besseli(1, kappalrl, 1) ./ expmabsrealkappalrl + ...
                                besselk(1, kappalrl, 1) ./ expkappalrl) ./ ...
                                (Ror1 .* besseli(0, kappalrl, 1) ./ expmabsrealkappalrl - ...
                                besselk(0, kappalrl, 1) ./ expkappalrl);
                        else
                            Rof = invkappamur .* besseli(1, kappalrl, 1) ./ besseli(0, kappalrl, 1);
                        end
                        
                        Ror1 = exp(log( (Rof .* besselk(0, kappalm1rl, 1) + besselk(1, kappalm1rl, 1) ) ./ ...
                            (Rof .* besseli(0, kappalm1rl, 1) - besseli(1, kappalm1rl, 1) ) ) - abs(real(kappalm1rl) ) - kappalm1rl);
                        
                    end        
                    
                    Ror = diag(Ror1);
    
    invW = UC - Ror * UD;
    E = ZKe + Ror * ZIe;
    
    FPsi(i_f) = -1j * omega(i_f) * 2 * pi * CT * mu0 / h * ...
        FI(i_f) * (ITZp * UC + KTZp * UD) / invW * E;
end

FPsi((ceil(Nf/2) + 1):end) = fliplr(conj(FPsi(2:ceil(Nf/2))));
% FPsi = 2 * pi * CT * mu0 / h * FPsi;

%% inverse Fourier transform and curve smoothing
Psiraw = ifft(FPsi);

% Psismooth = spline(t, real(Psiraw) );
% dPsismoothoverdt = fnder(Psismooth);
% EMF = -ppval(dPsismoothoverdt, t) * 1000;

% Post processing after inverse Fourier Transform
% Curve smoothing
% MTD V3 R4, April 17, 2015
Nch = tx_dur + rx_delay + rx_time;
tcs = (0:(Nch - 1)) * ms;
Deltatc = ones(1, Nch) * ms;
tce = tcs + Deltatc;

EMFa = zeros(1, Nch);

% ePDT V1 R1, May 11th, 2017
EMFsmooth = spline(t, Psiraw);

for ic = 1:Nch
    tc = linspace(tcs(ic), tce(ic), 1000);
    % MTD V3 R4, April 17, 2015
    % EMFsmooth = spline(t, EMF);
    EMFc = ppval(EMFsmooth, tc);

    EMFa(ic) = trapz(tc, EMFc) / Deltatc(ic);
end

% Only output decay curve after pulse
emf = EMFa((tx_dur + rx_delay + 1):end);

end


function sl0 = StruveL0(x)
%       ================================================
%       Purpose: Compute modified Struve function L0(x)
%       Input :  x   --- Argument of L0(x,x ?0)
%       Output:  SL0 --- L0(x)
%       ================================================

%pi=3.141592653589793;
s = 1.0;
r = 1.0;

if (x <= 20.0)
    a0 = 2.0 * x ./ pi;
    
    for k = 1:60
        r = r .* (x ./ (2.0 .* k + 1.0) ) .^ 2;
        s = s + r;

        if (abs(r ./ s) < 1.0e-12)
            break
        end
    end
    
    sl0 = a0 .* s;
else
    km = fix(0.5 * (x + 1.0) );

    if (x >= 50.0)
        km = 25;
    end
    
    for k = 1:km
        r = r .* ( (2.0 * k - 1.0) ./ x) .^ 2;
        s = s + r;

        if (abs(r ./ s) < 1.0e-12)
            break
        end
    end
    
    a1 = exp(x) ./ sqrt(2.0 * pi .* x);
    r = 1.0;
    bi0 = 1.0;
    
    for k = 1:16
        r = 0.125 * r .* (2.0 * k - 1.0) .^ 2 ./ (k .* x);
        bi0 = bi0 + r;

        if (abs(r ./ bi0) < 1.0e-12)
            break
        end
    end
    
    bi0 = a1 .* bi0;
    sl0 = -2.0 ./ (pi .* x) .* s + bi0;
end

% return;

end


function sl1 = StruveL1(x)
%       ================================================
%       Purpose: Compute modified Struve function L1(x)
%       Input :  x   --- Argument of L1(x,x ?0)
%       Output:  SL1 --- L1(x)
%       ================================================
% pi=3.141592653589793d0;

r = 1.0;
if (x <= 20.0)
    s = 0.0;

    for k = 1:60
        r = r .* x .* x ./ (4.0 * k .* k - 1.0);
        s = s + r;

        if (abs(r) < abs(s) * 1.0e-12)
            break
        end
    end
    
    sl1 = 2.0 / pi * s;
else
    s = 1.0;
    km = fix(0.5 * x);

    if (x > 50)
        km = 25;
    end
    
    for k = 1:km
        r = r .* (2.0 * k + 3.0) .* (2.0 * k + 1.0) ./ (x .* x);
        s = s + r;
        
        if (abs(r ./ s) < 1.0e-12)
            break
        end
    end
    
    sl1 = 2.0 / pi * (-1.0 + 1.0 ./ (x .* x) + 3.0 * s ./ x .^ 4);
    a1 = exp(x) ./ sqrt(2.0 * pi * x);
    r = 1.0;
    bi1 = 1.0;
    
    for k = 1:16
        r = -0.125 * r .* (4.0 - (2.0 * k - 1.0) .^ 2) ./ (k .* x);
        bi1 = bi1 + r;

        if (abs(r ./ bi1) < 1.0e-12)
            break
        end
    end
    
    sl1 = sl1 + a1 .* bi1;
end

% return

end


function result = IntxK1(x2, x1)

result = (pi / 2) * (x2 * (StruveL0(x2) * besselk(1, x2) + StruveL1(x2) * besselk(0, x2) ) - ...
                     x1 * (StruveL0(x1) * besselk(1, x1) + StruveL1(x1) * besselk(0, x1) ) );
                 
end


function result = IntxI1(x2, x1)

result = (pi / 2) * (x2 * (StruveL0(x2) * besseli(1, x2) - StruveL1(x2) * besseli(0, x2) ) - ...
                     x1 * (StruveL0(x1) * besseli(1, x1) - StruveL1(x1) * besseli(0, x1) ) );
                 
end


function Ror = CalculateRatio(L, kappa, k, r, mur)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

mur0 = 1;

kappalm1rl = kappa(:, L - 1) * r(L);
kappalrl = kappa(:, L) * r(L);
invkappamur = mur(L) ./ kappa(:, L) .* kappa(:, L - 1) / mur(L - 1);


Rof = -invkappamur .* (besselk(1, kappalrl, 1) ./ besselk(0, kappalrl, 1));
Rofl12 = Rof .* besselk(0, kappalm1rl, 1) + besselk(1, kappalm1rl, 1);
Rofl22 = Rof .* besseli(0, kappalm1rl, 1) - besseli(1, kappalm1rl, 1);
Ror = exp(log(Rofl12 ./ Rofl22) - abs(real(kappalm1rl) ) - kappalm1rl);

for el = (L - 1):-1:1
    
    if el > 1
        kappalm1rl = kappa(:, el - 1) * r(el);
        kappalrl = kappa(:, el) * r(el);
        invkappamur = mur(el) ./ kappa(:, el) .* kappa(:, el - 1) / mur(el - 1);
    else 
        kappalm1rl = k * r(el);
        kappalrl = kappa(:, el) * r(el);
        invkappamur = mur(el) ./ kappa(:, el) .* k / mur0;
    end
    
    if kappalrl <= 700
        expmabsrealkappalrl = exp(-abs(real(kappalrl) ) );
        expkappalrl = exp(kappalrl);
        
        Rof = invkappamur .* (Ror .* besseli(1, kappalrl, 1) ./ expmabsrealkappalrl + ...
                                     besselk(1, kappalrl, 1) ./ expkappalrl) ./ ...
                             (Ror .* besseli(0, kappalrl, 1) ./ expmabsrealkappalrl - ...
                                     besselk(0, kappalrl, 1) ./ expkappalrl);
    else
        Rof = invkappamur .* besseli(1, kappalrl, 1) ./ besseli(0, kappalrl, 1);
    end
    
   Ror = exp(log( (Rof .* besselk(0, kappalm1rl, 1) + besselk(1, kappalm1rl, 1) ) ./ ...
                  (Rof .* besseli(0, kappalm1rl, 1) - besseli(1, kappalm1rl, 1) ) ) - abs(real(kappalm1rl) ) - kappalm1rl);
              
end

end


function kc = Calculatekc(murcm, h, hc)
% Calculate root kc of transendental equation with optimization function
% All value calculated are real value 

func = @(kc) murcm * cot(kc * (h - hc) ) - tan(kc * hc);

n_h = 0:499;
kc_h_init = (2 * n_h + 1) * pi / (2 * (h - hc));
kc_h_out = SearchNearZero(func, kc_h_init);

% solve zero points of hc term solution
n_hc = ceil(length(kc_h_init) * hc / (h - hc));
kc_hc_init = (2 * (0:(n_hc - 1)) + 1) * pi / (2 * hc);
kc_hc_out = SearchDiffRange(func, kc_hc_init, kc_h_out);

kc_poss = sort([kc_h_out, kc_hc_out]);
kc = kc_poss(1:500);

end


function [ m ] = BinarySearchZero(func, l, r, epsilon)
%BinarySearchZero Binary search zero in a given range.
% Input:
%  func     monotonically decreasing function
%  l        left boundary
%  r        right boundary
%  epsilon  zero condition
%
% Output:
%  m        mid point that satisfy

% assert(l <= r);       % left not larger than right
% assert(func(l) > 0);  % opposite sign and monotonically decreasing
% assert(func(r) < 0);  % opposite sign and monotonically decreasing

m = (l + r) / 2;    % mid point
y_m = func(m);

% recursive: commented due to bad MATLAB suppport
% if abs(y_m) < epsilon
%     % do nothing
%     % return m;
% elseif y_m > 0
%     m = BinarySearchZero(func, m, r, epsilon);
% else  % y_m < 0
%     m = BinarySearchZero(func, l, m, epsilon);
% end

% iterative
while abs(y_m) > epsilon
    if (m == l) || (m == r)
        break;
    end
    
    if y_m > 0
        l = m;
        m = (l + r) / 2;
        y_m = func(m);
    else  % y_m < 0
        r = m;
        m = (l + r) / 2;
        y_m = func(m);
    end
end

end


function [ z ] = SearchNearZero(func, x)
%SearchNearZero Find the near zero in a monotonically decreasing function.
% Input:
%  func     monotonically decreasing function
%  x        random point
%
% Output:
%  z        zero point

epsilon = 1e-8;

z = zeros(size(x));
for i = 1:length(x)
    y = func(x(i));
    if abs(y) < epsilon    % satisfy zero condition
        z(i) = x(i);
        continue;
    end
    
    if y < 0
        inc = -0.1 * (x(2) - x(1)); % move to neg x to find pos func(x)
    else
        inc = 0.1 * (x(2) - x(1));  % move to pos x to find neg func(x)
    end
    
    x_inc = x(i) + inc;
    y_inc = func(x_inc);
    while (y * y_inc > 0)   % same sign
        x_inc = x_inc + inc;
        y_inc = func(x_inc);
    end
    
    if abs(y_inc) < epsilon    % satisfy zero condition
        z(i) = x_inc;
        continue;
    end
    
    % binary search
    if x_inc < x(i)
        l = x_inc;
        r = x(i);
    else
        l = x(i);
        r = x_inc;
    end
    z(i) = BinarySearchZero(func, l, r, epsilon);
end

end


function [ z ] = SearchDiffRange(func, x, grid)
%SearchDiffRange Search derivative range to find zero points.
% Input:
%  func     monotonically decreasing function
%  x        random point
%  grid     grid to locate x
%
% Output:
%  z        zero point

epsilon = 1e-8;

z = zeros(size(x));
for i = 1:length(x)
    % find the range between two plus sign zero points
    i_low = find(grid < x(i), 1, 'last');
    if i_low == length(grid)
        % break if find the end of grid
        n_z = i - 1;
        break;
    end
    
    interval = 1e-3;
    if isempty(i_low)
        % the zero point falls before the first value in grid
        kc_grid = 1e-3:1e-3:grid(1);
        f_grid = func(kc_grid);
        f_diff = diff(f_grid);
        i_r = find(f_diff > 1e3);
        i_l = i_r - 1;
    else
        % take derivative to find the two tangent sigularities
        i_limit = [];
        while (length(i_limit) < 2)
            kc_grid = grid(i_low):interval:grid(i_low + 1);
            f_grid = func(kc_grid);
            f_diff = diff(f_grid);
            i_limit = find(f_diff > 1e3);
            interval = interval * 1e-1;
        end

        i_l = i_limit(1) + 1;
        i_r = i_limit(2);
    end
    
    % set boundary
    l = kc_grid(i_l);
    r = kc_grid(i_r);
    
    while func(l) < 0  % maintain left side positive sign
        r = l;  % set previous left boudary to right
        l = kc_grid(i_l - 1);
        kc_grid = l:1e-3 * (r - l):r;
        f_grid = func(kc_grid);
        f_diff = diff(f_grid);
        i_l = find(f_diff > 1e3) + 1;
        l = kc_grid(i_l);
    end
    
    while func(r) > 0  % maintain right side negative sign
        l = r;
        r = kc_grid(i_r + 1);
        kc_grid = l:1e-3 * (r - l):r;
        f_grid = func(kc_grid);
        f_diff = diff(f_grid);
        i_r = find(f_diff > 1e3);
        r = kc_grid(i_r);
    end
    
    % binary search
    z(i) = BinarySearchZero(func, l, r, epsilon);
    
    % update valid solution counter
    n_z = i;
end

% output valid solution
z = z(1:n_z);

end

