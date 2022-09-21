function emf = CalcePDTBSegMuAvgEMF(n_p, od_p, tk_p, coil, sig_p, mur_p, ...
    tx_cur, tx_dur, rx_delay, rx_time)
%CalcePDTBAvgEMF  Forward Modeling capsulation for ePDT-B.
%Simulate hybrid core with segment mu.
% Input:
%  n_p          number of pipes
%  od_p         outer diameter of each pipe, unit: mm
%  tk_p         thickness of each pipe, unit: mm
%  coil         excitation coil configuration: 'S', 'M' or 'L'
%  sig_p        conductivity of each pipe, unit: MS/m
%  mur_p        relative permeability of each pipe
%  tx_cur       tx pulse amplitude, unit: mA
%  tx_dur       tx pulse duration, unit: ms
%  rx_delay     acq delay time, unit: ms
%  rx_time      acq total length, unit: ms
%
% Output:
%  emf          electromotive force (voltage in millivolt), 1 x rx_time

%% units
m = 1.0;
mm = 1.0e-3;
inch = 25.4 * mm;
MS_per_m = 1.0e6;

%% build tool params, get correction factor accordingly
% ePDT-B geometry
tool_params = struct;
tool_params.hc = 28 * inch / 2;
tool_params.rc = 12 * mm;
tool_params.hp = 26 * inch / 2;
tool_params.rp = 14.125 * mm;
tool_params.Np = 1000;
tool_params.re = 16.25 * mm;
tool_params.rih = 19 * mm;
tool_params.roh = 25.4 * mm;

% coil specific
if coil == 'S'
    tool_params.he = 1.44 * inch / 2;
    tool_params.Ne = 200;
    if n_p > 2   % three to five pipes
        tool_params.h = 1.8 * m;
        tool_params.No = 20;
    elseif n_p > 1   % two pipes
        tool_params.h = 1.5 * m;
        tool_params.No = 30;
    else   % one pipes
        tool_params.h = 1.2 * m;
        tool_params.No = 50;
    end
    
    t_seg = {11:15, 16:25, 26:40, 41:60};
    mur_core = [180, 240, 300, 420];
    correction = [0.35, 0.005];
elseif coil == 'M'
    tool_params.he = 5.76 * inch / 2;
    tool_params.Ne = 800;
    if n_p > 2   % three to five pipes
        tool_params.h = 1.8 * m;
        tool_params.No = 20;
    else   % two pipes
        tool_params.h = 1.5 * m;
        tool_params.No = 30;
    end
    
    t_seg = {51:76, 77:155, 156:300};
    mur_core = [600, 1200, 1800];
    
    if od_p(end) < 16 * 25.4
        correction = [1.25, 0.003];
    else
        correction = [0.7, 0];
    end
elseif coil == 'L'
    tool_params.he = 23.04 * inch / 2;
    tool_params.Ne = 3200;
    tool_params.h = 1.8 * m;
    tool_params.No = 20;
    
    t_seg = {201:400, 401:1008};
    mur_core = [6e4, 6e5];
    correction = [1, 1.5e-4];
else
    emf = [];
    return;
end

%% build pipe params
pipe_params = struct;
for i = 1:n_p
    pipe_params(i).od = od_p(i) * mm;
    pipe_params(i).tk = tk_p(i) * mm;
    pipe_params(i).sig = sig_p(i) * MS_per_m;
    pipe_params(i).mur = mur_p(i);
end

%% build cfg params
cfg_params = struct;
cfg_params.tx_cur = tx_cur;
cfg_params.tx_dur = tx_dur;
cfg_params.rx_delay = rx_delay;

%% calculate emf and combine segment data
% calculate segment mu
n_seg = length(t_seg);
v_seg = cell(1, n_seg);
for i = 1:n_seg
    t = t_seg{i};
    tool_params.mur_cm = mur_core(i);
    cfg_params.rx_time = t(end);
    v_seg{i} = EvalAvgEMF(tool_params, pipe_params, cfg_params);
end

% combine segments
v_c = CombineSegments(t_seg, v_seg, rx_time);

%% apply correction
emf = 1e3 * polyval(correction, real(v_c));

end