function [error, Psi] = optimFunLLErun(param)
global dispParam pulseProfile_init fitness F0

dispParam.paramTable.polyCoefs = param(1:end-1).';
D = -1 - 1i * dispParam.computeDispersion(); % in normalized units
D = ifftshift(D);
% initialize the seed (in one direction)
zeta = param(end);
try
    psi0 = pulseProfile_init;
catch
    psi0 = initialize_pulse(dispParam, 1, F0, startBranch, excite);
end

% Define the LLE type and initialization
nf = 0; %NO NOISE
lle = @(psi, t_evol, h, Nstore) LLE_Propagate_para(psi, F0, zeta, D, nf, t_evol, h, Nstore);

% actually run
% first, propagate for a bit
t_evol = 40;
Nstore = 1e3;
h = 2^-9;
Psi_tmp = lle(psi0, t_evol, h, Nstore);
% refine the solution by decreasing step size
t_evol = 10;
Nstore = 1e3;
h_large = h;
h_fine = 2^-14;
h = @(t) (h_large-h_fine) * (1+tanh(15*(.4 - t/t_evol)))/2 + h_fine;
Psi = lle(Psi_tmp, t_evol, h, Nstore);

% estimate the error
error = fitness(Psi);

end