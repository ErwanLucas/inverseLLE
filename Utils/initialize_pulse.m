function psi0 = initialize_pulse(dispParam, i_p, F0, startBranch, excite)
% defines the pulse initialization in the cavity

type = dispParam.ParameterFitType;
mu = dispParam.mu;
zeta = dispParam.paramTable{i_p, 'detuning'};

Dint = dispParam.computeDispersion(i_p);

switch type
    case {'Dictionary', 'Dict_multiShifts'}
        d2 = abs(diff(Dint, 2));
        d2 = d2(mu == 0) / 2;
        %d2 = min(N_modes/50, abs(median(diff(Dint,2))));
        d2 = min(length(mu)/50, d2);
    case 'poly'
        pc = dispParam.getPolyCoeffs(i_p);
        d2 = pc(end-2);
end
psi0 = init_cavity_field(mu, zeta, F0, startBranch, excite, abs(d2), 0);

if isrow(psi0)
    psi0 = psi0.';
end
