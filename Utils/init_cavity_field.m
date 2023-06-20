function psi0 = init_cavity_field(mu, Zeta, F0, startBranch, excite, d2, nf)

N_modes = length(mu);

verbose = false;

switch startBranch
    case 'top'
        psi0 = init_hom(Zeta,F0^2,2,verbose); % top branch
    case 'bottom'
        psi0 = init_hom(Zeta,F0^2,1,verbose); % bottom branch
    otherwise
        error('Unknown CW branch setting')
end

psi0 = psi0 * ones(N_modes, 1);

switch excite
    case 'soliton'
        T = linspace(-pi,pi,N_modes+1);
        T = T(1:end-1)';
        psi0 = psi0 + sqrt(2*Zeta)*...
            sech(sqrt(Zeta/d2)*T)*exp(0.5*1i);
    case 'platicon'
        % Normal dispersion case
        psi0 = 1 - double(abs(mu)<N_modes/8)';
        psi0 = init_hom(Zeta,F0^2, 2,verbose) .* psi0 + init_hom(Zeta,F0^2, 1,verbose);
    case 'rect'
        mask = abs(mu)<N_modes/8;
        psi0(mask) = F0*exp(1i*angle(psi0(mask)));
        [a,b]=butter(8,.1);
        psi0 = filtfilt(a,b,psi0);
    case 'bang'
        % big excitation, whatever you like
        psi0 = psi0 + sqrt(10*Zeta)*sech(sqrt(Zeta)/10*T)*exp(0.5*1i);
    case 'mod'
        psi0 = psi0.*(1-.5*cos(2*pi*(0:length(psi0)-1)'/length(psi0)));
    case 'gaussian'
        %psi0 = psi0 + sqrt(2*Zeta) * gausswin(N_modes, 2*pi*sqrt(Zeta/d2)) * exp(0.5*1i);
        psi0 = psi0 + init_hom(Zeta,F0^2, 2,verbose) * gausswin(N_modes, N_modes/20) * exp(0.5*1i);
    otherwise
        
end

if nf~=0
    psi0 = psi0 + nf * h *randn(1,N_modes) .*exp(1i*pi*rand(1,N_modes));
end

end

