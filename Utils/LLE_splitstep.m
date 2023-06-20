function psi_out = LLE_splitstep(psi, D_Hat, C_hat, S, yL, h, withXPMBiDir)
% psi is input field (in the time domain)
% D_Hat is the linear operator
% S is the pump function
% yL is the nonlinear factor (usually normalized to 1)
% h is the step size

% doesn't have self-steepening / Raman

% N is the nonlinear operator
%N = 1i * yL * abs(psi).^2 + S./(psi-eps);% solution to first order diff eq for nonlinear step
psi = exp(h * 1i * yL * abs(psi).^2).* (psi + S*h); % approx.

% going to the frequency domain (careful with the array dimensions)
psi = fft(psi);

% coupling step (if defined)
if ~isempty(C_hat)
    error('Bidir. coupling not implemented')
end

% solution to first order diff eq for linear step
psi = D_Hat .* psi;

% going back to time domain
psi = ifft(psi);

% XPM between cw/ccw if coupling and enabled
if ~isempty(C_hat) && withXPMBiDir
    error('Bidir. coupling not implemented')
end

% output
psi_out = psi;

return