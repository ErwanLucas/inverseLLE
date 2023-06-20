function [score, components] = fitness_maxLinePower(Psi, dispParam,  target_powerdB, target_widthModes)

todB = @(X) 10*log10(X);

% in terms of power level
spectr = spectrumF(Psi);

roi = abs(dispParam.mu) <= target_widthModes & dispParam.mu~=0;
soft_metric = @(x) (abs(x)-x).^2; % this clamps to zero the positive overshoots so that only power below is penalized
components.powerMetric = sum( soft_metric( todB(spectr(roi,:)) - target_powerdB ), 1)';


% Fit the target as close as possible where meaningful
u = abs(Psi);
components.pulse = ( 5 * (median(u) + mad(u))./max(u) )';
components = struct2table(components);

score = components.powerMetric; %sum(table2array(components),2); % for now, disregard the pulse metric

end
