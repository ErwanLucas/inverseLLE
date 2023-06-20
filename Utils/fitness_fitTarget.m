function [score, components] = fitness_fitTarget(Psi, dispParam, combTarget)

% in terms of power level
% Target_power_dB = -20;
% combTarget_width = 32;
% roi = abs(mu) <= combTarget_width & mu~=0;
% soft_metric = @(x) (abs(x)-x).^2; % this clamps to zero the positive overshoots so that only power below is penalized
% fitness = @(spectr) sum( soft_metric( todB(spectr(roi,:)) - Target_power_dB ), 1);

%fitness = @(spectr) sum((todB(spectr(roi,:)) - todB(combTarget(roi))).^2, 1);
%fitness = @(spectr) sum(( spectr(roi,:) - 1e-2 ).^2, 1); %
%fitness = @(spectr) sum(( todB(spectr(roi,:)) - Target_power_dB ).^2, 1); %

% mu_DW_target = 967;
% roi = abs(mu - mu_DW_target) < 5;
% fitness = @(spectr) sum((log10(spectr(roi,:)) + 3 ).^2, 1)/(9*sum(roi));

Target_Max = max(combTarget);

spectr = spectrumF(Psi);
u = abs(Psi);

% Fit the target as close as possible where meaningful
fitness_roi = dispParam.mu~=0 & combTarget.' >= Target_Max/1e4;
components.fit = 1e3/Target_Max * sum(( spectr(fitness_roi, :) - combTarget(fitness_roi) ).^2, 1).';
components.pulse = ( 5 * (median(u) + mad(u))./max(u) )';
components = struct2table(components);

score = components.fit; %sum(table2array(components),2);

end
