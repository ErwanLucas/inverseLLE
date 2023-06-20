function [score, components] = fitness_superOctave(Psi, dispParam, varargin)
% Fit the target as close as possible where meaningful
%fitness_roi = mu~=0 & combTarget.' >= Target_Max/1e4;
%fitness = @(spectr) 1/Target_Max * sum((spectr(fitness_roi, :) - combTarget(fitness_roi)).^2, 1);

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

p = inputParser;
addRequired(p,'Psi');
addRequired(p,'dispParam');
addOptional(p, 'selected', 1:dispParam.N_pop)
addParameter(p, 'TargetDW', [], @isnumeric)
parse(p, Psi, dispParam, varargin{:})

selected = p.Results.selected;
targetDW = p.Results.TargetDW;

if isnumeric(selected)
    N_slect = length(unique(selected));
elseif islogical(selected)
    N_slect = sum(selected);
else
    error('Unrecognised selection option used index array or logical mask')
end
    

if size(Psi, 2) < N_slect
    error('Number of selected individuals is smaller than the size of Psi')
end

spectr = spectrumF(Psi);
u = abs(Psi);

if isempty(targetDW)
    muTarget = dispParam.getTargetMuDW(selected);
else
    muTarget = targetDW*ones(length(selected), 1);
end

II = knnsearch(dispParam.mu', muTarget);
if size(Psi,2)>1
    II = sub2ind(size(Psi), II(selected)', selected);
else
    II = II(selected);
end

components.power = 0.1 * ( max(0,-log(spectr(II))) )';
components.pulse = ( 10 * (median(u) + mad(u))./max(abs(u)) )';
components = struct2table(components);

score =  sum(table2array(components),2);

end
