function [omega, disp_curves] = coupledRingDispersion(beta, mu, varargin)



if numel(beta)>5
    varargin{1} = beta(1);
    beta(1) = [];
end
beta = num2cell(beta);
[D2_1, D2_2, G, DD0, PCS] = beta{:};
DD1 = 0;

branch = 'top';
if ~isempty(varargin)
    if isnumeric(varargin{1})
        if varargin{1}>=0
            branch = 'top';
        else
            branch = 'bottom';
        end
    else
        branch = varargin{1};
    end
end

omega_A = D2_1/2 * mu.^2;
omega_B = D2_2/2 * mu.^2 + DD1 * mu + DD0;

D = sqrt(G^2 + (omega_A-omega_B).^2/4);
omega_P = (omega_A + omega_B)/2 + D;
omega_M = (omega_A + omega_B)/2 - D;

switch branch
    case 'top'
        omega = omega_P;
    case 'bottom'
        omega = omega_M;
end

omega = omega - omega(mu==0);
omega(mu~=0) = omega(mu~=0) + PCS;

varNames = {'omega_P', 'omega_M', 'omega_A', 'omega_B'};
disp_curves = cell2table(cellfun(@eval, varNames, 'un',0), 'VariableNames', varNames);

end

