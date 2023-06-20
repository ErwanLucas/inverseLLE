%% Parameter settings (N modes, pump power & co)
N_modes = 2^9;% 3500; %   % number of points
startBranch = 'bottom'; % initialise the comb from the top branch. The spontaneous formation should make the desired pulse. (althought it will be stochastic here to some extent)
excite = 'soliton';         % hard excitation yes/no
nf = 0;%1e-4;          % noise factor
mu = (0:N_modes-1)-N_modes/2; % mode index

polyOrder = 6;
F0_Max = sqrt(20);
Zeta = F0_Max^2 / 3;
target_varPower = 5; % xx dB variation allowed around the target power
target_powerdB = -19 + target_varPower/2;
target_halfWidthModes = 32;

%% Define the Dispersion_Parametrization
%staticParamsID list the index ot the parameters that should not be varied
dispParam = Dispersion_Parametrization_Poly(mu, 'polyOrder', polyOrder, 'F_Max', F0_Max);

%% Setup a comb target

[combTarget, pulseProfile_target, comb_init, pulseProfile_init, Dint_kerr] = dispParam.inverse_comb_dispersion(Zeta, F0_Max, ...
    'TargetType', 'spectrum', ...
    'FunctionTarget', @(x) 1./(1+((x)/target_halfWidthModes).^8) + eps, ... %@(x) shapingFunctions(x, 'raisedCos', 80,0.8), ...  % %@(x) exp(-(x/20).^2), ... @(x) exp(-(10*x).^2), ... % , ... % @(x) shapingFunctions(x-30, 'raisedCos', 80,0.8), ...% 
    'Scaling', .1 ...
    );
combTarget = abs(combTarget).^2;

% do not vary the odd terms (symmetric enforcement)
dispParam.paramSettings{'polyCoefs', 'evolParamIndex'}{:}(end-1:-2:1)=false;
% Do no very the pump power

plot(mu, dispParam.computeDispersion(), '--')
hold on
legend({'Computed inv. Kerr shift', 'smoothed', 'PolyFit', "Final dispersion" + newline + "parametrization"})
%plot(mu,Dint_kerr.KerrShift)
xlabel('Mode #')
ylabel('Normalized deviation')
box on
axx = gca;

figure
hold on
stem(mu, max(-100,pow2db(combTarget)), 'Marker', 'none', 'BaseValue', -100, 'color', 'r')
axis tight
xlim(axx.XLim)
xlabel('Mode #')
ylabel('Normalized power (dB)')
box on

%% Define the fitness
%fitness = @(Psi) fitness_fitTarget(Psi, dispParam, combTarget);
%fitness = @(Psi) fitness_maxLinePower(Psi, dispParam, target_powerdB, target_halfWidthModes); % fit for the most power in the bandwidth
fitness = @(Psi) fitness_LinePowerTarget(Psi, dispParam, target_powerdB, target_varPower, target_halfWidthModes); % fit for power target and a min flatness