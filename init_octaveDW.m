%% Parameter settings (N modes, pump power & co)
F0_Max = sqrt(50);      % Set the pump here, and let the system accomodate. dimensionless forcing term [ sqrt(P/P_thresh) ]
N_modes = 3500; %   % number of points
startBranch = 'bottom'; % initialise the comb from the top branch. The spontaneous formation should make the desired pulse. (althought it will be stochastic here to some extent)
excite = 'soliton';         % hard excitation yes/no
nf = 0;%1e-4;          % noise factor
mu = (0:N_modes-1)-N_modes/2; % mode index

%% Initialize the dispersion parameters object
polyOrder = 4;
dispParam = Dispersion_Parametrization_Poly(mu, 'polyOrder', polyOrder, 'F_Max', F0_Max);
dispParam.polyScaling = [0,1];
dispParam.octaveMode = 969;

% detuning definition
Zeta = F0_Max^2 - 2;
PCS = 5;

%% Load the sample data and fit the init poly to a realistic dispersion
% Initialize based on roots
d2 = 5e-3;
p = conv([0,-d2/dispParam.octaveMode,0,0], [1,-dispParam.octaveMode]);

%% parameters init append the PCS and detuning
dispParam.paramTable.polyCoefs = p;
dispParam.paramTable.PCS = PCS;
dispParam.paramTable.detuning = Zeta;
dispParam.paramTable.pumpPow = F0_Max;
%dispParam.paramTable{1,:} = [p, PCS, Zeta];

% fix the offset of the poly to avoid redundancy w/ PCS
dispParam.paramTable.polyCoefs(1,end-1:end) = 0;
dispParam.paramSettings{'polyCoefs', 'evolParamIndex'}{:}(end-1:end) = false;
dispParam.N_pop = 1;

figure
hold on
plot(mu, dispParam.computeDispersion())
drawnow()
%plot(mu, polyval(PolyCoefs, mu) + PCS, '--' )

%% Define the fitness
fitness = @(Psi) fitness_superOctave(Psi, dispParam, 1:size(Psi,2), 'TargetDW', dispParam.octaveMode);
