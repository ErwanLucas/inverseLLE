function [psi, t, var_genL, genL] = LLE_Propagate_para(psi, fpump, zeta, D, nf, t_evol, t_step, Nstore, varargin)

%% Parse parameters
p = inputParser;
addParameter(p, 'Coupling', [], @isnumeric)
addParameter(p, 'XPMBiDir', true, @islogical)
addParameter(p, 'storeAllEvo', false, @islogical)
addParameter(p, 'reportOption', 'none')
parse(p, varargin{:})
coupling = p.Results.Coupling;
withXPMBiDir = p.Results.XPMBiDir;

% Setup the reporting options
reporting = true;
reportCond = max(1, round(Nstore / 200));
switch p.Results.reportOption
    case 'plot'
        reporter = llePlotter(psi, zeta, fpump, t_evol);
    case 'waitbar'
        reporter = lleProgress(t_evol);
    otherwise
        reporting = false;
end

% initialize operators variables
D_hat = [];
C_hat = [];

%% Manage errors
N_modes = size(psi, 1);
NDir = size(psi, 2);
if ~isempty(coupling) && NDir < 2
    error('Coupling condition is activated and the field array dimension is too small')
end

if not(isa(t_step, 'function_handle') || isnumeric(t_step) && length(t_step) == 1)
    error('Time step must be defined as a function of time or as a constant')
end

if not(isa(zeta, 'function_handle') || isnumeric(zeta) && length(zeta) == 1)
    error('Detuning must be defined as a function of time or as a constant')
end

if not(isa(fpump, 'function_handle') || isnumeric(fpump) && length(fpump) == NDir)
    error('The pump term must be defined as a function of time or as a constant')
end

%% Manage the evolution cases
isTstepVariable = isa(t_step, 'function_handle');
if ~isTstepVariable
    h = t_step;
    updateOperators()
end

% Variable detuning?
isZetaVariable = isa(zeta, 'function_handle');
if ~isZetaVariable
    Z = zeta;
end

% Variable pump power?
F_mask = zeros(N_modes, 1);
F_mask(1) = 1;
F_mask = fft(F_mask);
isFVariable = isa(fpump, 'function_handle');
if ~isFVariable
    F = fpump;
    % driving pump definition
    FVec =  F_mask.* F;
end

% storage timestamps
t_store_flags = linspace(0, t_evol, Nstore);

storeFinalOnly = ~(isZetaVariable || isFVariable || p.Results.storeAllEvo);

%% generated light formula
gen_light = @(X) sum(abs(X - mean(X)).^2);
genL = zeros(Nstore, NDir);
genL(1, :) = gen_light(psi);

if ~storeFinalOnly
    t_store = zeros(Nstore, 1);
    Psi = zeros(length(psi), NDir, Nstore);
    Psi(:, :, 1) = psi;
    disp('All the evolution will be stored, no tolerance stop criterion')
end

%% initialize
psi_old = psi;
t_current = 0;
j = 2;

% init time step
if isTstepVariable
    h = t_step(t_current);
    updateOperators()
end

%% Main split step loop
while ~stopCondition(t_current, t_evol, psi, psi_old)

    t_current = t_current + h; % loop increment

    if isZetaVariable
        Z = zeta(t_current);
    end

    if isFVariable
        F = fpump(t_current);
        FVec = F_mask .* F;
    end

    % store previous state
    psi_old = psi;

    % split step process
    psi = LLE_splitstep(psi, exp(-1i * Z * h)*D_hat, C_hat, FVec, 1, h, withXPMBiDir);

    % storage and update
    if t_current >= t_store_flags(j)
        %fprintf('%.3g %% complete\n',t_current/t_evol*100);
        genL(j, :) = gen_light(psi);

        if ~storeFinalOnly
            Psi(:, :, j) = psi;
            t_store(j) = t_current;
        end

        % Update the time step if active
        if isTstepVariable
            h = t_step(t_current);
            updateOperators() % refresh the operators with the new step
        end

        j = j + 1; % increment

        if reporting && mod(j, reportCond) == 0
            reporter.report(t_current, psi)
        end
    end

    % noise addition
    if nf ~= 0
        noise = nf * h * randn(size(psi)) .* exp(1i*pi*randn(size(psi)));
        psi = psi + noise;
    end

end

%% end main loop

% average the last 2 steps, reduces a bit the noise
psi = (psi + psi_old) / 2;
t = t_current;

% generated light stdev for stability metric (takes only half of the evolution)
var_genL = std(genL(end / 2:end, :)) / mean(genL(end / 2:end, :));

% In case full evolution stored
if ~storeFinalOnly
    psi = Psi;
    t = t_store;
end

%% Stop loop condition, the criterion is hard coded, should be variable
    function test = stopCondition(t_current, t_evol, psi, psi_old)
        if storeFinalOnly
            if t_current == 0
                test = false;
            else
                test = t_current >= t_evol || norm(psi-psi_old) < 1e-10;
            end
        else
            test = t_current >= t_evol;
        end
    end

%% THIS Function updates the pre-computed dispersion and coupling operators
% avoids computing the exp of array at each step for speedup
    function updateOperators()
        D_hat = exp(h*D); % update dispersion operator
        C_hat = [cos(h * coupling / 2), 1i * sin(h * coupling / 2)]; %Same: precompute the coupling operator
    end

end
