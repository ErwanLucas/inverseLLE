 % Genetic dispersion optimisation for a target comb profile
addpath(genpath('./Utils/'))
clearvars
close all
clc
warning on
rng('shuffle')

testingToggle = false; % just runs the initial solution computed from the inverse Kerr shift

%% Setup the dispersion parametrization
run init_polyCombTarget.m
% run init_octaveDW.m

%% Genetic algo parameters
dispParam.N_pop = 208; % population of each generation, should be even in order to try positive / negative dispersion
maxIter = 300; % max number of generation
fitnessTol = 1e-6; % fitness target
mutation_proba = 0.15;

%% Setup the initial population based on the seed
dispParam.initialize_population(mutation_proba);
%dispParam.paramArray(3,3:end) = dispParam.paramArray(3,3:end) .* (double(rand(1, dispParam.N_pop-2)>.5)*2-1); % random sign flip

% plot the initial dispersion
figure
hold all
set(gca, 'ColorOrder', hsv(dispParam.N_pop))
Dint = dispParam.computeDispersion();
plot(mu,Dint(:,1), 'LineWidth', 2)
plot(mu,Dint(:,2:end))
ylim([-1,1]*200)

%% Run
if testingToggle
    nWorker = 0;
    popRange = 1;
else
    nWorker = feature('numcores');
    poolObj = parpool('local', nWorker);
    popRange = 1:dispParam.N_pop;
end

% setting up storage
dispParam_Storage = cell(maxIter,1);

% Storage defined based on the type of LLE
Psi = zeros(N_modes, dispParam.lleType, dispParam.N_pop, maxIter);

%Psi = randn(N_modes, dispParam.N_pop, maxIter);
var_GL = zeros(dispParam.N_pop, maxIter);
conv_flag = zeros(dispParam.N_pop, maxIter);
fitnessArray = fitnessTol * ones(dispParam.N_pop, maxIter);

i_gen = 1; % generation number

tst = datetime('now'); tst.Format = 'yyyy_MM_dd-HH.mm.ss';
tmpFile = sprintf('./Data/tmp_%s.mat', tst);
tic
% Stop criterion: fitness good enough and number of generations
while i_gen <= maxIter && ~any(fitnessArray(:,max(1,i_gen-1))<fitnessTol)
    
    % Store the current coefficients
    dispParam_Storage{i_gen} = dispParam.paramTable;
    
    % Propagate the states
    tic
    parfor(i_p = popRange, nWorker)
        %for i_p = 1
        % Setup the dispersion for the LLE run
        Dint = dispParam.computeDispersion(i_p);
        D = -1 - 1i * Dint; % in normalized units
        D = ifftshift(D);
        % initialize the seed (in one direction)
        zeta = dispParam.paramTable{i_p, 'detuning'};
        F0 = dispParam.paramTable{i_p, 'pumpPow'};
        try
            psi0 = pulseProfile_init;
        catch
            psi0 = initialize_pulse(dispParam, i_p, F0, startBranch, excite);
        end
        
        % Define the LLE type and initialization
        lle = @(psi, t_evol, h, Nstore) LLE_Propagate_para(psi, F0, zeta, D, nf, t_evol, h, Nstore);
        
        % actually run
        % first, propagate for a bit
        t_evol = 40;
        Nstore = 1e3;
        h = 2^-9;
        [Psi_tmp, var_GL(i_p,i_gen), ~] = lle(psi0, t_evol, h, Nstore);
        % refine the solution by decreasing step size
        t_evol = 10;
        Nstore = 1e3;
        h_large = h;
        h_fine = 2^-14;
        h = @(t) (h_large-h_fine) * (1+tanh(15*(.4 - t/t_evol)))/2 + h_fine;
        Psi(:,:,i_p,i_gen) = lle(Psi_tmp, t_evol, h, Nstore);
        
        %% test plots
        if testingToggle
            %plot
            Psi_end = Psi(:,:,i_p,i_gen);
            disp(fitness(Psi_end))
            
            figure
            hold on
            plot(abs(Psi_end))
            plot(abs(Psi_tmp))
            plot(abs(psi0))
            try
                plot(abs(pulseProfile_target))
            catch
                disp('no comb shape target set')
            end
            
            figure
            hold on
            plot(mu, pow2db(spectrumF(squeeze(Psi_end))), 'r', 'DisplayName', 'final')
            plot(mu, pow2db(spectrumF(squeeze(Psi_tmp))), 'b', 'DisplayName', 'coarse evol')
            plot(mu, pow2db(spectrumF(psi0)), 'k', 'DisplayName', 'init')
            %plot(mu,10*log10(double(roi)+eps))
            if strcmpi(dispParam.ParameterFitType, 'Dictionary')
                try
                    vline(dispParam.getTargetMuDW())
                catch
                    vline(dispParam.octaveMode)
                end
            end
            try
                plot(mu, pow2db(combTarget), '.', 'DisplayName', 'target')
            catch
                disp('no comb shape target set')
            end
            legend show
            axis tight
        end
        
    end
    fprintf('End propagation gen %i in %g s\n', i_gen, toc)
    
    if testingToggle
        return
    end
    
    % Fitness evaluation
    tmp_fitn = zeros(dispParam.N_pop, size(Psi, 2));
    for i_dir = 1:size(Psi, 2)
        tmp_fitn(:, i_dir) = fitness(squeeze(Psi(:,i_dir,:,i_gen)));
    end
    fitnessArray(:,i_gen) = min(tmp_fitn, [], 2);
    
    % generate the new population based on fitness
    dispParam.mutatePop(fitnessArray(:,i_gen), max( 0.01, mutation_proba * tanh(3*(1-i_gen/maxIter)) ));
    
    i_gen = i_gen + 1; % increment
    
    save(tmpFile, '-v7.3')
    
end

delete(poolObj);
clear poolObj

%% save the data
tst = datetime('now'); tst.Format = 'yyyy_MM_dd-HH.mm.ss';
save(fullfile('./Data', sprintf('%s_Genetic_evolution.mat', tst)), '-v7.3')
delete(tmpFile);

disp('DONE')
close all
clearvars
return