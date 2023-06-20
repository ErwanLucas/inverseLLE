close all
clearvars
clc

addpath('./Utils/')
postLoad = '';

fName = '2021_03_30-10.40.06_Genetic_evolution'; % Enter the filename of the genetic evolution result file (without .mat extension)


load(['./Data/' fName '.mat'])
eval(postLoad)

%% Evaluate the convergence
figure
semilogy(min(fitnessArray(:,1:i_gen-1)))
hold on
semilogy(mean(fitnessArray(:,1:i_gen-1)))

grid on
box on
xlabel('Generation #')
ylabel('Error (fitness)')
title(sprintf('Population size = %d', dispParam.N_pop))
legend({'min', 'mean'})
set(gcf, 'color','w')

%% Initial dispersion
dispParam.paramTable = dispParam_Storage{1};
Dint_init = dispParam.computeDispersion(1);

%% Find the optimum
ftnArraySrch = fitnessArray(:, 1:i_gen-1);
fitVals = sort(ftnArraySrch(:));
II = find(ftnArraySrch(:)==fitVals(1), 1, 'last'); % search for the nth lowest fitness
[indiv_opt, gen_opt] = ind2sub(size(fitnessArray),II)

dispParam.paramTable = dispParam_Storage{gen_opt};

if ndims(Psi) == 4
    Psi_opt = squeeze(Psi(:,:,indiv_opt,gen_opt));
else
    Psi_opt = squeeze(Psi(:,indiv_opt,gen_opt));
end

zeta = dispParam.paramTable.detuning(indiv_opt);

if ismember('pumpPow', dispParam.paramTable.Properties.VariableNames)
    F0 = dispParam.paramTable.pumpPow(indiv_opt); % New version when pump power varied as well
    F0 = F0 * [1,0];
end


% plot the comb spectrum
fg(1) = figure('Tag', 'CombOptimum');
hold on
h_plot=stem(mu,pow2db(spectrumF(Psi_opt)), 'BaseValue', -200, 'Marker', 'none', 'DisplayName', 'Result');
plot(mu,pow2db(spectrumF(squeeze(Psi(:,1,1)))), 'k--', 'DisplayName', 'Initial')
% try
%     plot(mu, pow2db(combTarget+eps), 'r-', 'DisplayName', 'Target')
% catch
%     disp('No Comb target specified')
%     try
%         vline(dispParam.getTargetMuDW())
%     catch
%         vline(969)
%     end
% end
axis tight
xlim([-1,1]*65) ; ylim([-100, max(ylim)])

% Plot the target power + BW
vline([-1,1]*target_halfWidthModes)
hline(target_powerdB)
%hline(target_powerdB + [-1,1]*target_varPower/2)

xlabel('Mode #')
ylabel('Power (norm. unit)')

% plot the pulse profile
fg(2) = figure('Tag', 'PulseOptimum');
hold on
plot(mu,abs(Psi_opt), 'DisplayName', 'Result');
plot(mu,abs(squeeze(Psi(:,1,1))), 'k:', 'DisplayName', 'initial')
axis tight
ylim(ylim + [0,1]*.1*diff(ylim))
try
    if exist('pulseProfile_target','var')
        plot(mu, abs(pulseProfile_target), 'r', 'DisplayName', 'Target')
    else
        plot(mu, abs(pulse), 'r', 'DisplayName', 'Target')
    end
catch
    disp('No Comb target specified')
end

% Post process the plots (title and so on ...)
for fig = fg
    figure(fig)
    box on
    leg = legend('show');
    title(gca, sprintf('\\zeta_0 = %.3g, F^2 = %.3g', zeta, F0.^2));
    savefig(fig, sprintf('./Figures/%s_%s.fig', tst, fig.Tag))
end

figure,
hold on
imagesc(var_GL(:, 1:i_gen-1))
axis xy
xlabel('Generation #')
ylabel('Population #')
plot(gen_opt,indiv_opt, 'wx')
axis tight


%% plot the dispersion
figure
hold on
Dint = dispParam.computeDispersion(indiv_opt);
plot(mu, Dint, '.-')
plot(mu, Dint_init, '-.k')
axis tight
lim_clamp = 500;
ylim([max(-lim_clamp, min(ylim)), min(lim_clamp, max(ylim))])
ylim(ylim + [-1,1]*.1*diff(ylim))
hline(0, 'k--')

xlabel('Mode #')
ylabel('Normalized deviation D_{int}')
box on

savefig(gcf, sprintf('./Figures/%s_DintOptimum.fig', tst))

% export the dispersion
T = table(mu', Dint, 'VariableNames', {'Mode_Number', 'Normalized_Deviation'});
writetable(T, sprintf('./Data/%s_Optim_Dispersion.txt',tst))

% export the parameters
T = dispParam.parametersDescriptor(indiv_opt)
% Manipulate a bit
writetable(T, sprintf('./Data/%s_Optim_Dispersion_Coeffs.txt',tst))


%% Query for scan computation
disp(repmat('=', 1, 50))
disp('Program Paused, press key to initiate simulation')
pause

%% Test the dispersion profile in a scan to see if it can be reached
% Tune some split step properties
nf = 1e-4;
D = -1 - 1i * Dint; % in normalized units
D = ifftshift(D);

t_evol = 500;
h = 2^-10;
Nstore = 2e3;
z_start = -1;
z_stop = (F0^2 + 2);
zrate = abs(z_start - z_stop)/t_evol;
zeta = @(t) z_start + (z_stop - z_start)  * t/t_evol;

% initialize the seed
psi0 = init_hom(zeta(1), F0^2, 1) * ones(1,N_modes);

tic
Psi_evol_fwd = LLE_Propagate_para(psi0.', F0, zeta, D, nf, t_evol, h, Nstore);
toc
Psi_evol_fwd = squeeze(Psi_evol_fwd);

mnvl = std( max(abs(Psi_evol_fwd),[],2) );
pks=arrayfun(@(i)findpeaks(abs(Psi_evol_fwd(:,i)), 'MinPeakProminence', 3*mnvl), 1:size(Psi_evol_fwd,2), 'un',0);
pks = cellfun(@numel, pks);

%% Find the section with single pulse
zpos = find(~[0 pks==1 0]);
[~, grpidx] = max(diff(zpos));
alph = .9;
restart_idx = round( (1-alph)*zpos(grpidx) + alph * (zpos(grpidx+1)-2) );

figure
imagesc(abs(Psi_evol_fwd))
vline(restart_idx, 'w')

figure
spectr = pow2db(abs(spectrumF(Psi_evol_fwd)));
imagesc(spectr)

figure
p=panel();
p.pack('h',2)
p(1).select()
mag = abs(Psi_evol_fwd(:,restart_idx));
[~,iM] = max(mag);
plot(fftshift(circshift(mag, -iM)))
axis tight
p(2).select()
stem(mu, spectr(:,restart_idx), 'BaseValue', -140, 'Marker', 'none')
axis tight


%% Tune backward
r = round(t_evol/h);
z_start = zeta(restart_idx/Nstore*t_evol);
z_stop = -1;
t_evol = abs(z_start - z_stop)/zrate;
zeta = @(t) z_start + (z_stop - z_start)  * t/t_evol;

tic
Psi_evol_bckw = LLE_Propagate_para(Psi_evol_fwd(:,restart_idx), F0, zeta, D, nf, t_evol, h, Nstore);
toc
Psi_evol_bckw = squeeze(Psi_evol_bckw);

figure
imagesc(abs(Psi_evol_bckw))

figure
spectr = pow2db(abs(spectrumF(Psi_evol_bckw)));
imagesc(spectr)

%% Concatenate the 2 maps
Psi_FB = [Psi_evol_fwd(:,1:restart_idx), Psi_evol_bckw];

[minFit,best_fit_idx] = min(fitness(Psi_FB));

figure
imagesc(abs(Psi_FB))
colormap(morgenstemning)
vline(best_fit_idx, 'w')

figure
spectr = pow2db(abs(spectrumF(Psi_FB)));
plot(mu, spectr(:,best_fit_idx))


return

%% Now we can also explore the other states
% Alternate metric computation
u = zeros(size(fitnessArray));
components = {};
for i_gen = 1:maxIter
    dispParam.paramTable = dispParam_Storage{i_gen};
    [u(:,i_gen), components{i_gen}] = fitness_superOctave(squeeze(Psi(:,:,:,i_gen)), dispParam, 1:size(Psi,3),'TargetDW',dispParam.octaveMode);
end
%% plot
metricPow = cell2mat(cellfun(@(x) x.power, components, 'UniformOutput', false));
metricPulse = cell2mat(cellfun(@(x) x.pulse, components, 'UniformOutput', false));

[minM, II] = min(metricPow(:));
[indiv_opt, gen_opt] = ind2sub(size(u),II);
fprintf('Max power criterion: fitn %g, powMetric %g, pulseMetric %g\n', u(indiv_opt, gen_opt),  metricPow(indiv_opt, gen_opt), metricPulse(indiv_opt, gen_opt))
Psi_opt = squeeze(Psi(:,:,indiv_opt,gen_opt));
figure
stem(mu,pow2db(spectrumF(Psi_opt)), 'Basevalue', -150, 'Marker', 'none')
axis tight
ylim([-150, 0])
dispParam.paramTable = dispParam_Storage{gen_opt};
yyaxis right
plot(mu, dispParam.computeDispersion(indiv_opt)), ylim([-1,1]*1e3)
try
    vline(dispParam.getTargetMuDW(indiv_opt))
catch
    vline(969)
end

[minM, II] = min(metricPulse(:));
[indiv_opt, gen_opt] = ind2sub(size(u),II);
fprintf('Max pulse criterion: fitn %g, powMetric %g, pulseMetric %g\n', u(indiv_opt, gen_opt),  metricPow(indiv_opt, gen_opt), metricPulse(indiv_opt, gen_opt))
Psi_opt = squeeze(Psi(:,:,indiv_opt,gen_opt));
figure
stem(mu,pow2db(spectrumF(Psi_opt)), 'Basevalue', -150, 'Marker', 'none')
axis tight
ylim([-150, 0])
dispParam.paramTable = dispParam_Storage{gen_opt};
yyaxis right
plot(mu, dispParam.computeDispersion(indiv_opt)), ylim([-1,1]*1e3)
try
    vline(dispParam.getTargetMuDW(indiv_opt))
catch
    vline(969)
end

[minM, II] = min(u(:));
[indiv_opt, gen_opt] = ind2sub(size(u),II);
fprintf('Best fitness: fitn %g, powMetric %g, pulseMetric %g\n', u(indiv_opt, gen_opt),  metricPow(indiv_opt, gen_opt), metricPulse(indiv_opt, gen_opt))
Psi_opt = squeeze(Psi(:,:,indiv_opt,gen_opt));
figure
hold on
stem(mu,pow2db(spectrumF(Psi_opt)), 'Basevalue', -150, 'Marker', 'none')
plot(mu, pow2db(spectrumF(squeeze(Psi(:,1,1)))))
axis tight
ylim([-150, 0])
dispParam.paramTable = dispParam_Storage{gen_opt};
yyaxis right
plot(mu, dispParam.computeDispersion(indiv_opt)), ylim([-1,1]*1e3)
try
    vline(dispParam.getTargetMuDW(indiv_opt))
catch
    vline(969)
end

%% Plot the optimum of first round
[F,IM]=min(fitnessArray(:,1))
h_plot=stem(mu,pow2db(spectrumF(squeeze(Psi(:,IM,1)))), 'BaseValue', -150, 'Marker', 'none');

