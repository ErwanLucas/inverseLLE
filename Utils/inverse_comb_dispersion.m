function [dispParam, comb, pulseProfile, Dint_smooth] = inverse_comb_dispersion(mu, F, varargin)

p = inputParser;
validFitType = {'poly','dualRing','DintVect'};
checkFitType = @(x) any(validatestring(x,validFitType));
validTarget = {'pulse', 'spectrum'};
checkTarget = @(x) any(validatestring(x,validTarget));

addRequired(p,'mu',@isnumeric);
addRequired(p,'F',@isnumeric);
addParameter(p,'TargetType','spectrum',checkTarget)
addParameter(p,'FunctionTarget',[],@(f) isa(f, 'function_handle'))
addParameter(p,'Scaling', .8)
addParameter(p,'ParameterFitType','poly',checkFitType)
addParameter(p, 'Verbose', false)
parse(p,mu,F,varargin{:})


p = p.Results;


switch p.TargetType
    case 'pulse'
        T = linspace(-pi, pi, length(p.mu)+1);
        T = T(1:end-1);
        pulse = p.FunctionTarget(T);
        comb = ifftshift(ifft(pulse));
        
    case 'spectrum'
        comb = p.FunctionTarget(mu);
end

% rescale the amplitude
comb = p.F * sqrt(p.Scaling /( 4 * sum(abs(comb(mu~=0)).^2) )) * comb;

% Add a background
comb(mu==0) = exp(1i*pi/4) * p.F/10;

pulseProfile = fft(fftshift(comb));

[KerrShift,~] = computeKerrResponse(pulseProfile, 1);
Dint = -(KerrShift - KerrShift(mu==0)); % flip the Kerr shift to get the dispersion

% reject stuff too weak
roi = 20*log10(abs(comb)) > -80; % this is an arbitrary criterion just to reject the values too low...

figure
hold on
plot(mu(roi), Dint(roi), '.')

% smooth to reject fast deviations such as at the pump mode
Dint_smooth = medfilt1(Dint,3);

% compute the PCS at the pump mode
PCS = Dint_smooth(mu==0) - Dint(mu==0);

% shift to zero the smoothed version for fitting
Dint_smooth = Dint_smooth - Dint_smooth(mu==0);


%% Now, fit with a dual ring structure ...
% Parameter listing: D2_1, D2_2, G, DD1, DD0, PCS

switch p.ParameterFitType
    case 'poly'
        pp = polyfit(mu(roi), Dint_smooth(roi), 12); % TODO: make the poly order as parameter
        
        plot(mu, polyval(pp, mu)+PCS)
        
        dispParam = [pp, PCS];
        
    case 'dualRing'
        D2_1 = mean(diff(Dint_smooth(roi & abs(mu)>combTarget_width),2));
        G = combTarget_width;
        D2_2 = median(diff(Dint_smooth(roi & abs(mu)<combTarget_width),2));
        DD0 = 10;
        beta0 = [sign(D2_1), D2_1, D2_2, G, DD0];
        
        fit_fun = @(beta,x) coupledRingDispersion([beta,0], x);
        
        plot(mu,fit_fun(beta0, mu))
        
        opts = statset('nlinfit');
        %opts.RobustWgtFun = 'bisquare';
        opts.MaxFunEvals = 1e3;
        opts.MaxIter = 1e4;
        beta = nlinfit(mu(roi)',Dint_smooth(roi)', fit_fun, beta0, opts); %'Weights',W
        beta(1) = sign(beta(1)); % the branch selection
        
        plot(mu, coupledRingDispersion([beta, PCS], mu))
        
        dispParam = [beta, PCS];
        
    case 'DintVect'
        dispParam = Dint;
end

ylim([max(min(ylim),-200), min(200, max(ylim))])

end