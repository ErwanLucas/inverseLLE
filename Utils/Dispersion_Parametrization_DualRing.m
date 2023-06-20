classdef Dispersion_Parametrization_DualRing < Dispersion_Parametrization
    
    properties
        polyOrder
        polyScaling
    end
    
    methods
        %% ==================================================
        function obj = Dispersion_Parametrization_DualRing(mu, varargin)
            p = inputParser;
            p.KeepUnmatched=true;
            addRequired(p,'mu');
            parse(p, mu, varargin{:})
            
            % Param settings defined as such: {'paramName', 'staticParamIndex', 'mutateFun', 'postMutateFun', 'descriptor', 'descriptorScaling'}
            paramSettings = {
                'D2_1', true, [], [], 'D2 ring 1', 1
                'D2_2', true, [], [], 'D2 ring 2', 1
                'DD0', true, [], [], 'resonance offset 1-2 at pump', 1
                'G', true, [], [], 'Inter ring coupling', 1
                'PCS', true, [], [], 'PCS', 1
                };
            
            obj@Dispersion_Parametrization(mu, 'dualRing', paramSettings, varargin{:})
        end
        % ==================================================
        
        
    end
    
    methods  (Access = protected)
        %% ==================================================
        function fitInit(obj, Dint_smooth, PCS, comb_target, roi, m, m_dw)
            % Parameter listing: D2_1, D2_2, G, DD1, DD0, PCS
            D2_1 = mean(diff(Dint_smooth(roi),2));
            G = combTarget_width;
            D2_2 = median(diff(Dint_smooth(roi),2));
            DD0 = 10;
            beta0 = [sign(D2_1), D2_1, D2_2, G, DD0];
            
            fit_fun = @(beta,x) coupledRingDispersion([beta,0], x);
            
            plot(obj.mu,fit_fun(beta0, obj.mu))
            
            opts = statset('nlinfit');
            %opts.RobustWgtFun = 'bisquare';
            opts.MaxFunEvals = 1e3;
            opts.MaxIter = 1e4;
            beta = nlinfit(obj.mu(roi)',Dint_smooth(roi)', fit_fun, beta0, opts); %'Weights',W
            beta(1) = sign(beta(1)); % the branch selection
            
            plot(obj.mu, coupledRingDispersion([beta, PCS], obj.mu))
            
            obj.paramTable{1,:} = [beta, PCS];
        end
        % ==================================================
        
        
        %% ==================================================
        % TODO: make sure the dispersion calc can handle array and compute
        % for a single index with varagin
        function Dint = dispersionFunction(obj, k)
            disp_params = table2array(obj.paramTable).';
            Dint = coupledRingDispersion(disp_params(:,k), obj.mu);
        end
        % ==================================================
        
        
    end
    
end