classdef Dispersion_Parametrization_Poly < Dispersion_Parametrization
       
    properties
        polyOrder
        polyScaling
        octaveMode
    end
    
    methods
        %% ==================================================
        function obj = Dispersion_Parametrization_Poly(mu, varargin)
            p = inputParser;
            p.KeepUnmatched=true;
            addRequired(p,'mu');
            addParameter(p, 'polyOrder', 12, @isnumeric)
            addParameter(p, 'polyScaling', [], @isnumeric)
            parse(p, mu, varargin{:})
            
            % Param settings defined as such: {'paramName', 'evolParamIndex(true) or static(false)', 'mutateFun', 'postMutateFun', 'descriptor', 'descriptorScaling'}
            paramSettings = {
                'polyCoefs', true(1,p.Results.polyOrder+1), @mutatePolyCoefs, [], cellstr('D_' + string(p.Results.polyOrder:-1:0)')', factorial(p.Results.polyOrder:-1:0); ...
                'PCS', true, [], [], 'PCS', 1
                };
            
            obj@Dispersion_Parametrization(mu, 'poly', paramSettings, varargin{:})
            obj.polyOrder = p.Results.polyOrder;
            obj.polyScaling = p.Results.polyScaling;
            obj.octaveMode = 0;
        end
        % ==================================================
        
        %% Convert the scaled polynomial coefficients to standard polynomial coeffs (has an effect if polyScaling is different than [0,1])
        function retval = getPolyCoeffs(obj, varargin)
            % convert 3 return parameter polynomial fit vector to the type returned if
            % one return value is requested:
            % p1 = polyfit(x,y,n);
            % [p2,S,mu] = polyfit(x,y,n);
            % p3 = polyfit_convert(p2);
            %  =>  p1 == p3
            
            if isempty(varargin)
                selected = 1:obj.N_pop;
            else
                selected = varargin{1};
            end
            
            p = obj.paramTable{selected,'polyCoefs'};
            n = obj.polyOrder;
            m = obj.polyScaling(1);
            s = obj.polyScaling(2);
            
            retval = zeros(size(p));
            for h = 1:size(p,1)
                for i = 0:n
                    for j = 0:i
                        retval(h, n+1-j) = retval(h, n+1-j) + p(h, n+1-i) * nchoosek(i,j)*(-m)^(i-j)/s^i;
                    end
                end
            end
        end
        
    end
    
    methods  (Access = protected)
        %% ==================================================
        function fitInit(obj, Dint_smooth, PCS, comb_target, roi, m, m_dw)
            
            [p0, ~,obj.polyScaling] = polyfit(obj.mu(roi), Dint_smooth(roi), obj.polyOrder); % TODO: make the poly order as parameter
            
            polyeval = @(p) polyval(p,obj.mu(roi), [], obj.polyScaling);
            %rescale = @(x) (x - obj.polyScaling(1))/obj.polyScaling(2);
            %scale = @(x) x * obj.polyScaling(2) + obj.polyScaling(1);
            
            if abs(m)>2 % DW adding
                fprintf('Asymmetric problem, adding a DW at %g\n', -m_dw)
                p0(1) = []; % reduce the order of the free poly
                polyFunc = @(p) polyeval(p) .* ( obj.mu(roi) + m_dw );
                objectiveFunc = @(p) sum(abs(comb_target(roi).').^2 .* (Dint_smooth(roi) - polyFunc(p)).^2);
                options = optimset('MaxFunEvals', 1e4);
                pp = fminsearch( objectiveFunc, p0, options );
                
                % add a root to trigger a DW
                dwpp = [obj.polyScaling(2), m_dw + obj.polyScaling(1)];
                pp = conv(dwpp,pp);
            else
                objectiveFunc = @(p) sum(abs(comb_target(roi).').^2 .* (Dint_smooth(roi) - polyeval(p)).^2);
                options = optimset('MaxFunEvals', 1e4);
                pp = fminsearch( objectiveFunc, p0, options );
            end
            
            plot(obj.mu(roi), Dint_smooth(roi) + PCS)
            plot(obj.mu, polyval(pp, obj.mu, [], obj.polyScaling) + PCS)
            
            % initialize the paramters
            obj.paramTable{1,'polyCoefs'} = pp;
            obj.paramTable{1,'PCS'} = PCS;
        end
        % ==================================================
        
        
        %% ==================================================
        function out = mutatePolyCoefs(obj, varName, n_p, idx_ref, mutation_proba)
            
            %varName = 'polyCoefs';
            polyscale = 0:(obj.polyOrder);
            polyscale = fliplr( polyscale );
            normArray = 1./factorial(polyscale);
            
            % the rest, I add random on the non static parameters
            mask = obj.paramSettings{varName,'evolParamIndex'}{:}; % deal with the static parameters
            N_params = sum(mask);
            mutate_prop = randn(1, N_params);
            mutate_sgn = double(rand(1, N_params)>mutation_proba)*2 - 1;
            sigmoid = @(x) 1/(1+exp(-(x-obj.N_pop*4/5-1)/(obj.N_pop/10))); % rescaling proba
            mutate_offs = sigmoid(n_p) * randn(1, N_params);
            polyscale = polyscale(mask);
            normArray = normArray(mask);
            
            out = obj.paramTable{idx_ref, varName}; % get the reference in the existing population
            
            out(mask) = out(mask) .* mutate_sgn .* (1 + mutation_proba * mutate_prop) + ... % and modify it
                normArray .* mutate_offs .* (.5*mutation_proba).^polyscale;
            
        end
        % ==================================================
        
        %% ==================================================
        % TODO: make sure the dispersion calc can handle array and compute
        % for a single index with varagin
        function Dint = dispersionFunction(obj, k)
            Dint = polyval(obj.paramTable{k, 'polyCoefs'}, obj.mu, [], obj.polyScaling);
            sel = obj.mu==0 | obj.mu==obj.octaveMode;
            Dint(~sel) = Dint(~sel) + obj.paramTable{k, 'PCS'};
        end
        % ==================================================
        
        
    end
    
end