classdef Dispersion_Parametrization < handle
    
    % TODO: Make subclasses for the various cases
    
    properties
        mu
        F_max
        ParameterFitType
        paramTable % this is a be a table, where the columns contains the parameters useful for the LLE computation
        paramSettings % a table that contains the fields names, their static parameters in each fields (thus the size of the fields),  function of evolution, their post treatment function ...
        N_pop
        lleType
    end
    
    methods (Abstract, Access = protected)
        paramArray_init = fitInit(obj, Dint_smooth, PCS, comb_target, roi, m, m_dw);
        %paramsPopulate(obj, mutation_proba)
        Dint = dispersionFunction(obj, k);
        %[disp_params, scaling, descr_lst] = parametersDescriptor(obj, selected);
    end
    
    %     methods (Static, Access = protected)
    %         function t = tableCleanup(t)
    %             for i = 1:width(t)
    %                 if iscell(t.(i)) && any(cellfun(@(x)isnumeric(x) || islogical(x), t.(i)))
    %                     t.(i) = cell2mat(t.(i));
    %                 end
    %             end
    %         end
    %     end
    
    methods
        %% ==================================================
        function obj = Dispersion_Parametrization(mu, ParameterFitType, paramSettings, varargin)
            validFitType = {'poly','dualRing','DintVect', 'Dictionary', 'Dict_multiShifts'};
            checkFitType = @(x) any(validatestring(x,validFitType));
            valid_mu = @(x) isnumeric(mu) && any(mu==0);
            p = inputParser;
            p.KeepUnmatched=true;
            addRequired(p,'mu', valid_mu);
            addRequired(p,'ParameterFitType', checkFitType);
            addRequired(p,'paramSettings', @(x) ~isempty(x) && iscell(x));
            addParameter(p, 'N_pop', 1, @isnumeric)
            addParameter(p, 'F_max', [], @isnumeric)
            parse(p, mu, ParameterFitType, paramSettings, varargin{:})
            
            obj.mu = p.Results.mu;
            obj.ParameterFitType = p.Results.ParameterFitType;
            obj.F_max = p.Results.F_max;
            
            % Param settings defined as such: {'paramName', 'evolParamIndex', 'mutateFun', 'postMutateFun', 'descriptor', 'descriptorScaling'}
            % add the detuning and pump power parameters and parse as a table
            paramSettings = p.Results.paramSettings;
            paramSettings(end+1,:) = {'detuning', true, [], @detuningClamp, 'Detuning', 1}; % normalized units by default. TODO: implement to convert to real pysical units as well
            paramSettings(end+1,:) = {'pumpPow', true, [], @pumpPowClamp, 'Pump_Power', 1};
            obj.paramSettings = cell2table(paramSettings,...
                'VariableNames', {'paramName', 'evolParamIndex', 'mutateFun', 'postMutateFun', 'descriptor', 'descriptorScaling'}); % define the column names
            obj.paramSettings.Row = obj.paramSettings.paramName; % name the table rows with the paramName
            
            % initialize the parameter table based on the initialization
            obj.paramTable = [];
            for i = 1:height(obj.paramSettings)
                obj.paramTable = [obj.paramTable, table(zeros(size(obj.paramSettings.evolParamIndex{i})), 'VariableNames', obj.paramSettings.paramName(i))];
            end
            obj.N_pop = p.Results.N_pop;
            obj.lleType = 1; % by default 1 = single LLE, 2 = coupled LLE
        end
        % ==================================================
        
        %% ==================================================
        function [comb_target, pulseProfile_target,comb_init, pulseProfile_init, D_int_Kerr] = inverse_comb_dispersion(obj, Zeta, F, varargin)
            validTarget = {'pulse', 'spectrum'};
            checkTarget = @(x) any(validatestring(x,validTarget));
            
            p = inputParser;
            addRequired(p,'Zeta',@isnumeric);
            addRequired(p,'F',@isnumeric);
            addParameter(p,'TargetType','spectrum',checkTarget)
            addParameter(p,'FunctionTarget',[], @(f) isa(f, 'function_handle'))
            addParameter(p,'Scaling', .8)
            addParameter(p,'Verbose', false, @islogical) % TODO: handles the verbose
            parse(p,Zeta,F,varargin{:})
            p = p.Results;
            
            switch p.TargetType
                case 'pulse'
                    T = linspace(-pi, pi, length(obj.mu)+1);
                    T = T(1:end-1);
                    pulse = fftshift( p.FunctionTarget(T) );
                    comb_nopump = ifftshift(ifft(pulse));
                    
                case 'spectrum'
                    comb_nopump = p.FunctionTarget(obj.mu);
            end
            
            % rescale the amplitude
            comb_nopump = p.F * sqrt(p.Scaling /( 4 * sum(abs(comb_nopump(obj.mu~=0)).^2) )) * comb_nopump.';
            
            % Add a background
            pump = init_hom(p.Zeta, p.F^2, 1, false) * exp(-1i/2);
            
            % export the target comb
            comb_target = comb_nopump;
            comb_target(obj.mu==0) = pump;
            pulseProfile_target = fftshift(fft(fftshift(comb_target)));
            
            % compute the center of mass
            m = sum(obj.mu.' .* abs(comb_nopump).^2)/sum(abs(comb_nopump).^2);
            m_dw = sign(m) * max(obj.mu)/2;
            
            if abs(m)>2
                % define a 'Dispersive wave'
                s = .5; % arbitrary width
                DW = 1./(1+((obj.mu + m_dw)/s).^2);
                %DW = double(obj.mu == -round(ms*m));
                
                % scale the DW to move back the center of mass to the pump
                sc = fsolve(@(x) sum(obj.mu .* abs(comb_nopump.' + x * DW).^2), sum(abs(comb_nopump).^2));
                %sc = sum(obj.mu.' .* abs(comb_nopump).^2)./sum(obj.mu .* abs(lorentz).^2);
                %sc = roots([  ])
                comb_init = comb_nopump + sc * DW.';
                
            else
                comb_init = comb_nopump;
            end
            
            comb_init(obj.mu==0) = pump; % add pump
            
            % compute the pulse profile
            pulseProfile_init = fftshift(fft(fftshift(comb_init)));
            
            % Kerr shift computation
            [KerrShift,~] = computeKerrResponse(pulseProfile_target, 1);
            DintKerr = (KerrShift - KerrShift(obj.mu==0)); % flip the Kerr shift to get the dispersion
            
            % reject stuff too weak
            roi = 20*log10(abs(comb_target)) > -100; % this is an arbitrary criterion just to reject the values too low...
            
            figure
            hold on
            plot(obj.mu(roi), DintKerr(roi), '.')
            
            % smooth to reject fast deviations such as at the pump mode
            Dint_smooth = medfilt1(DintKerr,3);
            %Dint_smooth = Dint;
            %Dint_smooth(obj.mu==0) = median(Dint_smooth(abs(obj.mu)<5));
            
            % compute the PCS at the pump mode
            PCS = Dint_smooth(obj.mu==0) - DintKerr(obj.mu==0);
            
            % shift to zero the smoothed version for fitting
            Dint_smooth = Dint_smooth - Dint_smooth(obj.mu==0);
            
            % Now, fit and initialize the paramTable
            fitInit(obj, Dint_smooth, PCS, comb_target, roi, m, m_dw);
            
            ylim([max(min(ylim),-200), min(200, max(ylim))])
            
            D_int_Kerr.Dint_smooth = Dint_smooth;
            D_int_Kerr.KerrShift = DintKerr;
            
            % add the detuning and the pump power ...
            obj.paramTable{1,'detuning'} = Zeta;
            obj.paramTable{1,'pumpPow'} = F;
            % set population to 1
            obj.N_pop = 1;
            
        end
        % ==================================================
        
        
        %% ==================================================
        function [Dint_smooth] = fitTargetDintProfile(obj, Dint, Zeta, F, varargin)
            
            p = inputParser;
            p.KeepUnmatched=true;
            addRequired(p,'Dint',@isnumeric);
            addRequired(p,'Zeta',@isnumeric);
            addParameter(p,'roi',[],@(x) isnumeric(x) || islogical(x))
            parse(p,Dint,Zeta,varargin{:})
            p = p.Results;
            
            if isempty(p.roi)
                p.roi = 1:length(Dint);
            end
            
            % smooth to reject fast deviations such as at the pump mode
            Dint_smooth = medfilt1(Dint,3);
            %Dint_smooth = Dint;
            %Dint_smooth(obj.mu==0) = median(Dint_smooth(abs(obj.mu)<5));
            
            % compute the PCS at the pump mode
            PCS = Dint_smooth(obj.mu==0) - Dint(obj.mu==0);
            
            % shift to zero the smoothed version for fitting
            Dint_smooth = Dint_smooth - Dint_smooth(obj.mu==0);
            
            % Now, fit and initialize the paramTable
            fitInit(obj, Dint_smooth, PCS, ones(size(Dint_smooth)), p.roi, 0, 0, varargin{:});
            
            ylim([max(min(ylim),-200), min(200, max(ylim))])
            
            % add the detuning
            obj.paramTable{1,'detuning'} = Zeta;
            obj.paramTable{1,'pumpPow'} = F;
            % set population to 1
            obj.N_pop = 1;
            
        end
        % ==================================================
        
        
        %% ==================================================
        function initialize_population(obj, mutation_proba)
            initPopSize = height(obj.paramTable);
            
            % loop through the parameter variables in paramTable and init
            % the variables with their corresponding functions
            varNames = obj.paramSettings.paramName;
            
            for n_p = (initPopSize+1):obj.N_pop
                % add row in table
                obj.paramTable{n_p, :} = missing;
                
                for i = 1:length(varNames)
                    % reuse the initial population to populate the new
                    % variables
                    idx_ref =  mod(n_p-1, initPopSize) + 1;
                    funCall = obj.paramSettings{varNames{i},'mutateFun'}{:};
                    if isempty(funCall)
                        obj.paramTable{n_p, varNames{i}} = defaultMutateFun(obj, varNames{i}, n_p, idx_ref, mutation_proba);
                    else
                        obj.paramTable{n_p, varNames{i}} = funCall(obj, varNames{i}, n_p, idx_ref, mutation_proba);
                    end
                end
            end
            
            % post process to clamp if necessary
            obj.postMutationProcess()
            
        end
        
        %% default mutation function for pop init.
        function out = defaultMutateFun(obj, varName, n_p, idx_ref, mutation_proba)
            % the rest, I add random on the non static parameters
            mask = obj.paramSettings{varName,'evolParamIndex'}{:}; % deal with the static parameters
            N_params = sum(mask);
            mutate_prop = randn(1, N_params);
            mutate_sgn = double(rand(1, N_params)>mutation_proba)*2 - 1;
            sigmoid = @(x) 1/(1+exp(-(x-obj.N_pop*4/5-1)/(obj.N_pop/10))); % rescaling proba
            mutate_offs = sigmoid(n_p) * randn(1, N_params);
            
            out = obj.paramTable{idx_ref, varName}; % get and duplicate the reference in the existing population
            
            out(mask) = out(mask) .* mutate_sgn .* (1 + mutation_proba * mutate_prop) + ... % and modify it
                mutate_offs .* (.5*mutation_proba);
        end
        % ==================================================
        
        %% ==================================================
        % TODO: make sure the disperison calc can handle array and compute
        % for a single index with varagin
        function Dint = computeDispersion(obj, varargin)
            
            if ~isempty(varargin)
                idxLst = varargin{1};
            else
                idxLst = 1:obj.N_pop;
            end
            
            Dint = zeros(length(obj.mu), length(idxLst));
            
            for j = 1:length(idxLst)
                k = idxLst(j);
                Dint(:,j) = obj.dispersionFunction(k);
            end
            
        end
        % ==================================================
        
        %% ==================================================
        function obj = mutatePop(obj, fitnessArray, mutation_proba)
            % Cast the table into and array, where select only the varied
            % parameters
            paramArray = table2array(obj.paramTable);
            evolMask = obj.paramSettings.evolParamIndex;
            evolMask = horzcat(evolMask{:});
            
            paramArrayEvol = paramArray(:, evolMask); % keep the evolving data
            
            % evolve the population, crossover and mutate
            paramArrayEvol = create_new_pop(paramArrayEvol.', fitnessArray, mutation_proba);
            
            % Cast back into a table, taking the static params in account
            % parameters
            paramArray(:, evolMask) = paramArrayEvol.'; % place new values back into larger array
            
            %now... put things back into table format
            varNames = obj.paramTable.Properties.VariableNames;
            for i=1:length(varNames)
                nVar = length(obj.paramSettings{varNames{i}, 'evolParamIndex'}{:});
                idx = 1:nVar;
                obj.paramTable{:, varNames{i}} = paramArray(:, idx);
                paramArray(:, idx) = []; % remove from table so that the next one takes over
            end
            
            obj.postMutationProcess(); % takes care of applying the clamping functions ...
        end
        
        %% Post mutation adjustments
        function postMutationProcess(obj)
            % call the subfunctions that clamp the parameter / post-process
            
            varNames = obj.paramSettings.paramName;
            
            for i = 1:length(varNames) % loop over the parameters variables
                postFun = obj.paramSettings.postMutateFun{i};
                if ~isempty(postFun) % if empty, do nothing by default
                    obj.paramTable{:, varNames{i}} = postFun(obj); % Apply the function, hmmm looks like the other arguments are not really passed
                end
            end
            
        end
        
        % clamp the pump power to stay below the max
        function F = pumpPowClamp(obj)
            F = obj.paramTable{:,'pumpPow'};
            F(F>obj.F_max) = obj.F_max;
            F(F<1) = 1; % clamp the min power to 1
        end
        
        % clamp the detuning to stay within a reasonable range
        function Z = detuningClamp(obj)
            Z = obj.paramTable{:,'detuning'};
            for i=1:obj.N_pop
                F0 = obj.paramTable{i,'pumpPow'};
                UpperDetunBound = 1.5*F0^2;
                while Z(i) < 0 || Z(i) > UpperDetunBound
                    Z(i) = UpperDetunBound * rand(1); %UpperDetunBound/2 + randn();
                end
            end
        end
        % ==================================================
        
        %% output the parameters to a readable format
        function tbl = parametersDescriptor(obj, varargin)
            
            if isempty(varargin)
                selected = 1:obj.N_pop;
            else
                selected = varargin{1};
            end
            
            tmp = obj.paramTable(selected, :); % crop the table to the selected data
            
            % scale the parameters to their according weights
            for i=1:height(obj.paramSettings)
                varName = obj.paramSettings.paramName{i};
                scaling = cell2mat(obj.paramSettings{varName, 'descriptorScaling'});
                tmp{:, varName} = tmp{:, varName} .* scaling;
            end
            
            tmp = table2array(tmp);
            names = cellfun(@(S) obj.paramSettings{S, 'descriptor'}{:} , obj.paramTable.Properties.VariableNames, 'un',0);
            % Un-nest the cell array
            while any(cellfun(@iscell,names))
                names = [names{cellfun(@iscell,names)} names(~cellfun(@iscell,names))];
            end
            
            tbl = array2table(tmp, 'VariableNames', names);
        end
        
        
    end
    
end