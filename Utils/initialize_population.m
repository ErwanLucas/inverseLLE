function [paramArray,staticParams] = initialize_population(paramArray_init, staticParamsID, ParameterFitType)

switch ParameterFitType
    case 'poly'
        
        
        
    case 'dualRing'
        paramArray = paramArray_init;
        
        staticParams = paramArray(staticParamsID,:); % this is the branch selection <0 bottom, >=0 top
        staticParams = repmat(staticParams, 1, (N_pop)/2);
        paramArray(staticParamsID,:) = [];
        
        % detuning definition: replace the sign colunm with detuning and send it at the end
        paramArray(1,:) = Zeta;
        paramArray = circshift(paramArray, -1, 1); % send the detuning at the end
        
        % the rest, I divide the dispersion and add random
        if N_pop>2
            for n_p=3:N_pop
                mutate = randn(size(paramArray,1)-1, 1);
                idx_ref = mod(n_p-1,2)+1;
                paramArray(:,n_p) = [...
                    paramArray(1:end-1, idx_ref) .* (1 + 2*mutation_proba*mutate) + mutation_proba*mutate;...
                    Zeta + 2*randn()];
            end
        end
        
        
end
