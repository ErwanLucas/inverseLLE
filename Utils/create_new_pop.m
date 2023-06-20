function new_param = create_new_pop(paramArray,fitnessArray, mutation_proba)

N_pop = size(paramArray,2);

average_param_scale = sum(fitnessArray.' .* paramArray,2)./sum(fitnessArray);

%% rank the individues by their fitness
[~,order] = sort(fitnessArray, 'ascend'); % larger fitness gets a larger ranking
ranking = zeros(length(fitnessArray),1);
ranking(order) = (1:length(fitnessArray))';
score = 1./sqrt(ranking); % This has the effect of inverting the scaling of the fitness the lowest 'fitness' individuals are ranked higher
score = 2*N_pop/sum(score) * score;
%sortrows([fitnessArray,ranking,score],1)

make_cum_distr = @(sco) [0,cumsum(sco/sum(sco))'];

%% init new population
new_param = zeros(size(paramArray));

% keep the fittest individual untouched
NhandPicked = 1;
[~, iM] = max(score);
new_param(:,1) = paramArray(:,iM);

% the rest is mutated
N_parents = 3;
for i=1:ceil(N_pop/N_parents)
    
    % could implement a tournament for selecting parents
    
    % check https://www.tutorialspoint.com/genetic_algorithms/genetic_algorithms_parent_selection.htm
    % Selection of "parents" N_parents
    % this is Fitness proportionate selection
    sel_P = discretize(rand(1,N_parents), make_cum_distr(score));
    sel_parents = paramArray(:,sel_P);
    
    % crossover between parents parameters, based on the fitness
    prm_lngth = size(sel_parents,1);
    
    % select the parent for each gene source
    parent4gene = repmat(1:N_parents, prm_lngth, 1);
    sel1 = rand(prm_lngth,1);
    sel1 = discretize(sel1, make_cum_distr(score(sel_P))); % replace cp with [0, N_parents]/N_parents for uniform sampling
    parent4gene = cell2mat(arrayfun(@(k) circshift(parent4gene(k,:),1-sel1(k)),1:prm_lngth,'uni',0)');
    
    for j=1:N_parents
        
        k = (i-1)*N_parents + j + NhandPicked;
        if k>N_pop
            break
        end
        
        I_sel = sub2ind(size(sel_parents), (1:prm_lngth)', parent4gene(:,j));
        offspring = sel_parents(I_sel); % crossover
        
        % could otherwise implement a weighted average...
        %p = sum( score(sel_parents).' .* p, 2)./sum(score(sel_parents));
        
        % add random variations, either proportional or based on the
        % average
        if rand() >= .5
            mutation_prop = double(rand(prm_lngth,1) <= mutation_proba) * mutation_proba .* randn(prm_lngth,1);
            offspring = offspring.*(1 + mutation_prop) ;
        else
            mutationR_pos = (rand(prm_lngth,1) <= mutation_proba);
            mutationReplace = average_param_scale .* randn(prm_lngth,1);
            offspring(mutationR_pos) = mutationReplace(mutationR_pos);
        end
        
        % add in the new population
        new_param(:,k) = offspring;
        
    end
    
end
end
