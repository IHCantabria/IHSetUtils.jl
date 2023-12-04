using Distributions
using Optim
using LinearAlgebra
using Random
using Base.Threads
# using CUDA


###########################################################
###########################################################
#######################NGSA2###############################
###########################################################
###########################################################

function nsga2(obj_fun, npar, ngen, npop, ub, lb, nobj)
    
    # Initialize population of candidate solutions
    population = randn(npop, npar)
    
    objectives = zeros(npop,nobj)
    # Evaluate initial population
    @threads for i in 1:npop
        objectives[i,:] .= obj_fun(((population[i,:] .- (-1)) ./ (1 .- (-1))).*(ub.-lb).+lb)
    end
    
    # objectives = obj_fun.(population)
    
    # Initialize iteration counter
    t = 1
    
    while t <= ngen
        # Select parents for breeding
        parents = tournament_selection(population, objectives)
        
        # Create offspring using genetic operators
        offspring = genetic_operators(parents)
        noff = size(offspring,1)
        
        offspring_obj = zeros(noff, nobj)
        # Evaluate offspring
        # offspring_obj = obj_fun.(offspring)
        @threads for i in 1:noff
            offspring_obj[i,:] .= obj_fun(((offspring[i,:].- (-1)) ./ (1 .- (-1))).*(ub.-lb).+lb)
        end
        
        # Update population with offspring and original solutions
        population, objectives = non_dominated_sort(
            vcat(population, offspring),
            vcat(objectives, offspring_obj)
        )
        # Check length of sorted population and add random individuals if necessary
        pop_len = size(population, 1)
        if pop_len < npop && t != ngen
            num_missing = npop - pop_len
            new_pop = randn(num_missing, npar)
            new_obj = zeros(num_missing, nobj)
            @threads for i in 1:num_missing
                new_obj[i,:] .= obj_fun(((new_pop[i,:].- (-1)) ./ (1 .- (-1))).*(ub.-lb).+lb)
            end
            population = vcat(population, new_pop)
            objectives = vcat(objectives, new_obj)
        elseif pop_len > npop
            population = population[1:npop, :]
            objectives = objectives[1:npop, :]
        end
        
        # Update iteration counter

        if t % Int(ngen/20) == 0
            println("Progress = $(100 * t / ngen) %")
            println("Generation: $t / $ngen")
        end
        t += 1

    end
    for i in axes(population, 1)
        population[i,:] = ((population[i,:].- (-1)) ./ (1 .- (-1))).*(ub.-lb).+lb
    end

    return population, objectives
end
export nsga2

###########################################################
###########################################################
#######################SPEA2###############################
###########################################################
###########################################################

function spea2(obj_fun, npar, ngen, npop, ub, lb, nobj)
    
    # Initialize population of candidate solutions
    population = rand(npop, npar)
    
    # Evaluate initial population
    objectives = zeros(npop, nobj)
    @threads for i in 1:npop
        objectives[i,:] .= obj_fun(population[i,:] .*(ub.-lb).+lb)
    end
    
    # Initialize iteration counter
    t = 1
    ff_val = zeros(size(objectives,1))

    while t <= ngen

        # Environmental selection
        archive, archive_fitness, _ = environmental_selection(population, objectives)
        # println("*******************************************")
        # println(archive_fitness)
        # println("*******************************************")
        
        # Select parents for breeding
        parents  = tournament_selection2(archive, archive_fitness)
        # parents = copy(archive)
        
        # println(idxSel)
        # println(size(parents))
        # Create offspring using genetic operators
        offspring = genetic_operators2(parents)
        
        # Evaluate offspring
        offspring_obj = zeros(size(offspring,1), nobj)
        @threads for i in axes(offspring, 1)
            offspring_obj[i,:] .= obj_fun(offspring[i,:] .*(ub.-lb).+lb)
        end
        
        # Combine parents, offspring, and archive
        combined_pop = vcat(population, offspring)
        combined_obj = vcat(objectives, offspring_obj)
        
        

        # Environmental selection to form next generation population
        population, ff_val ,objectives = environmental_selection(combined_pop, combined_obj)
        
        # # Sorting by fitness
        # idxSort = sortperm(ff_val)
        # population, objectives = population[idxSort], objectives[idxSort]
        
        pop_len = size(population, 1)

        if pop_len < npop && t != ngen
            num_missing = npop - pop_len
            new_pop = rand(num_missing, npar)
            new_obj = zeros(num_missing, nobj)
            @threads for i in 1:num_missing
                new_obj[i,:] .= obj_fun(new_pop[i,:] .*(ub.-lb).+lb)
            end
            population = vcat(population, new_pop)
            objectives = vcat(objectives, new_obj)
        elseif pop_len > npop
            population = population[1:npop, :]
            objectives = objectives[1:npop, :]
        end
        
        # Update iteration counter
        if t % Int(ngen/20) == 0
            println("Progress = $(100 * t / ngen) %")
            println("Generation: $t / $ngen")
        end
        t += 1

    end

    for i in eachindex(population[:,1])
        population[i,:] .= population[i, :] .* (ub.-lb).+lb
    end

    return population, objectives, ff_val
end
export spea2

###########################################################
###########################################################
#################Genetic Operations########################
###########################################################
###########################################################


function tournament_selection(population, objectives)
    npop, npar = size(population)
    num_selected = convert(Int, npop / 2)
    selected = zeros(num_selected, npar)
    
    for i in 1:num_selected
        competitors = rand(1:npop, 2, 1)
        if dominates(objectives[competitors[1], :], objectives[competitors[2], :])
            selected[i, :] = population[competitors[1], :]
        else
            selected[i, :] = population[competitors[2], :]
        end
    end
    
    return selected
end

function genetic_operators(parents)
    num_offspring, npar = size(parents)
    offspring = zeros(num_offspring, npar)
    
    for i in 1:2:num_offspring
        # Select parents to breed
        p1 = parents[i, :]
        p2 = parents[i + 1, :]
        
        # Apply crossover operator
        o1, o2 = crossover(p1, p2)
        
        # Apply mutation operator
        o1 = mutation(o1)
        o2 = mutation(o2)
        
        offspring[i, :] = o1
        offspring[i + 1, :] = o2
    end
    
    return offspring
end

function crossover(p1, p2)
    alpha = rand()
    o1 = alpha * p1 + (1 - alpha) * p2
    o2 = (1 - alpha) * p1 + alpha * p2
    return o1, o2
end

function dominates(obj1, obj2)
    # Determine if obj1 dominates obj2
    # A solution dominates another if it is better in all objectives
    # and not worse in at least one objective.
    
    dominates = true
    worse = false
    
    for i in eachindex(obj1)
        if obj1[i] > obj2[i]
            dominates = false
            break
        elseif obj1[i] < obj2[i]
            worse = true
        end
    end
    
    return dominates && worse
end

function non_dominated_sort(population, objectives)
    # println("**********************************")
    # println("OBJS - RP, RMSE, MSS")
    # println(objectives)
    # println("PARAMS")
    # println(population)
    # println("**********************************")
    # Sort population into non-dominated fronts
    fronts = sort_fronts(objectives)
    # println(fronts)

    # Select individuals from each front
    selected = select_individuals(fronts)
    # println(selected)
    
    # Update population and objectives arrays
    NEW_pop = population[selected, :]
    NEW_obj = objectives[selected, :]
    return NEW_pop, NEW_obj
end

function sort_fronts(objectives)
    # Sort population into non-dominated fronts
    npop = size(objectives, 1)

    dominated_by = zeros(Int, npop)
    n_dominated = zeros(Int, npop)

    fronts = [[]]

    for p in 1:npop
        for q in 1:npop
            if p == q
                continue
            end

            if dominates(objectives[p, :], objectives[q, :])
                dominated_by[q] += 1
            elseif dominates(objectives[q, :], objectives[p, :])
                n_dominated[p] += 1
            end
        end

        if n_dominated[p] == 0
            push!(fronts[1], p)
        end
    end

    i_front = 1
    while !isempty(fronts[i_front])
        next_front = []
        for p in fronts[i_front]
            for q in findall(dominated_by .== p)
                n_dominated[q] -= 1
                if n_dominated[q] == 0
                    push!(next_front, q)
                end
            end
        end
        push!(fronts, next_front[:])
        i_front += 1
    end
    
    # frt = 1
    # for i in eachindex(fronts)
    #     frt = [frt; Int.(fronts[i])]
    # end

    # return frt[2:end]

    return fronts[1:end-1]
end

function select_individuals(fronts)
    # Select individuals from each front
    selected = []
    for front in fronts
        n = length(front)
        if n > 0
            if n <= length(selected)
                selected[1:n] .= front
            else
                selected = [selected; Int.(front)]
            end
        end
    end
    
    return selected
end

function mutation(o)
    # Apply mutation operator to offspring
    # Add Gaussian noise with mean 0 and standard deviation 1 to each gene
    mutated = o .+ randn(size(o)) .* 0.2
    return mutated
end

function pairwise(fn, A, B)
    [fn(a,b) for a in eachrow(A), b in eachrow(B)]
end

function Euclidean(a,b)
    norm(a-b)
end

function environmental_selection(population, objectives)
    npop = size(population, 1)
    fitness_values = zeros(npop)
    # Calculate distances between solutions in the objective space
    dist = pairwise(Euclidean, objectives, objectives)
    k = mean(dist)/5
    m=2
    dist[Diagonal(trues(npop,npop))] .= Inf
    # Calculate raw fitness values for each solution
    for i in 1:npop
        fitness_values[i] = sum(exp.(-1 .* dist[i,:].^2 ./ (2 .* k^2)))
    end
    # Calculate the density estimate for each solution
    density_estimate = zeros(npop)
    for i in 1:npop
        sorted_distances = sort(dist[i,:])
        density_estimate[i] = 1 / (2 * k) * (1 / sorted_distances[m+1])
    end
    # Calculate combined fitness value for each solution
    combined_fitness = fitness_values .+ density_estimate
    # Environmental selection to form archive
    archive = population[findall(combined_fitness .<= median(combined_fitness)), :]
    archive_obj = objectives[findall(combined_fitness .<= median(combined_fitness)), :]
    archive_fitness = combined_fitness[findall(combined_fitness .<= median(combined_fitness))]


    return archive, archive_fitness, archive_obj, combined_fitness
end

function tournament_selection2(population, objectives)
    npop, npar = size(population)
    selected = zeros(npop, npar)
    idxSel = zeros(npop) 

    for i in 1:npop
        competitors = sample(1:npop, 2, replace=false)
        if dominates(objectives[competitors[1], :], objectives[competitors[2], :])
            selected[i, :] = population[competitors[1], :]
            idxSel[i] = competitors[1]
        else
            selected[i, :] = population[competitors[2], :]
            idxSel[i] = competitors[2]
        end
    end

    unique_indices = indexin(unique(selected[:,1]), selected[:,1])

    # _, unique_indices = unique(selected, dims=1)
    selected = selected[unique_indices, :]

    return selected
end

function tournament_selection3(population, fitness_values, tournament_size=2)
    npop, npar = size(population)
    selected = zeros(npop, npar)
    
    for i in 1:npop
        competitors = sample(1:npop, tournament_size, replace=false)
        best_index = argmin(fitness_values[competitors])
        selected[i, :] = population[competitors[best_index], :]
    end
    
    # remove individuals with duplicated values
    unique_indices = indexin(unique(selected[:,1]), selected[:,1])

    # _, unique_indices = unique(selected, dims=1)
    selected = selected[unique_indices, :]
    
    return selected
end

function genetic_operators2(parents)
    num_parents, npar = size(parents)
    num_offspring = num_parents - mod(num_parents, 2) # garante um número par de filhos
    offspring = zeros(num_offspring, npar)
    
    for i in 1:2:num_offspring
        # Select parents to breed
        p1 = parents[i, :]
        p2 = parents[i + 1, :]
        
        # Apply crossover operator
        o1, o2 = crossover(p1, p2)
        
        # Apply mutation operator
        o1 = mutation2(o1)
        o2 = mutation2(o2)
        
        offspring[i, :] = o1
        offspring[i + 1, :] = o2
    end
    
    # Se num_parents for ímpar, gerar um indivíduo adicional
    if num_offspring < num_parents
        p = parents[end, :]
        o = mutation2(p)
        offspring = vcat(offspring, o')
    end
    
    return offspring
end

function nsga2_lim(obj_fun, npar, ngen, npop, ub, lb, nobj)
    
    # Initialize population of candidate solutions
    population = rand(npop, npar)
    
    objectives = zeros(npop,nobj)
    # Evaluate initial population
    @threads for i in 1:npop
        objectives[i,:] .= obj_fun(((population[i,:] .*(ub.-lb).+lb)))
    end
    
    # objectives = obj_fun.(population)
    
    # Initialize iteration counter
    t = 1
    
    while t <= ngen
        # Select parents for breeding
        parents = tournament_selection(population, objectives)
        
        # Create offspring using genetic operators
        offspring = genetic_operators3(parents)
        noff = size(offspring,1)
        
        offspring_obj = zeros(noff, nobj)
        # Evaluate offspring
        # offspring_obj = obj_fun.(offspring)
        @threads for i in 1:noff
            offspring_obj[i,:] .= obj_fun(((offspring[i,:].*(ub.-lb).+lb)))
        end
        
        # Update population with offspring and original solutions
        population, objectives = non_dominated_sort(
            vcat(population, offspring),
            vcat(objectives, offspring_obj)
        )
        # Check length of sorted population and add random individuals if necessary
        pop_len = size(population, 1)
        if pop_len < npop && t != ngen
            num_missing = npop - pop_len
            new_pop = rand(num_missing, npar)
            new_obj = zeros(num_missing, nobj)
            @threads for i in 1:num_missing
                new_obj[i,:] .= obj_fun(((new_pop[i,:].*(ub.-lb).+lb)))
            end
            population = vcat(population, new_pop)
            objectives = vcat(objectives, new_obj)
        elseif pop_len > npop
            population = population[1:npop, :]
            objectives = objectives[1:npop, :]
        end
        
        # Update iteration counter

        if t % Int(ngen/20) == 0
            println("Progress = $(100 * t / ngen) %")
            println("Generation: $t / $ngen")
        end
        t += 1

    end
    for i in axes(population, 1)
        population[i,:] = population[i,:].*(ub.-lb).+lb
    end

    return population, objectives
end

function genetic_operators3(parents)
    num_offspring, npar = size(parents)
    offspring = zeros(num_offspring, npar)
    
    for i in 1:2:num_offspring
        # Select parents to breed
        p1 = parents[i, :]
        p2 = parents[i + 1, :]
        
        # Apply crossover operator
        o1, o2 = crossover(p1, p2)
        
        # Apply mutation operator
        o1 = mutation2(o1)
        o2 = mutation2(o2)
        
        offspring[i, :] = o1
        offspring[i + 1, :] = o2
    end
    
    return offspring
end

function mutation2(o)
    # Apply mutation operator to offspring
    mutated = o .+ randn(size(o)) .* 0.2
    clamp!(mutated, 0, 1)
    return mutated
end
###########################################################
###########################################################
########################PSO################################
###########################################################
###########################################################

function pso(f, lb, ub, n_particles, n_iterations, w, c1, c2)
    # Initialize particles
    particles = rand(Float64, n_particles, length(lb)) .* (ub .- lb) .+ lb
    velocities = zeros(Float64, n_particles, length(lb))
    pbest = particles
    pbest_fitness = [f(p) for p in pbest]
    gbest_index = argmin(pbest_fitness)
    gbest = copy(pbest[gbest_index, :])
    gbest_fitness = pbest_fitness[gbest_index]
    # Iterate
    for t = 1:n_iterations
        for i = 1:n_particles
            # Update velocity
            velocities[i, :] = w .* velocities[i, :] .+ 
                               c1 .* rand(Float64, length(lb)) .* (pbest[i, :] .- particles[i, :]) .+
                               c2 .* rand(Float64, length(lb)) .* (gbest .- particles[i, :])
            # Apply velocity limits
            velocities[i, :] = max.(velocities[i, :], lb .- particles[i, :])
            velocities[i, :] = min.(velocities[i, :], ub .- particles[i, :])
            # Update position
            particles[i, :] .+= velocities[i, :]
            # Apply position limits
            particles[i, :] = max.(particles[i, :], lb)
            particles[i, :] = min.(particles[i, :], ub)
            # Evaluate fitness
            fitness = f(particles[i, :])
            # Update personal best
            if fitness < pbest_fitness[i]
                pbest[i, :] = copy(particles[i, :])
                pbest_fitness[i] = fitness
            end
            # Update global best
            if fitness < gbest_fitness
                gbest = copy(particles[i, :])
                gbest_fitness = fitness
            end
        end
        # Print progress
        if t % Int(n_iterations/20) == 0
            println("Progress = $(100 * t / n_iterations) %")
        end
        
    end
    return gbest, gbest_fitness
end


###########################################################
###########################################################
#########################SA################################
###########################################################
###########################################################


function simulated_annealing(obj_func, lb, ub, args...; max_iter=1000, t0=1.0, alpha=0.95, verbose=false)
    # inicialização
    n_params = length(lb)
    x_curr = rand(lb:ub, n_params)
    f_curr = obj_func(x_curr, args...)
    f_best = f_curr
    x_best = x_curr
    t = t0
    iter = 0

    while iter < max_iter
        # perturbação
        x_proposed = x_curr + randn(n_params) .* t
        x_proposed = max.(lb, min.(ub, x_proposed)) # garante que x_proposed está dentro do intervalo lb e ub
        f_proposed = obj_func(x_proposed, args...)

        # aceitação/rejeição
        delta_f = f_proposed - f_curr
        if delta_f < 0 || rand() < exp(-delta_f / t)
            x_curr = x_proposed
            f_curr = f_proposed
            if f_curr < f_best
                f_best = f_curr
                x_best = x_curr
            end
        end

        # atualiza a temperatura e a contagem de iterações
        t *= alpha
        iter += 1

        # imprime informações
        if verbose && iter % 100 == 0
            println("iter=%d, f_best=%g\n", iter, f_best)
        end
    end

    return x_best, f_best
end

function DE(f::Function, lb::Vector, ub::Vector, npop::Int, ngen::Int;
        F::Float64=0.8, CR::Float64=0.9)
    dim = length(lb)
    pop = [rand(lb[i]:ub[i], npop) for i in 1:dim]
    best = pop[1]
    fbest = f(best)
    for i in 1:ngen
        for j in 1:npop
            r1, r2, r3 = sample([k for k in 1:npop if k ≠ j], 3, replace=false)
            x = pop[r1] + F * (pop[r2] - pop[r3])
            x = min.(max.(x, lb), ub)
            u = rand(dim) .< CR
            y = [u[k] ? x[k] : pop[j][k] for k in 1:dim]
            fy = f(y)
            if fy < fbest
                best = y
                fbest = fy
            end
            pop[j] = y
        end
    end
    return best, fbest
end


function ga(fitness::Function, lb::AbstractVector, ub::AbstractVector, n_pop::Int, n_gen::Int,
    pc::Float64=0.8, pm::Float64=0.1, elitism::Bool=true)

n_var = length(lb)
pop = rand(lb:ub, n_pop, n_var)
fitness_values = [fitness(pop[i, :]) for i in 1:n_pop]

if elitism
 elite_idx = argmin(fitness_values)
 elite = pop[elite_idx, :]
end

for i in 1:n_gen
 # Seleção
 fitness_prob = fitness_values ./ sum(fitness_values)
 parents_idx = sample(1:n_pop, n_pop, weights=fitness_prob, replace=true)
 parents = pop[parents_idx, :]
 
 # Cruzamento
 masks = rand(Bool, n_pop, n_var)
 children = zeros(Int, n_pop, n_var)
 for j in 1:2:n_pop-1
     if rand() < pc
         for k in 1:n_var
             if masks[j, k]
                 children[j, k] = parents[j, k]
                 children[j+1, k] = parents[j+1, k]
             else
                 children[j, k] = parents[j+1, k]
                 children[j+1, k] = parents[j, k]
             end
         end
     else
         children[j, :] = parents[j, :]
         children[j+1, :] = parents[j+1, :]
     end
 end
 
 # Mutação
 for j in 1:n_pop
     for k in 1:n_var
         if rand() < pm
             children[j, k] = rand(lb[k]:ub[k])
         end
     end
 end
 
 # Substituição
 pop = vcat(parents, children)
 fitness_values = [fitness(pop[i, :]) for i in 1:size(pop, 1)]
 
 if elitism
     best_idx = argmin(fitness_values)
     if fitness_values[best_idx] < fitness(elite)
         pop[best_idx, :] = elite
         fitness_values[best_idx] = fitness(elite)
         elite = pop[best_idx, :]
     end
 end
end

if elitism
 return elite
else
 best_idx = argmin(fitness_values)
 return pop[best_idx, :]
end
end
