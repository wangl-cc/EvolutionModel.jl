export CoevoPrPdModel

"""
    CoevoPrPdModel{I<:Integer, R<:Real} <: AbstractModel

A model sumulates the coevolation between prey and predator.
"""
struct CoevoPrPdModel{I<:Integer, R<:Real, FA<:Function, FB<:Function} <: AbstractModel
    "sizes of prey populations"
    X::Vector{I}
    "sizes of predator populations"
    Y::Vector{I}
    "baseline birthrate of prey populations"
    bX::R
    "baseline deathrate of prey population"
    dX::R
    "baseline deathrate of prey populations"
    dY::R
    "coefficients scaling prey growing ability"
    g::Vector{R}
    "ratio of predator growth to predation"
    k::Vector{R}
    "baseline k value"
    K::R
    "growth-defense trade-off coefficient"
    m::R
    "resource competition coefficient matrix of prey populations"
    payoff::Matrix{R}
    "coefficient scaling the intensity of competition"
    M::R
    "coefficient of the predation rate"
    p::R
    "mutation rate of prey populations"
    X_mutrate::R
    "mutation rate of predator populations"
    Y_mutrate::R
    "a function descirbing how mutations change the prey population"
    X_mutfunc::FA
    "a function descirbing how mutations change the pred population"
    Y_mutfunc::FB
end

function CoevoPrPdModel(X::AbstractVector{<:Integer},
                        Y::AbstractVector{<:Integer},
                        bX::Real,
                        dX::Real,
                        dY::Real,
                        g::AbstractVector{<:Real},
                        k::AbstractVector{<:Real},
                        K::Real,
                        m::Real,
                        payoff::AbstractMatrix{<:Real},
                        M::Real,
                        p::Real,
                        X_mutrate::Real,
                        Y_mutrate::Real,
                        X_mutfunc::Function,
                        Y_mutfunc::Function)
    return CoevoPrPdModel(promote(X, Y)...,
    promote_(bX, dX, dY, g, k, K, m, payoff, M, p, X_mutrate, Y_mutrate)...,
    X_mutfunc, Y_mutfunc)
end

function CoevoCompModel(model::CoevoPrPdModel)
    @copyfields X, bX, dX, g, payoff, M, X_mutrate = model
    X_mutfunc = model.X_mutfunc
    return CoevoCompModel(X, bX, dX, g, payoff, M, X_mutrate, X_mutfunc)
end

function gillespie(model::CoevoPrPdModel{I, R, <:Function, <:Function}, T::Real)where {I<:Integer, R<:Real}
    t = zero(R) # initialize t

    # initialize populations
    X_populations_current = Population.(t, model.X)
    Y_populations_current = Population.(t, model.Y)
    X_populations_all = copy(X_populations_current)
    Y_populations_all = copy(Y_populations_current)

    # destuct model parameters
    @copyfields bX, dX, dY, g, k, K, payoff, M, p, m, X_mutrate, Y_mutrate = model
    g_history = copy(g)
    k_history = copy(k)
    X_nomutrate = 1 - X_mutrate
    Y_nomutrate = 1 - Y_mutrate
    repM = 1/M
    X_mutfunc = model.X_mutfunc
    Y_mutfunc = model.Y_mutfunc

    # main loop
    while t <= T
        # check extinction
        length(X_populations_current) <= 0 && length(Y_populations_current) <=0 && break # both prey and predator extinction
        length(Y_populations_current) <= 0 && break # predator extinction

        # get current population sizes
        X_num = [getfield(p, :n) for p in X_populations_current]
        Y_num = [getfield(p, :n) for p in Y_populations_current]

        # transpose need vector
        trans_X_num = transpose(X_num)
        trans_g = transpose(g)

        # calculate reactions
        ## prey
        X_birth_nomut = @. bX * X_num * X_nomutrate * g
        X_birth_mut = @. bX * X_num * X_mutrate * g
        X_intrinsic_death = @. dX * X_num
        X_competitiondeath = @. trans_X_num / payoff * X_num * repM

        ## predator
        Y_intrinsic_death = @. dY * Y_num

        ## predations
        predatarion_all = @. Y_num * p * trans_g ^ (m * k / K) * trans_X_num

        predatation_nobirth = @. (1 -k) * predatarion_all
        predatation_birth_nomut = @. k * Y_nomutrate * predatarion_all
        predatation_birth_mut = @. k * Y_mutrate * predatarion_all

        # calculate τ and choose reaction
        τ, i, index = choosereaction(X_birth_nomut, X_birth_mut, X_intrinsic_death, X_competitiondeath,
            Y_intrinsic_death, predatation_nobirth, predatation_birth_nomut, predatation_birth_mut)

        t += τ # time increase

        # solve reactions
        if i == 1 # prey birth no mutate
            birth!(X_populations_current[index[1]], t)
        elseif i == 2 # prey birth mutate
            index = index[1]
            newpop = Population(t, 1)
            push!(X_populations_current, newpop)
            push!(X_populations_all, newpop)
            g, payoff = X_mutfunc(g, payoff, index)
            push!(g_history, g[end])
        elseif i in (3, 4) # prey death (intrinsic and competition)
            index = index[1]
            death!(X_populations_current[index], t)
            if X_populations_current[index].n == 0 # check extinction
                deleteat!(X_populations_current, index)
                deleteat!(g, index)
                payoff = payoff[1:size(payoff, 1) .!= index, 1:size(payoff, 2) .!= index]
            end
        elseif i == 5 # predator intrinsic death
            index = index[1]
            death!(Y_populations_current[index], t)
            if Y_populations_current[index].n == 0 # check extinction
                deleteat!(Y_populations_current, index)
                deleteat!(k, index)
            end
        elseif i == 6 # predatation without birth
            # prey death
            index = index[2]
            death!(X_populations_current[index], t)
            if X_populations_current[index].n == 0 # check extinction
                deleteat!(X_populations_current, index)
                deleteat!(g, index)
                payoff = payoff[1:size(payoff, 1) .!= index, 1:size(payoff, 2) .!= index]
            end
        elseif i == 7 # predatation with birth no mutate
            # predator birth
            birth!(Y_populations_current[index[1]], t)

            # prey death
            death!(X_populations_current[index[2]], t)
            if X_populations_current[index[2]].n == 0 # check extinction
                deleteat!(X_populations_current, index[2])
                deleteat!(g, index[2])
                payoff = payoff[1:size(payoff, 1) .!= index[2], 1:size(payoff, 2) .!= index[2]]
            end
        elseif  i == 8 # predatation with birth mutate
            # predator mutate
            newpop = Population(t, 1)
            push!(Y_populations_current, newpop)
            push!(Y_populations_all, newpop)
            k = Y_mutfunc(k, index[1])
            push!(k_history, k[end])

            # prey death
            death!(X_populations_current[index[2]], t)
            if X_populations_current[index[2]].n == 0 # check extinction
                deleteat!(X_populations_current, index[2])
                deleteat!(g, index[2])
                payoff = payoff[1:size(payoff, 1) .!= index[2], 1:size(payoff, 2) .!= index[2]]
            end
        end
    end
    return X_populations_all, Y_populations_all, g_history, k_history
end
