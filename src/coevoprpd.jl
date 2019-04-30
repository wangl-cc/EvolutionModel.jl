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
    "birthrates of prey populations"
    bX::Vector{R}
    "deathrate of prey population"
    dX::Vector{R}
    "deathrate of prey populations"
    dY::Vector{R}
    "Growth-defense trade off parameter"
    g::Vector{R}
    "ratio of predator growth to predation"
    k::Vector{R}
    "the base k value"
    K::R
    "resource competition coefficient matrix of prey populations"
    payoff::Matrix{R}
    "coefficient scaling the intensity of competition"
    M::R
    "coefficient of the predation rate"
    p::R
    "the growth-defense trade-off coefficient"
    m::R
    "mutation rate of prey populations"
    X_mutrate::R
    "mutation rate of predator populations"
    Y_mutrate::R
    "a function descirbing how mutations change the prey population"
    X_mutfunc::FA
    "a function descirbing how mutations change the pred population"
    Y_mutfunc::FB
end

function gillespie(model::CoevoPrPdModel{I, R, <:Function, <:Function}, T::Real)where {I<:Integer, R<:Real}
    t = zero(R)
    X_populations_current = Population.(t, model.X)
    Y_populations_current = Population.(t, model.Y)
    X_populations_all = copy(X_populations_current)
    Y_populations_all = copy(Y_populations_current)
    @copyfields bX, dX, dY, g, k, K, payoff, M, p, m, X_mutrate, Y_mutrate = model
    g_history = copy(g)
    k_history = copy(k)
    X_nomutrate = 1 - X_mutrate
    Y_nomutrate = 1 - Y_mutrate
    repM = 1/M
    X_mutfunc = model.X_mutfunc
    Y_mutfunc = model.Y_mutfunc
    while t <= T
        X_num = [getfield(p, :n) for p in X_populations_current]
        Y_num = [getfield(p, :n) for p in Y_populations_current]
        sum(X_num)<=0 && sum(Y_num)<=0 && break # both prey and predator extinct
        sum(Y_num)<=0 && break # predator extinct

        trans_X_num = transpose(X_num)
        trans_g = transpose(g)

        X_birth_nomut = @. bX * X_num * X_nomutrate
        X_birth_mut = @. bX * X_num * X_mutrate
        X_intrinsicdeath = @. dX * X_num
        X_competitiondeath = @. trans_X_num / payoff * X_num * repM

        Y_intrinsic_death = @. dY * Y_num

        predatarion_all = @. Y_num * p * trans_g ^ (m * k / K) * trans_X_num

        predatation_nobirth = @. (1 -k) * predatarion_all
        predatation_birth_nomut = @. k * Y_nomutrate * predatarion_all
        predatation_birth_mut = @. k * Y_mutrate * predatarion_all

        τ, i, index = findreaction(X_birth_nomut, X_birth_mut, X_intrinsicdeath, X_competitiondeath,
            Y_intrinsic_death, predatation_nobirth, predatation_birth_nomut, predatation_birth_mut)
        t += τ
        if i == 1
            index = index[1]
            X_populations_current[index].n += 1
            push!(X_populations_current[index].history, (t, X_populations_current[index].n))
        elseif i == 2
            index = index[1]
            newpop = Population(t, 1)
            push!(X_populations_current, newpop)
            push!(X_populations_all, newpop)
            bX, dX, g, payoff = X_mutfunc(bX, dX, g, payoff, index)
            push!(g_history, g[end])
        elseif i in (3, 4)
            index = index[1]
            X_populations_current[index].n -= 1
            push!(X_populations_current[index].history, (t, X_populations_current[index].n))
            if X_populations_current[index].n == 0
                deleteat!(X_populations_current, index)
                deleteat!(bX, index)
                deleteat!(dX, index)
                deleteat!(g, index)
                payoff = payoff[1:size(payoff, 1) .!= index, 1:size(payoff, 2) .!= index]
            end
        elseif i == 5
            index = index[1]
            Y_populations_current[index].n -= 1
            if Y_populations_current[index].n == 0
                deleteat!(Y_populations_current, index)
                deleteat!(dY, index)
                deleteat!(k, index)
            end
        elseif i == 6
            index = index[2]
            X_populations_current[index].n -= 1
            push!(X_populations_current[index].history, (t, X_populations_current[index].n))
            if X_populations_current[index].n == 0
                deleteat!(X_populations_current, index)
                deleteat!(bX, index)
                deleteat!(dX, index)
                deleteat!(g, index)
                payoff = payoff[1:size(payoff, 1) .!= index, 1:size(payoff, 2) .!= index]
            end
        elseif i == 7
            Y_populations_current[index[1]].n += 1
            push!(Y_populations_current[index[1]].history, (t, Y_populations_current[index[1]].n))

            X_populations_current[index[2]].n -= 1
            push!(X_populations_current[index[2]].history, (t, X_populations_current[index[2]].n))
            if X_populations_current[index[2]].n == 0
                deleteat!(X_populations_current, index[2])
                deleteat!(bX, index[2])
                deleteat!(dX, index[2])
                deleteat!(g, index[2])
                payoff = payoff[1:size(payoff, 1) .!= index[2], 1:size(payoff, 2) .!= index[2]]
            end
        elseif  i == 8
            X_populations_current[index[2]].n -= 1
            push!(X_populations_current[index[2]].history, (t, X_populations_current[index[2]].n))
            if X_populations_current[index[2]].n == 0
                deleteat!(X_populations_current, index[2])
                deleteat!(bX, index[2])
                deleteat!(dX, index[2])
                deleteat!(g, index[2])
                payoff = payoff[1:size(payoff, 1) .!= index[2], 1:size(payoff, 2) .!= index[2]]
            end

            newpop = Population(t, 1)
            push!(Y_populations_current, newpop)
            push!(Y_populations_all, newpop)
            dY, k = Y_mutfunc(dY, k, index[1])
            push!(k_history, k[end])
        end
    end
    return X_populations_all, Y_populations_all, g_history, k_history
end
