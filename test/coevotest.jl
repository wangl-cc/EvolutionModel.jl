using Distributions

# CoevoCompModel Test

function mutpayoff(payoff::AbstractMatrix{R}, i::Integer, σ::Real)where R <: Real
    l = size(payoff, 1)
    newpayoff = zeros(R, l + 1, l + 1)
    newpayoff[1:l, 1:l] = payoff
    @simd for k in 1:l
        @inbounds newpayoff[l + 1, k] = rand(TruncatedNormal(payoff[i, k], σ, 0, Inf))
        @inbounds newpayoff[k, l + 1] = rand(TruncatedNormal(payoff[k, i], σ, 0, Inf))
    end
    newpayoff[l + 1, l + 1] = rand(TruncatedNormal(payoff[i, i], σ, 0, Inf))
    return newpayoff
end

@inline function mutfunc(g::Vector, c::Matrix, i::Integer)
    push!(g, g[i])
    c = mutpayoff(c, i, 0.01)
    return g, c
end

model = CoevoCompModel([100], 1, 0.1, [1], ones(1,1), 100, 0.05, mutfunc)

ps, g = gillespie(model, 10)

# CoevoCompModel Test

const X_σ = 0.1
const Y_σ = 0.03

function X_mutfunc(g, payoff, i)
    push!(g, rand(TruncatedNormal(g[i], X_σ, 0.0, 1.0)))
    l = size(payoff, 1)
    payoff = ones(l+1, l+1)
    return g, payoff
end

function Y_mutfunc(k, i)
    k = push!(k, rand(TruncatedNormal(k[i], Y_σ, 0.0, 0.3)))
    return k
end

begin
    X = [1000]
    Y = [100]
    bX = 1
    dX = 0.1
    dY = 0.5
    g = [1.0]
    k = [0.3]
    K = 0.3
    m = 1.
    payoff = ones(1, 1)
    M = 1/0.00005
    p = 0.005
    X_mutrate = 0.001
    Y_mutrate = 0.01
end

model = CoevoPrPdModel(X, Y, bX, dX, dY, g, k ,K, m, payoff, M, p, X_mutrate, Y_mutrate, X_mutfunc, Y_mutfunc)

X, Y, g, k = gillespie(model, 500)

convertmodel = CoevoCompModel(model)

X, g = gillespie(convertmodel, 100)
