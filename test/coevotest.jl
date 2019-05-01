using Distributions

# CoevoCompModel Test

function mutcompmat(compmat::AbstractMatrix{R}, i::Integer, σ::Real)where R <: Real
    l = size(compmat, 1)
    newmat = zeros(R, l + 1, l + 1)
    newmat[1:l, 1:l] = compmat
    @simd for k in 1:l
        @inbounds newmat[l + 1, k] = rand(TruncatedNormal(compmat[i, k], σ, 0, Inf))
        @inbounds newmat[k, l + 1] = rand(TruncatedNormal(compmat[k, i], σ, 0, Inf))
    end
    newmat[l + 1, l + 1] = rand(TruncatedNormal(compmat[i, i], σ, 0, Inf))
    return newmat
end

@inline function mutfunc(b, d, c, i)
    b = push!(b, b[i])
    d = push!(d, d[i])
    c = mutcompmat(c, i, 0.01)
    return b, d, c
end

model = CoevoCompModel([100], [0.6], [0.1], ones(1,1), 0.05, 100, mutfunc)

ps = gillespie(model, 10)

# CoevoCompModel Test

const X_σ = 0.1
const Y_σ = 0.03

function X_mutfunc(bX, dX, g, payoff, i)
    push!(bX, bX[i])
    push!(dX, dX[i])
    push!(g, rand(TruncatedNormal(g[i], X_σ, 0.0, 1.0)))
    l = size(payoff, 1)
    payoff = ones(l+1, l+1)
    return bX, dX, g, payoff
end

function Y_mutfunc(dY, k, i)
    push!(dY, dY[i])
    k = push!(k, rand(TruncatedNormal(k[i], Y_σ, 0.0, 0.3)))
    return dY, k
end


X = [1000]
Y = [100]
bX = [1.0]
dX = [0.1]
dY = [0.5]
g = [1.0]
k = [0.3]
K = 0.3
payoff = ones(1, 1)
M = 1/0.00005
p = 0.005
m = 1.
X_mutrate = 0.001
Y_mutrate = 0.01

model = CoevoPrPdModel(X, Y, bX, dX, dY, g, k ,K, payoff, M, p, m, X_mutrate, Y_mutrate, X_mutfunc, Y_mutfunc)

X, Y, g, k = gillespie(model, 500)
