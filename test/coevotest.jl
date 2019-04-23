using Distributions

function mutcompmat(compmat::AbstractMatrix{R}, i::Integer, var::Real)where R <: Real
    l = size(compmat, 1)
    newmat = zeros(R, l + 1, l + 1)
    newmat[1:l, 1:l] = compmat
    @inbounds for k in 1:l
        newmat[l + 1, k] = rand(TruncatedNormal(compmat[i, k], var, 0, Inf))
        newmat[k, l + 1] = rand(TruncatedNormal(compmat[k, i], var, 0, Inf))
    end
    newmat[l + 1, l + 1] = rand(TruncatedNormal(compmat[i, i], var, 0, Inf))
    return newmat
end

function mutfunc(b, d, c, i)
    b = vcat(b, b[i])
    d = vcat(d, d[i])
    c = mutcompmat(c, i, 0.01)
    return b, d, c
end

m = CoevoCompModel([100], [0.6], [0.1], ones(1,1), 0.05, 100, mutfunc)

ps = gillespie(m, 100)
