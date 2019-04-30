# sample

using Random

# Copy from pkg StatsBase
# Link: https://github.com/JuliaStats/StatsBase.jl/blob/master/src/sampling.jl#L413
# [MIT License](https://github.com/JuliaStats/StatsBase.jl/blob/master/LICENSE.md)
function sample(rng::AbstractRNG, wv::AbstractVector)
    t = rand(rng) * sum(wv)
    n = length(wv)
    i = 1
    cw = wv[1]
    while cw < t && i < n
        i += 1
        @inbounds cw += wv[i]
    end
    return i
end

sample(wv::AbstractVector) = sample(Random.GLOBAL_RNG, wv)
sample(rng::AbstractRNG, a::AbstractVector, wv::AbstractVector) = a[sample(rng, wv)]
sample(a::AbstractArray, wv::AbstractVector) = sample(Random.GLOBAL_RNG, a, wv)

@inline function randchoice(w::AbstractArray{R})where R <: Real
    indices = vec(CartesianIndices(w))
    wv = vec(w)
    sample(indices, wv)
end

# sample end

# macro

macro copyfields(ex)
    fields, var = ex.args
    callcopy = [:(copy($var.$field)) for field in fields.args]
    ex.args[2] = Expr(:tuple, callcopy...)
    return esc(ex)
end

# macro end
