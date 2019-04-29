export gillespie, Population, AbstractModel

"""
    AbstractModel

Supertype of all model.
"""
abstract type AbstractModel end

"""
    Population{R<:Real, I<:Integer}
"""
mutable struct Population{R <: Real,I <: Integer}
    "current population size"
    n::I
    "population change history"
    history::Vector{Tuple{R,I}}
    function Population(t::R, n::I)where {R <: Real,I <: Integer}
        return new{R,I}(n, [(t, n)])
    end
end

"""
    gillespie(m::AbstractModel,T::Real)

A uniform interface to sumulate model `m` by gillespie algorithm.
"""
function gillespie(m::AbstractModel, T::Real)
    modeltype = typeof(m)
    println("There is no method to sumulate ", modeltype, "!")
end

"""
    findreaction(reactionrates::AbstractArray{R}...)where R<:Real

Calculate the reaction time `τ` and randomly choose an reaction from reactions with weights `reactionrates`.
"""
function findreaction(reactionrates::AbstractArray{R}...)where R <: Real
    sum_rs = [sum(i) for i in reactionrates]
    sum_all = sum(sum_rs)
    τ = -log(rand(R)) / sum_all
    i = randchoice(sum_rs)
    index = randchoice(reactionrates[i])
    return τ, i, index
end

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
    sample(indices, wv)[1]
end
