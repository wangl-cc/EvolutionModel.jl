export gillespie

abstract type AbstractModel end

mutable struct Population{R<:Real, I<:Integer}
    n::I
    history::Vector{Tuple{R, I}}
    function Population(t::R, n::I)where {R<:Real, I<:Integer}
        return new{R, I}(n, [(t, n)])
    end
end

# function findreaction(reactions::AbstractArray{R}...)where {R}
#     min_τ = Inf
#     indexof_min_τ = (0, 0)
#     for i in eachindex(reactions)
#         reaction = reactions[i]
#         τs = - log.(rand(R, size(reaction))) ./ reaction
#         τ, index = findmin(τs)
#         if τ <= min_τ
#             min_τ = τ
#             indexof_min_τ = (i, index)
#         end
#     end
#     return min_τ, indexof_min_τ
# end

function gillespie(m::AbstractModel, T::Number) end

function findreaction(reactionrates::AbstractArray{R}...)where R<:Real
    sum_rs = [sum(i) for i in reactionrates]
    sum_all = sum(sum_rs)
    τ = -log(rand())/sum_all
    i = randchoice(sum_rs)[1]
    index = randchoice(reactionrates[i])
    return τ, i, index
end

using Random

# Copy from pkg StatsBase
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

function randchoice(w::AbstractArray{R})where R<:Real
    indices = vec(CartesianIndices(w))
    wv = vec(w)
    sample(indices, wv)
end
