export gillespie

abstract type AbstractModel end

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
    i = randchoice(sum_rs)
    index = randchoice(reactionrates[i[1]])
    return τ, i, index
end

using StatsBase

function randchoice(l::AbstractArray{R})where R<:Real
    indices = vec(CartesianIndices(l))
    w = Weights(vec(l))
    sample(indices, w)
end
