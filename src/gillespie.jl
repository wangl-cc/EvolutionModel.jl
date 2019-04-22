export gillespie

abstract type AbstractModel end

function findreaction(reactions::AbstractArray{R}...)where {R}
    min_τ = Inf
    indexof_min_τ = (0, 0)
    for i in eachindex(reactions)
        reaction = reactions[i]
        τs = - log.(rand(R, size(reaction))) ./ reaction
        τ, index = findmin(τs)
        if τ <= min_τ
            min_τ = τ
            indexof_min_τ = (i, index)
        end
    end
    return min_τ, indexof_min_τ
end

function gillespie(m::AbstractModel, T::Number) end
