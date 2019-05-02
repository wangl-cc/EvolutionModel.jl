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
    "record population change history"
    history::Vector{Tuple{R,I}}
end

function Population(t::R, n::I)where {R <: Real,I <: Integer}
    return Population{R,I}(n, [(t, n)])
end

function birth!(p::Population, t::Real)
    p.n += 1
    push!(p.history, (t, p.n))
end

function death!(p::Population, t::Real)
    p.n -= 1
    push!(p.history, (t, p.n))
end

"""
    gillespie(m::AbstractModel,T::Real)

A uniform interface to sumulate model `m` by gillespie algorithm.
"""
function gillespie(m::AbstractModel, T::Real)
    println("There is no method to sumulate ", typeof(m), "!")
end

"""
    choosereaction(reactionrates::AbstractArray{R}...)where R<:Real

Calculate the reaction time `τ` and randomly choose an reaction from reactions with weights `reactionrates`.
"""
function choosereaction(reactionrates::AbstractArray{R}...)where R <: Real
    sum_rs = [sum(i) for i in reactionrates]
    sum_all = sum(sum_rs)
    τ = -log(rand(R)) / sum_all
    i = randchoice(sum_rs)[1]
    index = randchoice(reactionrates[i])
    return τ, i, index
end
