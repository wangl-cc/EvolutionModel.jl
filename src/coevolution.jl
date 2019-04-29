export CoevoCompModel

"""
    CoevoCompModel <: AbstractModel

A model sumulates the coevolation during competition.

# Example
```jldoctest
```
"""
struct CoevoCompModel{I <: Integer,R <: Real,F <: Function} <: AbstractModel
    "population sizes"
    populations::Vector{I}
    "birthrates of given populations"
    birthrates::Vector{R}
    "deathrates of given populations"
    deathrates::Vector{R}
    "payoff mutrate of given populations"
    payoff::Matrix{R}
    "mutation rate of given populations"
    mutrate::R
    "coefficient scaling the intensity of competition"
    M::R
    "a function descirbing how mutations change the population"
    mutfunc::F
end

function CoevoCompModel(p::Vector{<:Integer},
               b::Vector{RA},
               d::Vector{RB},
               c::Matrix{RC},
               mutrate::RD,
               M::RE,
               f::Function)where {RA <: Real, RB <: Real, RC <: Real, RD <: Real, RE <: Real}
    T = promote_type(RA, RB, RC, RD, RE)
    b = T.(b)
    d = T.(d)
    c = T.(c)
    mutrate = T(mutrate)
    M = T(M)
    return CoevoCompModel(p, b, d, c, mutrate, M, f)
end

function gillespie(m::CoevoCompModel{I,R,F}, T::Real)where {I <: Integer,R <: Real,F <: Function}
    t = zero(R)
    populations = Population.(t, m.populations)
    birthrates = copy(m.birthrates)
    deathrates = copy(m.deathrates)
    payoff = copy(m.payoff)
    populations_history = copy(populations)
    mutrate = m.mutrate
    nomutrate = 1 - mutrate
    mutfunc = m.mutfunc
    repM = 1/m.M
    while t <= T
        populations_num = [getfield(p, :n) for p in populations]
        sum(populations_num)<=0 && break
        birth_nomut = birthrates .* populations_num .* nomutrate
        birth_mut = birthrates .* populations_num .* mutrate
        death = deathrates .* populations_num
        comp  = reshape(populations_num, (1, length(populations_num))) ./ payoff .* populations_num .* repM
        τ, i, index = findreaction(birth_nomut, birth_mut, death, comp)
        t += τ
        if i == 1
            populations[index].n += 1
            push!(populations[index].history, (t, populations[index].n))
        elseif i == 2
            newpop = Population(t, 1)
            push!(populations, newpop)
            push!(populations_history, newpop)
            birthrates, deathrates, payoff = mutfunc(birthrates, deathrates, payoff, index)
        elseif i in (3, 4)
            populations[index].n -= 1
            push!(populations[index].history, (t, populations[index].n))
            if populations[index].n == 0
                deleteat!(populations, index)
                deleteat!(birthrates, index)
                deleteat!(deathrates, index)
                payoff = payoff[1:size(payoff, 1) .!= index, 1:size(payoff, 2) .!= index]
            end
        end
    end
    return populations_history
end
