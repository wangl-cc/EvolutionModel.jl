export CoevoCompModel

"""
    CoevoCompModel <: AbstractModel

A model sumulates the coevolation during competition.
"""
struct CoevoCompModel{I <: Integer,R <: Real,F <: Function} <: AbstractModel
    "population sizes"
    populations::Vector{I}
    "baseline birthrate of given populations"
    b::R
    "deathrate of given populations"
    d::R
    "coefficients scaling growing ability"
    g::Vector{R}
    "payoff mutrate of given populations"
    payoff::Matrix{R}
    "coefficient scaling the intensity of competition"
    M::R
    "mutation rate of given populations"
    mutrate::R
    "a function descirbing how mutations change the population"
    mutfunc::F
end

# promote args
function CoevoCompModel(p::AbstractVector{<:Integer},
               b::Real,
               d::Real,
               g::AbstractVector{<:Real},
               payoff::AbstractMatrix{<:Real},
               M::Real,
               mutrate::Real,
               f::Function)
    return CoevoCompModel(p, promote_(b, d, g, payoff, M, mutrate)..., f)
end

function gillespie(m::CoevoCompModel{I,R,<:Function}, T::Real)where {I <: Integer,R <: Real}
    t = zero(R) # initialize t

    # initialize populations
    populations_current = Population.(t, m.populations)
    populations_all = copy(populations_current)

    # destuct model parameters
    @copyfields b, d, g, payoff, mutrate = m
    g_history = copy(g)
    nomutrate = 1 - mutrate
    repM = 1 / m.M
    mutfunc = m.mutfunc

    # main loop
    while t <= T
        length(populations_current) <= 0 && break # check extinction

        # get current population sizes
        populations_num = [getfield(p, :n) for p in populations_current]
        trans_p_num = transpose(populations_num)

        # calculate reactions
        birth_nomut = @. b * g * populations_num * nomutrate
        birth_mut = @. b * g * populations_num * mutrate
        death = @. d * populations_num
        comp  = @. trans_p_num / payoff * populations_num * repM

        # calculate τ and choose reaction
        τ, i, index = choosereaction(birth_nomut, birth_mut, death, comp)
        index = index[1]

        t += τ # time increase

        # solve reactions
        if i == 1 # birth no mutate
            birth!(populations_current[index], t)
        elseif i == 2 # birh mutate
            newpop = Population(t, 1)
            push!(populations_current, newpop)
            push!(populations_all, newpop)
            g, payoff = mutfunc(g, payoff, index)
            push!(g_history, g[end])
        elseif i in (3, 4) # death (intrinsic and competition)
            death!(populations_current[index], t)
            if populations_current[index].n == 0 # check extinction
                deleteat!(populations_current, index)
                deleteat!(g, index)
                payoff = payoff[1:size(payoff, 1) .!= index, 1:size(payoff, 2) .!= index]
            end
        end
    end
    return populations_all, g_history
end
