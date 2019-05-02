export CoevoCompModel

"""
    CoevoCompModel <: AbstractModel

A model sumulates the coevolation during competition.
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

# promote args
function CoevoCompModel(p::Vector{<:Integer},
               b::Vector{RA},
               d::Vector{RB},
               c::Matrix{RC},
               mutrate::RD,
               M::RE,
               f::Function)where {RA <: Real,RB <: Real,RC <: Real,RD <: Real,RE <: Real}
    T = promote_type(RA, RB, RC, RD, RE)
    b = T.(b)
    d = T.(d)
    c = T.(c)
    mutrate = T(mutrate)
    M = T(M)
    return CoevoCompModel(p, b, d, c, mutrate, M, f)
end

function gillespie(m::CoevoCompModel{I,R,<:Function}, T::Real)where {I <: Integer,R <: Real}
    t = zero(R) # initialize t

    # initialize populations
    populations_current = Population.(t, m.populations)
    populations_all = copy(populations_current)

    # destuct model parameters
    @copyfields birthrates, deathrates, payoff, mutrate = m
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
        birth_nomut = @. birthrates * populations_num * nomutrate
        birth_mut = @. birthrates * populations_num * mutrate
        death = @. deathrates * populations_num
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
            birthrates, deathrates, payoff = mutfunc(birthrates, deathrates, payoff, index)
        elseif i in (3, 4) # death (intrinsic and competition)
            death!(populations_current[index], t)
            if populations_current[index].n == 0
                deleteat!(populations_current, index)
                deleteat!(birthrates, index)
                deleteat!(deathrates, index)
                payoff = payoff[1:size(payoff, 1) .!= index, 1:size(payoff, 2) .!= index]
            end
        end
    end
    return populations_all
end
