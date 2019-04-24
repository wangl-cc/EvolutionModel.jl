export CoevoCompModel

struct CoevoCompModel{I <: Integer,R <: Real,F <: Function} <: AbstractModel
    populations::Vector{I}
    birthrates::Vector{R}
    deathrates::Vector{R}
    compmat::Matrix{R}
    mutrate::R
    M::R
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

function gillespie(m::CoevoCompModel{I,R,F}, T::Number)where {I <: Integer,R <: Real,F <: Function}
    t = R(0.0)
    populations = [Population(t, i) for i in m.populations]
    birthrates = deepcopy(m.birthrates)
    deathrates = deepcopy(m.deathrates)
    compmat = deepcopy(m.compmat)
    populations_history = [populations...]
    mutrate = m.mutrate
    mutfunc = m.mutfunc
    repM = 1/m.M
    while t <= T
        populations_num = [i.n for i in populations]
        sum(populations_num)<=0 && break
        birth_nomut = birthrates .* populations_num * (1 - mutrate)
        birth_mut = birthrates .* populations_num * mutrate
        death = deathrates .* populations_num
        comp  = reshape(populations_num, (1, length(populations_num))) ./ compmat .* populations_num .* repM
        τ, i, index = findreaction(birth_nomut, birth_mut, death, comp)
        i = i[1]
        index = index[1]
        t += τ
        if i == 1
            populations[index].n += 1
            push!(populations[index].history, (t, populations[index].n))
        elseif i == 2
            newpop = Population(t, 1)
            populations = vcat(populations, newpop)
            populations_history = vcat(populations_history, newpop)
            birthrates, deathrates, compmat = mutfunc(birthrates, deathrates, compmat, index)
        elseif i in (3, 4)
            populations[index].n -= 1
            push!(populations[index].history, (t, populations[index].n))
            if populations[index] == 0
                deleteat!(populations, index)
                deleteat!(birthrates, index)
                deleteat!(deathrates, index)
                compmat = compmat[1:size(compmat, 1) .!= index, 1:size(compmat, 2) .!= index]
            end
        end
    end
    return populations_history
end
