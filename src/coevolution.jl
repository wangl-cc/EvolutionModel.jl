using Distributions

struct CoevoCompModel{I <: Integer,R <: Real,F <: Function} <: AbstractModel
    populations::Vector{I}
    birthrates::Vector{R}
    deathrates::Vector{R}
    compmat::Matrix{R}
    mutrate::R
    mutfunc::F
end

function gillespie(m::CoevoCompModel{I,R,F}, T::Number)where {I <: Integer,R <: Real,F <: Function}
    t = R(0.0)
    populations = deepcopy(m.populations)
    birthrates = deepcopy(m.birthrates)
    deathrates = deepcopy(m.deathrates)
    compmat = deepcopy(m.compmat)
    mutrate = m.mutrate
    mutfunc = m.mutfunc
    data_p = Vector{I}[deepcopy(populations)]
    data_t = R[t]
    while t <= T
        birth_nomut = birthrates .* populations * (1 - mutrate)
        birth_mut = birthrates .* populations * mutrate
        death = deathrates .* populations
        comp  = reshape(populations, (1, length(populations))) .* compmat .* populations
        τ, (index, i) = findreaction(birth_nomut, birth_mut, death, comp)
        # τ == Inf && println(birth_nomut, birth_mut, death, comp); break
        t += τ
        if index == 1
            populations[i] += 1
        elseif index == 2
            populations = vcat(populations, I(1))
            birthrates, deathrates, compmat = mutfunc(birthrates, deathrates, compmat, i)
        elseif index == 3
            populations[i] -= 1
            #  if populations[i] == 0
            #      deleteat!(populations, i)
            #      deleteat!(birthrates, i)
            #      deleteat!(deathrates, i)
            #      compmat = compmat[1:size(compmat, 1) .!= i, 1:size(compmat, 2) .!= i]
            #  end
        elseif index == 4
            populations[i[1]] -= 1
            #  if populations[i] == 0
            #      deleteat!(populations, i)
            #      deleteat!(birthrates, i)
            #      deleteat!(deathrates, i)
            #      compmat = compmat[1:size(compmat, 1) .!= i, 1:size(compmat, 2) .!= i]
            #  end
        end
        push!(data_p, deepcopy(populations))
        push!(data_t, t)
    end
    return data_p, data_t
end

function mutcompmat(compmat::AbstractMatrix{R}, i::Integer)where R <: Real
    l = size(compmat, 1)
    newmat = zeros(R, l + 1, l + 1)
    newmat[1:l, 1:l] = compmat
    for k in 1:l
        newmat[l + 1, k] = rand(TruncatedNormal(compmat[i, k], 1, 0, Inf))
        newmat[k, l + 1] = rand(TruncatedNormal(compmat[k, i], 1, 0, Inf))
    end
    newmat[l + 1, l + 1] = rand(TruncatedNormal(compmat[i, i], 1, 0, Inf))
    return newmat
end
