export ExpModel, LogistModel

struct ExpModel{I <: Integer,R <: Real} <: AbstractModel
    n::I
    b::R
    d::R
end

ExpModel(n::Integer, b::Real, d::Real) = ExpModel(n, promote(b, d)...)

function gillespie(m::ExpModel{I, R}, T::Real)where {I <: Integer,R <: Real}
    n = m.n
    b = m.b
    d = m.d
    t = R(0.0)
    history = [(t, n)]
    while t <= T
        reactions = [n*b, n*d]
        τ, (_, i) = findreaction(reactions)
        t += τ
        if i == 1
            n += 1
        elseif i == 2
            n -= 1
        end
        push!(history, (t, n))
    end
    return history
end

struct LogistModel{I <: Integer, R <: Real} <: AbstractModel
    n::I
    b::R
    d::R
    c::R
end

LogistModel(n::Integer, b::Real, d::Real, c::Real) = LogistModel(n, promote(b, d, c)...)
LogistModel(n::Integer, b::Real, d::Real; K::Real) = LogistModel(n, promote(b, d, (b-d)/K)...)

function gillespie(m::LogistModel{I, R}, T::Integer)where {I <: Integer, R <: Real}
    n = m.n
    b = m.b
    d = m.d
    c = m.c
    t = R(0.0)
    history = [(t, n)]
    while t <= T
        reactions = [n*b, n*d, n*c*n]
        τ, (_, i) = findreaction(reactions)
        t += τ
        if i == 1
            n += 1
        elseif i == 2
            n -= 1
        elseif i == 3
            n -= 1
        end
        push!(history, (t, n))
    end
    return history
end
