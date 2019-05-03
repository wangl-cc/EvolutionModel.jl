# sample

using Random

# Copy from pkg StatsBase
# Link: https://github.com/JuliaStats/StatsBase.jl/blob/master/src/sampling.jl#L413
# [MIT License](https://github.com/JuliaStats/StatsBase.jl/blob/master/LICENSE.md)
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

sample(wv) = sample(Random.GLOBAL_RNG, wv)

"""
    randchoice(w::AbstractArray{R})where R<:Real

Randomly choose a element of given weight array `w` return CartesianIndex.

# Example
```jldoctest
julia> v = [randchoice([1, 2, 3, 4]) for _ in 1:1000]
1000-element Array{CartesianIndex{1},1}:
 CartesianIndex(4,)
 ⋮
 CartesianIndex(3,)

julia> count(x->x==CartesianIndex(4,), v)/1000
0.412
```
"""
@inline function randchoice(w::AbstractArray{R})where R <: Real
    indices = vec(CartesianIndices(w))
    wv = vec(w)
    return indices[sample(wv)]
end

# sample end

function promote_(as...)
    T = Bool
    for a in as
        ta = typeof(a)
        if ta <: Number
            T = promote_type(T, ta)
        elseif ta <: AbstractArray
            T = promote_type(T, ta.parameters[1])
        else
            error("Unsupport type!")
        end
   end
   return [T.(a) for a in as]
end

# macro

"""
    copyfields(ex)

Copy fields of given struct.

# Example
```jldoctest
julia> struct M
           a
           b
           c
           d
       end

julia> m = M(1, 2, 3, 4)
M(1, 2, 3, 4)

julia> @copyfields a, b, c, d = m
(1, 2, 3, 4)

julia> a
1

julia> b
2

julia> c
3

julia> d
4
```
"""
macro copyfields(ex)
    fields, var = ex.args
    callcopy = [:(copy($var.$field)) for field in fields.args]
    ex.args[2] = Expr(:tuple, callcopy...)
    return esc(ex)
end

# macro end
