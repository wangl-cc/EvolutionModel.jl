using EvolutionModel
using Plots

plotlyjs()

function mutfunc(b, d, c, i)
    b = vcat(b, b[i])
    d = vcat(d, d[i])
    c = mutcompmat(c, i)
    return b, d, c
end

m = CoevoCompModel([100], [0.6], [0.1], ones(1,1), 0.05, 100., mutfunc)


ps = gillespie(m, 100)

plt = plot(legend=false)
for p in ps
    plot!(p.history)
end

i = rand(Int)
png(plt, joinpath(@__DIR__, "plt-$i.png"))
