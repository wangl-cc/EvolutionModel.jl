using EvolutionModel
using Plots

plotlyjs()

function mutfunc(b, d, c, i)
    b = vcat(b, b[i])
    d = vcat(d, d[i])
    c = EvolutionModel.mutcompmat(c, i)
    return b, d, c
end

m = CoevoCompModel([100], [0.6], [0.1], ones(1,1), 0.05, 500., mutfunc)


p, t = gillespie(m, 100)

plt = plot(t, [get(i, 1, 0) for i in p], legend=false)
for j in 2:length(p[end])
    plot!(plt, t, [get(i, j, 0) for i in p], )
end

i = rand(Int)
png(plt, joinpath(@__DIR__, "plt-$i.png"))
