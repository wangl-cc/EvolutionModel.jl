using HypothesisTests

# Parameters

n = 10
b = 0.5
d = 0.1
T = 10

repeatimes = 100

# ExpModel Tests

m = ExpModel(n, b, Float32(d))

test_results = Bool[]

for i in 1:3
    v = [gillespie(m, T).n for _ in 1:repeatimes]
    minv, maxv = confint(OneSampleTTest(v))
    expect = n*exp((b-d)*T)
    push!(test_results, minv < expect < maxv)
end

@test any(test_results)

# LogistModel Tests

c = 0.005

K = (b - d)/c

m = LogistModel(n, Float32(b), d, c)

test_results = Bool[]

for i in 1:10
    v = [gillespie(m, T).n for _ in 1:repeatimes]
    minv, maxv = confint(OneSampleTTest(v))
    expect = K*n/(n+(K-n)*exp(-(b-d)T))
    push!(test_results, minv < expect < maxv)
end

@test any(test_results)

m = LogistModel(n, Float32(b), d, K=K)

test_results = Bool[]

for i in 1:10
    v = [gillespie(m, T).n for _ in 1:repeatimes]
    minv, maxv = confint(OneSampleTTest(v))
    expect = K*n/(n+(K-n)*exp(-(b-d)T))
    push!(test_results, minv < expect < maxv)
end

@test any(test_results)
