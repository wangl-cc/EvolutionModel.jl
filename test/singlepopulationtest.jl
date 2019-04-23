# using HypothesisTests

# Parameters

n = 10
b = 0.5
d = 0.1
T = 10

repeatimes = 1000

# ExpModel Tests

m = ExpModel(n, b, Float32(d))

v = [gillespie(m, T)[end][2] for _ in 1:repeatimes]

# minv, maxv = confint(OneSampleTTest(v))
#
# expect = n*exp((b-d)*T)
#
# @test minv < expect < maxv

# LogistModel Tests

c = 0.005

K = (b - d)/c

m = LogistModel(n, Float32(b), d, c)

v = [gillespie(m, T)[end][2] for _ in 1:repeatimes]

# minv, maxv = confint(OneSampleTTest(v))

# expect = K*n/(n+(K-n)*exp(-(b-d)T))

# @test minv < expect < maxv

m = LogistModel(n, Float32(b), d, K=K)

v = [gillespie(m, T)[end][2] for _ in 1:repeatimes]

# minv, maxv = confint(OneSampleTTest(v))

# expect = K*n/(n+(K-n)*exp(-(b-d)T))

# @test minv < expect < maxv
