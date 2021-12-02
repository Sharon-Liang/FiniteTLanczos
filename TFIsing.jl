using FiniteTLanczos
using Printf

const z = pauli(:z)
const x = pauli(:x)
const i2 = eye(2)

g = 1.0
len = 3
sample = 30

h = TFIsing(1., g, L = len)
ti = @sprintf "Γ/J = %.1f, system size = %i" g len
s1 = @sprintf "fidelity of %i random samples of %i dimensional v0 = " sample 2^h.L
println(s1, fidelity(2^h.L, sample))


A = FED(h)
B = FTLM(h, R=sample)
C = OFTLM(h, R=sample)

mf1 = a -> FTLM(a, R = sample)
cf1 = a -> correlation_2time(20, 20, z, z, mf1(a))

c_ave = c_average(100, h, mf = mf1, cf = cf1)

beta = [i for i in range(0.1, 20, length = 50)]
tau = [i for i in range(0.1/200,20-0.1/200, length = 50)]
#f0 = [critical_zz_cor(t, 0.1) for t in tau]
f1 = [correlation_2time(t, 20, z, z, A) for t in tau]
f2 = [correlation_2time(t, 20, z, z, B) for t in tau]
#f3 = [correlation_2time(t, 0.1, z, z, C) for t in tau]



using PyPlot
figure()
#plot(tau, f0, "k", label="exact")
plot(tau, f1, "k--", label="FED")
plot(tau, f2, "--", label="FTLM")
#plot(tau, f3, "--", label="OFTLM")

xlabel("tau")
ylabel("correlation2time")
#xscale("log")
title(ti)
legend()
PyPlot.display_figs()
