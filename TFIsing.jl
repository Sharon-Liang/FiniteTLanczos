using FiniteTLanczos
using PyPlot
using Printf

z = pauli('z')
x = pauli('x')
g = 1.0
ti = @sprintf "Γ/J = %.1f" g
h = TFIsing(1.,g,L = 3)
A = FED(h)
B = FTLM(h)
C = OFTLM(h)

beta = [i for i in range(0.1, 20, length = 50)]
tau = [i for i in range(0.1/200,(20-0.1)/200, length = 20)]
#f0 = [critical_zz_cor(t, 0.1) for t in tau]
f1 = [correlation2time(t, 0.1, z, z, A) for t in tau]
f2 = [correlation2time(t, 0.1, z, z, B) for t in tau]
f3 = [correlation2time(t, 0.1, z, z, C) for t in tau]




figure()
#plot(tau, f0, "k", label="exact")
plot(tau, f1, "k--", label="FED")
plot(tau, f2, "--", label="FTLM")
plot(tau, f3, "--", label="OFTLM")

xlabel("tau")
ylabel("correlation2time")
#xscale("log")
title(ti)
legend()
PyPlot.display_figs()
