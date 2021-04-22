using FiniteTLanczos
using PyPlot
using Printf

z = pauli('z')

g = 1.0
ti = @sprintf "Γ/J = %.1f" g
h = TFIsing(1.,g,L = 10)
A = FED(h)
B = FTLM(h)
C = OFTLM(h);

beta = [i for i in range(1, 20, length = 50)]

f0 = [free_energy(1., g, β) for β in beta]
f1 = [free_energy(β, A) for β in beta ]
f2 = [free_energy(β, B) for β in beta ];
f3 = [free_energy(β, C) for β in beta ];




figure()
plot(beta, f0, "k", label="exact")
plot(beta, f1, "--", label="FED")
plot(beta, f2, "--", label="FTLM")
plot(beta, f3, "--", label="OFTLM")

xlabel("β")
ylabel("free_energy")
#xscale("log")
title(ti)
legend()
PyPlot.display_figs()
