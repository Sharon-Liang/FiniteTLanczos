using myED
using PyPlot
using Printf

z = pauli('z')

g = 0.1
ti = @sprintf "Γ/J = %.1f" g
h = TFIsing(1.,g,L = 10)
A = FED(h)

beta = [i for i in range(1, 20, length = 50)]
T = 1 ./ beta
#f1 = [1/ partitian(β, A) for β in beta ]
f2 = [structure_factor(0, β, z,z,A) for β in beta]
f3 = [structure_factor(0, β, z,z,A, η = 1.e-5) for β in beta]


figure()
#plot(T, f1, "k", label="theoretical")
plot(T, f2, "o--", markerfacecolor="none", label="FED, η = 0.05")
plot(T, f3, "o--", markerfacecolor="none", label="FED, η = 1.e-5")
xlabel("T")
ylabel(" S(ω - 0) ")
#xscale("log")
title(ti)
legend()
PyPlot.display_figs()
