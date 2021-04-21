using myED
using PyPlot

gamma = [i for i in range(0,2,length = 50)]
sx = zeros(50)
x = zeros(50)
for i = 1:50
    g = gamma[i]
    h = TFIsing(1.,g,L = 10)
    A = FED(h)
    x[i] = ave_sx(1., g, 20)
    sx[i] = thermal_average(20,pauli('x'),A)
end

figure()
plot(gamma, x, label="exact")
plot(gamma, sx, label="FED")
legend()
PyPlot.display_figs()

figure()
plot(beta, f2 .- f1)
PyPlot.display_figs()
