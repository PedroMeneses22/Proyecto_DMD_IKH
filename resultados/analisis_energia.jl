using Plots
using DelimitedFiles

casos = [
    ("resultados/delta0.2_eps0.01_nu0.01/diagnosticos.dat", "ν=1e-2", :red),
    ("resultados/delta0.2_eps0.01_nu0.001/diagnosticos.dat", "ν=1e-3", :orange),
    ("resultados/delta0.2_eps0.01_nu0.0001/diagnosticos_nu0.0001.dat", "ν=1e-4", :blue),
    ("resultados/delta0.2_eps0.01_nu1.0e-5/diagnosticos.dat", "ν=1e-5", :green),
]

p_energia   = plot(title="Energía cinética", xlabel="t",
                   ylabel="E(t)", yscale=:log10)
p_enstrofia = plot(title="Enstrofía",        xlabel="t",
                   ylabel="Z(t)")

for (path, label, color) in casos
    data = readdlm(path, comments=true)
    t = data[:, 1]
    E = data[:, 2]
    Z = data[:, 3]
    plot!(p_energia,   t, E, label=label, color=color, lw=2)
    plot!(p_enstrofia, t, Z, label=label, color=color, lw=2)
end

savefig(p_energia,   "comparacion_energia.png")
savefig(p_enstrofia, "comparacion_enstrofia.png")