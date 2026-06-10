# =============================================================================
# Simulación de Inestabilidad de Kelvin-Helmholtz
# Método pseudo-espectral 2D con vorticidad
# Por Pedro Damian Meneses Orozco
# Facultad de Ciencias 
# =============================================================================

using FFTW
using LinearAlgebra
using Printf
using Plots

# =============================================================================
# PARÁMETROS
# =============================================================================

const Nx  = 256          # Puntos en x
const Ny  = 256          # Puntos en y
const Lx  = 2π          # Longitud del dominio en x (periódico)
const Ly  = 2π          # Longitud del dominio en y (periódico)
const ν   = 1e-5        # Viscosidad cinemática (Re ~ 10^4)
const U0  = 1.0          # Velocidad de cizallamiento
const δ   = 0.2          # Espesor de la capa de mezcla
const ε   = 0.01         # Amplitud de la perturbación inicial
const dt  = 5e-4         # Paso temporal
const T   = 20.0         # Tiempo total de simulación
const Nt  = Int(T / dt)  # Número de pasos
const n_save = 200       # Guardar cada n_save pasos
const ouput_dir = "resultados/delta$(δ)_eps$(ε)_nu$(ν)"

# Crear directorio
mkpath(ouput_dir)

# =============================================================================
# GRILLA Y NÚMEROS DE ONDA
# =============================================================================

x = [Lx * i / Nx for i in 0:Nx-1]
y = [Ly * j / Ny - Ly/2 for j in 0:Ny-1]  # centrado en 0

# Números de onda (convención FFTW)
kx_vec = fftfreq(Nx, Nx / Lx) .* 2π
ky_vec = fftfreq(Ny, Ny / Ly) .* 2π

# Matrices 2D de números de onda: KX[j,i], KY[j,i]
KX = [kx_vec[i] for j in 1:Ny, i in 1:Nx]
KY = [ky_vec[j] for j in 1:Ny, i in 1:Nx]
K2 = KX.^2 .+ KY.^2

# Evitar división por cero en el modo k=0
K2_safe = copy(K2)
K2_safe[1, 1] = 1.0

# =============================================================================
# MÁSCARA DE DEALISING (regla de los 2/3)
# Elimina los modos de alto número de onda que generan aliasing
# en los productos no lineales.
# =============================================================================

kx_max = maximum(abs.(kx_vec))
ky_max = maximum(abs.(ky_vec))
mask = (abs.(KX) .<= 2/3 * kx_max) .& (abs.(KY) .<= 2/3 * ky_max)

# =============================================================================
# CONDICIÓN INICIAL
# Perfil de tangente hiperbólica + perturbación del modo más inestable
# =============================================================================

function condicion_inicial(x, y, U0, δ, ε, Lx)
    Nx = length(x)
    Ny = length(y)
    ω0 = zeros(Ny, Nx)
    for j in 1:Ny
        for i in 1:Nx
            # Vorticidad base: -dU/dy = -U0/δ * sech²(y/δ)
            base = -U0 / δ * sech(y[j] / δ)^2
            # Perturbación: modo fundamental en x
            perturb = ε * cos(2π * x[i] / Lx) * sech(y[j] / δ)^2
            ω0[j, i] = base + perturb
        end
    end
    return ω0
end

# =============================================================================
# DERIVADA TEMPORAL (lado derecho de ∂ω/∂t)
# Ecuación de vorticidad: ∂ω/∂t + u·∇ω = ν∇²ω
# En espectral: û = ∂ψ/∂y, v̂ = -∂ψ/∂x, ψ̂ = -ω̂/K²
# =============================================================================

function rhs!(dω, ω, KX, KY, K2_safe, K2, mask, ν, plan_fwd, plan_inv)
    # --- Espacio espectral ---
    ω_hat = plan_fwd * ω

    # Función de corriente: ψ̂ = -ω̂ / K²
    ψ_hat = -ω_hat ./ K2_safe

    # Velocidades en espectral
    u_hat =  im .* KY .* ψ_hat
    v_hat = -im .* KX .* ψ_hat

    # Gradiente de vorticidad en espectral
    dωdx_hat = im .* KX .* ω_hat
    dωdy_hat = im .* KY .* ω_hat

    # --- Espacio físico (para término no lineal) ---
    u    = real(plan_inv * u_hat)
    v    = real(plan_inv * v_hat)
    dωdx = real(plan_inv * dωdx_hat)
    dωdy = real(plan_inv * dωdy_hat)

    # Término advectivo (no lineal): u·∇ω
    advec = u .* dωdx .+ v .* dωdy

    # Volver al espectral y aplicar dealising
    advec_hat = plan_fwd * advec
    advec_hat .*= mask

    # Difusión: ν∇²ω = -ν·K²·ω̂
    difus_hat = -ν .* K2 .* ω_hat

    # RHS en espectral: -(advección) + difusión
    rhs_hat = -advec_hat .+ difus_hat

    # Volver al espacio físico
    dω .= real(plan_inv * rhs_hat)
    return nothing
end

# =============================================================================
# INTEGRADOR RK4
# =============================================================================

function rk4_step!(ω_new, ω, dt, KX, KY, K2_safe, K2, mask, ν, plan_fwd, plan_inv)
    dω1 = similar(ω)
    dω2 = similar(ω)
    dω3 = similar(ω)
    dω4 = similar(ω)
    tmp  = similar(ω)

    rhs!(dω1, ω,                   KX, KY, K2_safe, K2, mask, ν, plan_fwd, plan_inv)
    tmp .= ω .+ (dt/2) .* dω1
    rhs!(dω2, tmp,                  KX, KY, K2_safe, K2, mask, ν, plan_fwd, plan_inv)
    tmp .= ω .+ (dt/2) .* dω2
    rhs!(dω3, tmp,                  KX, KY, K2_safe, K2, mask, ν, plan_fwd, plan_inv)
    tmp .= ω .+ dt     .* dω3
    rhs!(dω4, tmp,                  KX, KY, K2_safe, K2, mask, ν, plan_fwd, plan_inv)

    ω_new .= ω .+ (dt/6) .* (dω1 .+ 2 .* dω2 .+ 2 .* dω3 .+ dω4)
    return nothing
end

# =============================================================================
# DIAGNÓSTICOS
# =============================================================================

function energia_cinetica(ω, K2_safe, plan_fwd)
    ω_hat = plan_fwd * ω
    ψ_hat = -ω_hat ./ K2_safe
    return 0.5 * sum(abs.(ψ_hat).^2) / length(ω)^2
end

function enstrofia(ω)
    return 0.5 * sum(ω.^2) / length(ω)
end

# =============================================================================
# MAIN
# =============================================================================

function main()
    println("="^60)
    println(" Kelvin-Helmholtz — Simulación Pseudo-Espectral 2D")
    println("="^60)
    println("  Grid:        $(Nx) × $(Ny)")
    println("  Re ~ $(Int(round(U0*Ly/ν)))")
    println("  dt = $dt,  T = $T,  pasos = $Nt")
    println("  Espesor ϵ = $(ε)")
    println("  Aplitud de la perturbación: δ = $(δ)")
    println("="^60)

    # Condición inicial
    ω = condicion_inicial(x, y, U0, δ, ε, Lx)

    # Planes FFTW (precompilados para máxima velocidad)
    plan_fwd = plan_fft(ω)
    plan_inv = plan_ifft(plan_fwd * ω)

    ω_new = similar(ω)

    # Listas para diagnósticos
    tiempos    = Float64[]
    energias   = Float64[]
    enstrofias = Float64[]

    # Colección de frames para animación
    frames = []

    # Guardar archivo .dat con los datos de energía y enstrofia
    data_file = open(joinpath(ouput_dir, "diagnosticos.dat"), "w")
    println(data_file, "# t     E(t)    Z(t)")
    println(data_file, "# nu = $(ν)  delta = $(δ)  eps = $(ε)  U0 = $(U0)")

    println("\nIniciando simulación...\n")
    t_inicio = time()

    for n in 1:Nt
        rk4_step!(ω_new, ω, dt, KX, KY, K2_safe, K2, mask, ν, plan_fwd, plan_inv)
        ω .= ω_new

        # Dealising adicional al campo
        ω_hat = plan_fwd * ω
        ω_hat .*= mask
        ω .= real(plan_inv * ω_hat)

        if n % n_save == 0
            t = n * dt
            E = energia_cinetica(ω, K2_safe, plan_fwd)
            Z = enstrofia(ω)
            push!(tiempos, t)
            push!(energias, E)
            push!(enstrofias, Z)

            # Frame de vorticidad
            push!(frames, copy(ω))

            # Escribir en el .dat
            @printf(data_file, "%.6f  %.10e  %.10e\n", t, E, Z)

            elapsed = time() - t_inicio
            @printf("  t = %6.2f  |  E = %.4e  |  Z = %.4e  |  %.1f s\n", t, E, Z, elapsed)
        end
    end
    #Cerrar .dat
    close(data_file)


    println("\nSimulación completada. Generando visualizaciones...")

    # --- Animación del campo de vorticidad ---
    clim_val = 5.0
    anim = @animate for (i, frame) in enumerate(frames)
        t = tiempos[i]
        heatmap(
            x, y, frame,
            clim=(-clim_val, clim_val),
            color=:RdBu,
            xlabel="x", ylabel="y",
            title=@sprintf("Vorticidad ω(x,y)   t = %.2f", t),
            aspect_ratio=:equal,
            size=(700, 600),
            colorbar_title="ω"
        )
    end
    gif(anim, joinpath(ouput_dir,"kh_delta$(δ)_eps$(ε).gif"), fps=10)
    println("  → kh_delta$(δ)_eps$(ε).gif guardado")

    # --- Energía cinética en el tiempo ---
    p1 = plot(tiempos, energias,
        xlabel="Tiempo", ylabel="Energía cinética E(t)",
        title="Evolución de la energía cinética",
        lw=2, color=:steelblue, legend=false,
        yscale=:log10, size=(700, 400)
    )
    savefig(p1, joinpath(ouput_dir,"en_delta$(δ)_eps$(ε).png"))
    println("  → en_delta$(δ)_eps$(ε).png guardado")

    # --- Enstrofía en el tiempo ---
    p2 = plot(tiempos, enstrofias,
        xlabel="Tiempo", ylabel="Enstrofía Z(t)",
        title="Evolución de la enstrofía",
        lw=2, color=:crimson, legend=false,
        size=(700, 400)
    )
    savefig(p2, joinpath(ouput_dir,"ens_delta$(δ)_eps$(ε).png"))
    println("  → ens_delta$(δ)_eps$(ε).png guardado")

    # --- Espectro de energía final ---
    ω_final = frames[end]
    ω_hat   = plan_fwd * ω_final
    ψ_hat   = -ω_hat ./ K2_safe
    E_hat   = 0.5 .* abs.(ψ_hat).^2

    k_mag = sqrt.(K2)
    k_bins = range(0, stop=maximum(k_mag), length=60)
    E_spectrum = zeros(length(k_bins)-1)
    k_centers  = zeros(length(k_bins)-1)

    for i in 1:length(k_bins)-1
        idx = (k_mag .>= k_bins[i]) .& (k_mag .< k_bins[i+1])
        E_spectrum[i] = sum(E_hat[idx])
        k_centers[i]  = 0.5 * (k_bins[i] + k_bins[i+1])
    end

    # Filtrar bins vacíos
    valid = E_spectrum .> 0
    p3 = plot(k_centers[valid], E_spectrum[valid],
        xlabel="Número de onda k",
        ylabel="E(k)",
        title="Espectro de energía (t final)",
        lw=2, color=:darkorange, legend=false,
        xscale=:log10, yscale=:log10,
        size=(700, 400)
    )
    # Referencia ley k^(-3) (turbulencia 2D)
    k_ref = k_centers[valid][5:end-5]
    E_ref = E_spectrum[valid][5] .* (k_ref ./ k_centers[valid][5]).^(-3)
    plot!(p3, k_ref, E_ref, lw=1.5, ls=:dash, color=:gray, label="k⁻³")

    savefig(p3, joinpath(ouput_dir,"espectro.png"))
    println("  → espectro.png guardado")

    println("\n¡Listo! Archivos generados:")
    println("  kh_delta$(δ)_eps$(ε).gif  — animación del campo de vorticidad")
    println("  en_delta$(δ)_eps$(ε).png           — energía cinética vs tiempo")
    println("  ens_delta$(δ)_eps$(ε).png         — enstrofía vs tiempo")
    println("  espectro.png          — espectro de energía E(k)")
end

main()