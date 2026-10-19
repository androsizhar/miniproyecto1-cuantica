using Plots, Roots

h = 1
L = 2
M = (pi^2)/8
U1 = 85
c = 2
U2 = c*U1
E = range(0.001, 84.999, length=3000)

Tangente = tan.((L/h)*sqrt.(2*M*E))
F2 = (sqrt.(E) .* (sqrt.(U1 .- E) .+ sqrt.(U2 .- E))) ./ (E .- sqrt.((U1 .- E) .* (U2 .- E)))

for i in 1:(length(Tangente)-1) #quiero quitar los valores que se rompen en la tangente
    if abs(Tangente[i+1] - Tangente[i]) > 160
        Tangente[i] = NaN
    end
end
for i in 1:(length(F2)-1) #quitamos lineas horizontales, ya que no representan algo
    if abs(F2[i+1] - F2[i]) > 20
        F2[i] = NaN
    end
end
g(e) = (L/h)*sqrt(2*M*e) - atan(sqrt(U1-e)/sqrt(e)) - atan(sqrt(U2-e)/sqrt(e)) #tomamos la fecuacion de 2.8
E_puntos = Float64[]
for n in 1:9 #encontramos que puntos corresponden
    f(e) = g(e) - (n - 1) * pi
    raiz = find_zero(f, (0.001, 84.999))
    push!(E_puntos, raiz)
end
Y_puntos = tan.((L/h) .* sqrt.(2*M .* E_puntos))# Evaluamos
grafica = plot(E, Tangente, ylims=(-160,160), label="tan(ωL)", xlabel="Energía E", ylabel="funciones(E)", title="Grafico de eigenenergias")
plot!(grafica,E,F2)
scatter!(E_puntos, Y_puntos, color=:red, markersize=5, label="Eigenenergías corte de funciones")
offsets_y = [25.0, -25.0, 70.0, -65.0, 25.0, -25.0, 25.0, 0.0, 25.0]

for i in 1:9
    e_val = round(E_puntos[i], digits=1)
    
    if Y_puntos[i] < -100
        annotate!(grafica, E_puntos[i] + 2.5, Y_puntos[i], text("E$i = $e_val", 7, :left, :bold))
    else
        annotate!(grafica, E_puntos[i], Y_puntos[i] + offsets_y[i], text("E$i = $e_val", 7, :center, :bold))
    end
end

display(grafica)
x = range(-1.5, L + 1.5, length=1000) #hagamos nuestros valores de x a graficar
region1 = x .< 0.0
region2 = (x .>= 0.0) .& (x .<= L)
region3 = x .> L
function phi_n(x, E_n, region1, region2, region3; h=1.0, L=2.0, M=(pi^2)/8.0, U1=85.0, U2=170.0)
    wn = sqrt(2 * M * E_n) / h
    k  = sqrt(2 * M * (U1 - E_n)) / h
    y  = sqrt(2 * M * (U2 - E_n)) / h
    beta = (sin(wn * L) * (k^2 + wn^2)) / (wn * (y + k))
    Calrevez = (1 / (2 * k)) + (L / 2) * (1 + (k^2 / wn^2)) + (sin(2 * wn * L) / (4 * wn)) * (1 - (k^2 / wn^2)) + (k / wn^2) * (sin(wn * L)^2) + (beta^2 / (2 * y))     
    C = 1 / sqrt(Calrevez)
    phi = zeros(length(x))
    phi[region1] .= C .* exp.(k .* x[region1])
    phi[region2] .= C .* (cos.(wn .* x[region2]) .+ (k / wn) .* sin.(wn .* x[region2]))
    phi[region3] .= C .* beta .* exp.(-y .* (x[region3] .- L))
    return phi
end
U = zeros(length(x))
U[region1] .= U1
U[region2] .= 0
U[region3] .= U2
graficas = plot(x,U,color=:black,label="U(x)",xlabel = "x",ylabel = "Energia",title = "Grafica de eigenestados",size(850,600))
for (n, E_n) in enumerate(E_puntos)
    phi_vec =phi_n(x, E_n, region1, region2, region3)
    plot!(graficas, x, E_n .+ 3* phi_vec, label="n = $n")
end
display(graficas)
b = 0.2
cancho = 0.4
v0 = 0.0
A = (2.0 / (pi * b^2))^(0.25)

phi_0 = A .* exp.(-((x .- cancho).^2) ./ (b^2)) .* exp.(1im .* (M * v0 / h) .* x)
dx = x[2] - x[1]
N_estados = length(E_puntos)
c_n = zeros(ComplexF64,N_estados)
for n in 1:length(E_puntos)
    phi_vec = phi_n(x, E_puntos[n], region1, region2, region3)
    c_n[n] = sum(phi_vec.*phi_0)*dx
end
for n in 1:length(c_n)
    println("c_$n = ", real.(c_n[n]))
end
suma = 0 
for or in 1:length(c_n)
    suma += (c_n[or]).^2
end
real(suma)
t_final = 5.0
t = range(0.0, t_final, length=200)
Nx = length(x)
Nt = length(t)
P_xt = zeros(Float64, Nx, Nt)
phi_m = zeros(Float64, Nx, N_estados)
for n in 1:N_estados
    phi_m[:, n] = phi_n(x, E_puntos[n], region1, region2, region3)
end
for j in 1:Nt
    tj = t[j]
    for k in 1:Nx
        xk = x[k]
        sumas = 0.0 + 0.0im
        
        for n in 1:N_estados
            phi_n_x = phi_m[k, n]
            cn = c_n[n]
            En = E_puntos[n]
            
            for m in 1:N_estados
                phi_m_x = phi_m[k, m]
                cm_conj = conj(c_n[m])
                Em = E_puntos[m]
                prob = cn * cm_conj * phi_n_x * phi_m_x * exp(-1im * (En - Em) * tj / h)
                sumas += prob
            end
        end
        P_xt[k, j] = real(sumas)
    end
end
grafica_2D = heatmap(x, t, P_xt',xlabel = "Posición x", ylabel = "Tiempo t",title = "Evolución de g(x,t)^2",color = :viridis,size = (800, 500))

display(grafica_2D)
valor_esperado = zeros(Float64, Nt)
for j in 1:Nt
    P_columna = P_xt[:, j]
    valor_esperado[j] = sum(x .* P_columna) * dx
end
valoresperadografica = plot(t, valor_esperado,linewidth = 2,color = :blue,xlabel = "Tiempo t",ylabel = "<x>(t)",title = "Valor esperado de x en el tiempo",size = (800, 400))
display(valoresperadografica)

Tope = maximum(P_xt)
animacion = @animate for j in 1:Nt
    tj = t[j]
    P_actual = P_xt[:, j]
    plot(
        x, P_actual,
        linewidth = 2,
        color = :blue,  label = "P(x,t) = |g(x,t)|²",xlabel = "x",ylabel = "P(x,t)",title = "Evolución del Pulso (t = $(round(tj, digits=2)))",xlims = (-0.3, L + 0.3),ylims = (0.0, Tope * 1.05),legend = :topright,size = (800, 450))
    vline!([0.0], color = :black, linestyle = :solid, linewidth = 1.5, label = "Paredes")
    vline!([L], color = :black, linestyle = :solid, linewidth = 1.5, label = "")
end
mp4(animacion, "pulso.mp4", fps = 15)
