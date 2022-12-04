using SpecialFunctions
using Plots

#Konstante
j = sqrt(complex(-1))
n = 2
m = 1
fnm = 3.5*10^9
a = 1.3 / 100
b = 2.3 / 100
c = 2.3 / 100
ϵr = 2.52
ϵ0 = 8.854 * 10^-12
µ0 = 4 * π * 10^-7
ω = 2 * π * fnm
k0 = ω * sqrt(µ0 * ϵ0)
knm = ω * sqrt(µ0 * ϵ0 * ϵr)
r = 1
V0 = 1

#Funkcija
Fnm(ρ) = (bessely(n-1,knm*a)-bessely(n+1,knm*a))/2 * besselj(n,knm*ρ) -
         (besselj(n-1,knm*a)-besselj(n+1,knm*a))/2 * bessely(n,knm*ρ)

#Funkcija 1
function Eθ(θ, ϕ = 0) 
    -(j)^n * (ℯ^(-j * k0 * r) / r) * (k0 * V0 / 2 * Fnm(c)) * cos(n * ϕ)
    *((b * Fnm(b) * (besselj(n-1, k0 * b * sin(θ)) - besselj(n+1, k0 * b * sin(θ))))
    -a * Fnm(a) * (besselj(n-1, k0 * a * sin(θ)) - besselj(n+1, k0 * a * sin(θ))))
    end

#Funkcija 2
function Eϕ(θ, ϕ = π/2) 
    -(j)^n * (ℯ^(-j * k0 * r) / r) * (k0 * V0 / 2 * Fnm(c)) * cos(θ) * sin(n * ϕ)
    *((b * Fnm(b) * (besselj(n-1, k0 * b * sin(θ)) + besselj(n+1, k0 * b * sin(θ))))
    -a * Fnm(a) * (besselj(n-1, k0 * a * sin(θ)) + besselj(n+1, k0 * a * sin(θ))))
    end

    x = range(start=0,stop=π/2,length=10^4)

    EθMAX = maximum(Eθ,x)
    EϕMAX = maximum(Eϕ,x)

    Eθn(θ) = 20 * log(abs(Eθ(θ))/abs(EθMAX))
    Eϕn(θ) = 20 * log(abs(Eϕ(θ))/abs(EϕMAX))

    xplot = range(start=0,stop=π/2,length=10^3)

begin
    #prvi primjer
    k1 = plot(Eθn , xplot, title="Prvi Kartezijev Eθn", size=(500,500), linewidth=1.3)
    k2 = plot(Eϕn , xplot, title="Prvi Kartezijev Eϕn", size=(500,500), linewidth=1.3)
    p1 = plot(Eθn , xplot, title="Prvi Polarni Eθn", size=(500,500), linewidth=1.3, proj = :polar)
    p2 = plot(Eϕn , xplot, title="Prvi Polarni Eϕn", size=(500,500), linewidth=1.3, proj = :polar)
end

begin
    #drugi primjer
    n = 3
    fnm = 5.1*10^9
    k3 = plot(Eθn , xplot, title="Drugi Kartezijev Eθn", size=(500,500), linewidth=1.3)
    k4 = plot(Eϕn , xplot, title="Drugi Kartezijev Eϕn", size=(500,500), linewidth=1.3)
    p3 = plot(Eθn , xplot, title="Drugi Polarni Eθn", size=(500,500), linewidth=1.3, proj = :polar)
    p4 = plot(Eϕn , xplot, title="Drugi Polarni Eϕn", size=(500,500), linewidth=1.3, proj = :polar)
end

begin
    #treci primjer
    n = 1
    m = 2
    fnm = 4.5*10^9
    a = 2 / 100
    b = 4.05 / 100
    c = 2.5 / 100
    k5 = plot(Eθn , xplot, title="Treći Kartezijev Eθn", size=(500,500), linewidth=1.3)
    k6 = plot(Eϕn , xplot, title="Treći Kartezijev Eϕn", size=(500,500), linewidth=1.3)
    p5 = plot(Eθn , xplot, title="Treći Polarni Eθn", size=(500,500), linewidth=1.3, proj = :polar)
    p6 = plot(Eϕn , xplot, title="Treći Polarni Eϕn", size=(500,500), linewidth=1.3, proj = :polar)
end

begin
    #cetvrti primjer
    #ovaj primjer je biran za usporedbu velikih brojeva
    n = 4
    m = 6
    fnm = 9*10^9
    a = 3 / 100
    b = 10 / 100
    c = 5 / 100
    k7 = plot(Eθn , xplot, title="Četvrti Kartezijev Eθn", size=(500,500), linewidth=1.3)
    k8 = plot(Eϕn , xplot, title="Četvrti Kartezijev Eϕn", size=(500,500), linewidth=1.3)
    p7 = plot(Eθn , xplot, title="Četvrti Polarni Eθn", size=(500,500), linewidth=1.3, proj = :polar)
    p8 = plot(Eϕn , xplot, title="Četvrti Polarni Eϕn", size=(500,500), linewidth=1.3, proj = :polar)
end

k12 = plot(k1, k2, layout=(2,1), legend = false)
p12 = plot(p1, p2, layout=(2,1), legend = false)
k34 = plot(k3, k4, layout=(2,1), legend = false)
p34 = plot(p3, p4, layout=(2,1), legend = false)
k56 = plot(k5, k6, layout=(2,1), legend = false)
p56 = plot(p5, p6, layout=(2,1), legend = false)
k78 = plot(k7, k8, layout=(2,1), legend = false)
p78 = plot(p7, p8, layout=(2,1), legend = false)

display(plot(k12,p12))
display(plot(k34,p34))
display(plot(k56,p56))
display(plot(k78,p78))

#=
png(k12,raw"C:\Users\gabri\Downloads\k12")
png(p12,raw"C:\Users\gabri\Downloads\p12")
png(k34,raw"C:\Users\gabri\Downloads\k34")
png(p34,raw"C:\Users\gabri\Downloads\p34")
png(k56,raw"C:\Users\gabri\Downloads\k56")
png(p56,raw"C:\Users\gabri\Downloads\p56")
png(k78,raw"C:\Users\gabri\Downloads\k78")
png(p78,raw"C:\Users\gabri\Downloads\p78")
=#
