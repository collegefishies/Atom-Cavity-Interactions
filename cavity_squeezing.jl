
#Written by Cristoph Hotter (Mar 3, 2023)
using Pkg
# Pkg.add(PackageSpec(name="QuantumCumulants", version="0.2.13"))
# Pkg.add(PackageSpec(name="ModelingToolkit", version="8.21.0"))
# Pkg.add(PackageSpec(name="OrdinaryDiffEq", version="6.11.2"))
# Pkg.add(PackageSpec(name="PyPlot", version="2.11.0"))

println("Adding QuantumCumulants")
using QuantumCumulants
println("Added QC")
using ModelingToolkit
println("Added ModelingToolkit")
using OrdinaryDiffEq
println("Added OrdinaryDiffEq")
using PyPlot; pygui(true)
println("Added PyPlot")
println("Adding Hilbert Space")
hc = FockSpace(:cavity)
ha = NLevelSpace(:atom, 3)
h = hc ⊗ ha

Transition(h,:σ,1,1)

@qnumbers a::Destroy(h)
σ(α,β,i) = IndexedOperator(Transition(h,:σ,α,β),i)

@cnumbers Γ κ Ω g λ N Δc Δ2 Δ3
# @syms t::Real
# @register_symbolic Ωt(t)
# @register_symbolic λt(t)

i = Index(h,:i,N,ha)
j = Index(h,:j,N,ha)

println("Adding Hamiltonian")
H0 = -Δc*a'a - Δ2*∑(σ(2,2,i),i) - Δ3*∑(σ(3,3,i),i)
Hda = Ω*∑( σ(2,1,i) + σ(1,2,i), i)
Hdc = λ*(a' + a)
Hint = g*∑( (a'*σ(2,3,i) + a*(σ(3,2,i))) ,i)
H = H0 + Hda + Hdc + Hint

J = [a, σ(2,3,i)]
rates = [κ, Γ]

println("Caclulating Equations")
ops = [a'a, σ(2,2,j)]
eqs = meanfield(ops, H, J;rates, order=2)

eqs_c = complete(eqs)
eqs_sc = scale(eqs_c)
println("Length of eqs_sc: ", length(eqs_sc))
@named sys = ODESystem(eqs_sc)

# Initial state
println("Setting up Initial State")
u0_ = zeros(ComplexF64, length(eqs_sc))
# System parameters
#Γ κ Ω g λ N Δc Δ2 Δ3
N_ = 3
Γ_ = 0.36
κ_ = 1.0
Ω_ = 3.28
g_ = 1.15
λ_ = 1 #10
Δc_ = -9.67
Δ3_ = Δc_
Δ2_ = 0

ps = [Γ, κ, Ω, g, λ, N, Δc, Δ2, Δ3]
p0_ = [Γ_, κ_, Ω_, 0g_, 0λ_, N_, Δc_, Δ2_, Δ3_]
#
tΠ2 = π/4Ω_
prob_ = ODEProblem(sys,u0_,(0.0, tΠ2), ps.=>p0_)

println("Solving the ODE.")
sol_ = solve(prob_,Tsit5();abstol=1e-10, reltol=1e-10, maxiters=1e7)
sol_[σ(2,2,1)][end]

###
u0 = sol_.u[end]
p0 = [Γ_, κ_, Ω_, g_, λ_, N_, Δc_, Δ2_, Δ3_]

T = 250
prob = ODEProblem(sys,u0,(0.0, T), ps.=>p0)

sol = solve(prob,Tsit5(),maxiters=1e7)
nt = sol[a'a]
t = sol.t

println("Plotting")
close("n(t")
figure("n(t)")
plot(t, nt)

