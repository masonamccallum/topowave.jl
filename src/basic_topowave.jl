using LinearAlgebra, Plots
using ProgressBars
using UnPack
using StartUpDG
using OrdinaryDiffEq

# ∫(∂ₜu+σ₁uₓ+σ₂uᵥ+imσ₃u)φₙdx = -∫̂n⋅(σ₁(u*-uᵏ)+σ₂(u*-uᵏ))φₙdΓ
# Stiffness Term: ∫(σ₁uₓ+σ₂uᵥ+imσ₃u)φₙdx  #NOTE:v=y wont let me type u\_y
# Boundary Term: -∫̂n⋅(σ₁(u*-uᵏ)+σ₂(u*-uᵏ))φₙdΓ

#TODO: intialize wave based on D. Relation
#TODO: define RHS! 

σ₁ = [0 1; 1 0];
σ₂ = [0 -im; im 0];
σ₃ = [1 0; 0 -1];
ħ = 1; 
initial(x,y) = ℯ^(-(1/200)*im*((x-50)^2+(y-50)^2));

begin
import_options = MeshImportOptions(true,true)
VXY, EToV, grouping = readGmsh2D_v4("data/mesh_no_pert_v4.msh", import_options) 
#VXY, EToV, grouping = readGmsh2D_v4("data/pert_mesh_v4.msh", import_options) 
rd = RefElemData(Tri(), 3);
md = MeshData(VXY, EToV, rd);
mp = MeshPlotter(rd, md)
md = make_periodic(md)
theme(:dracula)
plot(mp, size=(1000,1000),linecolor="white")
@unpack x,y = md;
scatter!(x,y, markersize=3, markerstrokewidth=0)
end

begin
@unpack x,y = md;
u = @. initial(x, y);
xflat = reshape(x, length(x), 1)
yflat = reshape(y, length(y), 1)
uflat = reshape(u, length(u), 1)
u_unique = unique(hcat(xflat, yflat, real.(uflat), imag.(uflat)), dims=1)
x = u_unique[:,1]
y = u_unique[:,2]
u_real = u_unique[:,3]
u_imag = u_unique[:,4]
Plots.scatter(x[1:end], y[1:end], u_real[1:end], leg=false, markersize=0.5)
end

function rhs!(du, u, parameters, t)
    @unpack rd, md = parameters
    @unpack Vf, Fmask, Dr, Ds, LIFT = rd
    @unpack rxJ, sxJ, ryJ, syJ, J, nxJ, nyJ, sJ = md
    @unpack uf, uP, flux_u, ur, us, dudx, lifted_flux = parameters

    #=
    mul!(uf, Vf, u)
    @. uP = uf[md.mapP]    
    @. flux_u = 0.5 * (uP - uf) * nxJ 

    mul!(ur, Dr, u)
    mul!(us, Ds, u)
    mul!(lifted_flux, LIFT, flux_u)
    @. dudx = rxJ * ur + sxJ * us
    @. du = dudx + lifted_flux
    du ./= -J
    =#
end

begin
tspan = (0.0, 1.0)
parameters = (; rd, md, uf=similar(md.xf), uP=similar(md.xf), flux_u=similar(md.xf),
    ur=similar(md.x), us=similar(md.x), dudx=similar(md.x), lifted_flux=similar(md.x))
prob = ODEProblem(rhs!, u, tspan, parameters)
sol = solve(prob, Tsit5(), reltol=1e-7, abstol=1e-7)
end