using LinearAlgebra, Plots
using ProgressBars
using UnPack
using StartUpDG
using OrdinaryDiffEq

#PDE: ħ∂ₜu = -σ₁uₓ-σ₂uᵥ-imσ₃u
#NOTE:v=y  VS code won't let me type u\_y
# ħ=1
# ∫(∂ₜu+σ₁uₓ+σ₂uᵥ+imσ₃u)φₙdx = -∫̂n⋅(σ₁(u*-uᵏ)+σ₂(u*-uᵏ))φₙdΓ (1 ≤ n ≤ Nₚ)
# uᵏ(x,t) = ∑uᵏ(xᵢ,t)ℓᵢᵏ(x)
# Mass Term: ∫(∂ₜu)φₙdx ≈ ħMᵏ∂ₜuᵏ
# Stiffness Term: ∫(σ₁uₓ+σ₂uᵥ+imσ₃u)φₙdx  
#                 ∫(φₙσ₁uₓ)dx + ∫(φₙσ₂uᵥ)dx + imᵏ∫(φₙσ₃u)dx ≈ -Dᵣσ₁u -Dₛσ₂u - i*m*M⁻¹M σ₃u
#                 ≈ -Dᵣσ₁u -Dₛσ₂u - i*m*σ₃u
# Boundary Term: -∫̂n⋅(σ₁(u*-uᵏ)+σ₂(u*-uᵏ))φₙdΓ ≈ 
# Putting it all together:
#   ∂ₜu = -M⁻¹sᵣσ₁u -M⁻¹Sₛσ₂u - M⁻¹⧆ - M⁻¹∫̂n⋅(σ₁(u*-uᵏ)+σ₂(u*-uᵏ))φₙdΓ
#   ∂ₜu = -Dᵣσ₁u -Dₛσ₂u - i*m*σ₃u - M⁻¹∫̂n⋅(σ₁(u*-uᵏ)+σ₂(u*-uᵏ))φₙdΓ

#TODO: intialize wave based on Dispersion Relation


#TODO: remove these global vars
σ₁ = [0 1; 1 0];
σ₂ = [0 -im; im 0];
σ₃ = [1 0; 0 -1];
ħ = 1; 

initial(x,y) = reshape([ 
    ℯ^(-(1/200)*im*((x-50)^2+(y-50)^2))
    ℯ^(-(1/200)*im*((x-50)^2+(y-50)^2))
],1,2);

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
u_temp = @. initial(x, y);
u = zeros(Complex,size(u_temp,1),size(u_temp,2),2)
for i in 1:size(u_temp,1)
    for j in 1:size(u_temp,2)
        u[i,j,1] = u_temp[i,j][1]
        u[i,j,2] = u_temp[i,j][2]
    end
end
Plots.scatter(x, y, real.(u[:,:,1]), leg=false, markersize=0.5)
end


#   ∂ₜu = -Dᵣσ₁u -Dₛσ₂u - i*m*σ₃u - M⁻¹∫̂n⋅(σ₁(u*-uᵏ)+σ₂(u*-uᵏ))φₙdΓ
function rhs!(du, u, parameters, t)
    @unpack rd, md = parameters
    @unpack Vf, Fmask, Dr, Ds, LIFT = rd
    @unpack rxJ, sxJ, ryJ, syJ, J, nxJ, nyJ, sJ = md
    @unpack uf, uP, flux_u, ur, us, dudx, lifted_flux = parameters

    σ₁ = [0 1; 1 0];
    σ₂ = [0 -im; im 0];
    σ₃ = [1 0; 0 -1];

    #mul!(uf, Vf, u)
    #@. uP = uf[md.mapP]    
    #@. flux_u = 0.5 * (uP - uf) * nxJ 

    uflat = transpose(reshape(u,size(u,1)*size(u,2),2))
    sig1_u = similar(uflat)
    sig2_u = similar(uflat)
    sig3_u = similar(uflat)
    mul!(sig1_u,σ₁,uflat)
    sig1_u = transpose(sig1_u)
    sig1_u = reshape(sig1_u,size(u,1),size(u,2),size(u,3))

    mul!(sig2_u,σ₂,uflat)
    sig2_u = transpose(sig2_u)
    sig2_u = reshape(sig2_u,size(u,1),size(u,2),size(u,3))

    Dr_test = reshape([Dr Dr],size(Dr,1),size(Dr,2),2)
    #ur =  Dr_test * sig1_u
    #=
    mul!(us, Ds, sig2_u)

    mul!(sig3_u,σ₃,u)
    sig3_u = im*sig3_u

    mul!(lifted_flux, LIFT, flux_u)
    @. dudx = rxJ * ur + sxJ * us - sig3_u
    @. du = dudx + lifted_flux
    du ./= -J
    =#
end;

begin
tspan = (0.0, 1.0)
parameters = (; rd, md, uf=similar(md.xf), uP=similar(md.xf), flux_u=similar(md.xf),
    ur=similar(zeros(size(md.x,1),size(md.x,2),2)), us=similar(md.x), dudx=similar(md.x), lifted_flux=similar(md.x))
end
prob = ODEProblem(rhs!, u, tspan, parameters)
sol = solve(prob, Tsit5(), reltol=1e-7, abstol=1e-7)