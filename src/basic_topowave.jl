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


initial(x,y) = reshape([ 
    ℯ^(-(1/200)*im*((x-50)^2+(y-50)^2))
    ℯ^(-(1/200)*im*((x-50)^2+(y-50)^2))
],1,2);

begin
import_options = MeshImportOptions(true,true)
file = "/home/mason/Documents/summer2022/topowave.jl/data/mesh_no_pert_v4.msh"
VXY, EToV, grouping = readGmsh2D_v4(file, import_options) 
#VXY, EToV, grouping = readGmsh2D_v4("data/pert_mesh_v4.msh", import_options) 
rd = RefElemData(Tri(), 3);
md = MeshData(VXY, EToV, rd);
md = make_periodic(md)
end


begin
 @info "Visualizing reference elements"
 @unpack r, s, rf, sf, rq, sq, Fmask = rd

 scatter(r, s, color=:blue, size=(800,600), markersize = 5,label="Mesh interior")
# scatter!(r[Fmask], s[Fmask], color=:red, markersize = 5, label="Mesh exterior")
# scatter!(rf, sf, color=:green, markersize = 5, label="Face mat")
# scatter(rq, sq, color=:black, markersize = 5, label="Quad mat")
# rp and sp are plotting nodes. I think these play nicer with plotting packages
#scatter!(rd.rp, rd.sp, color=:cyan, markersize = 5, label="equi nodes")
end

begin
# rd.Vf = vandermonde(...)/VDM interpolates from nodes to face nodes
# rd.Vf face quad interp map
# rd.Vq quad interp map
#  xf = rd.Vf * xyz
#  xq = rd.Vq * xyz

@unpack x,y = md; # Compute nodes: (x,y)=ψ(r,s) (equation 6.3 on page 172 of hesthaven)
@unpack Fmask = rd
@unpack xf,yf,xq,yq = md # 12 edge nodes and 12 quadrature notes per element

mp = MeshPlotter(rd, md)
plot(mp, size=(1500,1500),linecolor=:black)
#scatter!(x,y,markercolor="red",markersize=3,size=(1000,1000),markerstrokewidth=0,leg=false)
#scatter!(x[Fmask[:],:],y[Fmask[:],:],markercolor="white",markersize=3,size=(1000,1000),markerstrokewidth=0,leg=false)
#scatter!(xf,yf,markersize=2)
scatter!(xq,yq,markersize=2)
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

#   Jᵏ∂ₜu = -Dᵣσ₁u -Dₛσ₂u - i*m*σ₃u - M⁻¹∫̂n⋅(σ₁(u*-uᵏ)+σ₂(u*-uᵏ))φₙdΓ
#         TODO: Equation above missing geometric factors
#         TODO: understand mesh variety
#         TODO: add boundary conditions
#         TODO: add flux term
#         TODO: intialize wave based on Dispersion Relation

function rhs!(du, u, parameters, t)
    @unpack rd, md = parameters
    @unpack Vf, Fmask, Dr, Ds, LIFT = rd
    @unpack rxJ, sxJ, ryJ, syJ, J, nxJ, nyJ, sJ = md
    @unpack uf, uP, flux_u, ur, us, dudx, lifted_flux = parameters

    σ₁ = [0 1; 1 0];
    σ₂ = [0 -im; im 0];
    σ₃ = [1 0; 0 -1];

    #Vf: face quad interp map
    mul!(uf, complex.(Vf), u[:,:,1])
    @. uP = uf[md.mapP]
    @. flux_u = 0.5 * (uP - uf) * nxJ

    uflat = transpose(reshape(u,size(u,1)*size(u,2),2))
    sig1_u = similar(uflat)
    sig2_u = similar(uflat)
    sig3_u = similar(uflat)
    dr_sig1_u = similar(u)

    #=
    mul!(sig1_u,σ₁,uflat)
    sig1_u = transpose(sig1_u)
    sig1_u = reshape(sig1_u,size(u,1),size(u,2),size(u,3))
    dr_sig1_u = reshape([Dr*sig1_u[:,:,1] Dr*sig1_u[:,:,2]],size(u,1),size(u,2),2)

    dr_sig2_u = similar(u)
    mul!(sig2_u,σ₂,uflat)
    sig2_u = transpose(sig2_u)
    sig2_u = reshape(sig2_u,size(u,1),size(u,2),size(u,3))
    dr_sig2_u = reshape([Dr*sig2_u[:,:,1] Dr*sig2_u[:,:,2]],size(u,1),size(u,2),2)

    mul!(sig3_u,σ₃ ,uflat)
    sig3_u = transpose(sig3_u)
    sig3_u = reshape(sig3_u,size(u,1),size(u,2),size(u,3))
    return -dr_sig1_u-dr_sig2_u-im*0.2*sig3_u
    =#
end

begin
tspan = (0.0, 1.0)
parameters = (; rd, md, uf=similar(md.xf,Complex), uP=similar(md.xf,Complex), flux_u=similar(md.xf,Complex),
    ur=similar(zeros(size(md.x,1),size(md.x,2),2),Complex), us=similar(md.x,Complex), dudx=similar(md.x,Complex), lifted_flux=similar(md.x,Complex))
end

#prob = ODEProblem(rhs!, u, tspan, parameters)
#sol = solve(prob, Tsit5(), reltol=1e-7, abstol=1e-7)