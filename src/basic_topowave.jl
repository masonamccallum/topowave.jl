using LinearAlgebra, Plots
using ProgressBars
using UnPack
using StartUpDG

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
file = "data/mesh_no_pert_v4.msh"
VXY, EToV, grouping = readGmsh2D_v4(file, import_options) 
#VXY, EToV, grouping = readGmsh2D_v4("data/pert_mesh_v4.msh", import_options) 
rd = RefElemData(Tri(), 3);
md = MeshData(VXY, EToV, rd);
md = make_periodic(md)
end

begin
@unpack x,y = md; # Compute nodes: (x,y)=ψ(r,s) (equation 6.3 on page 172 of hesthaven)
@unpack Fmask = rd
@unpack xf,yf,xq,yq = md # 12 edge nodes and 12 quadrature notes per element

mp = MeshPlotter(rd, md)
plot(mp, size=(3500,3500),linecolor=:black)
scatter!(x[Fmask[:],:],y[Fmask[:],:],markercolor="black",markersize=3,size=(1000,1000),markerstrokewidth=0,leg=false)
end


#   Jᵏ∂ₜu = -Dᵣσ₁u -Dₛσ₂u - i*m*σ₃u - M⁻¹∫̂n⋅(σ₁(u*-uᵏ)+σ₂(u*-uᵏ))φₙdΓ
#         TODO: Equation above missing geometric factors
#         TODO: add boundary conditions
#         TODO: intialize wave based on Dispersion Relation

function rhs!(du, u, parameters) #TODO: combine terms and optimize memory usage 
    @unpack rd, md = parameters
    @unpack Vf, Fmask, Dr, Ds, LIFT = rd
    @unpack rxJ, sxJ, ryJ, syJ, J, nxJ, nyJ, sJ = md
    @unpack uf, uP, flux_u, ur, us, dudx, lifted_flux = parameters
    @unpack dr_sig1_u, ds_sig2_u = parameters

    σ₁ = [0 1; 1 0];
    σ₂ = [0 -im; im 0];
    σ₃ = [1 0; 0 -1];

    #   -Jᵏ∂ₜu = Dᵣσ₁u + Dₛσ₂u + i*m*σ₃u + M⁻¹∫̂n⋅(σ₁(u*-uᵏ)+σ₂(u*-uᵏ))φₙdΓ
    #
    #   FLUX TERM
    #   M⁻¹∫̂n⋅(σ₁(uᵏ-u*)+σ₂(uᵏ-u*))φₙdΓ
    mul!(uf[:,:,1], complex.(Vf), u[:,:,1])
    mul!(uf[:,:,2], complex.(Vf), u[:,:,2])

    @. uP[:,:,1] = uf[:,:,1][md.mapP]
    @. uP[:,:,2] = uf[:,:,2][md.mapP]

    @. flux_u[:,:,1] =  0.5*(uP[:,:,1] - uf[:,:,1])
    @. flux_u[:,:,2] =  0.5*(uP[:,:,2] - uf[:,:,2])

    #@. flux_u[:,:,1] =  0.5*(uP[:,:,1]) + uf[:,:,1]
    #@. flux_u[:,:,2] =  0.5*(uP[:,:,2]) + uf[:,:,2]

    flat_flux_u = transpose(reshape(flux_u,size(flux_u,1)*size(flux_u,2),2))
    sig1_flux_u = reshape(σ₁ * flat_flux_u,size(flux_u,1),size(flux_u,2),2)  # [10*1160*2]
    sig2_flux_u = reshape(σ₂ * flat_flux_u,size(flux_u,1),size(flux_u,2),2)

    @. flux_u = sig1_flux_u + sig2_flux_u
    @. flux_u[:,:,1] = nxJ * (flux_u[:,:,1]) + nyJ * (flux_u[:,:,1])
    @. flux_u[:,:,2] = nxJ * (flux_u[:,:,2]) + nyJ * (flux_u[:,:,2])

    lifted_flux[:,:,1] .= LIFT * flux_u[:,:,1]
    lifted_flux[:,:,2] .= LIFT * flux_u[:,:,2]

    # dudx
    # Dᵣσ₁u + Dₛσ₂u + i*m*σ₃u
    uflat = transpose(reshape(u,size(u,1)*size(u,2),2)) # [2*11600]
    sig1_u = reshape(σ₁ * uflat,size(u,1),size(u,2),2)  # [10*1160*2]
    sig2_u = reshape(σ₂ * uflat,size(u,1),size(u,2),2)
    sig3_u = reshape(σ₃ * uflat,size(u,1),size(u,2),2)

    dr_sig1_u[:,:,1] = Dr*sig1_u[:,:,1]
    dr_sig1_u[:,:,2] = Dr*sig1_u[:,:,2]

    ds_sig2_u[:,:,1] = Ds*sig2_u[:,:,1]
    ds_sig2_u[:,:,2] = Ds*sig2_u[:,:,2]
     
    @. dudx[:,:,1] = rxJ * dr_sig1_u[:,:,1] + sxJ * ds_sig2_u[:,:,1] + im*0.2*sig3_u[:,:,1] 
    @. dudx[:,:,2] = rxJ * dr_sig1_u[:,:,2] + sxJ * ds_sig2_u[:,:,2] + im*0.2*sig3_u[:,:,2] 
    @. du = dudx - lifted_flux
    du ./= -J
    return du
end

# 610517
# Test my RHS! using function I know derivative of
# m2: amd queue gpu node

begin
tspan = (0.0, 1.0)
@unpack x,y = md;
    u_temp = @. initial(x, y);
    u = zeros(Complex,size(u_temp,1),size(u_temp,2),2)
    for i in 1:size(u_temp,1)
        for j in 1:size(u_temp,2)
            u[i,j,1] = u_temp[i,j][1]
            u[i,j,2] = u_temp[i,j][2]
        end
    end


parameters = (; rd, md, uf=zeros(Complex,12,size(u)[2],2), uP=zeros(Complex,12,size(u)[2],2), flux_u=zeros(Complex,12,size(u)[2],2),
    ur=similar(zeros(size(md.x,1),size(md.x,2),2),Complex),ds_sig2_u=similar(u),dr_sig1_u=similar(u), us=similar(md.x,Complex), dudx=similar(u), lifted_flux=similar(u))
end

function topo2d(u,parameters,FinalTime)
    rk4a = [            0.0,
            -567301805773.0/1357537059087.0,
            -2404267990393.0/2016746695238.0,
            -3550918686646.0/2091501179385.0,
            -1275806237668.0/842570457699.0];
    rk4b = [ 1432997174477.0/9575080441755.0,
            5161836677717.0/13612068292357.0,
            1720146321549.0/2090206949498.0,
            3134564353537.0/4481467310338.0,
            2277821191437.0/14882151754819.0];
            
    time=0
    resU = zeros(Complex,size(u)) # Fix
    du = zeros(Complex,size(u)) # Fix
    xmin = 1.0 
    CFL = 0.05;
    dt=CFL*xmin;
    Nsteps = ceil(FinalTime/dt);
    dt = FinalTime/Nsteps
    anim1 = @animate for tstep in ProgressBar(1:Nsteps)
        for INTRK=1:5
            rhsU = rhs!(du,u,parameters);
            resU = rk4a[INTRK]*resU + dt*rhsU;
            u = u+rk4b[INTRK]*resU;
            Plots.scatter(x, y, real.(u[:,:,1]), leg=false, markersize=0.5, zlims=(-2,2))
        end
        time = time+dt;
    end
    gif(anim1, "./topo.gif", fps=60)
    return u 
end
