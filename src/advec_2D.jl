using LinearAlgebra, Plots
using StartUpDG
using OrdinaryDiffEq


VXY, EToV,grouping = readGmsh2D_v4("data/mesh_no_pert_v4.msh",true);
#VXY, EToV = readGmsh2D_v4("data/no_group_v4.msh");
rd = RefElemData(Tri(), 3);
md = MeshData(VXY, EToV, rd);

rd = RefElemData(Tri(), Polynomial(), 3)
md = MeshData(uniform_mesh(Tri(), 30)..., rd)
mp = MeshPlotter(rd, md)
plot(mp)

md = make_periodic(md)

@unpack x, y = md;
u = @. exp(-10 * ((x+0.5)^2 + y^2));

# du/dt + du/dx = 0, periodic
function rhs!(du, u, parameters, t)
    @unpack rd, md = parameters
    @unpack Vf, Fmask, Dr, Ds, LIFT = rd
    @unpack rxJ, sxJ, ryJ, syJ, J, nxJ, nyJ, sJ = md
    @unpack uf, uP, flux_u, ur, us, dudx, lifted_flux = parameters

    mul!(uf, Vf, u)
    @. uP = uf[md.mapP]    
    @. flux_u = 0.5 * (uP - uf) * nxJ 

    mul!(ur, Dr, u)
    mul!(us, Ds, u)
    mul!(lifted_flux, LIFT, flux_u)
    @. dudx = rxJ * ur + sxJ * us
    @. du = dudx + lifted_flux
    du ./= -J
end

tspan = (0.0, 2.0)
parameters = (; rd, md, uf=similar(md.xf), uP=similar(md.xf), flux_u=similar(md.xf),
    ur=similar(md.x), us=similar(md.x), dudx=similar(md.x), lifted_flux=similar(md.x))
prob = ODEProblem(rhs!, u, tspan, parameters)
sol = solve(prob, Tsit5(), reltol=1e-7, abstol=1e-7)

Nsteps = length(sol.u)
using ProgressBars
@info "Starting animation creation"
anim1 = @animate for i in ProgressBar(600:800)
    plot(x, y, sol.u[i], leg=false)
end
gif(anim1, "./Advect_2D.gif", fps=60)