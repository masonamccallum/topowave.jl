using OrdinaryDiffEq
using LinearAlgebra, Plots
using StartUpDG
using Javis

VXY, EToV = readGmsh2D("data/tri_unstruc_rect.msh");
rd = RefElemData(Tri(), 5);
md = MeshData(VXY,EToV,rd::RefElemData);
plot(MeshPlotter(rd,md))
md = make_periodic(md);
plot(MeshPlotter(rd,md))

@unpack x, y = md;
u = @. exp(-100 * (x^2 + y^2));

# du/dt + du/dx = 0, periodic
function rhs!(du, u, parameters, t)
    @unpack rd, md = parameters
    @unpack Vf, Fmask, Dr, Ds, LIFT = rd
    @unpack rxJ, sxJ, ryJ, syJ, J, nxJ, nyJ, sJ = md
    @unpack uf, uP, flux_u, ur, us, dudx, lifted_flux = parameters

    mul!(uf, Vf, u)  # uf = Vf * u
    @. uP = uf[md.mapP]    
    @. flux_u = 0.5 * (uP - uf) * nxJ 

    mul!(ur, Dr, u)
    mul!(us, Ds, u)
    mul!(lifted_flux, LIFT, flux_u)
    @. dudx = rxJ * ur + sxJ * us
    @. du = dudx + lifted_flux
    du ./= -J
end

tspan = (0.0, 1.0)
parameters = (; rd, md, uf=similar(md.xf), uP=similar(md.xf), flux_u=similar(md.xf), 
                ur=similar(md.x), us=similar(md.x), dudx=similar(md.x), lifted_flux=similar(md.x))
prob = ODEProblem(rhs!, u, tspan, parameters)
sol = solve(prob, Tsit5(), save_everystep=false, reltol=1e-7, abstol=1e-7)

scatter(x, y, zcolor=sol.u[end], leg=false)
#scatter(map(x->vec(rd.Vp*x), (md.xyz..., getindex.(sol.u[end], 1))), 
         #zcolor=vec(rd.Vp*getindex.(sol.u[end], 1)), msw=0, leg=false, cam=(0,0))
