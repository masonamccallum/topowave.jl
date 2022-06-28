using StartUpDG
using UnPack
using Plots
using OrdinaryDiffEq
using LinearAlgebra, Plots

rd = RefElemData(Quad(), 4);
md = MeshData(uniform_mesh(Quad(), 16)..., rd);
md = make_periodic(md)

@unpack x,y = md;
E = @. sin(pi*x)*sin(pi*y) .*(x<0);
H = zeros(Np,K)

function rhs!(dE, E, dH, H, parameters, t)
    @unpack rd, md = parameters
    @unpack Vf, Fmask, Dr, Ds, LIFT = rd
    @unpack rxJ, sxJ, ryJ, syJ, J, nxJ, nyJ, sJ = md
    @unpack uf, uP, flux_u, ur, us, dudx, lifted_flux = parameters

end
tspan = (0.0, 1.0)
parameters = (; rd, md, uf=similar(md.xf), uP=similar(md.xf), flux_u=similar(md.xf), 
                ur=similar(md.x), us=similar(md.x), dudx=similar(md.x), lifted_flux=similar(md.x));
prob = ODEProblem(rhs!, u, tspan, parameters);
sol = solve(prob, Tsit5(), save_everystep=false, reltol=1e-7, abstol=1e-7);