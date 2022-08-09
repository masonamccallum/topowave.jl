using UnPack
using Plots
using LinearAlgebra, Plots
using FastGaussQuadrature
using SpecialPolynomials
using SpecialFunctions
using SparseArrays
using ProgressBars
using Parameters
include("DG_1D.jl")

#1d solver for maxwells
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
    mapI=0; map0=0; vmapI=0; vmap0=0;

function rhs!(E,H,eps,mu,Np,k,Nfp,Nfaces)
    zimp = sqrt.(mu./eps);
    dE    = zeros(Nfp*Nfaces,k);
    dH    = zeros(Nfp*Nfaces,k);
    zimpm = zeros(Nfp*Nfaces,k);
    zimpp = zeros(Nfp*Nfaces,k);
    yimpm = zeros(Nfp*Nfaces,k);
    yimpp = zeros(Nfp*Nfaces,k);

    dE[:] = E[vmapM] - E[vmapP]
    dH[:] = H[vmapM] - H[vmapP];
    zimpm[:] = zimp[vmapM]
    zimpp[:] = zimp[vmapP]
    yimpm[:] = 1 ./zimpm
    yimpp[:] = 1 ./zimpp

    Ebc = -E[vmapB]
    Hbc = H[vmapB]
    dE[mapB] = E[vmapB] - Ebc
    dH[mapB] = H[vmapB] - Hbc
    #upwind flux
    fluxE = 1 ./(zimpm+zimpp).*(nx.*zimpp.*dH-dE)
    fluxH = 1 ./(yimpm+yimpp).*(nx.*yimpp.*dE-dH)

    rhsE = (-rx.*(Dr*H) + LIFT*(Fscale.*fluxE))./eps
    rhsH = (-rx.*(Dr*E) + LIFT*(Fscale.*fluxH))./mu
    return rhsE, rhsH
end

function Maxwell1D(E,H,Np,k,Nfp,Nfaces,eps,mu,x,FinalTime)
    time=0
    resE = zeros(Np,k)
    resH = zeros(Np,k)
    xmin = minimum(abs.(x[1,:]-x[2,:]))
    CFL = 1.0;
    dt=CFL*xmin;
    Nsteps = ceil(FinalTime/dt);
    dt = FinalTime/Nsteps
    anim1 = @animate for tstep in ProgressBar(1:Nsteps)
        for INTRK=1:5
            rhsE,rhsH = rhs!(E,H,eps,mu,Np,k,Nfp,Nfaces);
            resE = rk4a[INTRK]*resE + dt*rhsE;
            resH = rk4a[INTRK]*resH + dt*rhsH;
            E = E+rk4b[INTRK]*resE;
            H = H+rk4b[INTRK]*resH;
        end
        time = time+dt;
        p2=plot(x[:],E[:],leg=false);
        p1=plot(x[:],H[:],leg=false);
        plot!(p1,p2)
    end
    gif(anim1,"./maxwell_1D.gif",fps=60)
    return E,H
end

mesh = DG.Mesh1D(-1,1,80) #xmin,xmax,numberOfGridNodes
scheme = DG.Scheme_1D(N=6) #order, mesh

r = DG.JacobiGL(0,0,scheme.N);

V = DG.vandermonde1D(scheme.N,r);
Vr = DG.GradVandermonde1D(scheme.N,r);
Dr = DG.Dmatrix1D(scheme.N,r,V);
LIFT = DG.lift1D(V,scheme.Np,mesh.Nfaces,scheme.Nfp);

# Coords of nodes
va = mesh.EToV[:,1]';
vb = mesh.EToV[:,2]';
x = (ones(scheme.Np,1)*mesh.Vx[va] + 0.5*(r.+1)*(mesh.Vx[vb]-mesh.Vx[va])); # GLmesh

rx,J = DG.GeometricFactors1D(x,Dr);

# Masks for edge nodes
fmask1 = findall(abs.(r.+1).<scheme.NODETOL)';
fmask2 = findall(abs.(r.-1).<scheme.NODETOL)';
Fmask = [collect(fmask1); collect(fmask2)]';
Fx = x[Fmask[:],:];

# surface normals and inverse metrix at surface
nx = DG.Normals1D(scheme.Nfp,mesh.Nfaces,mesh.k);
Fscale = 1 ./J[Fmask[:],:];
# Connectivity matrix
EToE,EToF = DG.connect1D(mesh);  #TODO: package
# Connectivity maps
vmapM,vmapP,vmapB,mapB = DG.BuildMaps1D(Fmask,EToE,EToF,x,mesh,scheme);

#TODO: Material Parameters
eps1 = vcat(ones(Int(mesh.k/2)), 2*ones(Int(mesh.k/2)))';
mu1 = ones(1,mesh.k);

# Initial conditions
E = @. sin(pi*x).*(x<0);
H = zeros(scheme.Np,mesh.k);

FinalTime = 1.0;
epsilon = ones(scheme.Np,1)*eps1;
mu = ones(scheme.Np,1)*mu1
E,H = Maxwell1D(E,H,scheme.Np,mesh.k,scheme.Nfp,mesh.Nfaces,epsilon,mu,x,FinalTime)