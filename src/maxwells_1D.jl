using StartUpDG
using UnPack
using Plots
using OrdinaryDiffEq
using LinearAlgebra, Plots
using FastGaussQuadrature
using SpecialPolynomials

#1D solver for maxwells 

function MeshGen1D(xmin,xmax,k)
    Nv = k+1 # Number of vectors
    elementWidth = xmax - xmin
    Vx = [elementWidth*(i-1)/(k)+xmin for i in range(1,Nv)]
    EToV = zeros(k,2)
    for i = 1:k
        EToV[i,1] = i;
        EToV[i,2] = i+1;
    end
    return [Nv, Vx, k, EToV]
end

function GradVandermonde1D(N,r)
"""
Initialize the gradient of the modal basis
"""
    DVr = zeros(length(r),N+1) 
    for i in 0:N 
       GJ = GradJacobiP(r,0,0,i)
       DVr[:,i+1] = reshape(GJ,(length(r),1)) 
    end
    return DVr
end

function GradJacobiP(r, alpha, beta, N)
"""
    Evaluate the derivative of the Jacobi polynomial 
    at points r for order N
"""
    dP = zeros(length(r),1)
    if N!=0
        a = sqrt(N*(N+alpha+beta+1))
        coeffs = ones(N-1)
        JacobiP = a * Jacobi{alpha+1,beta+1}(coeffs)
        dP = JacobiP.(r)
    end
    return dP
end

function Dmatrix1D(N,r,V) #TODO: V matrix error effects this function [break point]
    Vr = GradVandermonde1D(N,r);
    Dr = Vr/V; # Vr*inv(V)
end

function vandermonde1D(N,x) #TODO: Not matching matlab [break point] norm needs fixed
    alpha = 0
    beta = 0
    V1D = zeros(length(x),N+1)
    for j = 1:N+1
        coeffs = ones(j)
        JacobiP = Jacobi{alpha,beta}(coeffs)
        V1D[:,j] = JacobiP.(x)./norm(JacobiP)
    end
    return V1D
end

function lift1D(Np,Nfaces,Nfp)
    Emat = zeros(Np,Nfaces*Nfp)
    Emat[1,1] = 1.0;
    Emat[Np,2] = 1.0;
    LIFT = V*(V'*Emat)
end

function GeometrixFactors1D(x,Dr)
    xr = Dr*x;
    J = xr;
    rx = 1 ./J
    return [xr,J,rx]
end

function Normals1D(Nfp,Nfaces,k)
    nx = zeros(Nfp*Nfaces,k);
    nx[1,:] .= -1.0;
    nx[2,:] .= 1;
    return nx;
end

#Order for Poly approximation
N = 6;
xmin = -2;
xmax = -xmin;
numGridPoints = 80;
Nv, Vx,k,EToV = MeshGen1D(xmin,xmax,numGridPoints);
Np = N + 1;
Nfp = 1;
Nfaces=2;


r,w = gausslobatto(N+1);
V = vandermonde1D(N,r);
Vr = GradVandermonde1D(N,r);
Dr = Dmatrix1D(N,r,V);
LIFT = lift1D(Np,Nfaces,Nfp);

# Coords of nodes
va = Int.(EToV[:,1]');
vb = Int.(EToV[:,2]');
x = @. Vx[va] + 0.5(r+1)*(Vx[vb]+Vx[va])

[rx,J] = GeometrixFactors1D(x,Dr)


#TODO: Material Parameters
eps1 = vcat(ones(Int(k/2)), 2*ones(Int(k/2)))';
mu1 = ones(1,k);
epsilon = ones(Np,1)*eps1;
mu = ones(Np,1)*mu1;

# Initial conditions
E = @. sin(pi*x).*(x<0);
H = zeros(Np,k);

function rhs!(E,H,eps,mu)
    Zimp = sqrt.(mu./eps);

    dE = zeros(Nfp,Nfaces,k);
    dH = zeros(Nfp,Nfaces,k);
    zimpm = zeros(Nfp,Nfaces,k);
    zimpp = zeros(Nfp,Nfaces,k);
    Yimpm = zeros(Nfp,Nfaces,k)
    Yimpp = zeros(Nfp,Nfaces,k)
end

tspan = (0.0,1.0)
prob = ODEProblem(rhs!,u,tspan);
sol = solve(prob,Tsit5(),reltol=1e-7,abstol=1e-7);