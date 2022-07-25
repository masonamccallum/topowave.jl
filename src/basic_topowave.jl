using LinearAlgebra, Plots
using ProgressBars
using UnPack

σ₁ = [0 1; 1 0];
σ₂ = [0 -im; im 0];
σ₃ = [1 0; 0 -1];
ħ = 1; 
#initial(x,y) = ℯ^(π*im*(x+y));
initial(x,y) = ℯ^(-100*((x - 0.5)^2 + (y - 0.5)^2));

#Note: This version of readGmsh2D is a local change I have
#made to the readGmsh2D function in StartUpDG
VXY, EToV, groups = readGmsh2D("data/pert_mesh.msh",true) 

begin
VXY, EToV, groups = readGmsh2D("data/pert_mesh.msh",true); #TODO: add groups to multiDomain simple .msh
rd = RefElemData(Tri(), 3);
md = MeshData(VXY,EToV,rd);
mp = MeshPlotter(rd,md)
md = make_periodic(md)
end

begin
@unpack x, y = md;
u = @. initial(x,y);
Plots.plot(x,y,real(u),leg=false)
#Plots.plot(x,y,imag(u),leg=false)
end

#TODO: construct rhs!
#TODO: readGMsh2D version 4 read function
#TODO: fork StartUpDG for read update