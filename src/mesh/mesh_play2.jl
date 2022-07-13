#add Gmsh
using Parameters
import Gmsh: gmsh
x=1
begin
rm("data/t2.msh")
gmsh.initialize()
gmsh.model.add("t2")
lc = 6e-2;

x=0
for i=1:61
    gmsh.model.geo.addPoint(x, sin(x), 0, lc)
    x += pi/30
end
s1 = gmsh.model.geo.addSpline(1:61) ## for this

R = gmsh.model.geo.addPoint(2*pi,0,0,lc);
RT = gmsh.model.geo.addPoint(2*pi,2,0,lc);
LT = gmsh.model.geo.addPoint(0,2,0,lc);
L = gmsh.model.geo.addPoint(0,0,0,lc);
gmsh.model.geo.synchronize()
line_R = gmsh.model.geo.addLine(R,RT)
line_T = gmsh.model.geo.addLine(RT,LT)
line_L = gmsh.model.geo.addLine(LT,L)

gmsh.model.geo.addCurveLoop([line_R,line_T,line_L,s1], 1)
#gmsh.model.geo.addPlaneSurface([c2], 9) # curve [2]
gmsh.model.geo.synchronize()
#gmsh.model.addPhysicalGroup(2, [8], 10)
#gmsh.model.addPhysicalGroup(2, [9], 11)
#gmsh.model.mesh.generate(2)
gmsh.write("data/t2.msh")
gmsh.clear()
gmsh.finalize()
end