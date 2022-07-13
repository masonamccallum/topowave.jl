#add Gmsh
using Parameters
import Gmsh: gmsh

gmsh.initialize()
gmsh.model.add("t1")
lc = 6e-2;
sp = spline_params(10,0.1,0.02);
d = domain(0,1,0,1,sp);
c1,c2,s1 = ConnectDomain(d,sp)
#ConnectDomain(d)

gmsh.model.geo.addPlaneSurface([c1, s1], 8) # [1] denotes which curve represents surface 8
gmsh.model.geo.addPlaneSurface([c2], 9) # curve [2]
gmsh.model.geo.synchronize()
gmsh.model.addPhysicalGroup(2, [8], 10)
gmsh.model.addPhysicalGroup(2, [9], 11)
gmsh.model.mesh.generate(2)
gmsh.write("data/t2.msh")
gmsh.clear()
gmsh.finalize()

@with_kw struct domain
    xmin::Float64
    xmax::Float64
    ymin::Float64
    ymax::Float64
    points::Vector{Int}
    knots::Vector{Int}

    function domain(xmin,xmax,ymin,ymax)
        dy = ymax - ymin;
        p1 = gmsh.model.geo.addPoint(xmin, ymin, 0, lc)
        p2 = gmsh.model.geo.addPoint(xmax, ymin, 0, lc)
        p3 = gmsh.model.geo.addPoint(xmax, ymax, 0, lc)
        p4 = gmsh.model.geo.addPoint(xmin, ymax, 0, lc)
        p5 = gmsh.model.geo.addPoint(xmin, ymin+dy/2, 0, lc)
        p6 = gmsh.model.geo.addPoint(xmax, ymax-dy/2, 0, lc)
        points = [p1, p2, p3, p4, p5, p6]
        knots = []
        return new(xmin,xmax,ymin,ymax,points, knots)
    end

    function domain(xmin,xmax,ymin,ymax,s_params::spline_params)
        @unpack num_spline_knots = s_params;
        @unpack irregularity_size_x = s_params;
        dy = ymax - ymin;
        irr_x = xmax -  irregularity_size_x * (xmax - xmin);
        irr_y = ymin+dy/2;
        p1 = gmsh.model.geo.addPoint(xmin, ymin, 0, lc);
        p2 = gmsh.model.geo.addPoint(xmax, ymin, 0, lc);
        p3 = gmsh.model.geo.addPoint(xmax, ymax, 0, lc);
        p4 = gmsh.model.geo.addPoint(xmin, ymax, 0, lc);
        p5 = gmsh.model.geo.addPoint(xmin, ymin+dy/2, 0, lc);
        p6 = gmsh.model.geo.addPoint(xmax, ymax-dy/2, 0, lc);
        p_irr = gmsh.model.geo.addPoint(irr_x,irr_y,0,lc);

        s_knots = spline_knots(xmin,xmax,ymin,ymax,sp);
        knots = zeros(Int,size(s_knots,1));
        for i=1:size(s_knots,1)
            knots[i] = gmsh.model.geo.addPoint(s_knots[i,2],s_knots[i,3],0,lc)
        end
        knots = reverse(knots)

        points = [p1, p2, p3, p4, p5, p6, p_irr];
        return new(xmin,xmax,ymin,ymax,points,knots)
    end
end

function ConnectDomain(domain,sp)
    @unpack points, knots = domain;
    l1 = gmsh.model.geo.addLine(points[1], points[2])
    l2 = gmsh.model.geo.addLine(points[6], points[2])
    l3 = gmsh.model.geo.addLine(points[6], points[7]) ## Swap this
    l4 = gmsh.model.geo.addLine(points[7], points[5])
    l5 = gmsh.model.geo.addLine(points[5], points[1])

    l6 = gmsh.model.geo.addLine(points[6], points[3])
    l7 = gmsh.model.geo.addLine(points[3], points[4])
    l8 = gmsh.model.geo.addLine(points[4], points[5])
    l9 = gmsh.model.geo.addLine(points[5], points[7])
    l10 = gmsh.model.geo.addLine(points[7], points[6])

    l11 = gmsh.model.geo.addLine(points[4], points[6]) #remove this

    s1 = gmsh.model.geo.addSpline(knots) ## for this
    gmsh.model.mesh.setTransfiniteCurve()
    c1 = gmsh.model.geo.addCurveLoop([l4, l5, l1, l6], 1)
    c2 = gmsh.model.geo.addCurveLoop([l6, l7, l11], 2)
    return c1, c2, s1
end

function ConnectDomain(domain)
    @unpack points = domain;
    l1 = gmsh.model.geo.addLine(points[1], points[2])
    l2 = gmsh.model.geo.addLine(points[6], points[2])
    l3 = gmsh.model.geo.addLine(points[6], points[5])
    l4 = gmsh.model.geo.addLine(points[5], points[1])
    l5 = gmsh.model.geo.addLine(points[6], points[3])
    l6 = gmsh.model.geo.addLine(points[3], points[4])
    l7 = gmsh.model.geo.addLine(points[4], points[5])


    c1 = gmsh.model.geo.addCurveLoop([l1, -l2, l3, l4], 1)
    c2 = gmsh.model.geo.addCurveLoop([-l3, l5, l6, l7], 2)
    return [c1, c2]
end

@with_kw struct spline_params
    num_spline_knots::Int;
    irregularity_size_x::Float64; #percentage of x range that is irregular
    irregularity_size_y::Float64; #percentage range of y irregularity
end

function  spline_knots(xmin,xmax,ymin,ymax, s_params) 
    @unpack num_spline_knots = s_params
    @unpack irregularity_size_x,irregularity_size_y  = s_params

    spline_x_start = xmax-irregularity_size_x*xmax
    x = range(spline_x_start,xmax,length=num_spline_knots)

    y_dev = (ymax-ymin)*irregularity_size_y;
    y = [ymax/2 rand((ymax/2)-y_dev:y_dev/100:(ymax/2)+y_dev,num_spline_knots-2)... ymax/2]

    t = range(0,1,length=num_spline_knots)
    spline_data = [[t[i], x[i], y[i], 0.0] for i in 1:num_spline_knots]
    spline_data = reduce(vcat,transpose.(spline_data))
    return spline_data
end