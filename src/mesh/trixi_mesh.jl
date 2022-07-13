using GLMakie, HOHQMesh 
using Distributions

begin
box_proj = newProject("box_two_circles","out")
getPolynomialOrder(box_proj)
getMeshFileFormat(box_proj)
setPolynomialOrder!(box_proj,4)
setMeshFileFormat!(box_proj, "ABAQUS")
lower_left = [0.0, 0.0, 0.0]
spacing = [1.0, 1.0, 0.0]
num_intervals = [30, 15, 0]
xmax = spacing[1]*num_intervals[1]
ymax = spacing[2]*num_intervals[2]
addBackgroundGrid!(box_proj, lower_left, spacing, num_intervals)
plotProject!(box_proj, GRID)

num_spline_knots = 10;
irregularity_size_x = 0.1; #percentage of x range that is irregular
irregularity_size_y = 0.03; #percentage range of y irregularity
spline_x_start = xmax-irregularity_size_x*xmax
x = range(spline_x_start,xmax,length=num_spline_knots)

y_dev = ymax*irregularity_size_y;
y = [ymax/2 rand((ymax/2)-y_dev:0.1:(ymax/2)+y_dev,num_spline_knots-2)... ymax/2]

t = range(0,1,length=num_spline_knots)
spline_data = [[t[i], x[i], y[i], 0.0] for i in 1:num_spline_knots]
spline_data = reduce(vcat,transpose.(spline_data))

inner_spline = newSplineCurve("Spline", num_spline_knots, spline_data)
inner_line1 = newEndPointsLineCurve("Line 1", [xmax, ymax/2,0.0],[xmax, ymax, 0.0])
inner_line4 = newEndPointsLineCurve("Line 4",  [xmax, ymax, 0.0],[0.0, ymax, 0.0])
inner_line3 = newEndPointsLineCurve("Line 3", [0.0, ymax, 0.0],[0.0, ymax/2, 0.0])
inner_line2 = newEndPointsLineCurve("Line 2", [0.0, ymax/2, 0.0], [spline_x_start, ymax/2, 0.0])
#inner_line = newEndPointsLineCurve("Line", [xmax, ymax/2, 0.0], [0.0, ymax/2, 0.0])
end
inner_line = newEndPointsLineCurve("Line", [xmax, ymax/2, 0.0],[0.0, ymax/2, 0.0], )
inner_line1 = newEndPointsLineCurve("Line 1", [xmax, ymax/2,0.0],[xmax, ymax, 0.0])
inner_line4 = newEndPointsLineCurve("Line 4",  [xmax, ymax, 0.0],[0.0, ymax, 0.0])
inner_line3 = newEndPointsLineCurve("Line 3", [0.0, ymax, 0.0],[0.0, ymax/2, 0.0])
#inner_line2 = newEndPointsLineCurve("Line 2", [0.0, ymax/2, 0.0], [spline_x_start, ymax/2, 0.0])

addCurveToOuterBoundary!(box_proj, inner_line)
addCurveToOuterBoundary!(box_proj, inner_line1)
addCurveToOuterBoundary!(box_proj, inner_line4)
addCurveToOuterBoundary!(box_proj, inner_line3)


#addCurveToOuterBoundary!(box_proj, inner_spline)
#addCurveToOuterBoundary!(box_proj, inner_line1)
#addCurveToOuterBoundary!(box_proj, inner_line4)
#addCurveToOuterBoundary!(box_proj, inner_line3)
#addCurveToOuterBoundary!(box_proj, inner_line2)

generate_mesh(box_proj)