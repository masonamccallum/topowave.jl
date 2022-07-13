using RecipesBase
using StartUpDG

# Our user-defined data type
struct T end
struct Mesh 
    VXY::Tuple{Vector{Float64},Vector{Float64}}
    EToV::Matrix
    function Mesh()
        VXY,EToV = StartUpDG.uniform_mesh(Tri(),2)
        return new(VXY,EToV)
    end
end

# This is all we define.  It uses a familiar signature, but strips it apart
# in order to add a custom definition to the internal method `RecipesBase.apply_recipe`
@recipe function plot(::T, n = 1; customcolor = :green)
    markershape --> :auto        # if markershape is unset, make it :auto
    markercolor :=  customcolor  # force markercolor to be customcolor
    xrotation   --> 45           # if xrotation is unset, make it 45
    zrotation   --> 90           # if zrotation is unset, make it 90
    rand(10,n)                   # return the arguments (input data) for the next recipe
end

@recipe function plot(::Mesh,VXY)
    markershape --> :auto
    return VXY
end
# ----------------------------

# Plots will be the ultimate consumer of our recipe in this example
using Plots
gr()

# This call will implicitly call `RecipesBase.apply_recipe` as part of the Plots
# processing pipeline (see the Pipeline section of the Plots documentation).
#   It will plot 5 line plots (a 5-column matrix is returned from the recipe).
#   All will have black circles:
#       - user override for markershape: :c == :circle
#       - customcolor overridden to :black, and markercolor is forced to be customcolor
#   If markershape is an unsupported keyword, the call will error.
#   By default, a warning will be shown for an unsupported keyword.
#   This will be suppressed for zrotation (:quiet flag).
plot(M)
plot(T(), 5; customcolor = :black, shape=:c)