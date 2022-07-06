"""
    function readGmsh2D(filename)

#reads triangular GMSH 2D file format 2.2 0 8. returns (VX, VY), EToV
reads triangular GMSH 2D file format 4.1 0 8. returns (VX, VY), EToV

# Examples
```julia
VXY, EToV = readGmsh2D("eulerSquareCylinder2D.msh")
```
"""
function readGmsh2D(filename)
    f = open(filename)
    lines = readlines(f)

    function findline(name, lines)
        for (i, line) in enumerate(lines)
            if line == name
                return i
            end
        end
    end

    node_start = findline("\$Nodes", lines) + 1
    Nv = parse(Int64, lines[node_start])
    VX, VY, VZ = ntuple(x -> zeros(Float64, Nv), 3)
    for i = 1:Nv
        vals = [parse(Float64, c) for c in split(lines[i+node_start])]
        # first entry =
        VX[i] = vals[2]
        VY[i] = vals[3]
    end

    elem_start = findline("\$Elements", lines) + 1
    K_all = parse(Int64, lines[elem_start])
    K = 0
    for e = 1:K_all
        if length(split(lines[e+elem_start])) == 8
            K = K + 1
        end
    end
    EToV = zeros(Int64, K, 3)
    sk = 1
    for e = 1:K_all
        if length(split(lines[e+elem_start])) == 8
            vals = [parse(Int64, c) for c in split(lines[e+elem_start])]
            EToV[sk, :] .= vals[6:8]
            sk = sk + 1
        end
    end

    EToV = EToV[:, vec([1 3 2])] # permute for Gmsh ordering

    EToV = correct_negative_Jacobians!((VX, VY), EToV)

    return (VX, VY), EToV
end