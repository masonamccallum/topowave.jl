module DG 
using Parameters
using LinearAlgebra
using SparseArrays
using SpecialFunctions

    @with_kw struct Mesh1D
        xmin::Real;
        xmax::Real;
        k::Real; 
        EToV::Matrix;
        Vx::Array;
        Nv::Int64;
        Nfaces::Int64;
    end

    function Mesh1D(xmin,xmax,k)
        Nv = k+1 # Number of vectors
        elementWidth = xmax - xmin
        Nfaces = 2;
        Vx = [elementWidth*(i-1)/(k)+xmin for i in range(1,Nv)]
        EToV = zeros(Int,k,2)
        for i = 1:k
            EToV[i,1] = i;
            EToV[i,2] = i+1;
        end
        return Mesh1D(xmin,xmax,k,EToV,Vx,Nv,Nfaces)
    end

    @with_kw  struct Scheme_1D
        N::Int = 6;
        Np::Int = 7;# N+1
        Nfp::Int = 1;
        NODETOL::Float64 = 1e-10;
    end

    function Scheme_1D(n) #TODO: Must provide default vals for all member vars need fix
        return Scheme_1D(N=n,Np=n+1)
    end

    function JacobiGL(alpha,beta,N)
        x = zeros(N+1,1);
        if N==1
            x[1] = -1.0;
            x[2] = 1.0;
            return x
        else
            xint,w = JacobiGQ(alpha+1,beta+1,N-2);
            x = [-1,xint,1]';
            x = [Float64(i) for element in x for i in element]
        end
        return x
    end

    function sub2ind(dim::Tuple,rows,cols)
        ind = zeros(Int,length(rows)) 
        @assert length(rows)==length(cols)
        for i=1:length(rows)
            @assert cols[i]<=dim[2]
            @assert rows[i]<=dim[1]
            ind[i] = ((cols[i]-1)*dim[1])+rows[i]
        end
        return ind
    end

    function JacobiGQ(alpha, beta, N)
        if N==0
            x[1] = -(alpha-beta)/(alpha+beta+2);
            w[1] = 2;
            return
        else
            J = zeros(N+1,N+1);
            h1 = 2*(0:N).+(alpha+beta)
        end
        J = Tridiagonal(zeros(N),-1/2*(alpha^2-beta^2)./(h1.+2)./h1,2 ./
            (h1[1:N].+2).*sqrt.((1:N).*((1:N).+(alpha+beta)).*
            ((1:N).+alpha).*((1:N).+beta)./(h1[1:N].+1)./(h1[1:N].+3)))

        J = Matrix(J);
        if (alpha+beta<10) #eps?
            J[1,1]=0.0;
        end
        J = J + J';

        V = eigvals(J)
        D = eigvecs(J)
        x = V
        w = (V[1,:]').^2*2*(alpha+beta+1)/(alpha+beta+1)*gamma(alpha+1)*
            gamma(beta+1)/gamma(alpha+beta+1)
        return x,w
    end

    function JacobiP(x,alpha,beta,N)
        xp=x;
        xp = reshape(xp,length(xp),1)
        dims = size(xp);
        if dims[2]==1
            xp = xp'
        end
        P = zeros(length(xp))
        PL = zeros(N+1,length(xp));
        #Initial values P_0(x) and P_1(x)
        gamma0 = 2^(alpha+beta+1)/(alpha+beta+1)*gamma(alpha+1)*
                gamma(beta+1)/gamma(alpha+beta+1);
        PL[1,:] .= 1 / sqrt(gamma0);
        if N==0
            P=PL';
            return P
        end
        gamma1 = (alpha+1)*(beta+1)/(alpha+beta+3)*gamma0;
        tt = ((alpha+beta+2)*xp/2 .+ (alpha-beta)/2)/sqrt(gamma1)
        PL[2,:]=tt;
        if N==1
            P=PL[N+1,:]'
            return P
        end
        aold = 2/(2+alpha+beta)*sqrt((alpha+1)*(beta+1)/(alpha+beta+3));
        for i in 1:N-1
            h1 = 2*i+alpha+beta;
            anew = 2/(h1+2)*sqrt((i+1)*(i+1+alpha+beta)*(i+1+alpha)*(i+1+beta)/(h1+1)/(h1+3))
            bnew = -(alpha^2-beta^2)/h1/(h1+2);
            PL[i+2,:] = 1/anew *(-aold*PL[i,:]' .+ (xp.-bnew).*PL[i+1,:]');
            aold = anew;
        end
        P = PL[N+1,:]'
        return P
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
            #JacobiP = a * Jacobi{alpha+1,beta+1}(coeffs)
            #dP = JacobiP.(r)
            dP = a*JacobiP(r[:],alpha+1,beta+1,N-1);
        end
        return dP
    end

    function Dmatrix1D(N,r,V) #TODO: V matrix error effects this function [break point]
        Vr = GradVandermonde1D(N,r);
        Dr = Vr/V; # Vr*inv(V)
    end

    function GradVandermonde1D(N,r)
        """
        # Initialize the gradient of the modal basis
        """
        DVr = zeros(length(r),N+1) 
        for i in 0:N 
        GJ = GradJacobiP(r,0,0,i)
        DVr[:,i+1] = reshape(GJ,(length(r),1)) 
        end
        return DVr
    end

    function vandermonde1D(N,x) #TODO: Not matching matlab [break point] norm needs fixed
        a= 0
        b= 0
        V1D = zeros(length(x),N+1)
        for j=1:N+1
            coeffs = ones(j)
            #JacobiPoly = Jacobi{a,b}(coeffs)
            #V1D[:,j] = JacobiPoly.(x)./norm(JacobiPoly)
            V1D[:,j] = JacobiP(x[:],0,0,j-1);
        end
        return V1D
    end

    function lift1D(V,Np,Nfaces,Nfp)
        Emat = zeros(Np,Nfaces*Nfp)
        Emat[1,1] = 1.0;
        Emat[Np,2] = 1.0;
        LIFT = V*(V'*Emat)
    end

    function GeometricFactors1D(x,Dr)
        #xr = Dr*x;
        J = Dr*x; #xr;
        rx = 1 ./ J
        return rx,J
    end

    function find(x, threshold=0)
    return findall(vcat((x.<threshold)...))
    end

    function Normals1D(Nfp,Nfaces,k)
        nx = zeros(Nfp*Nfaces,k);
        nx[1,:] .= -1.0;
        nx[2,:] .= 1;
        return nx;
    end

    function GeometricFactors1D(x,Dr)
        #xr = Dr*x;
        J = Dr*x; #xr;
        rx = 1 ./ J
        return rx,J
    end

    function BuildMaps1D(Fmask,EToE,EToF,x,m::Mesh1D, s::Scheme_1D)
        @unpack EToV, Nfaces, k = m
        @unpack Np,Nfp,NODETOL = s 
        nodeids = reshape(collect(1:Np*k),Np,k);
        vmapM = zeros(Int,Nfp,Nfaces,k);
        vmapP = zeros(Int,Nfp,Nfaces,k);

        for k1=1:k
            for f1=1:Nfaces
                vmapM[:,f1,k1] = nodeids[Fmask[:,f1],k1]; #sort right and left faces #TODO: ERROR
            end
        end

        for k1=1:k
            for f1=1:Nfaces
                # find neighbor
                k2 = Int.(EToE[k1,f1]);
                f2 = Int.(EToF[k1,f1]);

                vidM = Int.(vmapM[:,f1,k1]);
                vidP = Int.(vmapM[:,f2,k2]);
                x1 = x[vidM];
                x2 = x[vidP];
                D = (x1-x2).^2
                if (D.<NODETOL)==ones(size(D))
                    vmapP[:,f1,k1]=vidP;
                end
            end
        end

        vmapP=vmapP[:];
        vmapM=vmapM[:];

        mapB = findall(vmapP.==vmapM);
        vmapB = vmapM[mapB];
        mapI=1;
        map0=k*Nfaces
        vmapI=1;
        vmap0=k*Np;
        return vmapM,vmapP,vmapB,mapB
    end

    function connect1D(m::DG.Mesh1D)
        @unpack EToV, Nfaces, k= m
        k = size(EToV,1)::Int
        TotalFaces = Nfaces*k
        Nv = k + 1
        vn = [1,2]::Vector{Int}
        FToV = spzeros(TotalFaces,Nv) 
        sk = 1::Int
        for i in 1:k
            for face in 1:Nfaces
                FToV[sk,EToV[i,vn[face]]] = 1;
                sk = sk + 1
            end
        end

        # explore problem with minmax of interval being -1 to 1 vs. -2 to 2
        FToF = FToV*FToV' - sparse(I,TotalFaces,TotalFaces)
        inds = findall(FToF.==1)
        faces1 = getindex.(inds,1)
        faces2 = getindex.(inds,2)

        element1 = @. floor((faces1-1)/Nfaces) + 1;
        face1 = @. mod((faces1-1),Nfaces) + 1;
        element2 = @. floor((faces2-1)/Nfaces) + 1;
        face2 = @. mod((faces2-1),Nfaces) + 1;

        # Rearrange into Nelements x Nfaces sized arrays
        EToE = (1:k)*ones(Int,1,Nfaces);
        EToF = ones(Int,k,1)*(1:Nfaces)'
        ind = sub2ind((k,Nfaces),element1,face1);
        EToE[ind] = element2;
        EToF[ind] = face2;
        return EToE,EToF
    end
end