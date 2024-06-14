module PPPModule 
include("IPEAModule.jl")
using LinearAlgebra, .IPEAModule
const ipea = IPEAModule

function importCUDA()
    @eval begin
        using CUDA
    end
end

# Define function to get gamma_ii, which is the one center coloumb repulsion.
function calGammaSameSite(atomZ, atomZe, ia_param)
    ion_pot, elec_aff = ipea.getIPAndEA(atomZ, atomZe, ia_param)
    return ion_pot - elec_aff
end

# Define function to get gamma_ij, which is the two center coloumb repulsion
# Beveridge-Hinze
function calGammaCrossSiteBH(atom1_Z, atom1_Ze, atom2_Z, atom2_Ze, dist_12, ia_param)
    gamma_11 = calGammaSameSite(atom1_Z, atom1_Ze, ia_param)
    gamma_22 = calGammaSameSite(atom2_Z, atom2_Ze, ia_param)
    a12 = 2.0/(gamma_11 + gamma_22)
    r12 = dist_12
    gamma_12 = 1.0/(a12*exp(-0.5*(r12^2)/(a12^2)) + r12)
    return gamma_12
end

# Mataga-Nishimoto
function calGammaCrossSiteMN(atom1_Z, atom1_Ze, atom2_Z, atom2_Ze, dist_12, ia_param)
    gamma_11 = calGammaSameSite(atom1_Z, atom1_Ze, ia_param)
    gamma_22 = calGammaSameSite(atom2_Z, atom2_Ze, ia_param)
    a12 = 2.0/(gamma_11 + gamma_22)
    r12 = dist_12
    gamma_12 = 1/(a12 + r12)
    return gamma_12
end

function buildGammaMatrix(atomdata, dist_mat, natomic, param, ia_param)
    gamma_mat = zeros(Float64, (natomic, natomic))
    for i in 1:natomic
        gamma_mat[i,i] = calGammaSameSite(atomdata[i,1], atomdata[i,2], ia_param)
        if (param == "BH")
            for j in 1:natomic
                if (i!=j)
                    gamma_mat[i,j] = calGammaCrossSiteBH(atomdata[i,1], atomdata[i,2], 
                                                 atomdata[j,1], atomdata[j,2], dist_mat[i,j], ia_param)
                end
            end
        elseif (param == "MN")
            for j in 1:natomic
                if (i!=j)
                    gamma_mat[i,j] = calGammaCrossSiteMN(atomdata[i,1], atomdata[i,2], 
                                                 atomdata[j,1], atomdata[j,2], dist_mat[i,j], ia_param)
                end
            end
        end
    end
    return gamma_mat
end
            
# Now a lot of helper functions that are needed for 2ppi-2ppi overlap 
# Getting the exponent from the same site coloumb repulsion term
function getExponent(gamma)
    exponent = (1280.0/501.0)*gamma # got it from Beveridge-Hinze
    return exponent
end

# Calculating p and t parameters from exponents and distance
function calPAndTParameters(gamma_1, gamma_2, dist_12)
    exponent_1 = getExponent(gamma_1)
    exponent_2 = getExponent(gamma_2)
    r12 = dist_12
    p = 0.5*r12*(exponent_1 + exponent_2)
    t = (exponent_1 - exponent_2)/(exponent_1 + exponent_2)
    return (p, t)
end

function calA(k, p)
    diff_sum = 0.0
    for i in 1:(k+1)
        diff_sum = diff_sum + factorial(k)/(p^i * factorial(k - i + 1))
    end
    return exp(-p)*diff_sum
end

function calB(k, p, t)
    diff_sum = 0.0
    for i in 1:k+1
        diff_sum = diff_sum + (1 + (-1)^(k - i)*exp(2*p*t))*(factorial(k)/((p*t)^i * factorial(k - i + 1)))
    end
    return (-exp(-p*t))*diff_sum
end

# Bring the helper functions together
function calTwoPiTwoPiOverlap(gamma_1, gamma_2, dist_12, nndist)
    if (abs(dist_12 - nndist) > 0.1)
        return 0
    end
    p, t = calPAndTParameters(gamma_1, gamma_2, dist_12)
    overlap = 0
    if (p == 0)
        overlap = (1 - t^2)^2.5
    elseif (t == 0)
        overlap = exp(-p)*(1 + p + 0.4*p^2 + (1.0/15.0)*p^3)
    else
        overlap = (1.0/32.0) * p^5 * (1 - t^2)^2.5 * (calA(4, p)* (calB(0, p, t) - calB(2, p, t)) +
        calA(2, p)*(calB(4, p, t) - calB(0, p, t)) + calA(0, p)*(calB(2, p, t) - calB(4, p, t)))
    end
    return overlap
end

# Now that we have overlap, we can calculate beta
# Beveridge-Hinze parametrization
function calBetaBH(atom1_Ze, atom2_Ze, gamma_1, gamma_2, gamma_12, dist_12, nndist)
    Cfit = 0.545
    r12 = dist_12
    overlap = calTwoPiTwoPiOverlap(gamma_1, gamma_2, r12, nndist)
    beta = 0.5*(atom1_Ze + atom2_Ze)*overlap*(gamma_12 - (2*Cfit/r12))
    return beta
end

# Mataga-Nishimoto parametrization
function calBetaMN(atom1_Z, atom2_Z, dist_12, nndist)
    if (abs(dist_12 - nndist) > 0.1)
        return 0
    end
    if ((atom1_Z == 7 && atom2_Z == 6) 
        || (atom1_Z == 6 && atom2_Z == 7))
        return -2.576/27.2  # Taken from the Pariser-Parr papers (Part 2)
    elseif (atom1_Z == 6 && atom2_Z == 6)
        return -2.388/27.2
    else
        return -6442*exp(-3.009*dist_12)/27.2
    end
end

function buildBetaMatrix(atomdata, natomic, dist_mat, mindist, gamma_mat, param)
    beta_mat = zeros(Float64, (natomic, natomic))
    if (param == "BH")
        for i in 1:natomic
            for j in 1:natomic
                if (i!=j)
                    beta_mat[i,j] = calBetaBH(atomdata[i,2], atomdata[j,2], gamma_mat[i,i], gamma_mat[j,j], 
                                               gamma_mat[i,j], dist_mat[i,j], mindist[i])
                end
            end
        end
    elseif (param == "MN")
        for i in 1:natomic
            for j in 1:natomic
                if (i!=j)
                    beta_mat[i,j] = calBetaMN(atomdata[i,1], atomdata[j,1], dist_mat[i,j], mindist[i])
                end
            end
        end
    end
    return beta_mat
end
    
# Building the density matrix
function buildDensityMatrixCPU(coeffs, natomic, norbitals)
    #density_gpu = zeros(Float64, (natomic, natomic))
    density = 2*transpose(coeffs[1:norbitals,1:natomic])*coeffs[1:norbitals,1:natomic]
    #coeffs_gpu = CuArray(coeffs[1:norbitals,1:natomic])
    #density_gpu = transpose(coeffs_gpu)*coeffs_gpu
    #density = 2*Array(density_gpu)
    #=
    Threads.@threads for j in 1:natomic
        for k in 1:natomic
            for i in 1:norbitals
                density[j, k] = density[j,k] + 2*coeffs[i,j]*coeffs[i,k]
            end
        end
    end
    =#
    return density
end

function buildDensityMatrixGPU(coeffs, natomic, norbitals)
    coeffs_gpu = CuArray(coeffs[1:norbitals,1:natomic])
    density_gpu = CuArray(zeros(Float64, (natomic, natomic)))
    CUDA.CUBLAS.gemm!('T', 'N', 2.0, coeffs_gpu, coeffs_gpu, 0, density_gpu)
    return Array(density_gpu)
end

# Normalize vector
function normalizeVector(in_vec)
    norm = 0
    for i in eachindex(in_vec)
        norm = norm + in_vec[i]^2
    end   
    out_vec = (1.0/sqrt(norm))*in_vec
    return out_vec
end

# Build the core hamiltonian matrix
function buildCoreMatrix(atomdata, natomic, gamma_mat, beta_mat, ia_param)
    core = zeros(Float64, (natomic, natomic))
    Threads.@threads for i in 1:natomic
        ip_i, ea_i = ipea.getIPAndEA(atomdata[i, 1], atomdata[i, 2], ia_param)
        core[i,i] = -ip_i
        for j in 1:natomic
            if (j != i)
                core[i,i] = core[i,i] - atomdata[j, 2]*gamma_mat[i,j]
                core[i,j] = beta_mat[i,j]
            end
        end
    end
    return core
end

# Building the fock matrix
function buildFockMatrix(core, density, natomic, gamma_mat)
    fock = zeros(Float64, (natomic, natomic))
    Threads.@threads for i in 1:natomic
        fock[i,i] = core[i,i] + 0.5*density[i,i]*gamma_mat[i,i]
        for j in 1:natomic
            if (j != i)
                fock[i,i] = fock[i,i] + density[j,j]*gamma_mat[i,j]
                fock[i,j] = core[i,j] - 0.5*density[i,j]*gamma_mat[i,j]
            end
        end
    end
    return fock
end

# Solving the eigenproblem on CPU
function solveFockMatrixCPU(fock, natomic)
    sys = eigen(Symmetric(fock))
    newcoeffs = zeros(Float64, (natomic, natomic))
    for i in 1:natomic
        newcoeffs[i,:] = normalizeVector(sys.vectors[:,i])
    end
    return (sys.values, newcoeffs)
end

#Solving the eigenproblem on GPU
function solveFockMatrixGPU(fock, natomic)
        fock_gpu = CuArray(fock)
        values_gpu, vectors_gpu = CUDA.CUSOLVER.Xsyevd!('V','U',fock_gpu)
        values = Array(values_gpu)
        vectors = Array(vectors_gpu)
    newcoeffs = zeros(Float64, (natomic, natomic))
    for i in 1:natomic
        newcoeffs[i,:] = normalizeVector(vectors[:,i])
    end
    return (values, newcoeffs)
end

#Calculate error
function calError(newcoeffs, oldcoeffs, nmo)
    err = 0
    for i in 1:nmo
        err = err + abs(newcoeffs[i] - oldcoeffs[i])
    end
    return err
end

#Calculate total energies
function calTotalEnergy(fock, core, density, natomic)
    oes = 0
    for j in 1:natomic
        for k in 1:natomic
            oes = oes + 0.5*density[j,k]*(core[j,k] + fock[j,k])
        end
    end
    return oes
end

function doSCF(niter, tol, atomdata, natomic, norbitals, dist_mat, mindist, param, ia_param, print_orbitals, use_gpu)
    C =  sqrt(1/natomic)*ones(Float64, (norbitals, natomic))
    F = zeros(Float64, (natomic, natomic))
    #P = zeros(Float64, (natomic, natomic))
    gamma_mat = buildGammaMatrix(atomdata, dist_mat, natomic, param, ia_param)
    beta_mat = buildBetaMatrix(atomdata, natomic, dist_mat, mindist, gamma_mat, param)
    core = buildCoreMatrix(atomdata, natomic, gamma_mat, beta_mat, ia_param)
    epsilon = zeros(Float64, natomic)
    exchange = 0
    if (use_gpu)
        importCUDA()
        solver_function = solveFockMatrixGPU
        density_function = buildDensityMatrixGPU
    else
        solver_function = solveFockMatrixCPU
        density_function = buildDensityMatrixCPU
    end
    for iter in 1:niter
        println("SCF Iteration: ",iter)
        P = density_function(C, natomic, norbitals)
        F = buildFockMatrix(core, P, natomic, gamma_mat)
        new_epsilon, C = solver_function(F, natomic)
        error = calError(new_epsilon, epsilon, norbitals)
        println("Error: ", error)
        if (error < tol)
            total_energy = calTotalEnergy(F, core, P, natomic)
            exchange = calculateExchangeSingletTriplet(atomdata, C, natomic, norbitals, gamma_mat)
            printConvergedVariables(epsilon, C, natomic, norbitals, total_energy, exchange, print_orbitals)
            epsilon = new_epsilon
            break
        end
        epsilon = new_epsilon
    end
    return F, C, epsilon
end

function printConvergedVariables(Econv, Cconv, natomic, norbitals, Ten, exchange, print_orbitals)
    println("---------------------------------")
    println("CONVERGED\n")
    println("Orbital Eigenvalues (eV)")
    for i in 1:norbitals
        println(Econv[i]*27.2)
    end
    println("Unoccupied")
    for i in norbitals+1:natomic
        println(Econv[i]*27.2)
    end
    println("Total energy (eV): ", Ten*27.2)
    println("Delta E_ST (eV): ", exchange*27.2)
    if (print_orbitals)
        println("\nEigenvectors")
        for i in 1:norbitals
            for j in 1:natomic
                print(Cconv[i,j],"   ")
            end
            print("\n\n")
        end
        println("Unoccupied")
        for i in norbitals+1:natomic
            for j in 1:natomic
                print(Cconv[i,j],"   ")
            end
            print("\n\n")
        end
    end
    println("---------------------------------")
end

function calculateExchangeSingletTriplet(atomdata, conv_vectors, natomic, norbitals, gamma_mat)
    K = 0
    for i in 1:natomic
        for j in 1:natomic
            if (i == j)
                K = K + (conv_vectors[norbitals,i]^2)*(conv_vectors[norbitals+1,i]^2)*gamma_mat[i,i]
            else
                K = K + conv_vectors[norbitals,i]*conv_vectors[norbitals,j]*conv_vectors[norbitals+1,i]*conv_vectors[norbitals+1,j]*gamma_mat[i,j]
            end
        end
    end
    return 2*K
end

end