
include("PPPModule.jl")
include("InputModule.jl")
include("SystemModule.jl")
include("PlottingModule.jl")
import .PPPModule
import .InputModule
import .SystemModule
import .PlottingModule
using BenchmarkTools
const ppp = PPPModule
const inp = InputModule
const sys_geo = SystemModule
const plt = PlottingModule

#--------------------------------------------------
# The actual running happens here

inp.readInputFile()

if (inp.useGPUOption())
    using CUDA
end

total_atoms, total_electrons, atoms = sys_geo.readPosData(inp.getPositionFilename())
mindist, dist = sys_geo.createDistMatrix(atoms, total_atoms)

# Determine number of atomic basis elements
natombasis = total_atoms

# Determine number of molecular orbitals. In the closed shell case, it will just be half of npi
if (inp.isCalculationRestricted())
    nmo = Int64(total_electrons/2)
end

# Solve
@btime begin
F, conv_vectors, conv_values = ppp.doSCF(inp.getNumberOfIteration(), inp.getTolerance(), atoms, 
                                        natombasis, nmo, dist, mindist, inp.getParametrization(), 
                                        inp.getIASetType(), inp.printOrbitalOption(), inp.useGPUOption())
end

if (inp.plotEnergyOption())
    plt.plotEnergyLevels(conv_values, nmo, inp.getSystemName())
end

if (inp.plotOrbitalOption())
    plt.plotOrbitalOccupation(total_atoms, mindist, dist, atoms, 
    conv_vectors[nmo,:], "HOMO for "*inp.getSystemName(), inp.getSystemName()*"_homo.pdf", true)
    
    plt.plotOrbitalOccupation(total_atoms, mindist, dist, atoms, conv_vectors[nmo+1,:], 
    "LUMO for "*inp.getSystemName(), inp.getSystemName()*"_lumo.pdf", true)
end
#--------------------------------------------------