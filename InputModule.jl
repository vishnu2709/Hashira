module InputModule

input_values=Dict()

function setDefaultValues()
    input_values["sysname"] = "Unspecified"
    input_values["filename"] = "position.dat"
    input_values["mode"] = "PPP"
    input_values["restricted"] = "yes"
    input_values["print-orbitals"] = "no"
    input_values["plot-energies"] = "no"
    input_values["plot-orbitals"] = "no"
    input_values["niter"] = 100
    input_values["tol"] = 1e-5
    input_values["parametrization"] = "BH"
    input_values["ia"] = "BH"
    input_values["use-gpu"] = "no"
end

function readInputFile()
    setDefaultValues()
    indata = readlines("in.file")
    numlines = length(indata)
    for i in 1:numlines
        if contains(indata[i], "sysname")
            input_values["sysname"] = split(indata[i])[2]
        elseif contains(indata[i],"filename")
            input_values["filename"] = split(indata[i])[2]
        elseif contains(indata[i],"mode")
            input_values["mode"] = split(indata[i])[2]
        elseif contains(indata[i],"restricted")
            input_values["restricted"] = split(indata[i])[2]
        elseif contains(indata[i], "print-orbitals")
            input_values["print-orbitals"] = split(indata[i])[2]
        elseif contains(indata[i],"plot-orbitals")
            input_values["plot-orbitals"] = split(indata[i])[2]
        elseif contains(indata[i],"plot-energies")
            input_values["plot-energies"] = split(indata[i])[2]
        elseif contains(indata[i],"niter")
            input_values["niter"] = parse(Int64, split(indata[i])[2])
        elseif contains(indata[i],"tol")
            input_values["tol"] = parse(Float32, split(indata[i])[2])
        elseif contains(indata[i],"parametrization")
            input_values["parametrization"] = split(indata[i])[2]
        elseif contains(indata[i],"ia")
            input_values["ia"] = split(indata[i])[2]
        elseif contains(indata[i], "use-gpu")
            input_values["use-gpu"] = split(indata[i])[2]
        end
    end
end

function getSystemName()
    return input_values["sysname"]
end

function getPositionFilename()
    return input_values["filename"]
end

function getMode()
    return input_values["mode"]
end

function isCalculationRestricted()
    if (input_values["restricted"] == "yes")
        return true
    else
        return false
    end
end

function printOrbitalOption()
    if (input_values["print-orbitals"] == "yes")
        return true
    else 
        return false
    end
end

function plotOrbitalOption()
    if (input_values["plot-orbitals"] == "yes")
        return true
    else
        return false
    end
end

function plotEnergyOption()
    if (input_values["plot-energies"] == "yes")
        return true
    else
        return false
    end
end

function getNumberOfIteration()
    return input_values["niter"]
end

function getTolerance()
    return input_values["tol"]
end

function getParametrization()
    return input_values["parametrization"]
end

function getIASetType()
    return input_values["ia"]
end

function useGPUOption()
    if (input_values["use-gpu"] == "yes")
        return true
    else
        return false
    end
end

end