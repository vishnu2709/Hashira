module IPEAModule

# Define function that returns IP and EA according to atomic number
function getIPAndEA(Z, Ze, ia_param)
    if (ia_param == "BH")
        return getIPAndEA_BH(Z, Ze)
    elseif (ia_param == "PS1")
        return getIPAndEA_PS1(Z, Ze)
    elseif (ia_param == "PS2")
        return getIPAndEA_PS2(Z, Ze)
    end
end

# Beveridge-Hinze values
function getIPAndEA_BH(Z, Ze)
    if (Z == 6)
        IP = Float32(11.16/27.2) # Converting to hartrees
        EA = Float32(0.03/27.2) # Converting to hartrees
    elseif (Z == 7 && Ze == 1)
        IP = Float32(14.12/27.2)
        EA = Float32(1.78/27.2)
    elseif (Z == 7 && Ze == 2)
        IP = Float32(28.71/27.2)
        EA = Float32(11.96/27.2)
    elseif (Z == 8 && Ze == 1)
        IP = Float32(17.70/27.2)
        EA = Float32(2.47/27.2)
    elseif (Z == 8 && Ze == 2)
        IP = Float32(34.08/27.2)
        EA = Float32(15.30/27.2)
    end
    return (IP, EA)
end

# Pritchard-Skinner
function getIPAndEA_PS1(Z, Ze)
    if (Z == 6)
        IP = Float32(10.94/27.2)
        EA = Float32(0.28/27.2)
    elseif (Z == 7 && Ze == 1)
        IP = Float32(13.83/27.2)
        EA = Float32(0.85/27.2)
    elseif (Z == 7 && Ze == 2)
        IP = Float32(27.5/27.2)
        EA = Float32(13.79/27.2)
    elseif (Z == 8 && Ze == 1)
        IP = Float32(17.28/27.2)
        EA = Float32(2.70/27.2)
    elseif (Z == 8 && Ze == 2)
        IP = Float32(35.30/27.2)
        EA = Float32(19.85/27.2)
    end
    return (IP, EA)
end

# Pritchard-Skinner (alternate values)
function getIPAndEA_PS2(Z, Ze)
    if (Z == 6 && Ze == 1)
        IP = Float32(11.42/27.2)
        EA = Float32(0.58/27.2)
    elseif (Z == 6 && Ze == 2)
        IP = Float32(21.43/27.2)
        EA = Float32(9.26/27.2)
    elseif (Z == 7 && Ze == 1)
        IP = Float32(14.49/27.2)
        EA = Float32(1.58/27.2)
    elseif (Z == 7 && Ze == 2)
        IP = Float32(27.5/27.2)
        EA = Float32(13.79/27.2)
    elseif (Z == 8 && Ze == 1)
        IP = Float32(17.76/27.2)
        EA = Float32(4.85/27.2)
    elseif (Z == 8 && Ze == 2)
        IP = Float32(35.30/27.2)
        EA = Float32(19.85/27.2)
    end
    return (IP, EA)
end

end