module SystemModule

# Read position data
function readPosData(filename)
    posdata = readlines(filename)
    numlines = length(posdata)
    total_atoms = parse(Int64, posdata[1])
    total_electrons = 0
    atoms = zeros(Float32, (total_atoms, 5))
    for i in 3:numlines
        atom = split(posdata[i])
        for j in 1:5
            atoms[i-2, j] = parse(Float32, atom[j])
        end
        total_electrons = total_electrons + atoms[i-2, 2]
        atoms[i-2, 3:5] = atoms[i-2,3:5]/0.529177210903
    end
    new_atoms = transformCoordinates(total_atoms, atoms)
    return total_atoms, total_electrons, new_atoms
end

# If 2D material, then system will be rotated onto x-y plane.
function transformCoordinates(total_atoms, atomdata)
    # Shift origin to the first atom
    red_atoms = zeros(Float32, (total_atoms, 5))
    for i in 1:total_atoms
        red_atoms[i,1:2] = atomdata[i,1:2]
        for j in 3:5
            red_atoms[i,j] = atomdata[i,j] - atomdata[1,j]
        end
    end

    # Define the data structure for the transformed system
    new_atoms = zeros(Float32,  (total_atoms, 5))

    # Check if system is already along xy/yz/xz planes
    total_z = 0
    total_x = 0
    total_y = 0
    for i in 1:total_atoms
        total_x = total_x + abs(red_atoms[i,3])
        total_y = total_y + abs(red_atoms[i,4])
        total_z = total_z + abs(red_atoms[i,5])
    end
    
    # Rotate to xy-plane
    if (total_x < 0.05)
        for i in 1:total_atoms
            new_atoms[i,1:2] = red_atoms[i,1:2]
            for j in 3:5
                new_atoms[i,3] = red_atoms[i,5]
                new_atoms[i,4] = red_atoms[i,4]
                new_atoms[i,5] = -red_atoms[i,3]
            end
        end
        return new_atoms
    elseif (total_y < 0.05)
        for i in 1:total_atoms
            new_atoms[i,1:2] = red_atoms[i,1:2]
            for j in 3:5
                new_atoms[i,3] = red_atoms[i,3]
                new_atoms[i,4] = red_atoms[i,5]
                new_atoms[i,5] = -red_atoms[i,4]
            end
        end
        return new_atoms
    elseif (total_z < 0.05)
        return red_atoms
    end

    # The in-plane vector
    norm = red_atoms[2,3]^2 + red_atoms[2,4]^2 + red_atoms[2,5]^2
    inplanevector = sqrt(1.0/norm)*[red_atoms[2,3], red_atoms[2,4], red_atoms[2,5]]
    
    # Identify the plane
    planeA = [red_atoms[2,4]  red_atoms[2,5] ; red_atoms[3,4]  red_atoms[3,5]]
    planeb = [-red_atoms[2,3], -red_atoms[3,3]]
    plane_constants = planeA \ planeb
    
    # Identify the perpendicular in plane vector
    inplaneperpA = [red_atoms[2,4] red_atoms[2,5];plane_constants[1] plane_constants[2]]
    inplaneperpb = [-red_atoms[2,3], -1]
    inplaneperp = inplaneperpA \ inplaneperpb
    norm = 1 + inplaneperp[1]^2 + inplaneperp[2]^2
    inplaneperpvector = sqrt(1.0/norm)*[1, inplaneperp[1], inplaneperp[2]]
    
    # Identify the vector perpendicular to the plane
    perpplaneA = [red_atoms[2,4] red_atoms[2,5]; inplaneperpvector[2] inplaneperpvector[3]]
    perpplaneb = [-red_atoms[2,3], -inplaneperpvector[1]]
    perpplane = perpplaneA \ perpplaneb
    norm = 1 + perpplane[1]^2 + perpplane[2]^2
    perpplanevector = sqrt(1.0/norm)*[1, perpplane[1], perpplane[2]]
    
    # Construct transformation matrix
    Rmat = [inplanevector[1] inplaneperpvector[1] perpplanevector[1]; 
    inplanevector[2] inplaneperpvector[2] perpplanevector[2]; 
    inplanevector[3] inplaneperpvector[3] perpplanevector[3]]
    for i in 1:total_atoms
        new_atoms[i,1:2] = red_atoms[i,1:2]
        new_atoms[i,3:5] = Rmat \ red_atoms[i,3:5]
    end

    return new_atoms
end

# Define function to calculate distance between atoms
function calDist(atomdata, site_1, site_2)
    pos1 = atomdata[site_1, 3:5]
    pos2 = atomdata[site_2, 3:5]
    dist = 0.0
    for i in 1:3
        dist = dist + (pos1[i] - pos2[i])^2
    end
    return sqrt(dist)
end

# Create distance matrix and identify nearest neighbor distance
function createDistMatrix(atomdata, total_atoms)
    dist = zeros(Float32, (total_atoms, total_atoms))
    mindist = 1e5*ones(Float32, total_atoms)
    Threads.@threads for i in 1:total_atoms
        for j in 1:total_atoms
            dist[i, j] = calDist(atomdata, i, j)
            if (dist[i, j] > 0 && dist[i, j] < mindist[i])
                mindist[i] = dist[i,j]
            end
        end
    end
    return mindist, dist
end

end