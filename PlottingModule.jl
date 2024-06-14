module PlottingModule
using PyPlot, PerceptualColourMaps

function plotEnergyLevels(energies, norbitals, sysname)
    num_energies = length(energies)
    ax = PyPlot.subplot(111)
    for i in 1:norbitals-1
        ax.axhline(energies[i]*27.2,color="blue")
    end
    ax.axhline(energies[norbitals]*27.2,color="blue",label="Occupied")
    for i in norbitals+1:num_energies-1
        ax.axhline(energies[i]*27.2,color="red")
    end
    ax.axhline(energies[num_energies]*27.2,color="red",label="Unoccuped")

    ax.spines["right"].set_visible(false)
    ax.spines["top"].set_visible(false)
    ax.spines["bottom"].set_visible(false)
    ax.get_xaxis().set_visible(false)
    ax.set_ylabel("Orbital Energies (eV)")
    ax.legend(loc="best")
    ax.set_title("Energy levels of "*sysname)
    PyPlot.savefig(sysname*"_energies.pdf")
    PyPlot.show()
end

function plotOrbitalOccupation(total_atoms, mindist, dist_mat, atomdata, conv_vector, plot_title, filename, show)
    # Get x-y coordinates for points
    xlist = atomdata[:,3]
    ylist = atomdata[:,4]

    # Determine the x and y limits for the plots
    maxnndist = maximum(mindist) 
    xmin = minimum(xlist) - maxnndist/4
    xmax = maximum(xlist) + maxnndist/4
    ymin = minimum(ylist) - maxnndist/4
    ymax = maximum(ylist) + maxnndist/4

    # Square the coefficients to get occupation probabilities
    probabilities = zeros(Float64, length(conv_vector))
    sizes = 1000*ones(Int64, length(conv_vector))
    for i in eachindex(conv_vector)
        probabilities[i] = conv_vector[i]^2
    end
    
    # Plot the contours
    nc = 500
    map = cmap("L20")
    tricontourf(xlist, ylist, probabilities, nc, alpha=0.9, antialiased=true, cmap=ColorMap(map))
    colorbar()
    
    # Plot the atoms
    scatter(xlist, ylist, sizes, color="black")
    for i in 1:total_atoms
        text(xlist[i], ylist[i], Int64(atomdata[i,1]), ha="center", va="center", fontsize="18", color="white")
    end
   
    # Plot the near-neighbor bonds
    for i in 1:total_atoms
        for j in 1:total_atoms
            if (i != j && abs(dist_mat[i,j] - mindist[i]) < 0.5)
                plot([atomdata[i,3], atomdata[j,3]], 
                [atomdata[i,4], atomdata[j,4]], color="black")
            end
        end
    end

    # Set title for the plot
    PyPlot.title(plot_title, fontsize="20")

    # Specify limits for the axis and hide the axis
    xlim([xmin, xmax])
    ylim([ymin, ymax])
    PyPlot.axis("off")

    # Save the plot
    PyPlot.savefig(filename)

    # Show the plot
    if (show)
        PyPlot.show()
    end
end

end