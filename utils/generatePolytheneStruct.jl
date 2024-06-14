a = 1.35
n = parse(Int64, ARGS[1])

open("polythene.dat", "w") do file
    write(file, string(n))
    write(file, "\n")
    write(file, "Polythene\n")
    for i in 1:n
        write(file, "6 1 0 "*string((i-1)*a)*" 0\n")
    end
end
