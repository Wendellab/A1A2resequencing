import sys

populations = {}
with open(sys.argv[2]) as population_file:
    for line in population_file:
        data = line.strip('\n').split('\t')
        populations[data[0]] = data[1]

with open(sys.argv[1]) as structure_file:
    header = structure_file.readline().strip("\n")
    distance = structure_file.readline().strip("\n")

    print(header)
    print(distance)
    for line in structure_file:
        data = line.strip("\n").split(" ")
        if data[0] in populations:
            data[1] = populations[data[0]]
        else:
            sys.stderr.write("ERROR: no population information available.\t")
            sys.stderr.write(data[0]+"\n")
        print(" ".join(data))
