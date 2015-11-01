import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join
import argparse

def main():
    parser = argparse.ArgumentParser(description = "Plot data from MarI\O")
    parser.add_argument("dir", type=str, help="Directory in which the\
        data is stored")
    parser.add_argument("testname", type=str, help="the suffix of the\
        files to parse")
    args = parser.parse_args()
    files = [f for f in listdir(args.dir) if isfile(join(args.dir, f))]
    files = [f for f in files if args.testname in f and "data" in f]
    for file_ in files:
        f = open(args.dir + file_)
        g = int(f.readline())
        m = float(f.readline())
        a = float(f.readline())
        data[g] = dict(avgfit=a, maxfit=m)
    gens = []
    avgs = []
    maxes = []
    for gen, fit in data.items():
        gens.append(gen)
        avgs.append(fit["avgfit"])
        maxes.append(fit["maxfit"])
    plt.plot(gens, avgs, "b", gens, maxes, "r")
    plt.show()

if __name__ == "__main__":
    main()
