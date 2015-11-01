import argparse
from os import listdir
from os.path import isfile, join

def parse_mario_file(f):
    data = dict()
    data["gen"] = int(f.readline())
    data["max"] = float(f.readline())
    species = int(f.readline())
    totalfit = 0
    totalgenomes = 0
    for j in range(species):
        f.readline()
        f.readline()
        genomes = int(f.readline())
        totalgenomes += genomes
        for i in range(genomes):
            totalfit += float(f.readline())
            temp = f.readline()
            while not("done" in temp):
                temp = f.readline()
            genes = int(f.readline())
            for k in range(genes):
                f.readline()
    data["avg"] = totalfit / genomes 
    return data

def main():
    parser = argparse.ArgumentParser(description="Parses backup and .pool files\
        from MarI/O into .data files")
    parser.add_argument("dir", type=str, help="Directory in which the\
        data is stored")
    parser.add_argument("testname", type=str, help="the suffix of the\
        files to parse")
    args = parser.parse_args()
    files = [f for f in listdir(args.dir) if isfile(join(args.dir, f))]
    files = [f for f in files if f.endswith(args.testname) and f.startswith("backup")]
    for file_ in files:
        print file_
        f = open(join(args.dir, file_))
        data = parse_mario_file(f)
        print data
        f = open(\
            join(args.dir, "{}.{}.data".format(args.testname, str(data["gen"]))),\
            'w')
        f.write("{}\n{}\n{}".format(data["gen"], data["max"], data["avg"]))

if __name__ == "__main__":
    main()
