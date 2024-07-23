import sys

lines = open(sys.argv[1]).readlines()
print(len([i for i in range(len(lines[1])) if lines[1][i] == lines[7][i] and lines[3][i] == lines[5][i] and lines[1][i] != lines[3][i]]),
    len([i for i in range(len(lines[1])) if lines[1][i] == lines[5][i] and lines[3][i] == lines[7][i] and lines[1][i] != lines[3][i]]),
    len([i for i in range(len(lines[1])) if lines[1][i] == lines[3][i] and lines[5][i] == lines[7][i] and lines[1][i] != lines[5][i]]))