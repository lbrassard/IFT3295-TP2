import numpy as np
import pandas as pd

with open('sequences.fasta', 'r') as seq:
    lines = seq.readlines()

sequences = []
s = ''
for line in lines:

    if line == '\n':
        sequences.append(s)
        s = ''
    elif line[0] != '>':
        s += line.strip('\n')
sequences.append(s)

with open('BLOSUM62.txt', 'r') as scores:
    lines = scores.readlines()

kv = dict()
scoreMatrix = []

for line in lines:
    # create a key:pair value for each possible AA that are indicated in the first row
    if line[0] == ' ':
        x = line.replace(' ', '')
        x = x.replace('\n','')
        for i in range(len(x)):
            kv.update({x[i]: i})

    #Create the score matrix from the blosum file without the header and the first colum
    #Square matrix starting at 0,0 for A,A
    else:
        values = line[2:].replace(' ', ',').replace(',,', ',').strip(',\n')
        if values[0] == ',':
            values = values[1:]
        row = values.split(',')
        scoreMatrix.append(row)

# alignment for sequences

s1 = "ACACT"
s2 = "AAT"
a = 10
b = 1

#maybe look for word inversion here
V = np.ones(shape=(len(s2)  +1, len(s1) + 1)) * np.NINF
G = np.ones(shape=(len(s2) + 1, len(s1) + 1)) * np.NINF
F = np.ones(shape=(len(s2) + 1, len(s1) + 1)) * np.NINF
E = np.ones(shape=(len(s2) + 1, len(s1) + 1)) * np.NINF

#Initial conditions
G[0][0] = 0

for i in range(len(s2)+1):
    V[i][0] = -(a + b * i)
    F[i][0] = -(a + b * i)

for j in range(len(s1)+1):
    V[0][j] = -(a + b * j)
    E[0][j] = -(a + b * j)


def score(char1, char2):
    key1 = kv.get(char1)
    key2 = kv.get(char2)

    result = int(scoreMatrix[key1][key2])
    return result


for r in range(1, len(s2)+1):
    for c in range(1, len(s1)+1):

        char1 = s1[c - 1]
        char2 = s2[r - 1]
        # according to online resources

        # # update G matrix

        # m = G[int(r) - 1][int(c) - 1] + score(char1, char2)
        # ix = F[int(r) - 1][int(c) - 1] + score(char1, char2)
        # iy = E[int(r) - 1][int(c) - 1] + score(char1, char2)
        #
        # G[r][c] = max(m, ix, iy)
        #
        # # update Ix Matrix
        # m = G[int(r) - 1][int(j)] - a - b  # open a new gap in x
        # ix = F[int(r) - 1][int(j)] - b    # extend gap in x
        #
        # F[r][c] = max(m, ix)
        #
        # #update E Matrix
        # m = G[int(r)][int(j) - 1] - a - b  # open a new gap in y
        # iy = E[int(r)][int(j) - 1] - b    # extend gap in y
        #
        # E[r][c] = max(m, iy)

        # according to the slides

        # update g matrix
        G[r][c] = V[r - 1][c - 1] + score(char1, char2)

        # update E matrix
        v1 = E[r][c-1]
        v2 = V[r][c-1] - a
        E[r][c] = max(v1, v2) - b

        # update F matrix
        v1 = F[r-1][c]
        v2 = V[r-1][c] - a
        F[r][c] = max(v1, v2) - b

        # update V matrix
        V[r][c] = max(G[r][c], E[r][c], F[r][c])


print(V)
print(G)
print(E)
print(F)
