#!/bin/python3

'''
This is a python program that takes in the scores array at each iteration
of a run of pagerank. It then determines at what iteration each node 
reached a stable state.
'''

import sys

class iterlist:

    def __init__(self):
        self.lst = []   # list of list of (degree, termination iterator)

    def output(self, numnodes=10):
        self.lst = sorted(self.lst, key=lambda x : x[1])
        for i in range(min(numnodes, len(self.lst))):
            print("degree, term iter:", self.lst[i][0], ":", self.lst[i][1])

    def updatetermiter(self, d, j):
        self.lst[d][1] = j
        return(self.lst[d])

class scoresmatrix:

    def __init__(self):
        self.data = []
        self.il = iterlist()

    def populateiterlist(self, line):
        self.il.lst = []
        for degree in (line[1:]).split():
            self.il.lst.append([int(degree), 0])
        

    def computeconvergence(self, tolerance=0.00001):
        result = [0] * self.numnodes
        newres = []
        fs = self.finalscores
        for i in range(0, self.numnodes):
            total = 0
            for j in range(0, self.numiters):
                if abs(float(self.data[j][i]) - float(fs[i])) < tolerance:
                    newres.append((self.il).updatetermiter(i, j+1))
                    result[i] = j+1
                    break
        return sorted(newres, key=lambda x : -x[1])


    @property
    def numnodes(self):
        return len(self.data[0])

    @property
    def finalscores(self):
        return self.data[-1]

    @property
    def numiters(self):
        return len(self.data)

    

def main():
    assert len(sys.argv) >= 2
    sm = scoresmatrix()
    outnum = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    infile = open(sys.argv[1], "r")
    first = True
    for line in infile.readlines():
        if (line[0] == '$'):
            sm.populateiterlist(line)
        elif (line[0] != '>'):
            sm.data.append(line.split())

    print("numnodes: ", sm.numnodes)    
    print("numiters: ", sm.numiters)
    #print(sm.il.output())
    results = sm.computeconvergence()
    for r in range(min(len(results), outnum)):
        print(results[r])






if __name__ == '__main__':
    main()
