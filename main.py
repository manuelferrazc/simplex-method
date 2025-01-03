import util
import numpy as np
import argparse
from enum import Enum

class Status(Enum):
    FREE = 0
    GEZ = 1
    LEZ = 2

def zero(x,e=1e-8):
    return 0 if abs(x)<e else x

def printSolution():
    pass

def getConstraints(file):
    numvar = int(file.readline())
    numret = int(file.readline())
    var_map = [(0,None)]*numvar
    m = 0
    rest_var = file.readline().split(' ')
    rest_var = [int(i) for i in rest_var]
    for i in range(numvar):
        if rest_var[i]==0:
            var_map[i] = (Status.FREE,(m,m+1))
            m+=1
        elif rest_var[i]==1:var_map[i] = (Status.GEZ,m)
        else: var_map[i] = (Status.LEZ,m)
        m+=1
    return (numvar,numret,var_map,m)

def getObjectiveFunction(numVar,varMap,file):
    coef = file.readline().split(' ')
    c = []
    if coef[0][0]=='m':
        for i in range(numVar):
            f = float(coef[i+1])
            if varMap[i][0]==Status.GEZ: c.append(f)
            elif varMap[i][0]==Status.LEZ: c.append(-f)
            else:
                c.append(f)
                c.append(-f)
        if coef[0]=='min':
            c = [-x for x in c]
            return c,True
    else:
        for i in range(numVar):
            f = float(coef[i])
            if varMap[i][0]==Status.GEZ: c.append(f)
            elif varMap[i][0]==Status.LEZ: c.append(-f)
            else:
                c.append(f)
                c.append(-f)
    return c,False

def getMatrixFPI(numVar,numRet,varMap,numVarFPI,file):
    A = [[]]*numRet
    for i in range(numRet): A[i] = [0]*(numVarFPI)
    needsExtra = set()
    b = [0]*numRet
    for i in range(numRet):
        inequation = file.readline().split(' ')
        v=0
        if inequation[numVar]!='==':
            needsExtra.add(i)
        for j in range(numVar):
            if varMap[j][0]==Status.FREE:
                A[i][v] = float(inequation[j])
                v+=1
                A[i][v] = -float(inequation[j])
            elif varMap[j][0]==Status.LEZ:
                A[i][v] = -float(inequation[j])
            else: A[i][v] = float(inequation[j])
            v+=1

        b[i] = float(inequation[numVar+1])
        if inequation[numVar]=='>=':
            A[i] = [-x for x in A[i]]
            b[i] = -b[i]

    r = numVarFPI
    for i in range(numRet):
        for j in range(numVarFPI,numVarFPI+len(needsExtra)):
            A[i].append(0)
        if i in needsExtra:
            A[i][r]=1
            r+=1
    for i in range(numRet): A[i].append(b[i])
    return util.makeMatrixFullRank(np.array(A))[0]

def parseAndGetInput():
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', help='Name of the file to process')
    parser.add_argument('--decimals', type=int, default=4, help='Number of decimal digits to print numeric values')
    parser.add_argument('--digits', type=int, default=7, help='Total number of digits to print numeric values')
    parser.add_argument('--policy', type=str, default='largest', help='policy used to choose the pivot element')

    args = parser.parse_args()
    file = open(args.filename)
    numVar,numRet,varMap,numVarFPI = getConstraints(file)
    c,isMin = getObjectiveFunction(numVar,varMap,file)
    A = getMatrixFPI(numVar,numRet,varMap,numVarFPI,file)
    file.close()
    return numVar,numRet,varMap,numVarFPI,c,isMin,A,args

def hasBase(A):
    lines = set()
    baseMap = {}
    m = len(A[0])
    n = len(A)
    for column in range(m):
        z=0
        u=0
        i=0
        for line in range(n):
            A[line][column] = zero(A[line][column])
            if A[line][column]==0: z+=1
            elif zero(A[line][column]-1.0)==0 and u==0: u,i = u+1,line
            else: break
        if u==1 and z+u==n: # basic column
            if i not in lines:
                lines.add(i)
                baseMap[column] = i # maybe invert this ?? IDK

    return lines,baseMap

# def getTableau(): pass

def main():
    numVar,numRet,varMap,numVarFPI,c,isMin,A,args = parseAndGetInput()
    c = c + [0]*(len(A[0])-len(c)-1)

    x = hasBase(A)
    # if len(x[0])==len(A), then we already have a base

if __name__ == '__main__': main()