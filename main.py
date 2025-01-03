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

def getConstraints(askinput=False):
    numvar = int(input('Write the number of variables\n' if askinput else ''))
    numret = int(input('Write the number of constraints\n' if askinput else ''))
    var_map = [(0,None)]*numvar
    m = 0
    rest_var = input('Write the constraints on the variables (0:free, 1:>=0, -1:<=0)\n' if askinput else '').split(' ')
    rest_var = [int(i) for i in rest_var]
    for i in range(numvar):
        if rest_var[i]==0:
            var_map[i] = (Status.FREE,(m,m+1))
            m+=1
        elif rest_var[i]==1:var_map[i] = (Status.GEZ,m)
        else: var_map[i] = (Status.LEZ,m)
        m+=1
    return (numvar,numret,var_map,m)

def objectiveFunction(numVar,varMap,askinput=False):
    coef = input('Write the coefficients of the variables in the objective function\n\
If neither max or min expression is given, it\'s assumed that it is a maximization problem\n' if askinput else '').split(' ')
    # coef = [float(i) for i in coef]
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

def getMatrixFPI(numVar,numRet,varMap,numVarFPI,askinput=False):
    A = [[]]*numRet
    for i in range(numRet): A[i] = [0]*(numVarFPI)
    needsExtra = set()
    c = objectiveFunction(numVar,varMap,askinput)
    b = [0]*numRet
    for i in range(numRet):
        inequation = input(f'Write the inequation {i+1}\n' if askinput else '').split(' ')
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
    print(np.array(A))
    # next steps: remove dependences (use util) and return it
    # the next and final transformation will be make by another function i think



def main():
    askinput = True
    numVar,numRet,varMap,numVarFPI = getConstraints(askinput)
    getMatrixFPI(numVar,numRet,varMap,numVarFPI,askinput)

if __name__ == '__main__': main()