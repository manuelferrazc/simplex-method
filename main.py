import util
import numpy as np
import argparse
from enum import Enum

class Status(Enum):
    FREE = 0
    GEZ = 1
    LEZ = 2

def zero(x,e=1e-8):
    return 0.0 if abs(x)<e else x

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
    coef = input('Write the coefficients of the variables in the objective function\n' if askinput else '').split(' ')
    coef = [float(i) for i in coef]
    c = []
    for i in range(numVar):
        if varMap[i][0]==Status.GEZ: c.append(coef[i])
        elif varMap[i][0]==Status.LEZ: c.append(-coef[i])
        else:
            c.append(coef[i])
            c.append(-coef[i])
    return c

def getMatrixFPI(numVar,numRet,varMap,numVarFPI,askinput=False):
    A = [[]]*numRet
    for i in range(numRet): A[i] = [0]*(numVarFPI+1)
    needsExtra = []
    c = objectiveFunction(numVar,varMap,askinput)

    for i in range(numRet):
        inequation = input(f'Write the inequation {i+1}\n' if askinput else '').split(' ')
        v=0
        if inequation[numVar]!='==':
            needsExtra.append(i)
        for j in range(numVar):
            if varMap[j][0]==Status.FREE:
                A[i][v] = float(inequation[j])
                v+=1
                A[i][v] = -float(inequation[j])
            elif varMap[j][0]==Status.LEZ:
                A[i][v] = -float(inequation[j])
            else: A[i][v] = float(inequation[j])
            v+=1

        A[i][numVarFPI] = float(inequation[numVar+1])
        if inequation[numVar]=='>=': A[i] = [-x for x in A[i]]

    



def main():
    askinput = True
    numVar,numRet,varMap,numVarFPI = getConstraints(askinput)
    getMatrixFPI(numVar,numRet,varMap,numVarFPI,askinput)

if __name__ == '__main__': main()