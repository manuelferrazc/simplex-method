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

def printVarMap(varMap):
    print('Each variable Xi in the original problem is mapped to another variable Yj')
    print('If Xi is non-negative, then Xi = Yj, for some j')
    print('If Xi is non-positive, then Xi = -Yj, for some j')
    print('Finally, if Xi is free, then Xi = Yj - Yk, for some j and k')
    print('Yj is non-negative, for all j')
    
    for i in range(len(varMap)):
        if varMap[i][0]==Status.GEZ: print(f'X{i} = Y{varMap[i][1]}')
        elif varMap[i][1]==Status.LEZ: print(f'X{i} = -Y{varMap[i][1]}')
        else: print(f'X{i} = Y{varMap[i][1][0]} - Y{varMap[i][1][1]}')

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

    printVarMap(var_map)
    return (numvar,numret,var_map,m)

def printMatrixFPI(A,c,args):
    print('\nProblem given in the input after being placed in the SLE:')
    print('max ',end='')
    for i in range(len(c)): print(f"{'%+*.*f*' if i>0 else '%*.*f*'}Y{i} " % (args.digits, args.decimals, c[i]) ,end='')
    print('\n\nSubject to: ')
    for line in range(0,len(A)):
        for col in range(0,len(A[0])-1):
            print(f"{'%+*.*f*' if col>0 else '%*.*f*'}Y{col} " % (args.digits, args.decimals, abs(A[line][col])),end='')
        print(f'= %*.*f' % (args.digits,args.decimals,A[line][col-1]))

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
    negs = set()
    for i in range(numRet):
        A[i].append(b[i])
        if b[i]<0:
            negs.add(i)
            A[i] = [-x for x in A[i]]

    return util.makeMatrixFullRank(np.array(A))[0],negs

def parseAndGetInput():
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', help='Name of the file to process')
    parser.add_argument('--decimals', type=int, default=4,
                        help='Number of decimal digits to print numeric values')
    parser.add_argument('--digits', type=int, default=7,
                        help='Total number of digits to print numeric values')
    parser.add_argument('--policy',
                        type=str,
                        choices=['largest','smallest','bland'],
                        default='largest',
                        help='policy used to choose the pivot element')

    args = parser.parse_args()
    if args.decimals < 0:
        print('The number of decimal digits to print numeric values must be, at least, 0')
        exit(0)
    if args.digits < 5:
        print('The total number of digits to print numeric values must be greater than 4')
        exit(0)
    file = open(args.filename)
    numVar,numRet,varMap,numVarFPI = getConstraints(file)
    c,isMin = getObjectiveFunction(numVar,varMap,file)
    A,negs = getMatrixFPI(numVar,numRet,varMap,numVarFPI,file)
    file.close()
    return varMap,c,isMin,A,args,negs

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
        if u==1 and z+u==n:
            if i+1 not in lines:
                lines.add(i+1)
                baseMap[column+len(A)] = i+1
    return lines,baseMap

def pivot(A,i,j):
    A[i] = A[i]/A[i][j]

    for l in range(i):
        A[l] = A[l] - A[i]*A[l][j]
    for l in range(i+1,len(A)):
        A[l] = A[l] - A[i]*A[l][j]

def getTableau(A,c,lines,baseMap,negs):
    numrows,numcolumns = len(A), len(A[0]) 
    Al = np.zeros((1+numrows,numcolumns+2*numrows-len(lines)))
    for i in range(numcolumns-1): Al[0][numrows+i] = -c[i]
    for i in range(numrows):
        Al[i+1][i] = -1 if i in negs else 1
        Al[i+1][len(Al[0])-1] = A[i][numcolumns-1]
        for j in range(numcolumns-1):
            Al[i+1][j+numrows] = A[i][j]
    c=0
    
    for i in range(1,numrows+1):
        if i in lines: continue
        lines.add(i)
        baseMap[numrows+numcolumns+c-1] = i
        Al[i][numrows+numcolumns+c-1] = 1
        Al[0][numrows+numcolumns+c-1] = 1
        c+=1
    if c!=0:
        for i in range(numrows+numcolumns-1):
            Al[0][i] = 0
    for i in baseMap:
        if zero(Al[0][i])!=0:
            Al[0] = Al[0]-Al[baseMap[i]]*Al[0][i]
    return Al,c

def changeTableau(A,c,extra):
    numrows,numcolumns = len(A), len(A[0])
    Al = np.zeros((numrows,numcolumns-extra))
    for i in range(1,numrows):
        for j in range(numcolumns-extra-1): Al[i][j] = A[i][j]
        Al[i][numcolumns-extra-1] = A[i][numcolumns-1]
    
    for i in range(len(c)): Al[0][numrows+i-1] = -c[i]
    return Al

def printTableau(A,args):
    numrows,numcolumns = len(A), len(A[0])
    s = 6
    for i in range(numrows-1):
        print('%*.*f ' % (args.digits, args.decimals, A[0][i]), end='')
        s+=args.digits+1
    print('|| ', end='')
    for i in range(numrows-1,numcolumns-1):
        print('%*.*f ' % (args.digits, args.decimals, A[0][i]), end='')
        s+=args.digits+1
    print('|| %*.*f\n'% (args.digits,args.decimals,A[0][numcolumns-1]), end='')
    s+=args.digits
    print('-'*s)

    for i in range(1,numrows):
        for j in range(numrows-1):
            print('%*.*f ' % (args.digits, args.decimals,A[i][j]), end='')
        print('|| ', end='')
        for j in range(numrows-1,numcolumns-1):
            print('%*.*f ' % (args.digits, args.decimals,A[i][j]), end='')
        print('|| ', end='')
        print('%*.*f ' % (args.digits, args.decimals,A[i][numcolumns-1]))

def bland(A,cols): 
    c = []
    numlines,numcolumns = len(A),len(A[0])
    for i in range(numlines-1,numcolumns-1):
        A[0][i] = zero(A[0][i])
        if A[0][i]<0:
            m = [-1,-1]
            for j in range(1,len(A)):
                A[j][i] = zero(A[j][i])
                if A[j][i]<=0: continue
                if m[0]==-1 or A[j][len(A[0])-1]/A[j][i]<m[0]: m = [A[j][len(A[0])-1]/A[j][i],j]
            if m[0]== -1:
                print('Status: ilimitado')
                exit(0)
            else: return m[1],i,False
        elif A[0][i]==0 and i not in cols: c.append(i)
    for i in c:
        m = [-1,-1]
        for j in range(1,len(A)):
            if A[j][i]<=0: continue
            if m[0]==-1 or A[j][len(A[0])-1]/A[j][i]<m[0]: m = [A[j][len(A[0])-1]/A[j][i],j]
        if m[0]==-1: continue
        else: return m[1],i,True
    return -1,-1,True

def largest(A,cols):
    c = []
    numlines,numcolumns = len(A),len(A[0])
    b = -1
    for col in range(numlines-1,numcolumns-1):
        A[0][col] = zero(A[0][col])
        if A[0][col]>=0:
            if A[0][col]==0 and col not in cols: c.append(col)
            continue
        if b==-1 or A[0][b]>A[0][col]: b = col
    if b!=-1:
        r = -1
        for line in range(1,len(A)):
            if A[line][b]<=0: continue
            elif r==-1 or A[line][numcolumns-1]/A[line][b]<A[r][numcolumns-1]/A[r][b]:
                r = line
        if r==-1:
            print('Status: ilimitado')
            exit(0)
        else: return r,b,False
    for col in c:
        m = [-1,-1]
        for line in range(1,numlines):
            if A[line][col]<=0: continue
            if m[0]==-1 or A[line][numcolumns-1]/A[line][col]<m[0]: m = [A[line][numcolumns-1]/A[line][col],line]
        if m[0]==-1: continue
        else: return m[1],col,True
    return -1,-1,True

def smallest(A,cols):
    c = []
    numlines,numcolumns = len(A),len(A[0])
    b = -1
    for col in range(numlines-1,numcolumns-1):
        A[0][col] = zero(A[0][col])
        if A[0][col]>=0:
            if A[0][col]==0 and col not in cols: c.append(col)
            continue
        if b==-1 or A[0][b]<A[0][col]: b = col
    if b!=-1:
        r = -1
        for line in range(1,len(A)):
            if A[line][b]<=0: continue
            elif r==-1 or A[line][numcolumns-1]/A[line][b]<A[r][numcolumns-1]/A[r][b]:
                r = line
        if r==-1:
            print('Status: ilimitado')
            exit(0)
        else: return r,b,False
    for col in c:
        m = [-1,-1]
        for line in range(1,numlines):
            if A[line][col]<=0: continue
            if m[0]==-1 or A[line][numcolumns-1]/A[line][col]<m[0]: m = [A[line][numcolumns-1]/A[line][col],line]
        if m[0]==-1: continue
        else: return m[1],col,True
    return -1,-1,True

def select(A,cols,policy):
    if policy=='bland': return bland(A,cols)
    elif policy=='largest': return largest(A,cols)
    else: return smallest(A,cols)

def getSolution(A,cols):
    x = np.zeros(len(A[0])-len(A))
    for col in cols: x[col-len(A)+1] = zero(A[cols[col]][len(A[0])-1])
    return x

def getSolutionOriginal(varMap,xa):
    x = np.zeros(len(varMap))
    for i in range(len(varMap)):
        if varMap[i][0]==Status.FREE: x[i] = xa[varMap[i][1][0]]-xa[varMap[i][1][1]]
        elif varMap[i][0]==Status.GEZ: x[i] = xa[varMap[i][1]]
        else: x[i] = -xa[varMap[i][1]]
    return x

def getSolutionDual(A):
    x = np.zeros(len(A)-1)
    for i in range(len(A)-1):
        x[i] = zero(A[0][i])
    return x

def printArray(v,args):
    print('[',end='')
    for i in range(len(v)):
        print(f"{'    ' if i>0 else ''}%*.*f" % (args.digits,args.decimals, v[i]),end='')
    print(']')

def printSolution(A,varMap,base,args):
    if not hasattr(printSolution, "counter"):
        printSolution.counter = 0
    printSolution.counter += 1
    print(f'\n\nIteration {printSolution.counter}\nTableau:')
    printTableau(A,args)

    print('\nBase:', str([x for x in base])[1:-1])
    x = getSolution(A,base)
    print('Y (primal) = ',end='')
    printArray(x,args)
    x = getSolutionOriginal(varMap,x)
    print('X (originals variables) = ',end='')
    printArray(x,args)
    print('Dual\'s Solution:',end='')
    printArray(getSolutionDual(A),args)


def find(A,i,cols):
    for col in range(len(A[0])):
        if zero(A[i][col]-1)==0 and col in cols: return col

def selectAux(A,cols,extra,policy):
    c = []
    numlines,numcolumns = len(A),len(A[0])
    for col in range(numlines-1,numcolumns-1):
        A[0][col] = zero(A[0][col])
        if A[0][col]>=0:
            if A[0][col]==0 and col-numlines+1 not in cols: c.append(col)
            continue
        else:
            b = -1,[]
            for line in range(1,len(A)):
                if A[line][col]<=0: continue
                if b[0]==-1 or A[line][len(A[0])-1]/A[line][col]<b[0]:
                    b = A[line][len(A[0])-1]/A[line][col],[line]
                elif abs(A[line][len(A[0])-1]/A[line][col]-b[0])<0.00000001:
                    b[1].append(line)
            for line in b[1]:
                if find(A,line,cols)>=extra:
                    return line,col,False
    return select(A,cols,policy)
    
def findBase(A,varMap,c,cols,args,extra):
    print('Finding a new (and valid) base without extra variables (if it exists):')

    x,y,z = selectAux(A,cols,len(A[0])-1-extra,args.policy)
    e = extra
    while not z:
        d = find(A,x,cols)
        del cols[d]
        pivot(A,x,y)
        cols[y]=x
        x,y,z = selectAux(A,cols,len(A[0])-1-extra,args.policy)
    A[0][len(A[0])-1] = zero(A[0][len(A[0])-1])
    if A[0][len(A[0])-1]!=0:
        print('Status: inviavel')
        exit(0)
    
    print('\nNew valid base found:',str([x for x in cols])[1:-1])
    print('Tableau: ')
    printTableau(A,args)

    Al = changeTableau(A,c,extra)
    print('\nChanging the tableau to the original objective function:')
    printTableau(Al,args)
    b=True
    for i in cols:
        if zero(Al[0][i])!=0:
            if b:
                b=False
                print('\nFixing the tableau to make it\'s base:')
            Al[0] = Al[0]-Al[cols[i]]*Al[0][i]
    if not b:
        printTableau(Al,args)
        print()
    return Al

def simplex(A,varMap,c,isMin,args,extra,baseMap):
    if extra!=0: A = findBase(A,varMap,c,baseMap,args,extra)

    x,y,z = select(A,baseMap,args.policy)
    while not z:
        del baseMap[find(A,x,baseMap)]
        pivot(A,x,y)
        baseMap[y]=x
        printSolution(A,varMap,baseMap,args)
        print()
        x,y,z = select(A,baseMap,args.policy)

    print('\nOptimal value reached! Tableau: ')
    printTableau(A,args)
    print('\nOptimal value: %*.*f' % (args.digits,args.decimals,-A[0][len(A[0])-1] if isMin else A[0][len(A[0])-1]))
    print('Solutions:')
    sol = getSolution(A,baseMap)
    print('X = ',end='')
    printArray(getSolutionOriginal(varMap,sol),args)
    print('Dual\'s solution: ',end='')
    printArray(getSolutionDual(A),args,)

    if x!=-1:
        print('\n\nThis instance has infinite solutions! Another one:')
        del baseMap[find(A,x,baseMap)]
        baseMap[y]=x
        pivot(A,x,y)
        printTableau(A,args)
        sol = getSolution(A,baseMap)
        print('\nOptimal value: %*.*f (same)' % (args.digits,args.decimals,-A[0][len(A[0])-1] if isMin else A[0][len(A[0])-1]))
        print('Another solution to X:\nX = ',end='')
        printArray(getSolutionOriginal(varMap,sol),args)


def main():
    varMap,c,isMin,A,args,negs = parseAndGetInput()
    c = c + [0]*(len(A[0])-len(c)-1)
    printMatrixFPI(A,c,args)
    lines,baseMap = hasBase(A)
    A,ex = getTableau(A,c,lines,baseMap,negs)
    print('\nCreating the tableau for the matrix A', ':' if ex==0 else f' (with {ex} extra columns to create a base):')
    printTableau(A,args)
    print('\nColumns in the base:', str([x for x in baseMap])[1:-1])
    simplex(A,varMap,c,isMin,args,ex,baseMap)

if __name__ == '__main__': main()