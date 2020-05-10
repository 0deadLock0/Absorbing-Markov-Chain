def getMatrixInverse(matrix):
	l=len(matrix)
	inverseMatrix=getIdentity(l)
	for i in range(l):
		if matrix[i][i]==0:
			k=i+1
			while k<l and matrix[k][k]==0:
				k+=1
			matrix[i],matrix[k]=matrix[k],matrix[i]
			inverseMatrix[i],inverseMatrix[k]=inverseMatrix[k],inverseMatrix[i]	
		div=matrix[i][i]
		for j in range(l):
			matrix[i][j]/=div
			inverseMatrix[i][j]/=div
		for j in range(l):
			if i==j:
				continue
			sub=matrix[j][i]/matrix[i][i]
			for k in range(l):
				matrix[j][k]-=sub*matrix[i][k]
				inverseMatrix[j][k]-=sub*inverseMatrix[i][k]
	return inverseMatrix

def multipy(matrix1,matrix2):
	n=len(matrix1[0])
	m=len(matrix1)
	o=len(matrix2[0])
	matrix=[]
	for i in range(m):
		matrix.append([])
		for j in range(o):
			sum=0
			for k in range(n):
				sum+=matrix1[i][k]*matrix2[k][j]
			matrix[i].append(sum)

	return matrix

def subtract(matrix1,matrix2):
	n=len(matrix1)
	m=len(matrix1[0])
	matrix=[]
	for i in range(n):
		matrix.append([])
		for j in range(m):
			matrix[i].append(matrix1[i][j]-matrix2[i][j])

	return matrix

def getIdentity(size):
	identity=[[0 for _ in range(size)] for _ in range(size)]
	for i in range(size):
		identity[i][i]=1

	return identity

def getN(matrixQ):#The (i, j) entry of matrix N is the expected number of times the chain is in state j, given that the chain started in state i.
	I=getIdentity(len(matrixQ))
	temp=subtract(I,matrixQ)
	N=getMatrixInverse(temp)

	return N

def getQ(transientMatrix):#Q describes the probability of transitioning from some transient state.
	n=len(transientMatrix)
	matrixQ=[]
	for i in range(n):
		matrixQ.append([])
		for j in range(n):
			matrixQ[i].append(transientMatrix[i][j])

	return matrixQ

def getR(transientMatrix):# R describes the probability of transitioning from some transient state to some absorbing state.
	n=len(transientMatrix)
	m=len(transientMatrix[0])
	matrixR=[]
	k=0
	for i in range(n):
		matrixR.append([])
		for j in range(n,m):
			matrixR[k].append(transientMatrix[i][j])
		k+=1

	return matrixR

def getTransientStates(matrix):
	states=[]
	l=len(matrix)
	for i in range(l):
		if matrix[i][i]!=1:
			states.append(matrix[i])
	return states

def getAbsorbingStates(matrix):
	states=[]
	l=len(matrix)
	for i in range(l):
		if matrix[i][i]==1:
			states.append(matrix[i])
	return states

def getTransientStatesNumber(matrix):
	l=len(matrix)
	transientStates=[]
	for i in range(l):
		if matrix[i][i]!=1:
			transientStates.append(i)
	
	return transientStates

def getAbsorbingStatesNumber(matrix):
	l=len(matrix)
	absorbingStates=[]
	for i in range(l):
		if matrix[i][i]==1:
			absorbingStates.append(i)
	
	return absorbingStates

def getCanonicalMatrix(transitionMatrix):
	'''
	[Q R]
	[0 I]
	'''
	l=len(transitionMatrix)
	states=getTransientStatesNumber(transitionMatrix)+getAbsorbingStatesNumber(transitionMatrix)
	canonicalMatrix=[]
	k=0
	for i in states:
		canonicalMatrix.append([])
		for j in states:
			canonicalMatrix[k].append(transitionMatrix[i][j])
		k+=1

	return canonicalMatrix

def getCanonicalMatrixAndNewState(transitionMatrix,state):
	canonicalMatrix=getCanonicalMatrix(transitionMatrix)
	states=getTransientStatesNumber(transitionMatrix)
	newState=states.index(state)

	return canonicalMatrix,newState

def probabilitiesOfAbsorbing(transitionMatrix,state):
	if state>len(transitionMatrix):
		raise ValueError("State doesn't exist")
	if transitionMatrix==[[1]]:
		return [1]
	state-=1
	if transitionMatrix[state][state]==1:
		return [1]
	canonicalMatrix,state=getCanonicalMatrixAndNewState(transitionMatrix,state)
	transientStates=getTransientStates(canonicalMatrix)
	if len(transientStates)==1 or len(transientStates)==len(transitionMatrix)-1:
		return [1]
	Q=getQ(transientStates)
	R=getR(transientStates)
	N=getN(Q)
	B=multipy(N,R)
	absorbingStatesNumber=getAbsorbingStatesNumber(transitionMatrix)
	probabilities=[(absorbingStatesNumber[i]+1,B[state][i])  for i in range(len(absorbingStatesNumber))]
	
	return probabilities

# Sample-
# m=[[0,0.66,0.33,0,0],[0,0,0,0.42,0.57],[0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,1]]
# s=1
# p=probabilitiesOfAbsorbing(m,s)
# print(p)