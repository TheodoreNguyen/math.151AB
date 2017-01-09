import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt

#part A

def argand(w, m):		# takes in an array of eigenvalues w and plots on real/imag plane
	import matplotlib.pyplot as plt
	import numpy as np
	for x in range(len(w)):
		plt.plot([0,w[x].real],[0,w[x].imag],'ro',label='python', color='g')
	limit=np.max(np.ceil(np.absolute(w))) # set limits for axis
	plt.xlim((-limit,limit))
	plt.ylim((-limit,limit))
	plt.ylabel('Imaginary')
	plt.xlabel('Real')
	plt.title('Eigenvalues of 100 m x m matrices with m = %i'%(m))
	plt.show()

def randMatrixEigPlot(m, max):
	while m <= max:	# for square matrices of size m, up to 256
		matrixArr = []
		eigArr = []
		for i in range(100):	# doing this for 100 of m-sized matrices
			matrixArr.append(pow(m, -0.5) * np.random.randn(m, m) + 0)#generate a random mxm matrix
			w, v = LA.eig(matrixArr[i])
			for i in range(len(w)):
				eigArr.append(w[i])
		argand(eigArr, m)
		m = m * 2
		
	

def spectRadPlot(m, max):	#spectral radius rho(A) is the largest magnitude eigenvalue of A
	spectArr = []
	dimArr = []
	while m <= max:
		matrixArr = []
		#preSpectArr = []
		for i in range(100):	# doing this for 100 of m-sized matrices
			matrixArr.append(pow(m, -0.5) * np.random.randn(m, m) + 0)#generate a random mxm matrix
			w, v = LA.eig(matrixArr[i])
			spectArr.append(np.amax(np.absolute(w)))
			#preSpectArr.append(np.amax(np.absolute(w)))
			dimArr.append(m)
		#spectArr.append(np.amax(np.absolute(preSpectArr)))
		#dimArr.append(m)
		m = m * 2
	plt.plot(dimArr, spectArr, marker='o', linestyle='None', color='r')
	limit=np.max(np.ceil(np.absolute(spectArr)))
	llimit=np.min(np.floor(np.absolute(spectArr)))
	plt.xlim((0, max))
	plt.ylim((llimit, limit))
	plt.xlabel('m - Matrix Dimension')
	plt.ylabel('Spectral Radius')
	plt.title('Convergence of spectral radius as m approaches infinity')
	plt.show()

#----------------------------------------------------------------------------------------------------

#part B
	

#spectral radius(A) = largest absolute value of an eigenvalue of A
#2-norm(A) = induced spectral norm when p = 2 (Euclidean norm) = spectral norm(A) = 
# 		= SQRT(maximum_eigenvalue(A_conj_trans * A)) = maximum_singularvalue(A)			
def twoNormSpecTwo(m, max):
	spectArr = []
	dimArr = []
	singArr = []
	while m <= max:
		matrixArr = []
		for i in range(100):
			matrixArr.append(pow(m, -0.5) * np.random.randn(m, m) + 0)
			w, v = LA.eig(matrixArr[i])
			spectArr.append(np.amax(np.absolute(w)))
			U, s, V = LA.svd(matrixArr[i], full_matrices=True)
			singArr.append(np.amax(s))
			dimArr.append(m)
		m = m * 2
	#plot induced 2-norm against matrix size
	plt.plot(dimArr, singArr, marker='o', linestyle='None', color='c')
	plt.xlim((0, max))
	plt.ylim((np.amin(singArr), np.amax(singArr)))
	plt.xlabel('m - Matrix Dimension')
	plt.ylabel('2-norm')
	plt.title('Convergence of 2-norm as m approaches infinity')
	plt.show()
	#plot the ratio of induced 2-norm/spectral radius against matrix size
	# this is b/c spectral radius <= (2)-norm; ie, rho(A) <= ||A||
	ratioArr = []
	for i in range(len(dimArr)):
		ratioArr.append(singArr[i]/spectArr[i])
	plt.plot(dimArr, ratioArr, marker='o', linestyle='None', color='c')
	plt.xlim((0, max))
	plt.ylim((np.amin(ratioArr), np.amax(ratioArr)))
	plt.xlabel('m - matrix dimension')
	plt.ylabel('Ratio of matrix 2-norm to matrix spectral radius')
	plt.title('Spectral radius approaches approx 2 times 2-norm for increasing m')
	plt.show()
	
#--------------------------------------------------------------------------------
# part C
	
#the condition number of a matrix = the smallest singular value of the matrix
#in other words, cond_num(A) == sigma_min(A)	
#find proportion of matrices with sigma_min <= 2^-1, 4^-1,... 
#	plot probability distribution of smallest singular values


def numless(list, alpha):
	count = 0
	for value in list:
		if value <= alpha:
			count = count + 1
	return count

def condNums(m, max, alpha, alphamax):
	struct = {}
	dimArr = []
	condArr = []
	M = m
	while m <= max:
		matrixArr = []
		struct[m] = []
		for i in range(100):
			matrixArr.append(pow(m, -0.5) * np.random.randn(m, m) + 0)
			U, s, V = LA.svd(matrixArr[i], full_matrices=True)
			struct[m].append(np.amin(s))
			condArr.append(np.amin(s))
			dimArr.append(m)
		m = m * 2
	xarr = []
	yarr = []
	m = M
	while m <= max:
		xarr.append(m)
		yarr.append(np.mean(struct[m]))
		m = m * 2
	plt.plot(xarr, yarr, marker='x', linestyle=":", color='m')
	plt.xlim((0, max))
	plt.ylim((np.amin(yarr), np.amax(yarr)))
	plt.xlabel('m - matrix dimension')
	plt.ylabel('Average minimum singular value')
	plt.title('Correlation b/w conditional number and m')
	plt.show()
	'''return struct'''
	
	m = M
	ALPHA = alpha
	while m <= max:
		xarr = []
		yarr = []
		alpha = ALPHA
		while alpha <= alphamax:
			xarr.append(alpha)
			yarr.append(numless(struct[m], alpha))
			alpha = alpha + 0.02
		
		plt.plot(xarr, yarr, marker='x', linestyle=':', color='r')
		plt.xlim((0, alphamax))
		plt.ylim((0, 100))
		plt.xlabel('Threshold of singular value')
		plt.ylabel('Count of matrices w/ singular value less than threshold out of 100')
		plt.title('Proportion of Min singular-value < x for m = %i'%(m))
		plt.show()
		m = m * 2
		
#-----------------------------------------------------------------------------------------
#homework solved in order:
#PART A 
	#eigenvalues of 100 superimposed matrices for dimensions 8, 16, 32, 64, 128, 256
randMatrixEigPlot(8, 256)
	#spectral radius of a random matrix as n-> inf
spectRadPlot(8, 128)
#PART B
	#2-norm of random matrix as n-> inf
twoNormSpecTwo(8, 128)
#PART C
	# AVERAGE trend in condition numbers as dim increases ANNNNDDDDDDDDDDDDDDDDDDDd probability distribution of smallest singular values
condNums(8, 64, 0.02, 0.14)
