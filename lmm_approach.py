'''
lmm_approach.py

MEW 12/10/2014

Includes scripts for making a sparse covariance matrix from a dictionary
of trios, and for performing LMM with that sparse matrix and a list of
known phenotypes.
'''

try:
	import sys, os, operator

	''' For working with sparse adjacency matrices '''
	import scipy.sparse as sps
	import scipy.sparse.linalg as ssl
	import numpy as np

	''' For using a temporary directory and cleaning up when finished '''
	import shutil, tempfile

	''' For executing idcoefs and tsort from within Python '''
	import subprocess

except ImportError as e:
	print >> sys.stderr, 'Error importing modules:\n%s' % e
	sys.exit(1)
	
def construct_cov_matrix(trio_dict, gedcom_ids):
	'''
	Given a trio_dictionary (key: child GEDCOM id; value: tuple with mother and
	father's ids) containing n entries, constructs an nxn symmetric sparse
	adjacency matrix. Also returns a list giving the correspondence between
	indices of the adjacency matrix and GEDCOM ids.
	'''
	try:
		tmp_path = tempfile.mkdtemp(prefix='idcoefs_')
		
		''' prepare the input for idcoefs '''
		gedcom_ids.insert(0, '0') # a placeholder needed for idcoefs
		write_dag(trio_dict, tmp_path, gedcom_ids)
		tsort_and_reformat_for_idcoefs(trio_dict, gedcom_ids, tmp_path)
		gedcom_ids.pop(0)
		
		''' run idcoefs '''
		print '--- Begin idcoefs output ---'
		command = 'idcoefs -p %s -s %s -o %s' % (os.path.join(tmp_path,
			'tmp.ped'), os.path.join(tmp_path, 'tmp.sample'),
	    	os.path.join(tmp_path, 'idcoefs.out'))
		subprocess.call(command, shell=True)
		print '--- End idcoefs output ---'

		''' make a sparse covariance matrix from the output '''
		gen_cov_matrix = process_idcoefs_output(len(trio_dict.keys()), tmp_path)
	finally:
		shutil.rmtree(tmp_path)

	return gen_cov_matrix

def write_dag(trio_dict, tmp_path, gedcom_ids):
	try:
		dag_file = open(os.path.join(tmp_path, 'tmp.dag'), 'w', 1)

		for child, parents in trio_dict.items():
			mother, father = parents
			if mother is None:
				mother = '0'
			if father is None:
				father = '0'
			dag_file.write('%s\t%s\n' % (gedcom_ids.index(mother),
				gedcom_ids.index(child)))
			dag_file.write('%s\t%s\n' % (gedcom_ids.index(father),
				gedcom_ids.index(child)))
	except Exception as e:
		raise Exception('Error while writing directed acyclic graph:\n%s') % e
	finally:
		dag_file.close()

	return

def tsort_and_reformat_for_idcoefs(trio_dict, gedcom_ids, tmp_path):
	command = 'tsort %s > %s' % (os.path.join(tmp_path, 'tmp.dag'), 
			os.path.join(tmp_path, 'tmp.dag.sorted'))
	subprocess.call(command, shell=True)

	''' Make the input files needed for idcoefs using the result '''
	try:
		sorted_file = open(os.path.join(tmp_path, 'tmp.dag.sorted'), 'r')
		ped_file = open(os.path.join(tmp_path, 'tmp.ped'), 'w', 1)
		sample_file = open(os.path.join(tmp_path, 'tmp.sample'), 'w', 1)

		for line in sorted_file:
			id = int(line.strip())
			if id == 0:
				continue
			(mother, father) = trio_dict[gedcom_ids[id]]
			if mother is None:
				mother = '0'
			if father is None:
				father = '0'

			ped_file.write('%d\t%d\t%d\n' % (id, gedcom_ids.index(mother), gedcom_ids.index(father)))
			sample_file.write('%d\n' % id)
	finally:
		sorted_file.close()
		ped_file.close()
		sample_file.close()

	return

def process_idcoefs_output(n, tmp_path):
	try:
		idcoefs_file = open(os.path.join(tmp_path, 'idcoefs.out'), 'r')
		row_num = []
		col_num = []
		vals = []

		for line in idcoefs_file:
			'''
			Get the numerical identifier (which is the same as the covariance
			matrix index for each of the two people compared on this line). We
			subtract one because idcoefs used '0' to indicate a founder, but we
			need to start our indices from there.
			'''
			fields = line.strip().split('\t')
			i = int(fields[0]) - 1
			j = int(fields[1]) - 1

			''' Wright's coefficient of relationship is 2 * theta_ij '''
			wrights_coeff = 2 * float(fields[2]) + float(fields[4]) + \
			    float(fields[6]) + float(fields[8]) + 0.5 * float(fields[9])
			if wrights_coeff == 0:
				continue

			'''
			The covariance matrix is symmetric, so we add two entries unless
			we're on the diagonal
			'''
			row_num.append(i)
			col_num.append(j)
			vals.append(wrights_coeff)
			if i != j:
				row_num.append(j)
				col_num.append(i)
				vals.append(wrights_coeff)
		gen_cov_matrix = sps.coo_matrix((vals, (row_num, col_num)),
			shape=(n,n))
	finally:
		idcoefs_file.close()

	return gen_cov_matrix

'''
The following subroutines are for simulating phenotypes based on a genetic
covariance matrix and heritability h^2, and for revealing the appropriate
subset of phenotypes prior to LMM.
'''
def simulate_phenotypes(gen_cov_matrix, h2, gedcom_ids, **kwargs):
	'''
	Does not simulate the genetic effects: instead, simulates the phenotypes
	directly by updating the covariance matrix to reflect noise.
	'''
	n = len(gedcom_ids)
	if 'mean' in kwargs:
		mean_array = np.ones(n) * kwargs['mean']
	else:
		mean_array = np.ones(n) * 0.0

	cov_matrix = h2 * gen_cov_matrix + (1 - h2) * np.eye(n)

	phenotypes = np.random.multivariate_normal(mean_array, cov_matrix, 1)
	phenotypes = np.array(phenotypes).reshape(-1,).tolist()
	phenotype_dict = dict(zip(gedcom_ids, phenotypes))

	return phenotype_dict

def simulate_genetic_effects_and_noise(gen_cov_matrix, h2, gedcom_ids, **kwargs):
	'''
	Simulates the genetic effects using the genetic covariance matrix. Then
	adds an appropriate amount of noise on top to get the phenotypes.
	'''
	n = len(gedcom_ids)
	if 'mean' in kwargs:
		mean_array = np.ones(n) * kwargs['mean']
	else:
		mean_array = np.ones(n) * 0.0

	genetic_effects = np.random.multivariate_normal(mean_array, h2*gen_cov_matrix, 1)
	genetic_effects = np.array(genetic_effects).reshape(-1,).tolist()
	genetic_effect_dict = dict(zip(gedcom_ids, genetic_effects))

	noise = np.random.normal(0.0, (1-h2)**0.5, n)
	noise = np.array(noise).reshape(-1,).tolist()
	phenotypes = [i+j for i,j in zip(genetic_effects, noise)]
	phenotype_dict = dict(zip(gedcom_ids, phenotypes))

	return genetic_effect_dict, phenotype_dict

'''
The following subroutines are for calculating the LMM once we have the
sparse genetic covariance matrix in hand.

Our model:

Y = X \Beta + Z \Gamma + \Epsilon

(Actually in this script we will not include fixed effects, but I'll describe
how it would be implemented for f fixed effects.) Assume we have a covariance
matrix for n relatives and that we know the phenotypes for a subset k of them.
Let's also assume that I've ordered the individuals so that the first k rows/
columns of the covariance matrix correspond to the individuals of known
(and, let's say, quantitative) phenotype. Then:

* Y is a k x 1 column vector of observed phenotypes

* X is a k x f matrix. The first column is all ones; the rest correspond
  to the f-1 non-constant fixed effects traits.

* \Beta is a f x 1 column vector. The first row is the mean value of the
  trait; the rest are the weights of the other fixed effects.

* Z is k x n matrix, all zeros except for ones along the first diagonal [up
  to position (k,k)].

* \Gamma is an n x 1 column vector of latent individual genetic effects
  distributied with mean 0 and covariance = the genetic covariance matrix.

To estimate \Beta and \Gamma, we will:

* Calculate the inverse of the (sparse) covariance matrix via Cholesky
  decomposition: A = LDL', A^-1 = (D^-1 * L^-1)'(D^-1 * L^-1)

* From the heritability, estimate the variance ratio alpha.

* Define a "left-hand-side" matrix:
   ________________________________________________
  |              |                                |
  |   X' * X     |            X' * Z              |
  |   f x f      |            f x n               |
  |______________|________________________________|
  |              |                                |
  |   Z' * X     | Z' * Z + \alpha * inv(gen_cov) |
  |   k * f      |             n x n              |
  |______________|________________________________|

  ...and calculate its inverse

* Define a "right-hand-side" matrix:
 __________
|          |
|  X' * Y  |
|  f x 1   |
|__________|
|          |
|  Z' x Y  |
|  n x 1   |
|__________|

* Get the best estimates of the \beta_i and \gamma_i by multiplying:
  inv(LHS) * RHS

* Get the standard error in these estimates from:
  diag(inv(LHS)) * \sigma_e^2 (if known)

* Get the "reliabilities" ("accuracies" squared) from:
  1 - diag(inv(LHS)) * \alpha
  ...and then we'll be done!
'''

def predict_genetic_effects_only_fixed_effect_is_mean(h2, gen_cov_matrix, gedcom_ids, evidence_dict):
	'''
	It is easier to work with alpha than with the heritability h^2.
	alpha = sigma_e^2 / sigma_a^2 = (1 - h^2)/h^2
	'''
	alpha = (1 - h2)/h2

	''' Get the matrices we need for LMM. gen_cov_matrix is assumed sparse. '''
	(n, n) = gen_cov_matrix.shape

	'''
	There is a modest increase in speed by calculating the inverse by
	decomposition. However, scikits.sparse.cholmod is installed incorrectly
	on my laptop, so I have suppressed that for now.
	gen_cov_factor = cholesky(gen_cov)
	gen_cov_inverse = gen_cov_factor.inv()
	'''
	gen_cov_inverse = ssl.inv(sps.csc_matrix(gen_cov_matrix))

	X, Y, Z = make_XYZ_only_fixed_effect_is_mean(n, gedcom_ids, evidence_dict)
	Ztranspose = Z.transpose()
	Ztranspose_Z = Ztranspose.dot(Z)
	Xtranspose = X.transpose()
	Xtranspose_X = Xtranspose.dot(X) 
	Xtranspose_Z = Xtranspose.dot(Z) 
	Ztranspose_X = Ztranspose.dot(X)

	LHS_fourth_quadrant = Ztranspose_Z + alpha * gen_cov_inverse
	LHS = sps.vstack((sps.hstack((Xtranspose_X, Xtranspose_Z),
		format='csc'), sps.hstack((Ztranspose_X,LHS_fourth_quadrant),
		format='csc')), format='csc')

	'''
	See note above about inversion using decomposition. The speed increase
	was more modest here since the LHS matrix is not as sparse.
	LHS_factor = cholesky(LHS)
	LHS_inverse = LHS_factor.inv()
	'''
	LHS_inverse = ssl.inv(LHS)
	
	RHS = sps.vstack((Xtranspose.dot(Y), Ztranspose.dot(Y)), format='csc')
	
	solution = np.array(((LHS_inverse.dot(RHS)).todense())).reshape(-1,).tolist()
	'''
	solution is an (n + f + 1) x 1 column vector. Since the number
	of fixed effects f is currently zero, solution[0] is the fitted
	population mean and solution[1] through solution[n] are the
	fitted genetic effects. Assuming we don't care about the population
	mean, so I remove that.
	'''
	mean = solution.pop(0)

	return mean, solution


def make_XYZ_only_fixed_effect_is_mean(n, gedcom_ids, evidence_dict):
	k = len(evidence_dict.keys())
	X = sps.csc_matrix(np.ones((k,1)))

	''' Build Z as a sparse matrix '''
	row_num = range(0,k)
	vals = [1] * k
	col_num = []
	Y = np.zeros((k,1))
	revealed_phenotype_ids = []
	counter = 0
	for i, gedcom_id in enumerate(gedcom_ids):
		if gedcom_id in evidence_dict.keys():
			col_num.append(i)
			Y[counter,0] = evidence_dict[gedcom_id]
			counter += 1

	'''
	for gedcom_id in evidence_dict.keys():
		col_num.append(gedcom_ids.index(gedcom_id))
		Y[counter,0] = evidence_dict[gedcom_id]
		counter += 1
	'''
	Z = sps.coo_matrix((vals, (row_num, col_num)), shape=(k,n))
	Z.tocsc()
	Y = sps.csc_matrix(Y)

	return X, Y, Z

