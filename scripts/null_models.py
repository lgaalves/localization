#!/usr/bin/python
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd 

from multiprocessing import Pool
from functools import partial


countries=pd.read_csv("../data/countries.csv",delimiter=",",index_col=None)
sectors=pd.read_csv("../data/sectors.csv",delimiter="\t",index_col=None)

def read_edge_df(year=2010):
	data_path="../data/{}/multilayer/edges-{}.csv".format(year,year)
	print(data_path)
	edge_df=pd.read_csv(data_path,index_col=0)
	return edge_df

def to_multiplex(edge_df):
    df=edge_df[['source_node','source_layer','destination_node','weight']]
    df=df.groupby(['source_node','source_layer','destination_node']).sum().reset_index()
    df['destination_layer']=df['source_layer']
    return df

def to_simplex(edge_df):
    df=edge_df[['source_node','destination_node','weight']]
    df=df.groupby(['source_node','destination_node']).sum().reset_index()
    return df

def edges_to_supramatrix(edge_df):
	number_of_countries=43
	number_of_layers=56
	m=np.zeros((number_of_layers*number_of_countries,number_of_layers*number_of_countries))
	for alpha in range(1,number_of_layers+1,1):
		for beta in range(1,number_of_layers+1,1):
			edges=edge_df[(edge_df.source_layer==alpha)& (edge_df.destination_layer==beta)][["source_node","destination_node",'weight']]
			edges=np.array(edges)
			for edge in edges:
				m[number_of_countries*(alpha-1):number_of_countries*alpha,
				  number_of_countries*(beta-1):number_of_countries*beta][int(edge[0]-1)][int(edge[1]-1)]=edge[2]
	return m

def multiplex_supra_matrix(edge_df):
	number_of_countries=43
	number_of_layers=56
	m=np.zeros((number_of_layers*number_of_countries,number_of_layers*number_of_countries))
	for beta in range(1,number_of_layers+1,1):
		for alpha in range(1,number_of_layers+1,1):
			edges=edge_df[(edge_df.source_layer==alpha) & (edge_df.destination_layer==beta)][["source_node","destination_node",'weight']]
			edges=np.array(edges)
			for edge in edges:
				m[number_of_countries*(alpha-1):number_of_countries*alpha,
				  number_of_countries*(beta-1):number_of_countries*beta][int(edge[0]-1)][int(edge[1]-1)]=edge[2]
	return m

def to_simplex_matrix(edge_list):
	import scipy.sparse as sparse

	arr = np.array(edge_list)
	m,n = tuple(arr.max(axis=0)[:2]+1)
	shape=(int(m),int(n))
	print(shape)
	coo = sparse.coo_matrix((arr[:, 2], (arr[:, 0], arr[:, 1])), shape=shape,
							dtype=arr.dtype)
	return coo.toarray()

def convert_to_matrices(edge_df,year):
	print('Convert edge list to multilayer matrix: {}'.format(year))
	ml=edges_to_supramatrix(edge_df)
	return ml

def null_model(supramatrix,n_nodes=43,n_layers=56):
	""" 
	This model reshuffle the links and weights by blocks/layer of the supra-matrix. 
	Thus, the reshuffled matrix has the same degree, strength, and weight distributon of both layer and node.
	"""
	m=supramatrix.copy()
	for i in range(0,n_layers):
		for j in range(0,n_layers):
			for k in range(0,n_nodes):
				np.random.shuffle(m[i*n_nodes:(i+1)*n_nodes,j*n_nodes:(j+1)*n_nodes][k])
	return m

def mnet_eigen_centrality(M,node_labels,n_layers,n_nodes,max_iter=50):
	"""
	Calculates the eingenvector centrality for a multilayer network. 

	
	De Domenico, M. et al. Ranking in interconnected multilayer
	networks reveals versatile nodes. Nat. Commun. 6:6868 doi: 10.1038/ncomms7868
	(2015).
	
	input: 
		- M is an numpy array with the supra matrix of the multilayer network.
		- n_layers is the number of layers
		- n_nodes is the number of nodes in each layer
	output: 
		- centralitymnet is the centrality of each node using the whole multilayer structure.
	
	"""
	import scipy as sp
	from scipy.sparse import linalg 
	eigenvalue, eigenvector = linalg.eigs(M.T, k=1, which='LR' ,maxiter=max_iter )
	largest = eigenvector.flatten().real
	norm = sp.sign(largest.sum())*sp.linalg.norm(largest)
	centrality = largest/norm
	centrality=largest.reshape((n_layers,n_nodes))
	u=np.array([1 for i in range(0,n_layers)])
	centralitymnet=u.dot(centrality)
	norm = sp.sign(centralitymnet.sum())*sp.linalg.norm(centralitymnet)
	centralitymnet = centralitymnet/norm
	if node_labels==False:
		 centralitymnet=dict(zip(np.arange(0,len(centralitymnet)),centralitymnet))
	else:
		centralitymnet=dict(zip(node_labels,centralitymnet))
	return centralitymnet

def core_ipr_ensemble(m,n_layers,i):
	null_m=null_model(supramatrix=m,n_nodes=43,n_layers=n_layers)
	centrality_null=mnet_eigen_centrality(null_m,node_labels=list(countries.Country),n_layers=n_layers ,n_nodes=43,max_iter=1000000)
	bc=np.array(list(centrality_null.values()))
	ipr=np.sum(bc**4)
	return ipr

def calculate_samples(m_list,nsamples,which,simplex=True,processes=2):

	if simplex:
		n_layers=1
	else:
		n_layers=56

	results=[]
	ensembles=[]
	for m in m_list:
		if which!='import':
			m=m.T

		#Calculate centrality for real network
		centrality=mnet_eigen_centrality(m,node_labels=list(countries.Country),n_layers=n_layers ,n_nodes=43,max_iter=1000000)
		bc=np.array(list(centrality.values()))
		ipr_real=np.sum(bc**4)


		pool=Pool(processes=processes) #choose process as the number of cores available in your machine
		func=partial(core_ipr_ensemble,m,n_layers)
		ipr_ensemble=pool.map(func, np.arange(nsamples))
		pvalue=(np.sum(np.array(ipr_ensemble)>=ipr_real)/len(ipr_ensemble))
		results.append([ipr_real,np.mean(ipr_ensemble),np.std(ipr_ensemble),np.percentile(ipr_ensemble,2.5),np.percentile(ipr_ensemble,97.5),pvalue])
		ensembles.append(ipr_ensemble)
	return results,ensembles

def main():
	# Set parameters
	NSAMPLES=1000

	# Read data
	print('Reading data \n')
	m_list_multilayer=[]
	for year in range(2000,2015):
		print('Reading edge list {} \n'.format(year))
		edge_df=read_edge_df(year=year)
		m=convert_to_matrices(edge_df,year)

		m_list_multilayer.append(m)

	# Multilayer
	print('Calculating CI for multilayer')
	print('Imports')
	values,ensembles=calculate_samples(m_list_multilayer,NSAMPLES,'import',simplex=False)
	df=pd.DataFrame(values,columns=['ipr_real','mean','std','lb','ub','pvalue'])
	df.to_csv('../results/null_model/multilayer_import.csv')
	print('Exports')
	values,ensembles=calculate_samples(m_list_multilayer,NSAMPLES,'export',simplex=False)
	df=pd.DataFrame(values,columns=['ipr_real','mean','std','lb','ub','pvalue'])
	df.to_csv('../results/null_model/multilayer_export.csv')

if __name__ == "__main__":
	main()