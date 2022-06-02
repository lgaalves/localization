import warnings
warnings.filterwarnings('ignore')

#Standard libraries
import pathlib
import os
import random
from scipy import stats

#Third part libraries
import pandas as pd
import numpy as np
import networkx as nx

import seaborn as sns
import matplotlib.mlab as ml
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib.ticker as plticker

def stdfigsize(scale=1, nx=1, ny=1, ratio=1.3):
    """
    Returns a tuple to be used as figure size.
    -------
    returns (7*ratio*scale*nx, 7.*scale*ny)
    By default: ratio=1.3
    If ratio<0 them ratio = golden ratio
    """
    if ratio < 0:
        ratio = 1.61803398875
    return((7*ratio*scale*nx, 7*scale*ny))

def stdrcparams(usetex=False):
    """
    Set several mpl.rcParams and sns.set_style for my taste.
    ----
    usetex = True
    ----
    """
    sns.set_style("white")
    sns.set_style({"xtick.direction": "in",
                 "ytick.direction": "in"})
    rcparams = {
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica'],
    'axes.labelsize': 28,
    'axes.titlesize': 28,
    'legend.fontsize': 20,
    'ytick.right': 'off',
    'xtick.top': 'off',
    'ytick.left': 'on',
    'xtick.bottom': 'on',
    'xtick.labelsize': '25',
    'ytick.labelsize': '25',
    'axes.linewidth': 2.5,
    'xtick.major.width': 1.8,
    'xtick.minor.width': 1.8,
    'xtick.major.size': 14,
    'xtick.minor.size': 7,
    'xtick.major.pad': 10,
    'xtick.minor.pad': 10,
    'ytick.major.width': 1.8,
    'ytick.minor.width': 1.8,
    'ytick.major.size': 14,
    'ytick.minor.size': 7,
    'ytick.major.pad': 10,
    'ytick.minor.pad': 10,
    'axes.labelpad': 15,
    'axes.titlepad': 15,
    'axes.spines.right': False,
    'axes.spines.top': False
}
    mpl.rcParams.update(rcparams) 

mpl.rcParams['lines.linewidth'] = 5
mpl.rcParams['lines.color'] = '#3690c0'

stdrcparams(usetex=True)
figsize=stdfigsize(ratio=-1)
xs,ys=figsize

def custom_frame(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.tick_params(axis='x',length=10,direction='out')
    ax.tick_params(axis='x',which='minor',direction='out')
    ax.tick_params(axis='y',length=10,direction='out')
    ax.tick_params(axis='y',which='minor',direction='out')

plt.rcParams['legend.title_fontsize'] = '20'
plt.rcParams['pdf.fonttype'] = 42 
def eigenvector_plot(x,y,df):
    fig, ax = plt.subplots(figsize=(xs,ys))
    norm = mpl.colors.Normalize(vmin=0, vmax=43)
    cmap = plt.cm.viridis_r
    mapcolor = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cmaplist = [mapcolor.to_rgba(i) for i in range(0,43)]
    ax.bar(x,y,color=cmaplist)
    ax.set_xlim(-1,43)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
    labels=["",""]
    labels.extend(list(df.Countries))
    ax.set_xticklabels(labels,rotation=90,fontsize=15,ha='center')
    return ax

def to_single_layer(edges_list):
    single_layer_net = {}

    for e in edges_list:
        src_node, src_layer, dst_node, dst_layer, weight = e
        if not (src_node in single_layer_net):
            single_layer_net[src_node] = {}
        if not (dst_node in single_layer_net[src_node]):
            single_layer_net[src_node][dst_node] = 0.0
        single_layer_net[src_node][dst_node] += weight

    edges_list=[]
    for src_node in single_layer_net:
        for dst_node in single_layer_net[src_node]:
            edges_list.append([src_node, dst_node,single_layer_net[src_node][dst_node]])

    return edges_list

def read_edge_list(year=2010):
    data_path="../data/{}/multilayer/WIOT{}_Nov16_ROW.edges".format(year,year)
    edges_list, num_layers, num_countries=load_network(data_path)
    edge_df=pd.DataFrame(edges_list,columns=["source_node","source_layer","destination_node","destination_layer","weight"])
    return edge_df

def read_edge_df(year=2010):
    data_path="../data/{}/multilayer/edges-{}.csv".format(year,year)
    edge_df=pd.read_csv(data_path)[['source_node', 'source_layer', 'destination_node','destination_layer', 'weight']]
    return edge_df

def read_supra_matrix(year=2010):
  edge_df=read_edge_df(year=year)
  number_of_countries=43
  number_of_layers=56
  m=np.zeros((number_of_layers*number_of_countries,number_of_layers*number_of_countries))
  for alpha in range(1,number_of_layers+1,1):
      for beta in range(1,number_of_layers+1,1):
          edges=edge_df[(edge_df.source_layer==alpha) & (edge_df.destination_layer==beta)][["source_node","destination_node",'weight']]
          edges=np.array(edges)
          for edge in edges:
              m[number_of_countries*(alpha-1):number_of_countries*alpha,
                number_of_countries*(beta-1):number_of_countries*beta][int(edge[0]-1)][int(edge[1]-1)]=edge[2]
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

def f(x,a,b):
    return a*x+b

def conf_intervals_pearsonr(x,y,samples=100,a=0.05):
    data=pd.DataFrame(np.transpose([x,y]),columns=["x","y"])
    
    x,y=np.array(data.x),np.array(data.y)
    rho,pvalues=stats.pearsonr(x,y)
    rhos=[]
    for i in range(0,samples):
        ndata=data.sample(frac=1,replace=True)
        x,y=np.array(ndata.x),np.array(ndata.y)
        rho,pvalues=stats.pearsonr(x,y)
        rhos.append(rho)
    rho=np.average(rho)
    rhosinf,rhossup=pd.Series.quantile(pd.Series(rhos),a/2),pd.Series.quantile(pd.Series(rhos),1-a/2)
    return rho,rhosinf,rhossup