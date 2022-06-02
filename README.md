# The rise and fall of countries in the global value chains

##### What is this repository for?
This repository contains the code for the analysis reported in [Scientific Reports 12, 9086 (2022)](https://www.nature.com/articles/s41598-022-12067-x.pdf).

# Summary

<image src='featured.png' />

Countries become global leaders by controlling international and domestic transactions connecting geographically dispersed production stages. We model global trade as a multi-layer network and study its power structure by investigating the tendency of eigenvector centrality to concentrate on a small fraction of countries, a phenomenon called localization transition. We show that the market underwent a significant drop in power concentration precisely in 2007 just before the global financial crisis. That year marked an inflection point at which new winners and losers emerged and a remarkable reversal of leading role took place between the two major economies, the US and China. We uncover the hierarchical structure of global trade and the contribution of individual industries to variations in countries’ economic dominance. We also examine the crucial role that domestic trade played in leading China to overtake the US as the world’s dominant trading nation. There is an important lesson that countries can draw on how to turn early signals of upcoming downturns into opportunities for growth. Our study shows that, despite the hardships they inflict, shocks to the economy can also be seen as strategic windows countries can seize to become leading nations and leapfrog other economies in a changing geopolitical landscape.


### Disclaimer 

Please note that this repository is not intended for wide-spread distribution. We are only making the code available so that other researchers may reproduce the results published in our manuscript. 

### Directory Structure 

* [data](data) -- WIOD data;
* [notebooks](notebooks) -- Jupyter notebooks;
* [src](src) -- Python modules;
* [results](results) -- Results of simulations;
* [figures](figures) -- Figures are saved here;

# Data 

The WIOD data is freely available for download [here](http://www.wiod.org/database/wiots16).
Download and store the dataset at the folder called data Then run the script `raw_to_edges.py` to generate the multilayer edge data. This is step is necessary before running the Jupyter notebooks. 

# Jupyter notebooks' index

1. [Main data analysis notebook](notebooks/01_data_analysis.ipynb)
2. [Synthetic multilayer network models](notebooks/02_null_model.ipynb)
3. [Radial cluster plot](notebooks/03_radial_cluster.ipynb)

## Requirements

#### Python
```
Python version       : 3.8.12
IPython version      : 7.18.1
matplotlib: 3.4.3
seaborn   : 0.11.2
sys       : 3.8.12  
pandas    : 1.3.4
numpy     : 1.21.4
palettable: 3.3.0
networkx  : 2.4
graph_tool: 2.43 
scipy     : 1.7.3
cartopy   : 0.20.2
watermark : 2.3.0
```
#### R (for radial cluster plot only)

```
ggplot2
gridExtra
dplyr
grid
RColorBrewer
ggrepel
ggthemes
viridisLite
factoextra
```


# Reference

Luiz G. A. Alves, Giuseppe Mangioni, Francisco A. Rodrigues, Pietro Panzarasa, and Yamir Moreno, The rise and fall of countries in the global value chains. [Scientific Reports 12, 9086 (2022)3](https://www.nature.com/articles/s41598-022-12067-x.pdf)