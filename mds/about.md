<!--- 
  NOTE: this file is modified by running `redocument`, via `README.Rmd`
  only the Version line is modified.
 -->

## NetPredictor: Analyzing and Predicting Links in Bipartite Network data

* **Version:** 1.0.1
* ![Status](http://img.shields.io/badge/status-In_development_%28UNSTABLE%29-red.svg?style=flat)
* [![License](http://img.shields.io/badge/license-AGPL--3-orange.svg?style=flat)](https://www.gnu.org/licenses/agpl-3.0.html)
* **Source:** [GitHub](https://github.com/abhik1368/Shiny_NetPredictor)
* **Problems:** [Bug reports and feature reqests](https://github.com/abhik1368/Shiny_NetPredictor/issues)
* **Authors:** [Abhik Seal](https://www.linkedin.com/in/abseal)

# Start Prediction

This tab allows you to select your data source for use with the Netpredictor
system.  Every other part of Netpredictor requires you having dataset loaded first
from the start prediction tab. Until you select a custom data or example dataset
the following tabs won't work.

#### Custom Data

* **Drug-Target Biparite Network :**
Click the 'Browse' button to select your local Drug Target binary matrix csv file you wish to
upload . The matrix rows should be protein names and columns should be drugnames. The file should end in the extension `.csv`.

* **Drug Similarity Matrix :**.
Click the 'Browse' button to select your local Drug Similarity Matrix csv file you wish to
upload . The file should end in the extension `.csv`.

* **Target Similarity Matrix :**
Click the 'Browse' button to select your local Target Similarity Matrix csv file you wish to
upload .  The file should end in the extension `.csv`.

The count of drugs and targets should match all the given three matrices. Also the order in which drugnames and proteins are given in the Drug-Target bipartite network this should be maintained drug-similairty matrix and Target similarity matrix, otherwise the results will be misleading.

#### Example Data

Data have already been pre-processed for use with the Netpredictor System, so they should load very quickly. Four types of datasets are provided namely Enzyme, GPCR,Ion Channel and Nuclear Receptor from [Yamanishi's website](http://cbio.ensmp.fr/~yyamanishi/pharmaco/).

#### Network Prediction Algorithms

Four Types of network algorithms are being implemented in the current Netpredictor system which are described below.

* **HeatS :**
The algorithm based on the recommendation techniques developed by [Zhou et al](http://www.pnas.org/content/107/10/4511.full).It is analogous to heat diffusion across the user-object network. The algorithm works on bipartite network projection technique implementing the concept of resources transfer within the network.

* **Network Based Inference(NBI) :**
This method is also based on the two phase resource transfer as the HeatS algorithm. In addition this method extends extends the recommendation model by similarity between small molecules and sequence similarity between targets. This is basedo on the DT-Hybrid algorithm discussed by [Alaimo et al](http://bioinformatics.oxfordjournals.org/content/29/16/2004.long). **alpha** and **lamda** are  tuning parameters (value between 0 and 1) to adjust the performance of the algorithm.


* **Random walk with Restart(RWR) :**
This method based on the concept of widely known pagerank algorithm which computes random walk-based “distance” from a node to every other nodes in a network. More Information can be found by [Chen et al](http://pubs.rsc.org/en/Content/ArticleLanding/2012/MB/c2mb00002d#!divAbstract) and [Seal et al.](http://www.jcheminf.com/content/7/1/40).Parameter optimization is done in this implementation we optimized **eta** to 0.99 and **restart** paramter is provided with ranges between 0-1.


* **Netcombo :**
NetCombo algorithm fuses the results of NBI and RWR. Its takes the parameters for both the algorithms and averages the results. 

#### Bipartite Network Modules
It searches for bipartite modules in the network using label propgation BRIM algorithm The algorithm consists of two stages.First during the LP phase, neighbouring nodes (i.e. those
which share links) exchange their labels representing the community they belong to, with each node receiving the most common label amongst its neighbours. We iterate this process until
densely connected groups of nodes reach a consensus of what is the most representative label, as indicated by the fact that the modularity is not increased by additional exchanges. Secondly,the BRIM algorithm(2) refines the partitions found with label propagation.

(1) Liu, X. & Murata, T. (2010) Community detection in large-scale bipartite networks.Information
and Media Technologies, 5, 184–192.

(2) Barber, M. (2007) Modularity and community detection in bipartite networks. Physical Review E, 76, 066102.

(3)[Optimization of bipartite modularity using LP-BRIM](https://github.com/PoisotLab/lpbrim)

# Statistical Analysis

This tab allows you to perform statistical analysis on your loaded dataset. In advanced analysis we remove links from the network randomly based frequency of drug target associations. The method then re-predicts the links performs statistical analysis how well those remove links can be predicted and which algorithm shows better performance. Five of the statistical metrics used to show the performance **AUAC**(Area under the accumulation curve),**auc**(Area under curve), **auctop**(10%) ,**bedroc** and **enrichment factor**. The results dynamically adds up the table and at each run to the table. One can select three algorithms for calculation of the statistics

* **Choose Random links to be removed :**
Remove links from the network randomly

* **Frequency of associations between Biparite Nodes :**
Remove the nodes from network which are having more than certain frequency of drug target associations.

# Permutation Testing

In the advanced analysis tab perform permutation testing can be done on your loaded dataset. From the given network random permutations are computed using the given number of permutations and only the 
significant links with the user given pvalues are kept.

* **Choose number of random permutations :**
Perform random permutations of the network. As the number of random permutations increases the
computing time increases 

* **Keep Significant links of pvalue :**
set the threshold for links to keep in the results table.

# Search Drugbank 

This tab allows users to search predicted drugbank interactions using either Network based inference(NBI) or Random walk with restart (RWR). The tab allows users to search by drugbank id and proteins by hugo gene names. The data is accessed using RSQLite package. The dat table list the drugbank id, uniprot id,protein names,significant values, ATC Codes, Drug Categories, predicted and true interactions. 

# Ontology and pathway search

This tab allows users to identify functional characteritics of the predicted and true interactions using gene ontology(BP,CC,MF) features.Users can select gene ontology with the depth of the level or use pathway to get enrich pathways for the given set of genes. The textbox accepts (,) separated hugo gene names. 

### Acknowledgements

This project's rapid creation would not have been possible without
the numerous excellent R packages available to us.  We wish to
acknowledge these and thank the authors for their work:

* [shiny](http://cran.r-project.org/web/packages/shiny/index.html), [markdown](http://cran.r-project.org/web/packages/markdown/index.html), and 
* [RRO Open](https://mran.revolutionanalytics.com/open/),[Netpredictor](https://github.com/abhik1368/netpredicter),[igraph](https://cran.r-project.org/web/packages/igraph/index.html)
* Some amazing cool shiny packages like [shinysky](https://github.com/AnalytixWare/ShinySky),[shinyBS](https://github.com/ebailey78/shinyBS),[shinythemes](https://github.com/rstudio/shinythemes)
* [Netpredictor](https://github.com/abhik1368/netpredictor)
* Some javascript packages such as [visNetwork](http://dataknowledge.github.io/visNetwork/),[data.tables](https://cran.r-project.org/web/packages/data.table/index.html)




