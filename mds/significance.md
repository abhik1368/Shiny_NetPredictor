## load data

This tab allows you to select your data source for use with the Netpredictor
system.  Every other part of Netpredictor requires you having dataset loaded first
from the start prediction tab. Until you select a custom data or example dataset
the following tabs won't work.

#### Custom Data

* **Drug-Target Biparite Network**
Click the 'Browse' button to select your local Drug Target binary matrix csv file you wish to
upload .  The file should end in the extension `.csv`.

* **Drug Similarity Matrix**.
Click the 'Browse' button to select your local Drug Similarity Matrix csv file you wish to
upload .  The file should end in the extension `.csv`.

* **Target Similarity Matrixs**
Click the 'Browse' button to select your local Target Similarity Matrix csv file you wish to
upload .  The file should end in the extension `.csv`.



#### Example Data

Data have already been pre-processed for use with the Netpredictor System, so they should load very quickly. Four types of datasets are provided namely Enzyme, GPCR,Ion Channel and Nuclear Receptor from [Yamanishi's website](http://cbio.ensmp.fr/~yyamanishi/pharmaco/).

#### Network Prediction Algorithms

Four Types of network algorithms are being implemented in the current Netpredictor system which are described below.
* **HeatS**
    
    The algorithm based on the recommendation techniques developed by [Zhou et al](http://www.pnas.org/content/107/10/4511.full).It is analogous to heat diffusion across the user-object network. The algorithm works on bipartite network projection technique implementing the concept of resources transfer within the network.

* **Network Based Inference(NBI)**
    
    This method is also based on the two phase resource transfer as the HeatS algorithm. In addition this method extends extends the recommendation model by similarity between small molecules and sequence similarity between targets. This is basedo on the DT-Hybrid algorithm discussed by [Alaimo et al](http://bioinformatics.oxfordjournals.org/content/29/16/2004.long). **alpha** and **lamda** are  tuning parameters (value between 0 and 1) to adjust the performance of the algorithm.


* **Random walk with Restart(RWR)**
    
    This method based on the concept of widely known pagerank algorithm which computes random walk-based “distance” from a node to every other nodes in a network. More Information can be found by [Chen et al](http://pubs.rsc.org/en/Content/ArticleLanding/2012/MB/c2mb00002d#!divAbstract) and [Seal et al.](http://www.jcheminf.com/content/7/1/40).Parameter optimization is done in this implementation we optimized **eta** to 0.99 and **restart** paramter is provided with ranges between 0-1.

* **Netcombo**

NetCombo algorithm fuses the results of NBI and RWR. Its takes the parameters for both the algorithms and averages the results. 
                                                                                                                                                                                                                 
                                                                                                                                                                                                                 