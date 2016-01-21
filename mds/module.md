## Bipartite Network Modules

It searches for bipartite modules in the network using label propgation BRIM algorithm

#### LP&BRIM Algorithm

The algorithm consists of two stages.First during the LP phase, neighbouring nodes (i.e. those
which share links) exchange their labels representing the community they belong to, with each node receiving the most common label amongst its neighbours. We iterate this process until
densely connected groups of nodes reach a consensus of what is the most representative label, as indicated by the fact that the modularity is not increased by additional exchanges. Secondly,the BRIM algorithm(2) refines the partitions found with label propagation.

(1) Liu, X. & Murata, T. (2010) Community detection in large-scale bipartite networks.Information
and Media Technologies, 5, 184–192.

(2) Barber, M. (2007) Modularity and community detection in bipartite networks. Physical Review E, 76, 066102.


