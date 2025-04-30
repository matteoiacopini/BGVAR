# BGVAR

This archive contains data and codes for replicating the analysis in
Casarin, R., Iacopini, M., Molina, G. Ter Horst, E., Espinasa, R., Sucre, C., and Rigobon, R. (2020), _Multilayer Network Analysis of Oil Linkages_, Econometrics Journal (23):2, 269-296 (https://doi.org/10.1093/ectj/utaa003)


### Main codes description

* OilNetworkRep.m
  <br> This code allows for the extraction of the contemporaneous and lagged multi-layer networks in Section 3 of the paper. The pictures of the networks presented in this paper have been generated with Gephi, available at: https://gephi.org/

* OilNetworkSequential.m
  <br> This code allows for the sequential extraction of the contemporaneous and lagged multi-layer networks in Section 4 of the paper. Figures of Section 4 can be generated with this code.

* OilNetworkSequentialParallelRep.m
  <br> This code is a parallel implementation of the sequential estimation code OilNetworkSequential.m


### Subfolders description

* data
  <br> This folder contains the dataset used (AllData.mat) and the labels of the series (ListVar.mat).

* functions
  <br> This folder contains all functions used by OilNetworkSequential.m and OilNetwork.m
  * adj2edgeL.m                <br> converts an adjacency matrix into an edge list
  * CONVERGENCE.m              <br> implements convergence diagnostics and statistics
  * Estimate_BIC.m             <br> evaluates the BIC
  * Gibbs.m                    <br> runs the Gibbs (used only by OilNetworkSequentialParallelRep.m)
  * LOG_SCORE.m                <br> evaluates the log-score of the nework
  * PROC_DATA.m                <br> applies data transformation
  * SAMPLE_BGMAR_DAG.m         <br> samplse the lagged DAG          (global sampler)
  * SAMPLE_BGMIN_DAG.m         <br> samples the contemporeanous DAG (global sampler)
  * SAMPLE_BGMAR_DAGblocks.m   <br> samples the lagged DAG          (blocked sampler, Section 2 of the paper)
  * SAMPLE_BGMIN_DAGblocks.m   <br> samples the contemporeanous DAG (blocked sampler, Section 2 of the paper)
  * smth.m                     <br> smooths a time series (used by OilNetworkSequentialRep.m)
  * transform.m                <br> transforms raw time into stationary series

* results
  <br> This folder contains some results obtained with the codes described in this document, including the .mat files with the MCMC output and the .csv files to be used as inputs by Gephi.
