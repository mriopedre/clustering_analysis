# Clustering of fragments in MD trajectory

Python code to calculate clustering of fragments (chains, molecules...) in MD trajectories.

```python

#Initialize the class with the info of the trajectory
#Use a topology file with bond information (e.g., `.tpr`) and a trajectory file.
fol = 'example_folder/'
cluster = Clustering_analysis(fol+'X.tpr', fol+'X.xtc')

#Make the selection of the fragments - can be done by doing a selection of the residues and using the select_fragments method
sel = cluster.u.select_atoms('resname A B')
sel = cluster.select_fragments(sel)

#Perform the clustering analysis
#The cutoff is the distance in Angstroms to consider a contact
#The output is a matrix with the time, size of the largest cluster and number of clusters
cluster_data = cluster.clustering_analysis(sel, cutoff = 3.5)

#Stores the data in a xvg file, ready to plot with xmgrace or matplotlib
cluster.save_clusters_to_xvg(cluster_data, output_file=fol+'fragments_data.xvg')
```
