#Initialize the class with the info of the trajectory
#Use tpr (or any topology with bond information) and trajectory files
fol = '/wrk/HS-GI-S6/'
cluster = Clustering_analysis(fol+'04-HS-GI-S6-md.tpr', fol+'04-HS-GI-S6-md.xtc')

#Make the selection of the fragments - can be done by doing a selection of the residues and using the select_fragments method
sel = cluster.u.select_atoms('resname AGLCN AIDOA')
sel = cluster.select_fragments(sel)

#Perform the clustering analysis
#The cutoff is the distance in Angstroms to consider a contact
#The output is a matrix with the time, size of the largest cluster, number of clusters and the clusters themselves
cluster_data = cluster.clustering_analysis(sel, cutoff = 3.5)

#Stores the data in a xvg file, ready to plot with xmgrace or matplotlib
#Only the first 3 columns are saved, the clusters are not saved. They can be used for further analysis if needed
cluster.save_clusters_to_xvg(cluster_data, output_file=fol+'fragments_data.xvg')
