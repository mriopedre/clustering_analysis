# %%
#Example Usage
fol = '/wrk/HS-GI-S6/'

cluster = Clustering_analysis(fol+'04-HS-GI-S6-md.tpr', fol+'04-HS-GI-S6-md.xtc')

sel = cluster.u.select_atoms('resname AGLCN AIDOA')
sel = cluster.select_fragments(sel)

cluster_data = cluster.clustering_analysis(sel, cutoff = 3.5)

cluster.save_clusters_to_xvg(cluster_data, output_file=fol+'fragments_data.xvg')
