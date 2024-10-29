#%%
import MDAnalysis as md
import numpy as np
from MDAnalysis.analysis import contacts
from tqdm import tqdm

class Clustering_analysis():
    def __init__(self, topology, trajectory):
        # Load the topology and trajectory
        # The topology should be in tpr form (or any format that includes bond information)

        self.u = md.Universe(topology, trajectory)

    def select_fragments(self, sel):
        """Helps to select the fragments from a previous selection.

        Args:
            sel (AtomGroup): An MDAnalysis selection object.

        Returns:
            list: A list of AtomGroups, each containing all atoms in a fragment.
        """

        return [f.select_atoms('all') for f in sel.atoms.fragments]

    def _calc_contacts(self, sel1, sel2, radius = 3.5):
        """ Calculate contacts between 2 selections.
        Radius in Angstrom
        """
        dist = contacts.distance_array(sel1.positions, sel2.positions, self.u.dimensions)
        # determine which distances <= radius
        n_contacts = contacts.contact_matrix(dist, radius).sum()
        return n_contacts

    def _check_adjacent(self, matrix, row, visited):
        """Recursively finds all connected elements in a contact matrix row.

        This function identifies all elements connected to the starting row in 
        the contact matrix. It marks each visited row in the `visited` set, 
        avoiding repeats. This is used to calculate the size of a cluster.

        Args:
            matrix (np.array): A square matrix where a value of 1 at (i, j) 
                means row i is connected to row j. 0 otherwise.
            row (int): Index of the current row to check for connections.
            visited (set): A set of row indices that have been visited, 
                which is updated as the function runs.
        """

        indices = np.where(matrix[row] == 1)[0] #Indices where there is a contact in a given row
        if len(indices) > 0: # check if there are any contacts in that row
            for index in indices:
                if index not in visited: # It stays here until it finds a row with no contacts, or where all contacts have been visited
                    visited.add(index)
                    self._check_adjacent(matrix, index, visited)


    def _calculate_cluster(self, matrix):
        """Calls the check adjacent in a recursive way and stores the information obtained from it.

        Args:
            matrix (np.array): contact matrix

        Returns:
            int, list: n of molecules in the largest cluster, and list of clusters
        """
        largest_cluster = 0
        clusters = []
        for i in range(matrix.shape[0]):
            visited = set()
            self._check_adjacent(matrix, i, visited) #Visited is updated in place
            if len(visited) > largest_cluster: largest_cluster = len(visited)
            clusters.append(visited)
        clusters = list(set(tuple(s) for s in clusters))
        
        return largest_cluster, clusters

    def _gen_contact_matrix(self, sel, cutoff):
        """Generates a binary contact matrix for a given step indicating whether each pair of chains is in contact.
        Returns:
            np.array: A binary matrix (n_sel1 x n_sel2) where each element is 1 if 
                    two mols are in contact within a cutoff radius, and 0 otherwise.
        """

        matrix = []
        for s1 in sel:
            cts = [self._calc_contacts(s1, s2, radius= cutoff) for s2 in sel] #Row of the contact matrix - contacts between s1 and all s2 in sel
            matrix.append(cts)
        matrix = np.stack(matrix, axis=0)
        matrix = (matrix > 0).astype(int)
        return matrix

    def clustering_analysis(self, sel, cutoff = 3.5):
        """ Main hub for the clustering analysis.
        It performs the clustering analysis for each frame of the trajectory and returns the data.

        Returns a matrix with:
        np.array(time[ns], size of the largest cluster [n_molecules], number of clusters [n] )"""

        return_matrix = []
        for ts in tqdm(self.u.trajectory):
            largest_cluster, clusters = self._calculate_cluster(self._gen_contact_matrix(sel, cutoff))
            return_matrix.append([ts.time/1e3, largest_cluster, len(clusters), clusters])
        return np.stack(return_matrix, axis = 0)
    
    def save_clusters_to_xvg(self, cluster_data, output_file='clusters_data.xvg'):
        """Save the data as svg file.

        Args:
            cluster_data (np.array): results from clustering_analysis
            output_file (str, optional): Name for saveing file. Defaults to 'clusters_data.xvg'.
        """
        with open(output_file, 'w') as f:
            # Write header comments for readability in xmgrace
            f.write("@    title \"Cluster Information\"\n")
            f.write("@    xaxis  label \"Time [ns]\"\n")
            f.write("@    s0 legend \"Largest cluster [n_mol]\"\n")
            f.write("@    s1 legend \"Number of clusters [n_mol]\"\n")
            
            # Write fragment data with clear formatting
            for i, row in enumerate(cluster_data):
                f.write(f"{row[0]:.3f} {row[1]} {row[2]} ")
                f.write("\n")  # Marks the end of a dataset in XVG

#Example Usage

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