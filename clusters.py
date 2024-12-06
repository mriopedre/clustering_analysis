import MDAnalysis as md
import numpy as np
from MDAnalysis.analysis import contacts
from tqdm import tqdm
from multiprocessing import Pool, cpu_count

class ClusteringAnalysis:
    def __init__(self, topology, trajectory):
        self.u = md.Universe(topology, trajectory)

    def select_fragments(self, sel):
        """Helps to select the fragments from a previous selection."""
        return [f.select_atoms('all') for f in sel.atoms.fragments]

    def _calc_contacts(self, sel1, sel2, radius=3.5):
        """Calculate contacts between 2 selections."""
        dist = contacts.distance_array(sel1.positions, sel2.positions, self.u.dimensions)
        n_contacts = contacts.contact_matrix(dist, radius).sum()
        return n_contacts

    def _check_adjacent(self, matrix, row, visited):
        """Recursively finds all connected elements in a contact matrix row."""
        indices = np.where(matrix[row] == 1)[0]
        for index in indices:
            if index not in visited:
                visited.add(index)
                self._check_adjacent(matrix, index, visited)

    def _calculate_cluster(self, matrix):
        """Calculate the largest cluster size and list of clusters."""
        largest_cluster = 0
        clusters = []
        for i in range(matrix.shape[0]):
            visited = set()
            self._check_adjacent(matrix, i, visited)
            if len(visited) > largest_cluster:
                largest_cluster = len(visited)
            clusters.append(visited)
        clusters = list(set(tuple(sorted(s)) for s in clusters))
        return largest_cluster, clusters

    def _gen_contact_matrix(self, sel, cutoff):
        """Generates a binary contact matrix for a given step indicating whether each pair of chains is in contact."""
        matrix = []
        for s1 in sel:
            cts = [self._calc_contacts(s1, s2, radius=cutoff) for s2 in sel]
            matrix.append(cts)

        if not matrix:  # Handle empty matrix
            return np.zeros((len(sel), len(sel)), dtype=int)

        matrix = np.stack(matrix, axis=0)
        return (matrix > 0).astype(int)

    def _process_frame(self, ts_index):
        """Process a single frame of the trajectory."""
        ts = self.u.trajectory[ts_index]
        matrix = self._gen_contact_matrix(self.sel, self.cutoff)
        largest_cluster, clusters = self._calculate_cluster(matrix)
        return [ts.time / 1e3, largest_cluster, len(clusters)]

    def clustering_analysis(self, sel, cutoff=3.5, n_jobs=None):
        """Perform clustering analysis using parallel processing."""
        self.sel = sel
        self.cutoff = cutoff

        n_jobs = n_jobs or cpu_count()
        frame_indices = list(range(len(self.u.trajectory)))

        with Pool(processes=n_jobs) as pool:
            results = list(tqdm(pool.imap(self._process_frame, frame_indices), total=len(frame_indices)))
        
        return np.array(results)

    def save_clusters_to_xvg(self, cluster_data, output_file='clusters_data.xvg'):
        """Save the data as an XVG file."""
        with open(output_file, 'w') as f:
            f.write("@    title \"Cluster Information\"\n")
            f.write("@    xaxis  label \"Time [ns]\"\n")
            f.write("@    s0 legend \"Largest cluster [n_mol]\"\n")
            f.write("@    s1 legend \"Number of clusters [n_mol]\"\n")
            for row in cluster_data:
                f.write(f"{row[0]:.3f} {row[1]} {row[2]}\n")


if __name__ == "__main__":
    # Example Usage
    #Initialize the class with the info of the trajectory
    #Use a topology file with bond information (e.g., `.tpr`) and a trajectory file.
    fol = 'example_folder/'
    cluster = ClusteringAnalysis(fol+'X.tpr', fol+'X.xtc')

    #Make the selection of the fragments - can be done by doing a selection of the residues and using the select_fragments method
    sel = cluster.u.select_atoms('resname A B')
    sel = cluster.select_fragments(sel)

    #Perform the clustering analysis
    #The cutoff is the distance in Angstroms to consider a contact
    #The output is a matrix with the time, size of the largest cluster, number of clusters and the clusters themselves
    cluster_data = cluster.clustering_analysis(sel, cutoff = 3.5)

    #Stores the data in a xvg file, ready to plot with xmgrace or matplotlib
    #Only the first 3 columns are saved, the clusters are not saved. They can be used for further analysis if needed
    cluster.save_clusters_to_xvg(cluster_data, output_file=fol+'fragments_data.xvg')