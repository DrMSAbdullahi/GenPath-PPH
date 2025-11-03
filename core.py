"""
Summary:
    Core Python module for computing persistent path homology (PPH) on networks
    and extracting KEGG pathway adjacency matrices as part of the GenPath-PPH framework.

Classes:
    - GenPathHomology: Computes persistent path homology (PPH) from a network and associated
        data, such as gene expression or other quantitative features.
    - PathwayDataProcessor: Constructs gene interaction networks from KEGG pathway files, 
        maps genes to KEGG identifiers, and generates adjacency matrices filtered to the 
        genes present in the expression dataset for subsequent PPH analysis.

Original Author:
    Dong Chen

Created and Modified:
    2022-04-12 : 2022-09-20

Modified by:
    M.S. Abdullahi, 2024-09-21

Dependencies:
    Python 3.12.6
    numpy 2.1.1
    scipy 1.14.1
    pandas
    copy
    argparse
    sys
    socket
    Bio.KEGG.REST
    Bio.KEGG.KGML.KGML_parser
    collections
    pickle
"""

# ----------------------------
# Standard library imports
# ----------------------------
import sys                   # System-specific parameters and functions
import copy                  # For deep copy of objects if needed
import argparse              # For parsing command-line arguments
import pickle                # Saving/loading Python objects
import socket                # Setting network timeouts for KEGG queries
from collections import defaultdict  # Dictionary with default factory for adjacency and mappings

# ----------------------------
# Third-party library imports
# ----------------------------
import numpy as np           # Numerical computations and arrays
import pandas as pd          # DataFrame operations for gene expression and adjacency matrices
from scipy.stats import pearsonr  # Compute correlations
from Bio.KEGG.REST import kegg_find   # Query KEGG database
from Bio.KEGG.KGML.KGML_parser import read  # Parse KGML pathway files


class PathwayDataProcessor:
    count = 0 # Static count
    def __init__(self, kgml_file, data_gene_symbols, pathway_gene_symbols, pathway_df_symbols, priority_list):
        self.kgml_file = kgml_file
        self.data_gene_symbols = data_gene_symbols
        self.pathway_gene_symbols = pathway_gene_symbols
        self.pathway_df_symbols = pathway_df_symbols
        self.priority_list = priority_list
        
        # Load or create mappings
        self.gene_to_kegg = self.load_or_create_gene_to_kegg_mapping()
        self.kegg_to_gene = self.load_or_create_kegg_to_gene_mapping()
        
    def load_or_create_gene_to_kegg_mapping(self):
        try:
            with open('entrez_symbols_to_kegg_ids.pkl', 'rb') as f:
                return pickle.load(f)
        except (FileNotFoundError, EOFError):
            gene_to_kegg = {}
            for gene_symbol in self.data_gene_symbols:
                PathwayDataProcessor.count += 1  # Increment the counter for each gene
                print(f"Processing gene {PathwayDataProcessor.count}/{len(self.data_gene_symbols)}: {gene_symbol}")
                gene_to_kegg[gene_symbol] = self.get_kegg_ids(gene_symbol)
            with open('entrez_symbols_to_kegg_ids.pkl', 'wb') as f:
                pickle.dump(gene_to_kegg, f)
            print("All genes processed.")
            return gene_to_kegg

    def load_or_create_kegg_to_gene_mapping(self):
        try:
            with open('kegg_ids_to_entrez_symbols.pkl', 'rb') as f:
                return pickle.load(f)
        except (FileNotFoundError, EOFError):
            kegg_to_gene = defaultdict(list)
            for gene_symbol, kegg_ids in self.gene_to_kegg.items():
                for kegg_id in kegg_ids:
                    kegg_to_gene[kegg_id].append(gene_symbol)
            with open('kegg_ids_to_entrez_symbols.pkl', 'wb') as f:
                pickle.dump(kegg_to_gene, f)
            return kegg_to_gene
    
    def get_kegg_ids(self, gene_symbol, failed_genes_file="failed_genes.npy"):
        failed_genes = []
        timeouts = [None, 10, 20, 30, 40, 50, 60]  # Default timeout, then increasing timeouts
        for timeout in timeouts:
            try:
                # Set the socket timeout
                socket.setdefaulttimeout(timeout)
            
                # Attempt to retrieve KEGG IDs
                response = kegg_find("hsa", gene_symbol).read()
                if response:
                    lines = response.strip().split('\n')
                    kegg_ids = [line.split('\t')[0] for line in lines if gene_symbol.upper() in line.split('\t')[1].upper()]
                    return kegg_ids
                return []
            except Exception as e:
                print(f"Attempt with timeout {timeout} failed for {gene_symbol}: {e}")
    
        # If all attempts fail, return an empty list or handle the error
        print(f"Error retrieving KEGG IDs for {gene_symbol} after multiple attempts.")
        failed_genes.append(gene_symbol)
    
        # Save the failed genes to a .npy file
        np.save(failed_genes_file, np.array(failed_genes))
        return []

    def get_gene_symbols_if_in_pathway(self, kegg_id):
        if kegg_id in self.priority_list:
            return self.priority_list[kegg_id]
        gene_symbols = self.kegg_to_gene.get(kegg_id, [])
        filtered_symbols = [symbol for symbol in gene_symbols if symbol in self.pathway_gene_symbols]
        if filtered_symbols:
            return filtered_symbols[0]
        return None

    def expand_genes(self, gene_entry):
        return gene_entry.split()

    def get_adjacency_matrix(self):
        pathway = read(open(self.kgml_file, 'r'))
        
        genes = []
        for entry in pathway.entries.values():
            if entry.type == 'gene':
                genes.append(entry.name)
        
        unique_nodes = list(set(genes))
        
        all_genes = set()
        for node in unique_nodes:
            genes = self.expand_genes(node)
            all_genes.update(genes)

        all_genes = list(all_genes)
        
        gene_to_gene_relations = defaultdict(list)

        for relation in pathway.relations:
            entry1 = relation.entry1
            entry2 = relation.entry2
        
            if entry1.type == 'gene' and entry2.type == 'gene':
                genes1 = self.expand_genes(entry1.name)
                genes2 = self.expand_genes(entry2.name)
            
                for gene1 in genes1:
                    for gene2 in genes2:
                        if gene1 != gene2:
                            subtypes = [subtype for subtype in relation.subtypes]
                            gene_to_gene_relations[(gene1, gene2)].extend(subtypes)
        
        all_genes = sorted(all_genes)
        gene_index = {gene: idx for idx, gene in enumerate(all_genes)}

        adjacency_matrix = np.zeros((len(all_genes), len(all_genes)), dtype=int)

        for (gene1, gene2) in gene_to_gene_relations.keys():
            i = gene_index[gene1]
            j = gene_index[gene2]
            adjacency_matrix[i, j] = 1

        adjacency_matrix_df = pd.DataFrame(adjacency_matrix, index=all_genes, columns=all_genes)

        new_labels = [self.get_gene_symbols_if_in_pathway(kegg_id) for kegg_id in adjacency_matrix_df.index]
        adjacency_matrix_df.index = new_labels
        adjacency_matrix_df.columns = new_labels

        new_labels_filtered = [label for label in new_labels if label in self.pathway_df_symbols]

        new_adjacency_matrix_df = adjacency_matrix_df.loc[
            adjacency_matrix_df.index.isin(new_labels_filtered),
            adjacency_matrix_df.columns.isin(new_labels_filtered)
        ]
 
        # Sort row and column labels alphabetically
        new_adjacency_matrix_df = new_adjacency_matrix_df.sort_index().sort_index(axis=1)
        
        return new_adjacency_matrix_df

class GenPathHomology(object):
    def __init__(self, initial_axes=None):
        self.initial_axes = initial_axes
        self.initial_vector_x = np.array([1, 0, 0])
        self.initial_vector_y = np.array([0, 1, 0])
        self.initial_vector_z = np.array([0, 0, 1])
        self.save_temp_result = False
        self.total_edges_num = None
        return None

    @staticmethod
    def vector_angle(v0, v1):
        """_summary_ Calculate the angle between vector v0 and vector v1 in degree.

            Args:
                v0 (array): n dimension vector, n >= 2
                v1 (array): n dimension vector, n >= 2

            Returns:
                angle (int): angle in degree.
        """
        v0_u = v0 / np.linalg.norm(v0)
        v1_u = v1 / np.linalg.norm(v1)
        angle = np.degrees(np.arccos(np.clip(np.dot(v0_u, v1_u), -1.0, 1.0)))
        return angle

    @staticmethod
    def remove_loops(edges):
        """_summary_Remove the loops of the digraph.
            Args:
                edges (array): shape = [n, 2]

            Returns:
                edges (array): shape = [n-m, 2], m is the number of the loops
        """
        loop_idx = []
        loop_nodes = []
        for i, e in enumerate(edges):
            if e[0] == e[1]:
                loop_idx.append(i)
                loop_nodes.append(e[0])
        if len(loop_nodes) > 0:
            print(f'Warning, loops on node {loop_nodes} were removed.')
        edges = np.delete(edges, loop_idx, axis=0)
        return edges

    @staticmethod
    def split_independent_component(edges, nodes):
        """_summary_ If the digraph is not fully connected, then splitting it into independent components.
                    Using the depth first search (DFS) algorithms to split the undirected graph.

            Args:
                edges (array): shape = [n, 2]
                nodes (array): shape = [k, ], k is the number of the whole graph

            Returns:
                all_components (list): the nodes' set of independent components
        """
        # convert into str
        node_map_idx = {node: idx for idx, node in enumerate(nodes)}

        # adjacency list of the graph
        graph = [[] for i in range(len(nodes))]
        for i, one_edge in enumerate(edges):
            u, v = one_edge
            # Assuming graph to be undirected.
            graph[node_map_idx[u]].append(v)
            graph[node_map_idx[v]].append(u)

        # components list
        all_components = []
        visited = [False for n in nodes]

        def depth_first_search(node, component):
            # marking node as visited.
            visited[node_map_idx[node]] = True

            # appending node in the component list
            component.append(node)
            # visiting neighbours of the current node
            for neighbour in graph[node_map_idx[node]]:
                # if the node is not visited then we call dfs on that node.
                if visited[node_map_idx[neighbour]] is False:
                    depth_first_search(neighbour, component)
            return None

        for i, one_node in enumerate(nodes):
            if visited[i] is False:
                component = []
                depth_first_search(one_node, component)
                all_components.append(component)

        return all_components

    @staticmethod
    def split_independent_digraph(all_components, edges):
        """_summary_: If the digraph is not fully connected, then splitting it into independent components.
                    Using the depth first search (DFS) algorithms to split the undirected graph.

            Args:
                all_components (list): the nodes' set of independent components
                edges (array): shape = [n, 2]

            Returns:
                all_digraphs (list): a list of digraphs, each digraph contains a list of edges.
        """
        all_digraphs = [[] for i in all_components]
        edges_visited = [False for i in edges]
        for i_c, component in enumerate(all_components):
            for i_e, edge in enumerate(edges):
                if (edges_visited[i_e] is False) and (edge[0] in component or edge[1] in component):
                    all_digraphs[i_c].append(edge)
                    edges_visited[i_e] = True
            if len(component) == 1 and np.shape(all_digraphs[i_c])[0] < 1:
                all_digraphs[i_c].append(component)

        return all_digraphs

    def utils_generate_allowed_paths(self, edges, max_path):
        """_summary_Generate the allowed paths of the digraph.
            Args:
                edges (array): shape = [n, 2], int or str

            Returns:
                allowed_path_str (dict): all paths from dimension 0 to dimension max_path
        """
        # digraph info
        nodes = np.unique(edges)
        nodes_num = len(nodes)
        nodes_idx_map = {node: idx for idx, node in enumerate(nodes)}

        # edges matrix, start->end = row->column
        edge_matrix = np.zeros([nodes_num, nodes_num])
        for i, edge in enumerate(edges):
            edge_matrix[nodes_idx_map[edge[0]], nodes_idx_map[edge[1]]] = 1

        # path_0 = vertex set
        allowed_path = {0: [np.array([n]) for n in nodes]}
        allowed_path_str = {0: [str(n) for n in nodes]}

        # path_(1 to max_path)
        for i in range(0, max_path+1):
            allowed_path[i+1] = []
            allowed_path_str[i+1] = []
            for path_previous in allowed_path[i]:
                for node in nodes:
                    if edge_matrix[nodes_idx_map[path_previous[-1]], nodes_idx_map[node]]:
                        new_path = np.append(path_previous, node)
                        allowed_path[i+1].append(new_path)
                        allowed_path_str[i+1].append('->'.join([str(one_node) for one_node in new_path]))
                    else:
                        continue
        return allowed_path_str

    def utils_unlimited_boundary_operator(self, allowed_path, max_path):
        """_summary_: Generate the n-th boundary matrix for mapping (n+1)-path to (n)-path.

            Args:
                path_n_1 (list): list of (n+1)-paths, each path stores in array
                n (int): n should >= 1, n is the dimension of the boundary matrix

            Returns:
                unlimited_boundary_mat (array): the matrix representation of the n-th boundary matrix
        """
        # For D_0, matrix is [0]*len(nodes)
        boundary_map_matrix = {0: np.zeros([len(allowed_path[0]), ])}
        boundary_mat_matrix_rank = {0: 0}
        allowed_path_idx_argument = {0: [1]*len(allowed_path[0])}

        for n in range(1, max_path+2):
            boundary_map_dict = {}
            boundary_operated_path_name_collect = []

            allowed_path_n_types = len(allowed_path[n])
            if allowed_path_n_types == 0:
                boundary_map_matrix[n] = np.zeros([1, len(allowed_path[n-1])])
                boundary_mat_matrix_rank[n] = 0
                allowed_path_idx_argument[n] = [1] * len(allowed_path[n-1])
                break

            for i_path, path in enumerate(allowed_path[n]):

                # split the path into nodes with idx
                path_node_idx = path.split('->')

                # record the result path after boundary operation
                boundary_operated_path_info = {}
                for i_kill in range(n+1):
                    # kill the  i_kill-th vertex
                    temp_path = np.delete(path_node_idx, i_kill)
                    temp_path_str = '->'.join([str(pp) for pp in temp_path])
                    boundary_operated_path_info[temp_path_str] = (-1)**(i_kill)

                    # record all possible n_path
                    boundary_operated_path_name_collect.append(temp_path_str)
                boundary_map_dict[path] = copy.deepcopy(boundary_operated_path_info)

            # generate the boundary matrix, D; row_p * column_p = n_1_path * n_path
            considered_operated_path_name = np.unique(boundary_operated_path_name_collect + allowed_path[n-1])
            unlimited_boundary_mat = np.zeros([allowed_path_n_types, len(considered_operated_path_name)])
            for i_path, (n_1_path_str, operated_n_path_dict) in enumerate(boundary_map_dict.items()):
                for j, n_path in enumerate(considered_operated_path_name):
                    if n_path in operated_n_path_dict:
                        unlimited_boundary_mat[i_path, j] = operated_n_path_dict[n_path]

            # collect informations
            boundary_map_matrix[n] = unlimited_boundary_mat
            boundary_mat_matrix_rank[n] = np.linalg.matrix_rank(unlimited_boundary_mat)
            allowed_path_idx_argument[n] = [1 if tpn in allowed_path[n-1] else 0 for tpn in considered_operated_path_name]

        return boundary_map_matrix, boundary_mat_matrix_rank, allowed_path_idx_argument

    def path_homology_for_connected_digraph(self, allowed_path, max_path):
        """Calculate the Betti numbers for path homology up to a given dimension.

        Args:
            allowed_path (dict): The dictionary of paths of different dimensions.
            max_path (int): Maximum dimension for homology calculation.

        Returns:
            np.array: The Betti numbers from dimension 0 to max_path.
        """
        betti_numbers = np.array([0] * (max_path + 1))

        boundary_map_matrix, boundary_mat_matrix_rank, allowed_path_idx_argument =\
            self.utils_unlimited_boundary_operator(allowed_path, max_path)

        betti_0 = len(allowed_path[0]) - boundary_mat_matrix_rank[1]
        betti_numbers[0] = betti_0

        for n in range(1, max_path + 1):
            if len(allowed_path[n]) == 0:
                break

            dim_0 = len(allowed_path[n]) - boundary_mat_matrix_rank[n]

            dim_An_Bn = np.linalg.matrix_rank(
                np.vstack([
                    np.eye(len(allowed_path_idx_argument[n+1])) * allowed_path_idx_argument[n+1],
                    boundary_map_matrix[n+1]
                ])
            )
            dim_1 = len(allowed_path[n]) + boundary_mat_matrix_rank[n+1] - dim_An_Bn
            betti_numbers[n] = dim_0 - dim_1

        return betti_numbers

    def path_homology(self, edges, nodes, max_path):
        """Calculate the Betti numbers for all connected components of a graph.

        Args:
            edges (array): Array of edges in the graph.
            nodes (array): Array of nodes in the graph.
            max_path (int): Maximum dimension for path homology.

        Returns:
            np.array: The Betti numbers summed over all connected components.
        """

        edges, nodes = edges.astype(str), nodes.astype(str)

        all_components = GenPathHomology.split_independent_component(edges, nodes)
        all_digraphs = GenPathHomology.split_independent_digraph(all_components, edges)

        betti_numbers = []
        for edges in all_digraphs:
            if np.shape(edges)[1] <= 1:
                betti_numbers.append(np.array([1] + [0] * max_path))
            else:
                edges = GenPathHomology.remove_loops(edges)
                if np.shape(edges)[0] == 0:
                    betti_numbers.append(np.array([1] + [0] * max_path))
                    continue
                allowed_path = self.utils_generate_allowed_paths(edges, max_path)
                betti_numbers.append(self.path_homology_for_connected_digraph(allowed_path, max_path))
        return np.sum(betti_numbers, axis=0)

    def persistent_path_homology(self, cloudpoints, points_weight, max_path, filtration=None, distance_type='1-abs-correlation'):
        """Calculate persistent path homology with distance-based filtration.

        Args:
            cloudpoints (array): Coordinates of the points.
            points_weight (list): Weights of the points.
            max_path (int): Maximum path length.
            filtration (array): Filtration values (optional).
            distance_type (str): Distance metric to use ('euclidean' or '1-correlation' or '1-abs-correlation').


        Returns:
            list: Betti numbers over filtration steps.
        """
        points_num = np.shape(cloudpoints)[0]

        # Compute the appropriate distance matrix based on the provided distance type
        if distance_type == 'euclidean':
            # Compute the Euclidean distance matrix
            distance_matrix = np.linalg.norm(cloudpoints[:, None, :] - cloudpoints[None, :, :], axis=2)
            # Compute the correlation-based distance matrix    
        elif distance_type in ['1-correlation', '1-abs-correlation']:
            distance_matrix = self.compute_correlation_distance_matrix(cloudpoints, dist_type=distance_type)
        else:
            raise ValueError(f"Unknown distance type: {distance_type}")
        
        # Initialize the fully connected map based on point weights
        fully_connected_map = np.zeros([points_num, points_num], dtype=int)
        for i in range(points_num):
            for j in range(points_num):
                if i == j:
                    continue
                if points_weight[i] <= points_weight[j]:
                    fully_connected_map[i, j] = 1
        
        max_distance = np.max(distance_matrix)
        self.total_edges_num = np.sum(np.abs(fully_connected_map))
        self.max_distance = max_distance

        if filtration is None:
            filtration = np.arange(0, np.max(distance_matrix) + 0.1, 0.1)

        all_betti_num = []
        snapshot_map_temp = np.ones([points_num] * 2, dtype=int)

        for snapshot_dis in filtration:
            snapshot_map = (distance_matrix <= snapshot_dis).astype(int) * fully_connected_map
            if np.array_equal(snapshot_map, snapshot_map_temp):
                all_betti_num.append(all_betti_num[-1])
                continue
            snapshot_map_temp = snapshot_map.copy()

            edges = np.column_stack(np.nonzero(snapshot_map))
            betti_numbers = self.path_homology(edges, np.arange(points_num), max_path)
            all_betti_num.append(betti_numbers)

        return all_betti_num

    def compute_correlation_distance_matrix(self, cloudpoints, dist_type='1-correlation'):
        """Compute the correlation distance matrix using 1 - corr(x, y) or 1 - |corr(x, y)|.
	    dist_type (str, optional): Correlation distance type can be either '1-correlation' or '1-abs-correlation'.
	"""
    
        points_num = cloudpoints.shape[0]
        distance_matrix = np.zeros((points_num, points_num))

        # Compute the pairwise correlation distance
        for i in range(points_num):
            for j in range(i + 1, points_num):
                corr, _ = pearsonr(cloudpoints[i], cloudpoints[j])
                if dist_type == '1-correlation':
                    corr_distance = 1 - corr  # 1 - corr(x, y)
                elif dist_type == '1-abs-correlation':
                    corr_distance = 1 - abs(corr)  # 1 - |corr(x, y)|
                else:
                    raise ValueError(f"Unknown distance type: {dist_type}")
                    
                distance_matrix[i, j] = distance_matrix[j, i] = corr_distance

        return distance_matrix

    def edges_at_specific_filtration(self, edges_at_each_filtration, filtration=None, specific_value=None):
        """Retrieve edges at a specific filtration value.

        Args:
            edges_at_each_filtration (List): List of edges over all filtration values.
            filtration (list): Filtration values used for obtaining the edges at each filtration.
            specific_value (float, optional): Filtration value to retrieve the digraph.

        Returns:
            list: Edges at the specific filtration value.
        """

        # Default filtration if not provided
        if filtration is None:
            filtration = np.arange(0, 10.1, 0.1)

        edges_at_specific_filtration = None  # To capture edges at the specific filtration value

        # Check if the current filtration value matches the specific_value
        if specific_value is not None:
            for filt_value, edge_at_filt in zip(filtration, edges_at_each_filtration):
                if np.isclose(filt_value, specific_value, atol=1e-8):
                    edges_at_specific_filtration = edge_at_filt
                    break

        return edges_at_specific_filtration

    def persistent_path_homology_from_digraph(self, cloudpoints, all_edges, target_dimension, filtration=None, distance_type='1-abs-correlation'):
        """Calculate persistent path homology for KEGG digraph with filtration.

        Args:
            cloudpoints (array): Coordinates of points.
            all_edges (array): Array of edges.
            target_dimension (int): Dimension to retrieve Betti numbers.
            filtration (array, optional): Filtration values.
            distance_type (str, optional): Distance function to be used ('euclidean' or '1-correlation' or '1-abs-correlation')

        Returns:
            list: Betti numbers at each filtration step.
            list: Edges at the each filtration value (over all filtration values).
        """
        points_num = cloudpoints.shape[0]

        # Compute the appropriate distance matrix based on the provided distance type
        if distance_type == 'euclidean':
            # Compute the Euclidean distance matrix
            distance_matrix = np.linalg.norm(cloudpoints[:, None, :] - cloudpoints[None, :, :], axis=2)
            # Compute the correlation-based distance matrix
        elif distance_type in ['1-correlation', '1-abs-correlation']:
            distance_matrix = self.compute_correlation_distance_matrix(cloudpoints, dist_type=distance_type)
        else:
            raise ValueError(f"Unknown distance type: {distance_type}")

        fully_connected_map = np.zeros([points_num, points_num], dtype=int)

        # Fill the fully_connected_map with the given edges
        for edge in all_edges:
            fully_connected_map[edge[0], edge[1]] = 1

        # Default filtration if not provided
        if filtration is None:
            filtration = np.arange(0, 10.1, 0.1)

        all_betti_num = []
        edges_at_each_filtration = []  # To hold edges at each filtration step

        # Iterate over the filtration steps
        for snapshot_dis in filtration:
            snapshot_map = (distance_matrix <= snapshot_dis).astype(int) * fully_connected_map

            # Extract the edges present in this snapshot
            current_edges = np.column_stack(np.nonzero(snapshot_map))

            # Compute Betti numbers for the current edge set
            betti_numbers = self.path_homology(np.array(current_edges), np.arange(points_num), target_dimension)
            all_betti_num.append(betti_numbers[target_dimension])

            # Store the edges for this filtration step
            edges_at_each_filtration.append(current_edges)

        return all_betti_num, edges_at_each_filtration
    
# ----------------------------------------
# Demo/test block to test PPH on a toy datasaet
# ----------------------------------------

if __name__ == "__main__":
    import argparse
    import numpy as np
    import pandas as pd
    from core import GenPathHomology
    
    parser = argparse.ArgumentParser(description="GenPath Persistent Path Homology core module")
    parser.add_argument("--test", action="store_true", help="Run a simple PPH test")
    args = parser.parse_args()
    
    if args.test:
        # Toy expression data
        # Define genes and samples
        genes = [f"g{i}" for i in range(1, 11)]
        samples = ["C1", "C2", "C3", "D1", "D2", "D3"]

        # Simulated gene expression data
        data = np.array([
            [3.224085, 4.238909, 3.573100, 2.281019, 2.904727, 3.252311],
            [2.801195, 3.629210, 3.212241, 2.832321, 3.747732, 2.882458],
            [4.568956, 3.268862, 3.689043, 2.781017, 3.834468, 3.507576],
            [3.531844, 2.599495, 2.618524, 3.107632, 3.402103, 4.700751],
            [4.695837, 4.608166, 4.677947, 3.766694, 1.826884, 4.893230],
            [4.001411, 4.770909, 4.668929, 2.451091, 4.178800, 2.988391],
            [4.279513, 4.427967, 5.506402, 2.991138, 3.061895, 2.712762],
            [4.223187, 3.225830, 4.907987, 2.766738, 2.861530, 2.355171],
            [3.765850, 3.609896, 4.185733, 1.280467, 1.476483, 2.873749],
            [4.121385, 4.167441, 2.948542, 1.734737, 3.004740, 3.177959]
            ])

        # Create DataFrame
        toy_expression_df = pd.DataFrame(data, index=genes, columns=samples)

        # Toy edges
        toy_edges_adj = np.array([
            [0,1],[1,2],[2,3],[3,0],
            [5,6],[5,7],[8,6],[8,7],
            [4,1]
        ])

        # Filtration
        toy_filtration_abs = np.arange(0, 1, 0.01)

        # Initialize GenPathHomology
        pph = GenPathHomology()

        # Compute Betti numbers for dimension 0 and 1
        betti_0, edges_0 = pph.persistent_path_homology_from_digraph(
            np.array(toy_expression_df),
            toy_edges_adj,
            target_dimension=0,
            filtration=toy_filtration_abs,
            distance_type='1-abs-correlation'
        )

        betti_1, edges_1 = pph.persistent_path_homology_from_digraph(
            np.array(toy_expression_df),
            toy_edges_adj,
            target_dimension=1,
            filtration=toy_filtration_abs,
            distance_type='1-abs-correlation'
        )

        print("Betti numbers for dimension 0 across filtration:")
        print(betti_0)
        print("\nBetti numbers for dimension 1 across filtration:")
        print(betti_1)

