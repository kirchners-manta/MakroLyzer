import networkx as nx
from scipy.spatial import cKDTree
import numpy as np
from IPython.display import display
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import csv

from PolyLyzer import dictionaries

class GraphManager(nx.Graph):
    def __init__(self, data=None, pbc=None):
        super().__init__()
        if data is None:
            # Handle the case where no data is provided
            pass
        elif isinstance(data, GraphManager):
            # Copy constructor from another GraphManager instance
            self.add_nodes_from(data.nodes(data=True))
            self.add_edges_from(data.edges(data=True))
        elif isinstance(data, nx.Graph):
            # Initialize from a NetworkX graph object
            self.add_nodes_from(data.nodes(data=True))
            self.add_edges_from(data.edges(data=True))
        else:
            # Handle other types of data, such as initialization from raw data
            self.create_graph(data, pbc)


    def create_graph(self, atomData, pbc=None):
        """
        Create a graph from the atom data - Part of the constructor.
        """
        covalentRadii = dictionaries.dictCovalent()
        bondDistances = {(e1, e2): (covalentRadii[e1] + covalentRadii[e2]) * 1.15 
                         for e1 in covalentRadii for e2 in covalentRadii}
        elements = atomData['atom']
        coordinates = atomData[['x', 'y', 'z']]
        tree = cKDTree(coordinates.values)

        # Add nodes with atom labels and indices
        for index, row in atomData.iterrows():
            self.add_node(index, index=index, element=row['atom'], x=row['x'], y=row['y'], z=row['z'])
 
        # Find pairs of atoms within the max possible bond distance
        maxBondDistance = max(bondDistances.values())
        possibleBonds = tree.query_pairs(maxBondDistance)
        for i, j in possibleBonds:
            element_pair = (elements[i], elements[j])
            maxDistance = bondDistances.get(element_pair, float('inf'))
            distance = np.linalg.norm(coordinates.iloc[i] - coordinates.iloc[j])
            # PBC
            if pbc is not None:
                xdist = coordinates.iloc[i]['x'] - coordinates.iloc[j]['x']
                ydist = coordinates.iloc[i]['y'] - coordinates.iloc[j]['y']
                zdist = coordinates.iloc[i]['z'] - coordinates.iloc[j]['z']
                # Apply periodic boundary conditions
                x_distance = xdist - pbc[0] * np.round((xdist) / pbc[0])
                y_distance = ydist - pbc[1] * np.round((ydist) / pbc[1])
                z_distance = zdist - pbc[2] * np.round((zdist) / pbc[2])
                distance = np.sqrt(x_distance**2 + y_distance**2 + z_distance**2)
                
            if distance <= maxDistance:
                self.add_edge(i, j)

        # Add bond order
        for node in self.nodes:
            self.nodes[node]['degree'] = self.degree[node]
    

    def remove_1order(self):
        """
        Remove 1-order nodes from the graph.

        Returns:
            GraphManager: A new GraphManager instance with 1-order nodes removed.
        """
        newGraph = GraphManager(self)
        for node in list(newGraph.nodes):
            if newGraph.nodes[node]['degree'] == 1:
                newGraph.remove_node(node)
        return newGraph
    
    def update_degree(self):
        """
        Update the degree of each node in the graph.

        Returns:
            GraphManager: The updated graph with new degrees.
        """
        for node in self.nodes:
            self.nodes[node]['degree'] = self.degree[node]
        return self
    
    
    def surrounding(self):
        """
        Find the surrounding atoms of the graph.

        Returns:
            Modification of the graph: Adding a new attribute 'surroundingAtoms' to the graph.
        """
        for node in self.nodes:
            # Get the original atom label
            atomLabel = self.nodes[node]['element']
            
            # Get the neighbors of the atom
            neighbors = list(self.neighbors(node))
            neighborLabels = [self.nodes[neighbor]['element'] for neighbor in neighbors]
            # sort the labels alphabetically
            neighborLabels.sort()
            
            # Additional label
            addLabel = atomLabel + '_' + ''.join(neighborLabels)
            
            # Add additional label
            self.nodes[node]['surroundingAtoms'] = addLabel    
            
        return self
    
    def get_subgraphs(self):
        """
        Get the connected components of the graph.

        Returns:
            list of GraphManager: A list of GraphManager instances representing the connected components.
        """
        subgraphs = []
        for component in nx.connected_components(self):
            subgraph = GraphManager(self.subgraph(component))
            subgraphs.append(subgraph)
        return subgraphs
    
    
    def find_longest_path(self, startAtom=None):
        """
        Find the longest path in the graph using Depth First Search (DFS).

        Args:
            GraphManager: The graph to search.

        Returns:
            list of nodes: The longest path found in the graph.
        """
        def dfs(current, path, visited):
            visited.add(current)
            max_path = path

            for neighbor in self.neighbors(current):
                if neighbor not in visited:
                    new_path = dfs(neighbor, path + [neighbor], visited)
                    if len(new_path) > len(max_path):
                        max_path = new_path

            visited.remove(current)
            return max_path

        longest = []
        # Do dfs for each node with degree == 1
        for node in self.nodes:
            if self.degree[node] != 1:
                continue
            # If the type of the start atom is provided, the longest path is
            # searched only from atoms of this type
            if startAtom is not None and self.nodes[node]['element'] != startAtom:
                continue
            # Perform DFS from the current node
            path = dfs(node, [node], set())
            if len(path) > len(longest):
                longest = path
                
        return longest
    
    def distance(self, node1, node2):
        """
        Calculate the distance between two nodes in the graph.

        Args:
            node1 (int): The index of the first node.
            node2 (int): The index of the second node.

        Returns:
            float: The distance between the two nodes.
        """
        coords1 = np.array([self.nodes[node1]['x'], self.nodes[node1]['y'], self.nodes[node1]['z']])
        coords2 = np.array([self.nodes[node2]['x'], self.nodes[node2]['y'], self.nodes[node2]['z']])
        return np.linalg.norm(coords1 - coords2)
    
    def vector(self, node1, node2):
        """
        Calculate the vector between two nodes in the graph.

        Args:
            node1 (int): The index of the first node.
            node2 (int): The index of the second node.

        Returns:
            np.ndarray: The vector from node1 to node2.
        """
        coords1 = np.array([self.nodes[node1]['x'], self.nodes[node1]['y'], self.nodes[node1]['z']])
        coords2 = np.array([self.nodes[node2]['x'], self.nodes[node2]['y'], self.nodes[node2]['z']])
        return coords2 - coords1
    
    def dihedral(self, node1, node2, node3, node4):
        """
        Calculate the dihedral angle between four nodes in the graph.

        Args:
            node1 (int): The index of the first node.
            node2 (int): The index of the second node.
            node3 (int): The index of the third node.
            node4 (int): The index of the fourth node.

        Returns:
            float: The dihedral angle in degrees.
        """
        
        # Get vectors
        b1 = self.vector(node1, node2)
        b2 = self.vector(node2, node3)
        b3 = self.vector(node3, node4)
        
        # Calculate normal vectors
        n1 = np.cross(-b1, b2)
        n2 = np.cross(-b2, b3)
        
        # Calculate the angle between the normal vectors which is 
        # identical to the angle between the two planes (dihedral angle)
        
        # absolute value of the angle
        cos_angle = np.dot(n1, n2) / (np.linalg.norm(n1) * np.linalg.norm(n2))
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        dihedral = np.degrees(np.arccos(cos_angle))
        
        return dihedral
    
    def get_dihedrals(self):
        # Remove 1-order nodes, find subgraphs and surrounding atoms
        GraphWithout1order = self.remove_1order()
        GraphWithout1order.surrounding()
        GraphWithout1order.update_degree()
        # Get the subgraphs of the graph
        subgraphs = GraphWithout1order.get_subgraphs()
        dihedrals = []
        # For each subgraph, find the longest path
        for subgraph in subgraphs:
            longestPath = subgraph.find_longest_path()
            # For each node in the longest path, find the dihedral angles
            for i in range(len(longestPath) - 3):
                node1 = longestPath[i]
                node2 = longestPath[i + 1]
                node3 = longestPath[i + 2]
                node4 = longestPath[i + 3]
                dihedral = subgraph.dihedral(node1, node2, node3, node4)
                # round dihedral to integers
                dihedral = round(dihedral)
                dihedrals.append((dihedral))
                
        # Group dihedrals 
        Dihedrals = dict(sorted({x: dihedrals.count(x) for x in set(dihedrals)}.items(), key=lambda item: item[0]))
        # Convert the counts to a list of tuples
        dihedrals = [(k, v) for k, v in Dihedrals.items()]
        # Sort the dihedrals by size
        dihedrals = sorted(dihedrals, key=lambda x: x[0])
        # Print dihedrals to .cvs file
        with open('dihedrals.csv', 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=';')
            writer.writerow(['Dihedral', 'Count'])
            for dihedral in dihedrals:
                writer.writerow([dihedral[0], dihedral[1]])
    
    def find_and_tag_patterns(self, patterns, startAtom=None):
        """
        Find and tag specific patterns in the graph.

        Args:
            patterns (list of tuples): List of tuples representing the patterns to find.

        Returns:
            list of tuples: List of tuples representing the found patterns.
        """
        # Remove 1-order nodes, find subgraphs and surrounding atoms
        GraphWithout1order = self.remove_1order()
        GraphWithout1order.surrounding()
        GraphWithout1order.update_degree()
        subgraphs = GraphWithout1order.get_subgraphs()
        
        # Dictionary to store the fragment IDs
        fragmentIDs = {}
        fragmentID = 0
        
        
        # Iterate through the subgraphs and find the patterns
        for subgraph in subgraphs:
            longestPath = subgraph.find_longest_path(startAtom)

            # Add side chains to the longest path
            for node in longestPath:
                for neighbor in GraphWithout1order.neighbors(node):
                    if neighbor not in longestPath:
                        longestPath.insert(longestPath.index(node) + 1, neighbor)
                    
            # Create a list of labels for the longest path
            labels = [subgraph.nodes[node]['surroundingAtoms'] for node in longestPath]
            
            # Slide over the longest path and check for patterns
            for pattern in patterns:
                i = 0
                while i <= len(labels) - len(pattern):
                    # Check if the current window matches the pattern
                    if labels[i:i + len(pattern)] == pattern:
                        # Check if one of the nodes in the window is already tagged
                        if any(longestPath[i + j] in fragmentIDs for j in range(len(pattern))):
                            # If so, skip this window
                            i += 1
                            continue
                        
                        # Tag the nodes with the current fragment ID 
                        for j in range(len(pattern)):                       
                            fragmentIDs[longestPath[i + j]] = fragmentID
                        fragmentID += 1
                        i += len(pattern)
                    else:
                        i += 1
                        
            # Add the fragment IDs to the nodes in self
            for node in longestPath:
                if node in fragmentIDs:
                    self.nodes[node]['fragmentID'] = fragmentIDs[node]
                    self.fragents_for_neighbor(longestPath, node)
                    # Neighbors get the same fragment ID
                    #for neighbor in self.neighbors(node):
                    #    if self.degree[neighbor] == 1:
                    #        self.nodes[neighbor]['fragmentID'] = fragmentIDs[node]
                else:
                    self.nodes[node]['fragmentID'] = -1
                    
                    
    def fragents_for_neighbor(self, longestPath, node):
        """
        Assign fragment IDs to the neighbors of the given node.

        Args:
            node: The index of the node to assign fragment IDs to.

        Returns:
            None
        """
        # Get the fragment ID of the node
        fragmentID = self.nodes[node]['fragmentID']
        
        # Iterate over the neighbors of the node
        for neighbor in self.neighbors(node):
            # If the neighbor has a fragment ID of -1, or no fragment ID,
            # assign the fragment ID of the node to the neighbor
            if  'fragmentID' not in self.nodes[neighbor] and neighbor not in longestPath:
                self.nodes[neighbor]['fragmentID'] = fragmentID
                # If neighbor is not degree 1, recursively call the function
                if self.degree[neighbor] != 1:
                    self.fragents_for_neighbor(longestPath, neighbor)
                    
    def saturate(self):
        """
        Saturate the Polymers by adding atoms at the ends of the longest path.

        Returns:
            GraphManager: A new GraphManager instance with saturated Polymers.
        """
        GraphWithout1order = self.remove_1order()
        GraphWithout1order.surrounding()
        GraphWithout1order.update_degree()
        subgraphs = GraphWithout1order.get_subgraphs()
        
        # Iterate through the subgraphs and find the patterns
        for subgraph in subgraphs:
            longestPath = subgraph.find_longest_path()
            
            # find nodes in longest path with degree 1
            start_nodes = [node for node in longestPath if GraphWithout1order.nodes[node]['degree'] == 1]
            
            # Add OH to node if C atom
            for startNode in start_nodes:
                if startNode in self.nodes and self.nodes[startNode]['element'] == 'C':
                    # Add OH to the start node
                    self.add_OH_to_carboxylic_acid(startNode)
                elif startNode in self.nodes and self.nodes[startNode]['element'] == 'N':
                    # Add H to the start node
                    self.add_H_to_amide(startNode)
        # return the modified graph
        return self
            
    def add_OH_to_carboxylic_acid(self, node):
        """
        Add OH to a carboxylic acid node (assumed to be the carbon in -COOH).
        """
        import numpy as np

        coords = np.array([self.nodes[node]['x'], self.nodes[node]['y'], self.nodes[node]['z']])
        
        neighbor_C = None
        neighbor_O = None

        # Find the neighbors of the C atom
        # -> C and O
        for neighbor in self.neighbors(node):
            element = self.nodes[neighbor]['element']
            pos = np.array([self.nodes[neighbor]['x'], self.nodes[neighbor]['y'], self.nodes[neighbor]['z']])
            if element == 'C':
                neighbor_C = pos
            elif element == 'O':
                neighbor_O = pos
        if neighbor_C is None or neighbor_O is None:
            raise ValueError("Carboxylic acid carbon must have one C and one O neighbor.")
        
        # calculate the coordinates of the new O atom
        coords_mid = (neighbor_O - neighbor_C) / 2 + neighbor_C 
        length_OC_single = 1.34
        vecOH = (coords_mid - coords)/np.linalg.norm(coords_mid - coords) * length_OC_single
        coords_OH = coords - vecOH
        
        # Add the new oxygen atom
        new_index = max(self.nodes) + 1
        self.add_node(new_index, index=new_index, element='O', x=coords_OH[0], y=coords_OH[1], z=coords_OH[2])
        self.add_edge(node, new_index)
        
        # calculate the coordinates of the new H atom
        length_OH = 0.96
        vec_CC = (coords - neighbor_C)/np.linalg.norm(coords - neighbor_C)
        coords_ = coords_OH + vec_CC * length_OH
        
        # Add the new hydrogen atom
        new_index = max(self.nodes) + 1
        self.add_node(new_index, index=new_index, element='H', x=coords_[0], y=coords_[1], z=coords_[2])
        self.add_edge(new_index, new_index - 1)
        
    def add_H_to_amide(self, node):
        """
        Add H to an amide node (assumed to be the nitrogen in -NH2).        
        """
        coords = np.array([self.nodes[node]['x'], self.nodes[node]['y'], self.nodes[node]['z']])
        neighbor_C = None
        
        # Find the neighbors of the N atom
        # -> C, H
        for neighbor in self.neighbors(node):
            element = self.nodes[neighbor]['element']
            pos = np.array([self.nodes[neighbor]['x'], self.nodes[neighbor]['y'], self.nodes[neighbor]['z']])
            if element == 'C':
                # remember the C atom
                Cnode = neighbor
                neighbor_C = pos
            elif element == 'H':
                neighbor_H = neighbor
        if neighbor_C is None:
            raise ValueError("Amide nitrogen must have one C neighbor.")
        
        # Get the two H atoms coordinated to the N atom
        C_Hnodes = []
        for neighbor in self.neighbors(Cnode):
            element = self.nodes[neighbor]['element']
            pos = np.array([self.nodes[neighbor]['x'], self.nodes[neighbor]['y'], self.nodes[neighbor]['z']])
            if element == 'H':
                # remember the H atom
                C_Hnodes.append(neighbor)
        if len(C_Hnodes) != 2:
            raise ValueError("Amide carbon must have two H neighbors.")
        
        # coords of the H atoms
        coords_H1 = np.array([self.nodes[C_Hnodes[0]]['x'], self.nodes[C_Hnodes[0]]['y'], self.nodes[C_Hnodes[0]]['z']])
        coords_H2 = np.array([self.nodes[C_Hnodes[1]]['x'], self.nodes[C_Hnodes[1]]['y'], self.nodes[C_Hnodes[1]]['z']])
        
        # get C-H vectors
        vec_C_H1 = (neighbor_C - coords_H1)/np.linalg.norm(neighbor_C - coords_H1)
        vec_C_H2 = (neighbor_C - coords_H2)/np.linalg.norm(neighbor_C - coords_H2)
        
        # calculate the coordinates of the new H atoms
        length_NH = 1.02
        coords_new_H1 = coords + vec_C_H1 * length_NH
        coords_new_H2 = coords + vec_C_H2 * length_NH
        
        # Update the coordinates of the neighbor H atom
        self.nodes[neighbor_H]['x'] = coords_new_H1[0]
        self.nodes[neighbor_H]['y'] = coords_new_H1[1]
        self.nodes[neighbor_H]['z'] = coords_new_H1[2]
        # Add the new hydrogen atom
        new_index = max(self.nodes) + 1
        self.add_node(new_index, index=new_index, element='H', x=coords_new_H2[0], y=coords_new_H2[1], z=coords_new_H2[2])
        self.add_edge(node, new_index)
        
        #if np.linalg.norm(coords_new_H1 - neighbor_H) > 0.8:
        #    new_index = max(self.nodes) + 1
        #    self.add_node(new_index, index=new_index, element='H', x=coords_new_H1[0], y=coords_new_H1[1], z=coords_new_H1[2])
        #    self.add_edge(node, new_index)
        #elif np.linalg.norm(coords_new_H2 - neighbor_H) > 0.8:
        #    new_index = max(self.nodes) + 1
        #    self.add_node(new_index, index=new_index, element='H', x=coords_new_H2[0], y=coords_new_H2[1], z=coords_new_H2[2])
        #    self.add_edge(node, new_index)
        #else:
        #    raise ValueError("No new H atom can be added to the amide nitrogen.")
        
    def get_chemicalFormulas(self):
        """
        Get the chemical formulas of the polymers.

        Returns:
            dict: A dictionary with the chemical formulas of the polymers and their counts.
        """
        # Get the subgraphs of the graph
        subgraphs = self.get_subgraphs()
        formulas = {}
        
        # Iterate through the subgraphs and count the elements
        for subgraph in subgraphs:
            formula = {}
            for node in subgraph.nodes:
                element = subgraph.nodes[node]['element']
                if element not in formula:
                    formula[element] = 0
                formula[element] += 1
            
            # sort alphabetically
            formula = dict(sorted(formula.items()))
            # Convert the formula to a string
            formula_str = ''.join([f"{k}{v}" for k, v in formula.items()])
            formulas[formula_str] = formulas.get(formula_str, 0) + 1
            
        # Sort the formulas by their counts
        formulas = dict(sorted(formulas.items(), key=lambda item: item[1], reverse=True))   
        # Convert the counts to a list of tuples
        formulas = [(k, v) for k, v in formulas.items()]
        
        return formulas
    
                
    def draw_graph(self):
        """Draws the largest connected component of the graph."""
        subgraphs = self.get_subgraphs()
        largestSubgraph = max(subgraphs, key=len)
        
        # 3D visualization
        pos = {node: (self.nodes[node]['x'], self.nodes[node]['y'], self.nodes[node]['z']) for node in largestSubgraph.nodes}
        
        # Scatter plot
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(*zip(*pos.values()), c='b', marker='o')
        
        # Draw edges
        for edge in largestSubgraph.edges:
            x = [self.nodes[edge[0]]['x'], self.nodes[edge[1]]['x']]
            y = [self.nodes[edge[0]]['y'], self.nodes[edge[1]]['y']]
            z = [self.nodes[edge[0]]['z'], self.nodes[edge[1]]['z']]
            ax.plot(x, y, z, color='gray', alpha=0.5)
        
        # Add labels and indices
        for node, (x, y, z) in pos.items():
            label = f"{self.nodes[node]['element']}_{node}"
            label += '\n' + str(self.nodes[node]['surroundingAtoms'])
            ax.text(x, y, z, label, color='black', fontsize=6)
                
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        
        plt.title('3D Scatter Plot of Graph Nodes')
        plt.show()
        
    def write_xyz(self, filename='saturatedGraph.xyz'):
        """
        Write the graph data to an XYZ file.
        """
        with open(filename, 'w') as f:
            f.write(f"{len(self.nodes)}\n")
            f.write("Graph data\n")
            for node in self.nodes:
                element = self.nodes[node]['element']
                x = self.nodes[node]['x']
                y = self.nodes[node]['y']
                z = self.nodes[node]['z']
                f.write(f"{element} {x} {y} {z}\n")


    def write_fragment_data_to_csv(self, filename):
        """
        Write the fragment data to a CSV file.

        Args:
            filename (str): The name of the output CSV file.
        """
        with open(filename, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=';')
            writer.writerow(['Index', 'Element', 'X', 'Y', 'Z', 'FragmentID'])

            for node in self.nodes:
                data = self.nodes[node]
                writer.writerow([
                    data['index'],
                    data['element'],
                    data['x'],
                    data['y'],
                    data['z'],
                    data['fragmentID']
                ])

