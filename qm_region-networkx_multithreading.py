#!/usr/bin/env -S uv run --script
#
# /// script
# dependencies = ["ase", "networkx"]
# ///

"""Finds an ideal QM region partition within an Atoms object, for QM models.

Ideal here is defined as having the lowest edge cut, that is the fewest
bonds between the QM region and the surrounding MM region.

An example script is given after the function definitions.

How the algorithm works
-----------------------
First, a trial QM region is selected by naively picking neighbours around a
seed selection until n atoms are selected. Then, atoms selected as QM are
swapped with MM atoms repeatedly. At each step, the swaps which result in the
lowest edge cuts are determined, and then one is taken randomly. This process
is repeatedly for a number of steps and then the best results are collected.
This is repeated across several threads. Finally, the results from each thread
are compared, and the best solution saved for analysis.

Further Work
------------
- Have criteria to prioritise QM regions over others with equal edge cuts.
ie. average step count from seed.
- Allow users to add atoms that are always excluded from the QM region.
- Hypothetically - bond weighting.
"""

#from ase.io import read, write
from ase.build import bulk
from ase.neighborlist import natural_cutoffs
from ase.neighborlist import NeighborList
import networkx as nx
import random
import threading

##CLASSES AND FUNCTIONS##

class Edge_minim:
    """Stores data for a network of atoms and a QM region within them.

    Parameters
    ----------
    atoms : ase.Atoms object
    seed_qm : set[int] 
        All atoms that must be in the QM region.
    n : int
        Number of intended atoms in QM region.
    max_steps : int
        Maximum number of steps in the edge minimisation loop.
    n_threads : int
        Number of threads to use in multithreading.
    G : Networkx Graph object
        A graph object equivalent of atoms.
    """
    def __init__(self, atoms, seed_qm={0,}, n=40, max_steps=5000, n_threads=7):
        self.atoms = atoms  ##PROBABLY NOT NECESSARY
        self.seed_qm = seed_qm
        self.n = n
        self.max_steps = max_steps
        ##Add atom tags for easier interpretation of view##
        count = 0
        for atom in self.atoms:
            atom.tag = count
            count += 1
        ##Generate the Networkx Graph object##
        cutoffs = natural_cutoffs(atoms, mult=0.65)  #0.65 Ti--O, 0.8
        Neighborlist = NeighborList(cutoffs,
                                    self_interaction=False,
                                    bothways=True)
        Neighborlist.update(atoms)
        Matrix = Neighborlist.get_connectivity_matrix(sparse=False)
        self.G = nx.Graph(incoming_graph_data=Matrix)
        ##Run the edge minimising algorithm##
        self.best_qm_regions = []
        ##Multithreading##
        threads = []
        for _ in range(n_threads):
            threads.append(threading.Thread(target=self.find_best_qm_regions))
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()
        ##Find the lowest edge cuts among multithreading results##
        self.best_edge_cut = -1
        for i in self.best_qm_regions:
            if self.best_edge_cut == -1 or i[0] < self.best_edge_cut:
                self.best_edge_cut = i[0]
                self.best_qm_region = i[1]

    def find_best_qm_regions(self):
        """Function to be multithreaded. Finds QM regions with low edge cuts"""
        M = __Minimiser__(self.G, self.seed_qm, self.n, self.max_steps)
        for region in M.best_qm_regions:
            self.best_qm_regions.append((M.best_edge_cut, region))
        return

class __Minimiser__():
    """Copies data from Edge_minim, to run a single random search.
    Uses Networkx to optimise the QM region for minimum QM-MM bonds.
    Note that the optimised QM region does not maintain atom ratios.
    
    Parameters
    ----------
    n_qm_regions : int
        How many variants of the ideal QM region are saved.
    qm : set[int]
        All atoms in the QM region.
    valen : dict
        keys : int
            Atom id.
        values : int
            Valency - Number of QM neighbours minus the number of MM neighbours
            of that atom.
    best_qm_regions : list[set[int]]
        A list of up to 3 QM regions that have the lowest found edge cut.
    """
    def __init__(self, G, seed_qm, n, max_steps, n_qm_regions=3):
        ##User defined variables: n_qm_regions##
        self.n_qm_regions = n_qm_regions
        ##Inherited variables: G, seed_qm, n, max_steps##
        self.G, self.seed_qm, self.n, self.max_steps = G, seed_qm, n, max_steps
        ##Self contained variables: qm, edge_cut, best_edge_cut, valen##
        self.qm = self.seed_qm.copy()
        while len(self.qm) < n:
            for neighbor in nx.node_boundary(self.G, self.qm):
                if len(self.qm) < n:
                    self.qm.add(neighbor)
                else:
                    break
        self.edge_cut = nx.cut_size(self.G, self.qm)
        self.best_edge_cut = self.edge_cut
        self.valen = {}
        for s in self.qm:
            self.valen[s] = self.valency(s)
        for s in nx.node_boundary(self.G, self.qm):
            self.valen[s] = self.valency(s)
        ##Run the minimisation algorithm##
        self.best_qm_regions = self.minimise_edge_boundary()

    def valency(self, s):
        """Returns the number of QM neighbours minus the number of
        MM neighbours of site s."""
        qm_bonds = nx.cut_size(self.G, {s,}, self.qm)
        mm_bonds = self.G.degree[s] - qm_bonds
        return qm_bonds - mm_bonds  #2*qm_bonds - self.G.degree[s]

    def minimise_edge_boundary(self):
        count = 0
        best_qm_regions = [self.qm.copy()]
        excluded = [-1] * 2  #2 per step of prevention
        a, h = self.next_swap(set(excluded))
        while a >= 0 and count < self.max_steps:  #a == -1 if no good steps
            count += 1
            self.swap(a, h)
            ##If the new QM region has the best edge_cut, save it##
            if self.edge_cut < self.best_edge_cut:
                best_qm_regions = [self.qm.copy()]
                self.best_edge_cut = self.edge_cut
            elif self.edge_cut == self.best_edge_cut:
                if (len(best_qm_regions)<self.n_qm_regions and
                    self.is_unique_graph(self.qm.copy(), best_qm_regions)):
                    best_qm_regions.append(self.qm.copy())
            ##Prevent a or h from being used in the next swap##
            for i in [a, h]:
                excluded.pop(0)
                excluded.append(i)
            ##Generate the next swap##
            a, h = self.next_swap(set(excluded))
        return best_qm_regions
    
    def next_swap(self, excluded):
        """Returns the best next swap atom and hole.
        Returns -1 if best swap is worse than maximum.

        excluded : set
            QM or MM atoms that are excluded from possible swaps. Either
            because they've been used in a recent swap, or they're not
            permitted to be in the QM region.
        """
        maximum = -2
        swap = []
        for a in self.qm - excluded - self.seed_qm:
            for h in nx.node_boundary(self.G, self.qm):
                if h not in excluded:
                    sq = self.swap_qual_calc(a,h)
                    if sq == maximum:
                        swap.append((a, h))
                    elif sq > maximum:
                        maximum = sq
                        swap = [(a, h)]
        if len(swap) == 0:
            return -1, -1
        ah = random.choice(swap)
        return ah[0], ah[1]
    
    def swap_qual_calc(self, a, h):
        """Returns the net change in valency moving a to h."""
        return self.valen[h] - self.valen[a] + self.adj_diff(a, h)
    
    def adj_diff(self, a, h):
        """Returns -2 if a and h are adjacent, 0 otherwise."""
        if (a, h) in self.G.edges():
            return -2
        return 0
    
    def swap(self, old, new):
        """Moves atom "old" to MM site "new" and modifies self.variables."""
        self.edge_cut -= self.swap_qual_calc(old, new)
        ##Generate neighbour sets##
        old_n = {n for n in self.G.neighbors(old)}
        new_n = {n for n in self.G.neighbors(new)}
        ##Modify which atoms are in the QM region##
        self.qm = (self.qm - {old,}) | {new,}
        ##Modify valencies of neighbour atoms##
        for s in old_n:  #Neighbours of old QM atom lose valency.
            self.valen[s] -= 2
        for s in new_n:  #Neighbours of new QM atom gain valency.
            if s in self.valen.keys():
                self.valen[s] += 2
            else:  #It's a new MM atom outside the previous range.
                self.valen[s] = self.valency(s)
        return
    
    def is_unique_graph(self, graph, list_of_graphs):
        """Returns True if graph does not have an isomorph in list_of_graphs"""
        graph |= set(x for x in nx.node_boundary(self.G, graph))
        for g in list_of_graphs:
            mm_boundary = set(x for x in nx.node_boundary(self.G, g))
            if nx.vf2pp_is_isomorphic(self.G.subgraph(graph),
                                      self.G.subgraph(g|mm_boundary)):
                return False
        return True

##USER DEFINED VARIABLES##

seed_region = {2190,
               2192, 1653, 1709, 1703, 2134, 2140,
               2196, 1651, 1645, 2184, 2185, 2191, 1699, 1705, 2131, 2137}
             #{1362,}

from ase import Atoms
#9x9x9 centre:{2190,}
a, b, c = 4.59, 4.59, 2.96
op = 0.305
positions = [[0, 0, 0], [0.5*a, 0.5*b, 0.5*c],
             [op*a, op*b, 0], [(1-op)*a, (1-op)*b, 0],
             [(0.5-op)*a, (0.5+op)*b, 0.5*c], [(0.5+op)*a, (0.5-op)*b, 0.5*c]]
bulk_atoms = Atoms('Ti2O4', positions=positions, cell=[a, b, c])
"""
#7x7x7 centre:{1362,}
a, b, c = 4.16, 4.16, 4.16
positions = [[0, 0, 0], [0.5*a, 0.5*b, 0], [0.5*a, 0, 0.5*c],
             [0, 0.5*b, 0.5*c], [0.5*a, 0.5*b, 0.5*c], [0.5*a, 0, 0],
             [0, 0.5*b, 0], [0, 0, 0.5*c]]
bulk_atoms = Atoms('Mg4O4', positions=positions, cell=[a, b, c])
"""
from ase.build import make_supercell
from ase.visualize import view
bulk_atoms = make_supercell(bulk_atoms, [[9,0,0],
                                         [0,9,0],
                                         [0,0,9]])

##MAIN BODY##

n_max = 66
for n_tot in range(n_max, n_max+1):
#for n_tot in range(len(seed_region), n_max+1):
    EM = Edge_minim(bulk_atoms, seed_qm=seed_region, n=n_tot, max_steps=500,
                    n_threads=5)
    qm = EM.best_qm_region
    print("Number of QM atoms:", n_tot)
    print("Number of QM--MM bonds:", nx.cut_size(EM.G, qm))
    print("Ratio:", nx.cut_size(EM.G, qm)/n_tot)
    ##Uncomment the following to view the final QM region using ase.##
    #del bulk_atoms[[a.index for a in bulk_atoms if a.index not in qm]]
    #view(bulk_atoms)
