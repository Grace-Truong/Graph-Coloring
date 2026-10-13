#!/usr/bin/env python
from sage.all_cmdline import *   # import sage library
import copy
import itertools
import sage.all
import geogebra_graph
import time
import multiprocessing
import argparse

class Configuration:
    """
    A class used to represent a specific configuration we are trying
    
    Attributes:
        num_colors: number of colors we are trying to use to extend the coloring
        file_name: the Geogebra file that contains our configuration
    """
    def __init__(self, num_colors, file_name):
        self.num_colors = num_colors
        self.S = set()
        self.U = set()

        self.G= geogebra_graph.geogebra_to_graph(file_name)
        
        #relabeling graph edges from the geogebra file
        frontier_color = (0,0,255,0.0)                                                      #identifying frontier edges with the color blue
        frontier_edges = [e for e in self.G.edges() if e[2]["color"] == frontier_color]     #list of frontier edges in the graph
        for e in self.G.edges():
            self.G.set_edge_label(e[0],e[1],None)                                           #erasing all the extra information from Geogebra
            
        end_points_frontiers = list(set(sum([list(e[:2])                                    #all vertices incident to the frontier edges
            for e in frontier_edges],[])))
        config_vertices = [v for v in end_points_frontiers if self.G.degree(v) > 1]         #all vertices incident to the frontier edges AND are in the boundary
        
        #establishing the cyclic order of the verticies in the configuration
        pos = self.G.get_pos()
        angle_config_vertices = sorted([(atan2(pos[v][1], pos[v][0]),v) for v in config_vertices])
        cyclic_vertices = [x[1] for x in angle_config_vertices]
        self.degrees = [-1] * len(cyclic_vertices)                                          # degree is the number of frontier edges incident to the boundary
        for i,v in enumerate(cyclic_vertices):
            frontiers_at_v = [e for e in frontier_edges if v in e]
            self.degrees[i] = len(frontiers_at_v)
            for j,e in enumerate(frontiers_at_v):
                self.G.set_edge_label(*(e[:2]), (i,j))

    def draw_graph(self):
        self.G.show(edge_labels=True)

    def does_extend(self, frontier_coloring):
        """
        Using a MIP to check if the frontier coloring can extend into a coloring for self.G
    
        Args:
            frontier_coloring: A specific coloring that we are testing -- a tuple of frozen sets ({1,},{1,2,3},{0,})
    
        Returns:
            True if the coloring extends and False otherwise
            frontier_coloring: the frontier coloring associated with the T/F value
        """
    
        frontier_coloring_tuple = tuple(tuple(i) for i in frontier_coloring)
        p = MixedIntegerLinearProgram()
        x= p.new_variable(binary=True)
        L= self.G.line_graph()
        for e in self.G.edges(): 
            p.add_constraint(p.sum(x[e,i] for i in range(self.num_colors)) == 1)                    #every edge must have one color
        
        for e,f in L.edges(labels=False):
            for i in range(self.num_colors):
                p.add_constraint(x[e,i] + x[f,i] <=1)                                               #adjacent edges must have different colors
        
        for e in self.G.edges():
            if e[2]: #is not none
                tuple_index,color_index = e[2]
                p.add_constraint(x[e, frontier_coloring_tuple[tuple_index][color_index] ] == 1)       #considers the established coloring of the frontier edges
        
        p.set_objective(0)
        try:
            p.solve()
            return frontier_coloring, True
        except Exception as e:
            if "no feasible" not in str(e).lower():
                print(f"An exception occurred: {type(e).__name__}")
                print(f"Exception message: {e}")
            return frontier_coloring, False
        
    def init_colorings(self):
        """
        Generates a list of all the possible frontier colorings and organizes them into ones that are extendable (S) or not (U)
        """    
        C=list(range(self.num_colors))                                  # list of colors, for convenience
        
        all_colorings = itertools.product(                     # making all the possible tuple of frozensets
            *(
              (frozenset(x) for x in itertools.combinations(C,d))
              for d in self.degrees
             )
            )
        with multiprocessing.Pool() as pool:
            results = pool.map(self.does_extend, all_colorings, chunksize = 100)
        
        for coloring, does_extend in results:
            if does_extend:
                self.S.add(coloring)
            else:
                self.U.add(coloring)
            
        print(len(self.S), len(self.U))
    
    def create_auxiliary(self, frontier_coloring, i, j):
        """
        Creating the auxiliary graph for a specific pair of colors for a specific frontier coloring
    
        Args:
            frontier_coloring: A specific coloring that we are testing -- a tuple of frozen sets ({1,},{1,2,3},{0,})
            i, j : pair of colors
    
        Returns:
            vertex_index: the set indices from the tuple of the frontier coloring that represents the vertices in the auxiliary graph 
            aux_edges: list of hyperedges in the auxiliary
        """
        aux_edges = list()
        vertex_index = set()
        couples = set()
        
        for k, colors_incident_k in enumerate(frontier_coloring):   #creates set of verticies that will be in aux and creates set of couples
            degree = len(colors_incident_k & {i,j})                 #how many times i and j are at vertex k
            if degree >= 1:
                vertex_index.add(k)
                if degree == 2:
                    couples.add(k)
        singles = vertex_index - couples
        
        for start, end in itertools.combinations_with_replacement(vertex_index, 2): #making all the Kempe chains and testing if they extend
            test_coloring_new = tuple((colors_incident_k ^ {i,j})
                if ((k==start or k == end) and (k in singles)) 
                else colors_incident_k
                for k, colors_incident_k in enumerate(frontier_coloring))
                
            if test_coloring_new not in self.S:
                if start != end:  
                    for middle in powerset(couples-{start,end}):
                        kempe = set(middle) | {start,end}
                        aux_edges.append(kempe)
                else:
                    aux_edges.append({start})
        return vertex_index, aux_edges
        
        
    def are_Crossing(self, edge1, edge2): 
        """
        Checks if two edges in the auxiliary are crossing or not
        
        Args:
            edge1: the first hyperedge/Kempe chain
            edge2: the second hyperedge/Kempe chain
            
        Returns:
            True if the edges are crossing and False if not
        """
        is_crossing = True
        if len(edge1) == 1 or len(edge2)==1:        #trailing off Kempe chains are not crossing
            is_crossing = False
        else:
            minimum = min(min(edge1), min(edge2))   #edges of Kempe are in order so just need to test the 1st values for the min
            if minimum in edge1:
                start_edge = copy.deepcopy(edge1)
                test_edge = copy.deepcopy(edge2)
            else:
                start_edge = edge2
                test_edge = edge1
            maximum = max(max(test_edge), max(start_edge))
            start_edge = start_edge|{maximum}
            start_range = list(start_edge)
            for k in range(len(start_range)-1):
                test_range = set(range(start_range[k], start_range[k+1]+1))
                intersection = test_range & test_edge
                if intersection == test_edge:
                    is_crossing = False
        return is_crossing
        
    def has_NCPM(self, aux_vertices, aux_edges):
        """
        Uses a MIP to check to see if there is a non-crossing perfect matching in the auxiliary graph
        
        Args: 
            aux_vertices: the vertices that are in the auxiliary
            aux_edges: the edges that are in the auxiliarty
            
        Return:
            True if there is a matching and false if not
        """
        if not aux_vertices:                    #if empty
            return True
        p = MixedIntegerLinearProgram()
        x= p.new_variable(binary=True)
        matching = set()
    
        for i,j in itertools.combinations(range(len(aux_edges)),2):         #checks if there are any crossing edges in the auxiliary, if there is then choose one 
            if self.are_Crossing(aux_edges[i], aux_edges[j]) == True:  
                p.add_constraint(x[i]+x[j] <= 1)
                        
        for v in aux_vertices:                                              #constraint to make sure all vertices are covered for the NCPM
            adj_edges = list()
            for j in range(len(aux_edges)):
                if v in aux_edges[j]:
                    adj_edges.append(j)
            if adj_edges:           #non-empty
                p.add_constraint(p.sum(x[j] for j in adj_edges) ==1)
            else:                   #if adj_edges is empty then there is no edges in aux that contain v and so no perfect matching
                return False
                
        p.set_objective(0)
        try:
            p.solve()
            return True
        except Exception as e:
            if "no feasible" not in str(e).lower():
                print(f"An exception occurred: {type(e).__name__}")
                print(f"Exception message: {e}")
            return False
            
    def is_good(self, frontier_coloring):
        """
        For a frontier coloring, we are creating the auxiliary and trying to find a NCPM for each possible pair of colors/
        
        Args:
            frontier_coloring: the frontier coloring that we are testing 
            
        Return:
            False if all pairs of colors has a NCPM in it's auxiliary graph
            True if a single pair of colors do not have a NCPM in the auxiliary graph
            frontier_coloring the coloring we inputted
        """
        CL = copy.deepcopy(list(frontier_coloring))
        for i, j in itertools.combinations(range(self.num_colors),2):
            aux_vertices, aux_edges= self.create_auxiliary(CL, i, j)
            if self.has_NCPM(aux_vertices,aux_edges) == False:
                return frontier_coloring, True
        return frontier_coloring, False
    
    def test_reducibility(self):
        """
        Takes the frontier colorings in U check if they are good by checking the auxiliary graph
        """
        while True:
            start = time.perf_counter()
            print(f'Starting while loop, {len(self.U)=}')
            remove_from_U = set()
            U_list = list(self.U)
            with multiprocessing.Pool() as pool:
                results = pool.map(self.is_good, U_list, chunksize = 100)
            results = map(self.is_good, U_list)

            for frontier_coloring, is_good in results:
                if is_good:
                    remove_from_U.add(frontier_coloring)
            if remove_from_U: #non-empty
                self.U = self.U - remove_from_U
                self.S = self.S.union(remove_from_U)
                if not self.U: #U is empty
                    print(f'U is empty, therefore it is reducible. YAY!')
                    end = time.perf_counter()
                    print(f"Time taken for test_reducibility loop: {end - start:.6f} seconds")
                    break
            else:
                print(f'Nothing removed from U, {len(self.U)=}')
                end = time.perf_counter()
                print(f"Time taken for test_reducibility loop: {end - start:.6f} seconds")
                break
            end = time.perf_counter()
            print(f"Time taken for test_reducibility loop: {end - start:.6f} seconds")
            
    def find_equiv_classes(self, UL):
        """
        When there are no more good frontier colorings that can be found in U we try to combine the frontier colorings that are actually equivalent
        """
        equivalence_classes=[]
        for colornum,coloring in enumerate(UL):
            for orbit in equivalence_classes:
                if coloring in orbit:
                    break
            else:
                print(f"Adding a new class, {len(equivalence_classes)=}, {coloring=}")
                orbit=set()          
                for i,pi in enumerate(itertools.permutations(range(6))):
                    image=tuple(
                           tuple(sorted([pi[z] for z in x]))
                           for x in coloring)
                    orbit.add(image)
                equivalence_classes.append(orbit)

    def is_reducible(self):
        """
        Testing if the the configuration we inputted is reducible or not.
        This is what we run 
        """
        start = time.perf_counter()
        self.init_colorings()
        end = time.perf_counter()
        print(f"Time taken for init_colorings: {end - start:.6f} seconds")
        self.test_reducibility()
        if self.U: #non-empty
            UL = sorted(list(self.U))
            #self.find_equiv_classes(UL)
        else:
            print("YAY! It's reducible")

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('filename')
    args = parser.parse_args()
    config = Configuration(6, args.filename)
    config.is_reducible()
    
    # Demonstration code
    #print(f"Total CPU cores available: {os.cpu_count()}")
    #config = Configuration(6, "chordal4cycle.ggb")
    #config.is_reducible()