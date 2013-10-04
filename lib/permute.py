
from copy import copy
from random import shuffle, sample
from random import shuffle

class NetBalancedPermuter:
    """

        Encapsulates the permutation logic for an input heat set. Permutes
        Node scores with other nodes of similar network degree by sorting
        all nodes by degree, and binning them into blocks of a set size. 
        Permutations are done only within blocks, so that the degree distribution
        of input nodes is preserved.
    """


    def __init__(self, network, up_set):
        """
            Input:
                network: net[source] = [(i, t)]
                up_set: up_set[node] = score
        """
    

        # store node-degrees for all in network
        self.degrees = {}

        # the set of initial nodes to permute within blocks: save it to the 
        # instantiated object here
        self.nodes = up_set.keys()
        # heuristic: block needs to be significantly larger than the input set size
        BLOCK_MULTIPLE = 10
        self.block_size = len(self.nodes)*BLOCK_MULTIPLE
        self.scores = {}
        for node in self.nodes:
            # save the scores as a set of tuples 
            self.scores[(node, str(up_set[node]))] = 1

        # Compute total degree for each node in the network
        for source in network:

            if source not in self.degrees:
                self.degrees[source] = 0

            for (i, target) in network[source]:    
                self.degrees[source] += 1

                if target not in self.degrees:
                    self.degrees[target] = 0
                # add a degree for the incoming edge
                self.degrees[target] += 1    

        # reverse-sort the degrees of all nodes in the network.
        self.sorted_degrees = sorted(self.degrees.items(), key=lambda x:x[1], reverse=True)

            
    def permuteBlock(self, block):
        """
        Take a block of nodes and randomly shuffle using python's random.shuffle method.      

        Input:

            An array of node labels
        
        Returns:
        
            A hash mapping the original nodes to the nodes to swap with each. 
        """
        # make a copy
        orig = copy(block)
        b = copy(block)
        map = {}
        # shuffle the copy in-place with the random.shuffle() method.
        shuffle(b)
        for i in range(0, len(b)):
            # build the mapping from original to new
            map[orig[i]] = b[i]

        return map

    def permuteOne(self):
        """
        Generate one permutation of scores for all nodes, and return a hash of { node : score }
        for each. 
        """    
        group_count = 0
        permuted_scores = {}
        # initialize a new block
        block = []
        for (node, degree) in self.sorted_degrees:
            block.append(node)    
            group_count += 1
            # reset every time we use the block size
            if group_count % self.block_size == 0:
                # permute the order of this <block_size> block
                map = self.permuteBlock(block)
                for (node, score) in self.scores:
                    if node in map:
                        permuted_scores[map[node]] = float(score)
                block = []    

        return permuted_scores    

    def permute(self, iterations):
        """
        Generate an array of random permutations of node scores.
        
        Input:
            iteration: the number of permutations to generate

        Returns:
            an array of hashes: each hash indexes the nodes to permuted scores
        """    
        permuted = []    
        for i in range(0, iterations):
            permuted.append(self.permuteOne())

        return permuted

<<<<<<< HEAD
class SupervisedPermuter:
    """
        Encapsulates the permutation logic for an input heat set. Permutes
        Node scores with a supplied background list of other nodes
        all nodes by degree, and binning them into blocks of a set size. 
        Permutations are done only within blocks, so that the degree distribution
        of input nodes is preserved.
    """
    
    def __init__(self, network, up_set, background_nodes):
        """
            Input:
                network: net[source] = [(i, t)]
                up_set: up_set[node] = score
                background_nodes: a set of nodes to swap with the input set
        """
    

        # the set of initial nodes to permute within blocks: save it to the 
        # instantiated object here
        self.nodes = up_set.keys()
        self.bg = background_nodes

        self.scores = {}
        for node in self.nodes:
            # save the scores as a set of tuples 
            self.scores[(node, str(up_set[node]))] = 1

        # build an index of all nodes in the network
        valid_network_nodes = set()
        for s in network:
            valid_network_nodes.add(s)
            for (i, t) in network[s]:
                valid_network_nodes.add(t)

        # set the background, validate 
        self.bg = background_nodes.intersection(valid_network_nodes)
        if len(self.bg) == 0:
            raise Exception("Error: none of the supplied background nodes are in the network!")
        elif len(self.bg) <= len(self.nodes):
            raise Exception("Error: background validation set is no larger than input set!")

        self.input_set_size = len(self.nodes)

    def permute(self, iterations):
        """
        Generate an array of random permutations of node scores.
        
        Input:
            iteration: the number of permutations to generate

        Returns:
            an array of hashes: each hash indexes the nodes to permuted scores
        """    
        permuted = []    
        for i in range(0, iterations):
            permuted.append(self.permuteOne())

        return permuted

    def permuteOne(self):
        """
        for each. 
        """    

        # container for scores
        swapped_scores = {}
        rand_sample = sample(self.bg, self.input_set_size)
        i = 0
        for (node, score) in self.scores:
            swapped_scores[rand_sample[i]] = float(score)
            i += 1

        return swapped_scores
