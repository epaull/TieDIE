#!/usr/bin/env	python2.7

# Requires NETWORKX 1.7+
import networkx as nx
import re

class Pathway:

	def __init__(self, network=None, validator=None):
		'''
		Input:
			network - Object in {'source':set( (interaction, target), (interaction, target)...}
			Data format.

		'''
		if network:
			self.setGraph(network)
		if validator:
			self.validator = validator
	

	def setGraph(self, network):
		'''
			Set from data
		'''
		# multi graph required for parallel edges (representing undirected links)
		self.G = nx.MultiDiGraph()
		for source in network:
			for (i, target) in network[source]:
				action, type = Pathway.classifyInteraction(i)		
				self.G.add_edge(source, target, type=type, action=action, i=i)
				# these are non-directed interactions: duplicate the graph
				if type == "INTERACTS":
					self.G.add_edge(target, source, type=type, action=action, i=i)
					

	def getActionPath(self, path):

		action = 1
		for i in range(0, len(path)-1):
			action *= self.G[path[i]][path[i+1]][0]['action']
	
		return action	

	def pathToEdges(self, path):

		edges = set()
		for i in range(0, len(path)-1):
			interaction = self.G[path[i]][path[i+1]][0]['i']
			edges.add( (path[i], interaction, path[i+1]) )

		return edges

	def addPathsFrom(self, partial_path, path_ptr, source, targets, max_depth):
		"""
			partial_path: list of genes required to get to this source node

		"""
		if not max_depth > 0:
			return

		# may not actually be in the graph: check
		if source not in self.G.nodes():
			return
	
		for nbr in self.G.neighbors(source):

			# skip paths that have already been covered
			if nbr in path_ptr['covered_nodes']:
				continue

			if nbr in targets:
				# append this target to the partial path
				# and add the full completed path to the set  
				prev_path = list(partial_path)
				prev_path.append(nbr)
				path_ptr['full_paths'].add(tuple(prev_path))
				for n in prev_path:
					path_ptr['covered_nodes'].add(n)
			else:
				# decrement the depth
				prev_path = list(partial_path)
				prev_path.append(nbr)
				self.addPathsFrom(prev_path, path_ptr, nbr, targets, max_depth-1)


	def allPaths(self, sources, targets, max_depth, valid_action=None):
		"""
		Return all paths connecting any of <sources> to <targets>
		up to the maximum specified depth. Recursive function. 
		"""

		# a paths is stored here as an ordered list of nodes
		paths = {}
		# complete source->target paths (lists)
		paths['full_paths'] = set()
		# keep an index of nodes contained in full paths
		paths['covered_nodes'] = set()
		
		for source in sources:
			self.addPathsFrom([source], paths, source, targets, max_depth)
	
		edges = set()	
		for path in paths['full_paths']:
			# validate paths
			if not self.validator.validate(path, self):
				continue
			edges = edges.union(self.pathToEdges(path))

		return edges

	def printEdgeList(self, edges):

		for edge in edges:
			print "\t".join(edge)

	def writeEdgeList(self, edges, fh):

		for edge in edges:
			fh.write("\t".join(edge)+"\n")

	@staticmethod
	def classifyInteraction(i):
		componentRE = re.compile("^-?component>$")
		activatingRE = re.compile("^-?(\S)>$")
		inactivatingRE = re.compile("^-?(\S)\|$")
		rewiredAC = re.compile("^-?REWIRED>$")
		rewiredIN = re.compile("^-?REWIRED\|$")
	
		if componentRE.match(i):
			return (0, "component")
		elif activatingRE.match(i):
			type = activatingRE.match(i)
			return (1, type.group(1))
		elif inactivatingRE.match(i):
			type = inactivatingRE.match(i)
			return (-1, type.group(1))
		elif rewiredAC.match(i):
			type = "REWIRED"
			return (1, type.group(1))
		elif rewiredIN.match(i):
			type = "REWIRED"
			return (-1, type.group(1))
		else:
			# default to activating links for HPRD or other protein
			# component links. These are bi-directional
			return (1, "INTERACTS")

	def parseNet(self, network):
		'''
			Parse network from sif file, and set the internal graph
		'''
		self.G = nx.MultiDiGraph()
		for line in open(network, 'r'):

			parts = line.rstrip().split("\t")
			source = parts[0]
			interaction = parts[1]
			target = parts[2]

			# zero actions: component links	
			action, type = Pathway.classifyInteraction(interaction)		
			# skip component links, and anything else that doesn't have an associated action value
			if action == 0:
				continue
			self.G.add_edge(source, target, type=type, action=action, i=interaction)
			if type == "INTERACTS":
				# these are bi-directional links. Simply add a link the other way
				self.G.add_edge(target, source, type=type, action=action, i=interaction)
					

class TrivialPathValidator:
	'''
		Takes a path, and a graph object, and validates it based on a given set of rules
		implements validate(), 
	'''
	def __init__(self):
		return None
	
	def validate(self, path, pathwayObj):
		return True

class BasicPathValidator:
	'''
		Takes a path, and a graph object, and validates it based on a given set of rules
		implements validate(), 
	'''
	def __init__(self, input_sets):
		'''
			Valid paths must go from source to target while passing through a 'signaling'
			intermediate set.
			self.sets['source'] = set()
			self.sets['signaling'] = set()
			self.sets['target'] = set()
		'''	
		self.sets = input_sets

	def validate(self, path, pathwayObj):

		source_set = self.sets['source']
		target_set = self.sets['target']

		# should never occur: this is an inline test of the code
		if path[0] not in source_set:
			raise Exception("Error: path start not a source node:"+path[0])
		if path[-1] not in target_set:
			raise Exception("Error: path end not a target node:"+path[-1])

		# check that the pathway action checks out with the source/target actions
		path_action = pathwayObj.getActionPath(path)
		valid_action = 1
		if self.sets['source_actions'][path[0]] == "-":
			valid_action = -1*valid_action
		if self.sets['target_actions'][path[-1]] == "-":
			valid_action = -1*valid_action
		
		if path_action != valid_action:
			return False	

		return True

	
class SHERPAValidator:
	'''
		Takes a path, and a graph object, and validates it based on a given set of rules
		implements validate(), 
	'''
	def __init__(self, input_sets):
		'''
			Valid paths must go from source to target while passing through a 'signaling'
			intermediate set.
			self.sets['source'] = set()
			self.sets['signaling'] = set()
			self.sets['target'] = set()
		'''	
		self.sets = input_sets

	def isValid(self, path, pathwayObj):
		'''
			Validate ordered set of nodes on a given path:
			must contain each of the required categories of genes
			
			Returns:
				True or False 
		'''

		source_set = self.sets['source']
		target_set = self.sets['target']
		signaling_set = self.sets['signaling']

		all = source_set.union(target_set).union(signaling_set)

		# should never occur: this is an inline test of the code
		if path[0] not in source_set:
			raise Exception("Error: path start not a source node:"+path[0])
		if path[-1] not in target_set:
			raise Exception("Error: path end not a target node:"+path[-1])

		# first gene must be a 'linker' gene
		if path[1] in all:
			return False
	
		# paths must be longer than 3 steps to include at least one of the signaling set
		if len(path) < 4:
			return False

		# all subsequent genes must be in the signaling set
		for i in range(2, len(path)-1):
			if path[i] not in signaling_set:
				return False

		#print "linker\t"+path[1]	
		# validated
		return True

	def validate(self, path, pathwayObj):

		path_action = pathwayObj.getActionPath(path)
		valid_action = 1
		if self.sets['source_actions'][path[0]] == "-":
			valid_action = -1*valid_action
		if self.sets['target_actions'][path[-1]] == "-":
			valid_action = -1*valid_action
		
		if path_action != valid_action:
			return False	

		return self.isValid(path, pathwayObj)


class TriplesValidator:
	'''
		Takes a path, and a graph object, and validates it based on a given set of rules
		implements validate(), 
	'''
	def __init__(self, input_sets):
		'''
			Valid paths must go from source to target while passing through a 'signaling'
			intermediate set.
			self.sets['source'] = set()
			self.sets['signaling'] = set()
			self.sets['target'] = set()
		'''	
		self.sets = input_sets

		self.consider_signs = True
		if input_sets['ignore_signs'] is True:
			self.consider_signs = False	
		else:
			self.consider_signs = True

	def isValid(self, path, pathwayObj):
		'''
			Validate ordered set of nodes on a given path:
			must contain each of the required categories of genes
			
			Returns:
				True or False 
		'''

		source_set = self.sets['source']
		target_set = self.sets['target']
		signaling_set = self.sets['signaling']

		source_action = float(self.sets['source_actions'][path[0]]+"1")
		all = source_set.union(target_set).union(signaling_set)
		source_target_sets = source_set.union(target_set)

		# should never occur: this is an inline test of the code
		if path[0] not in source_set:
			raise Exception("Error: path start not a source node:"+path[0])
		if path[-1] not in target_set:
			raise Exception("Error: path end not a target node:"+path[-1])

		# paths must be at least 3 steps (and include a 'signaling' gene)
		if len(path) < 3:
			return False

		# check signaling genes
		for i in range(1, len(path)-1):
			if not path[i] in signaling_set:
				return False
			elif self.consider_signs is True:
				# validate sign of signaling gene

				this_gene_action = self.sets['signaling_actions'][path[i]]
				subpath = path[0:i+1]
				# validate path action and intermediate signal are consistent
				if source_action*pathwayObj.getActionPath(subpath) != this_gene_action:
					return False

		# validated
		return True

	def validate(self, path, pathwayObj):

		path_action = pathwayObj.getActionPath(path)
		valid_action = 1
		if self.sets['source_actions'][path[0]] == "-":
			valid_action = -1*valid_action
		if self.sets['target_actions'][path[-1]] == "-":
			valid_action = -1*valid_action
		
		if path_action != valid_action and self.consider_signs:
			return False	

		return self.isValid(path, pathwayObj)

