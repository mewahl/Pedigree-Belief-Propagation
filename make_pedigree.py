# make_pedigree.py
#
# MEW 11/06/2014
#
# Python 2.7 used
# On laptop: export PYTHONPATH=~/Library/pygraphviz/

try:
	import sys, os, operator 
	reload(sys)
	sys.setdefaultencoding('utf8')
	reload(sys) 

	# For pedigree visualization
	import pygraphviz as pgv
	import networkx as nx

	# To remove slashes from GEDCOM file names
	import re

except ImportError as e:
	print >> sys.stderr, 'Error importing modules:\n%s' % e
	sys.exit(1)

def make_pedigree_with_traits(trio_dict, gedcom_dict, evidence_dict, scores, color, spline_type):
	# Use this function to draw a pedigree where the name, sex, and GEDCOM id
	# of each person is indicated as well as the true phenotype (if known) and
	# the ``score" in (0,1) -- which might be a scaled genetic effect or a
	# probability that a given allele drawn from that person is dominant, etc.
	# -- is used to color the person's square or circle.

	pedigree = nx.DiGraph()

	# Make a subcluster for each individual
	# At the same time, populate a dictionary of all families.
	marriage_dict = {}
	for gedcom_id, parents in trio_dict.items():
		# Make a node for the child in the trio
		entry = gedcom_dict[gedcom_id]
		draw_person(entry[0], gedcom_id, entry[1], evidence_dict.get(gedcom_id, None),
			scores[gedcom_id], color, pedigree)

		# If parents exist, give the child a hanger (which will connect to the parents
		# eventually) and add the child to the list of this couple's children.
		if parents == (None, None):
			continue
		draw_hanger(gedcom_id, pedigree)
		marriage_dict = update_marriage_dict(gedcom_id, parents, marriage_dict)

	# Currently "pedigree" is a directed graph, which is appropriate for keeping track
	# of which generation is older. However, we will eventually want to create some
	# subgraphs/clusters of equal rank, which will require converting to an agraph. This
	# variable will be a list of lists, each list corresponding to a subgraph we need to
	# make later.
	same_rank_subgraphs_to_make = []

	for marriage_node, children in marriage_dict.items():
		# Unite the couple, drawing a marriage node between them.
		pedigree, temp = draw_marriage(marriage_node, pedigree)
		same_rank_subgraphs_to_make.append(temp)

		# If the number of children is even, add a dongle from the marriage node so that the
		# children will be centered on either side. Otherwise, connect the hanger node of the
		# middle child directly to the marriage node.
		pedigree, temp = hang_children(marriage_node, children, pedigree)
		same_rank_subgraphs_to_make.append(temp)

	# Convert to agraph to build the equal rank subgraphs.
	pedigree = nx.to_agraph(pedigree)
	for item in same_rank_subgraphs_to_make:
		pedigree.add_subgraph(item,rank='same')

	# Format connecting lines (tried them all, nothing really looks good unfortunately)
	pedigree.graph_attr['splines']=spline_type

	return pedigree

def draw_person(name, gedcom_id, sex, evidence, score, color, pedigree):
	''' Convert the score to a printable string '''
	if isinstance(score, float):
		score_string = '%0.4f' % score
	elif isinstance(score, int):
		score_string = '%d' % score
	elif isinstance(score, list) and len(score) == 3 and isinstance(score[0],
		float):
		score_string = 'p(aa)=%0.4f\\np(Aa)=%0.4f\\np(AA)=%0.4f' % tuple(score)
	else:
		raise Exception('Did not recognize score format')

	'''
	Remove the slashes some programs place around last names in GEDCOM files.
	Add GEDCOM id to the node label.
	'''
	node_label = re.sub('/', '', name) + '\\n(%s)\\n%s' % (gedcom_id,
		score_string)

	'''
	If any evidence was given, add it to the node label and specify a thicker
	border around the node.
	'''
	node_style = 'filled'
	if isinstance(evidence, (int, long)):
		evidence_string = '%d' % evidence
	elif isinstance(evidence, float):
		evidence_string = '%.3f' % evidence
	elif isinstance(evidence, str):
		evidence_string = '%.3f' % evidence
	elif isinstance(evidence, list):
		if isinstance(evidence[0], str):
			evidence_string = ' or '.join(evidence)
		elif isinstance(evidence[0], (float)):
			evidence_string = ','.join(['%0.4f' % i for i in evidence])

	if evidence is not None:
		node_style += ',bold'
		node_label += '\\nEvidence: %s' % evidence_string

	# Handle gender
	if sex == 'M':
		node_shape = 'box'
		hex_color = '#FFFFAA'
	else:
		node_shape = 'circle'
		hex_color = '#AAFFFF'

	if color is not None:
		hex_color = score_to_hex_code(score, color)

	pedigree.add_node(gedcom_id, shape=node_shape, style=node_style,
		fillcolor=hex_color, height='2', width='2',
		fixedsize=True, label=node_label)
	return

def score_to_hex_code(score, color):
	'''
	If the score is given as a list, we need to convert it to a single number.
	The probability that an allele selected from the person at random is a
	recessive allele will do.
	'''
	if isinstance(score, list) and len(score) == 3 and isinstance(score[0],
		float):
		score = score[0] + 0.5 * score[1]

	''' Currently handles only a few color options '''
	score_in_hex = format(int(score*255), '02x')
	inverse_score_in_hex = format(int((1-score)*255), '02x')
	hex_color = '#'
	if color == 'red':
		hex_color += 'FF' + inverse_score_in_hex + inverse_score_in_hex
	elif color == 'blue':
		hex_color += inverse_score_in_hex + inverse_score_in_hex + 'FF'
	elif color == 'green':
		hex_color += inverse_score_in_hex + 'FF' + inverse_score_in_hex
	else:
		if color != 'gray':
			temp = 'Did not recognize color "%s" during pedigree construction: substituting gray instead.' % color
			hex_color += score_in_hex + score_in_hex + score_in_hex
	return hex_color

def draw_hanger(gedcom_id, pedigree):
	pedigree.add_node('hanger' + gedcom_id, shape='diamond', style='filled',
		fillcolor='#000000', label='', height='0.1', width='0.1')
	pedigree.add_edge('hanger' + gedcom_id, gedcom_id, dir='none',
		tailport='s', headport='n', weight=10)
	pedigree.subgraph(nbunch=['hanger' + gedcom_id, gedcom_id])
	return

def update_marriage_dict(gedcom_id, parents, marriage_dict):
	mother, father = parents
	marriage_node = '%s_to_%s' % (mother, father)
	if marriage_node in marriage_dict:
		''' Add this child to the list of previously known children '''
		children = marriage_dict[marriage_node]
		children.append(gedcom_id)
		marriage_dict[marriage_node] = children
	else:
		marriage_dict[marriage_node] = [gedcom_id]
	return marriage_dict

def draw_marriage(marriage_node, pedigree):
	mother, father = marriage_node.split('_to_')

	pedigree.add_node(marriage_node, shape='diamond', style='filled',
		fillcolor='#AAAAAA', label='', height='0.1', width='0.1')
	pedigree.add_edge(father, marriage_node, dir='none', arrowhead='none',
		weight=10)
	pedigree.add_edge(marriage_node, mother, dir='none', arrowhead='none',
		weight=10)

	return pedigree, [father, marriage_node, mother]

def hang_children(marriage_node, children, pedigree):
	nodes_to_group = ['hanger' + child for child in children]

	'''
	If the number of children is even, add a dongle from the marriage node so
	that the children will be centered on either side. Otherwise, connect the
	hanger node of the middle child directly to the marriage node.
	'''

	if len(children) % 2 == 0:
		# Draw a "dongle" that will connect horizontally to the "hangers" of the children
		pedigree.add_node('dongle' + marriage_node, shape='diamond', style='filled',
			fillcolor='#AAAAAA', label='', height='0.1', width='0.1')
		pedigree.add_edge(marriage_node, 'dongle' + marriage_node, arrowhead='none',
			tailport='s', headport='n')
		pedigree.subgraph([marriage_node, 'dongle' + marriage_node])
		nodes_to_group.insert(len(children)/2, 'dongle' + marriage_node)
	else:
		pedigree.add_edge(marriage_node, nodes_to_group[len(children)/2], arrowhead='none',
			tailport='s', headport='n')

	# Connect the hanger nodes of the children (and the dongle, if one exists)
	for i in xrange(0,len(nodes_to_group)-1):
		pedigree.add_edge(nodes_to_group[i], nodes_to_group[i+1], dir='none', weight=2)

	return pedigree, nodes_to_group

def make_pedigree_simple(trio_dict, gedcom_dict):
	# Use this function to draw a pedigree that indicates name, sex, and
	# GEDCOM id but does not provide any information on phenotype.

	pedigree = nx.DiGraph()

	# Make a subcluster for each individual
	# At the same time, populate a dictionary of all families.
	marriage_dict = {}
	for child, parents in trio_dict.items():
		# Make a subcluster for the child
		entry = gedcom_dict[child]
		name = re.sub('/', '', entry[0]) + '\n(%s)' % child
		if entry[1] == 'F':
			pedigree.add_node(child, shape='circle', style='filled', fillcolor='#FFAAAA', label=name, height='2', width='2', fixedsize=True)
		elif entry[1] == 'M':
			pedigree.add_node(child, shape='box', style='filled', fillcolor='#AAAAFF', label=name, height='2', width='2', fixedsize=True)
		else:
			pedigree.add_node(child, shape='diamond', style='filled', fillcolor='#AAAAAA', label=name, height='2', width='2', fixedsize=True)
		# All individuals should have either two parents or no parents after GEDCOM import
		# since dummy parents are added at that stage. Check only the mother.
		mother, father = parents
		if mother is None:
			continue

		# If a mother and father exist, we will need a "hanger" above the child's node to
		# attach to the "dongle" coming down from the marriage node.
		pedigree.add_node('hanger' +  child, shape='diamond', style='filled', fillcolor='#000000', label='', height='0.1', width='0.1')
		pedigree.add_edge('hanger' + child, child, dir='none', tailport='s', headport='n')
		pedigree.subgraph(nbunch=['hanger' + child, child])

		# Now let's make sure we keep track of which marriage this child descends from
		node = '%s_to_%s' % (mother, father)
		if node in marriage_dict:
			# Add this child to the list of previously known children
			children = marriage_dict[node]
			children.append(child)
			marriage_dict[node] = children
		else:
			marriage_dict[node] = [child]

	# When we are done drawing the nodes and edges, we'll need to draw some subgraphs where
	# all elements have equal rank. This will require switching from a DiGraph to an agraph.
	# Retain a list of lists of nodes to join in equal-rank subgraphs for later.
	same_rank_subgraphs_to_make = []
	for node, children in marriage_dict.items():
		mother, father = node.split('_to_')

		# Draw a subcluster for the union
		pedigree.add_node(node, shape='diamond', style='filled', fillcolor='#AAAAAA', label='', height='0.1', width='0.1')
		pedigree.add_edge(father, node, dir='none', weight=10)
		pedigree.add_edge(node, mother, dir='none', weight=10)
		same_rank_subgraphs_to_make.append([father,node,mother])

		# If the number of children is even, add a dongle from the marriage node so that the
		# children will be centered on either side. Otherwise, connect the hanger node of the
		# middle child directly to the marriage node.
		nodes_to_group = ['hanger' + child for child in children]
		if len(children) % 2 == 0:
			# Draw a "dongle" that will connect horizontally to the "hangers" of the children
			pedigree.add_node('dongle' + node, shape='diamond', style='filled', fillcolor='#AAAAAA', label='', height='0.1', width='0.1')
			pedigree.add_edge(node, 'dongle' + node, arrowhead='none', tailport='s', headport='n')
			pedigree.subgraph([node, 'dongle' + node])
			nodes_to_group.insert(len(children)/2, 'dongle' + node)
		else:
			pedigree.add_edge(node, nodes_to_group[len(children)/2], arrowhead='none', tailport='s', headport='n')
		same_rank_subgraphs_to_make.append(nodes_to_group)

		# Connect the hanger nodes of the children (and the dongle, if one exists)
		for i in xrange(0,len(nodes_to_group)-1):
			pedigree.add_edge(nodes_to_group[i], nodes_to_group[i+1], dir='none')

	# To build the equal rank subgraphs, we need to convert to an agraph
	pedigree = nx.to_agraph(pedigree)
	for item in same_rank_subgraphs_to_make:
		pedigree.add_subgraph(item,rank='same')

	#pedigree.graph_attr['splines']='polyline'
	pedigree.graph_attr['splines']=False
	return pedigree
