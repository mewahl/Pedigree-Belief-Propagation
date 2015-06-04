'''
Defines the following classes:
* MendelVariableNode
* PriorMendelFactorNode
* TrioMendelFactorNode
* MendelPedigree

Typical usage:
* Create a belief propagation network by creating a new MendelPedigree by
  providing:
  - The frequency of the 'a' allele to be assumed for founders when no other
    evidence is provided (see below)
  - A dictionary of child-parent trios (key: child's identifer; values: list of
  	both parents' identifiers)
  - A dictionary containing any evidence available on individuals' genotypes
    (key: person's identifer; values: either a list of possible genotypes, i.e.
    a subset of ['aa', 'Aa', 'AA'], or a list containing exact probabilities of
    aa, Aa, and AA genotypes)
* Run belief propagation until convergence (or a capped round number)
* Calculate marginal probabilities for each genotype for each individual
'''

import re

''' A message is considered different if its L2 distance is > epsilon '''
epsilon = 0.0000001

'''
Probability of a child's genotype (rows correspond to values in genotype_dict)
given the parents' genotypes (columns are parent1_key*3 + parent2_key)
 \ M:      aa      |      Aa      |      AA      |
  \F: aa | Aa | AA | aa | Aa | AA | aa | Aa | AA |
C: \ ____|____|____|____|____|____|____|____|____|
aa  | 1  | .5 | 0  | .5 | .25|  0 | 0  | 0  | 0  |
Aa  | 0  | .5 | 1  | .5 | .5 | .5 | 1  | .5 | 0  |
AA  | 0  | 0  | 0  | 0  | .25| .5 | 0  | .5 | 1  |
'''
child_given_parents =   [[1.0, 0.5, 0.0, 0.5, 0.25, 0.0, 0.0, 0.0, 0.0],
                       	 [0.0, 0.5, 1.0, 0.5,  0.5, 0.5, 1.0, 0.5, 0.0],
					   	 [0.0, 0.0, 0.0, 0.0, 0.25, 0.5, 0.0, 0.5, 1.0]]

genotype_dict = {'aa': 0, 'Aa': 1, 'AA': 2}

'''
Probability of a parent's genotype (rows correspond to values in genotype_dict)
given the child's and other parent's genotypes (columns are child_key*3+parent_key)
  \ C:       aa        |       Aa        |       AA        |
   \P1: aa | Aa  | AA  | aa  | Aa  | AA  | aa  | Aa  | AA  |
P2: \ _____|_____|_____|_____|_____|_____|_____|_____|_____|
 aa  | 2/3 | 2/3 |  0  |  0  | 1/3 | 2/3 |  0  |  0  |  0  |
 Aa  | 1/3 | 1/3 |  0  | 1/3 | 1/3 | 1/3 |  0  | 1/3 | 1/3 |
 AA  |  0  |  0  |  0  | 2/3 | 1/3 |  0  |  0  | 2/3 | 2/3 |
'''
parent_given_others = [[0.0]*9]*3
parent_given_others =	[[2/3.0, 2/3.0, 0.0,   0.0, 1/3.0, 2/3.0, 0.0,   0.0,   0.0],
						 [1/3.0, 1/3.0, 0.0, 1/3.0, 1/3.0, 1/3.0, 0.0, 1/3.0, 1/3.0],
						 [  0.0,   0.0, 0.0, 2/3.0, 1/3.0,   0.0, 0.0, 2/3.0, 2/3.0]]


class MendelPedigree:
	'Belief propagation network for Mendelian genotypes in pedigrees'

	def __init__(self, f, trio_dict, *args):
		'''
		f is the frequency of the 'a' allele in founders without other evidence.
		trio_dict: key is child's id, value is list of parents' ids (None if the
			child is a founder)
		optional argument is a dictionary of genotypic evidence (key: person's
			id, value: several options, see handling in PriorFactorNode class)
		'''
		if len(args) == 0:
			evidence_dict = {}
		elif len(args) == 1:
			evidence_dict = args[0]
		else:
			raise Exception('Too many optional arguments (0 or 1 expected):' + \
				'\n%s' % ' '.join([str(i) for i in args]))

		''' Make variable nodes for every person '''
		self.ids = trio_dict.keys()
		self.ids.sort(natcasecmp)
		n = len(self.ids)
		self.variable_nodes = []
		for i in xrange(0,n):
			self.variable_nodes.append(MendelVariableNode())

		'''
		We make two kinds of factor nodes:
		* PriorMendelFactorNode: specifies any evidence we have on the prior
		  probabilities of each genotype. Needed for every founder and for
		  every other person for whom evidence is available.
		* TrioMendelFactorNode: needed for every child-parent trio.
		'''
		self.factor_nodes = []
		for i, person in enumerate(self.ids):
			parent1, parent2 = trio_dict[person]
			if parent1 == None: # parent2 is also unknown
				''' This is a founder -- make only a PriorMendelFactorNode '''
				if person in evidence_dict:
					self.factor_nodes.append(PriorMendelFactorNode(self.variable_nodes[i],
						f=f, evidence=evidence_dict[person]))
				else:
					self.factor_nodes.append(PriorMendelFactorNode(self.variable_nodes[i],
						f=f))
			else:
				'''
				This is a child: need a TrioMendelFactorNode and possiby a
				PriorMendelFactorNode as well.
				'''
				if person in evidence_dict:
					self.factor_nodes.append(PriorMendelFactorNode(self.variable_nodes[i],
						evidence=evidence_dict[person]))
				parent1 = self.ids.index(parent1)
				parent2 = self.ids.index(parent2)
				self.factor_nodes.append(TrioMendelFactorNode(self.variable_nodes[i],
					self.variable_nodes[parent1], self.variable_nodes[parent2]))

		self.bp_run = False
		return

	def bp(self, max_round_number):
		rounds = 0
		while (True and rounds < max_round_number):
			messages_sent = 0
			for variable_node in self.variable_nodes:
				messages_sent += variable_node.send_pending_messages()
			for factor_node in self.factor_nodes:
				messages_sent += factor_node.send_pending_messages()
			if messages_sent == 0:
				break
			rounds += 1
		self.bp_run = True
		return rounds

	def marginalize_all(self):
		marginals = {}
		for person in self.ids:
			marginals[person] = self.marginalize_one(person)
		return marginals

	def marginalize_one(self, id_to_marginalize):
		if self.bp_run == False:
			raise Exception('Marginalization can only be performed after ' + \
				'belief propagation has been performed successfully.')
		try:
			i = self.ids.index(id_to_marginalize)
		except Exception as e:
			raise Exception('Error finding the index of the identifier for' + \
				' marginalization:\n%s' % e)
		messages = self.variable_nodes[i].own_messages
		p = [1.0, 1.0, 1.0]
		for message in messages:
			if message is None:
				raise Exception('Marginalization not possible for this ' + \
					'node: some messages not received.')
				break
			for i in range(0,3):
				p[i] *= message[i]
		if sum(p) == 0:
			raise Exception('Not possible to normalize marginalized ' + \
				'probabilities for this node since all probabilities ' + \
				'are zero.')
			return [0.0,0.0,0.0]
		return [i/sum(p) for i in p]

class MendelVariableNode:
	'Base class for all variable nodes'

	def __init__(self):
		self.factor_nodes = []
		self.own_messages = []
		self.pending_outgoing_messages = []
		return

	def add_factor_node(self, factor_node, *args):
		self.factor_nodes.append(factor_node)
		if len(args) == 1:
			self.own_messages.append(args[0])
		else:
			self.own_messages.append(None)

		'''
		When a new factor node is added, we need to reset which messages are
		pending. (We might have had enough information to send a message
		previously but not enough information now.) Since we are still building
		the network, the variable nodes have not sent any messages yet, so
		there will be a message pending if messages have been received from all
		other factor nodes. Note that if there are no other factor nodes, then
		this is a terminal node and there is a pending outgoing message.
		'''
		n = len(self.factor_nodes)
		self.pending_outgoing_messages = [True]*n
		for i in xrange(0, n):
			for j in xrange(0, n):
				if i == j:
					continue
				if self.own_messages[j] is None:
					self.pending_outgoing_messages[i] = None
		return

	def receive_message(self, factor_node, incoming_message):
		''' Find the node that sent this message; store the message appropriately '''
		sender = None
		for i, item in enumerate(self.factor_nodes):
			if item is factor_node:
				if did_message_change(self.own_messages[i], incoming_message) is False:
					self.own_messages[i] = incoming_message
					return
				self.own_messages[i] = incoming_message
				sender = i
		if sender is None:
			raise Exception('Could not find sender!')

		'''
		Now that we have received a message from the sender, we may need to
		send new messages to the other factor nodes. For each of the factor
		nodes besides the sender, check whether any message has been sent yet
		(in which case the value of pending_outgoing_messages[i] will be False
		instead of None); if so, we need to send a new message. Otherwise,
		check if we now have enough information to send the first outgoing
		message to this factor node (i.e., whether every other factor node has
		sent a message already) and if so, indicate there is a message pending.
		'''
		n = len(self.factor_nodes)
		for i in range(0,n):
			if i == sender:
				continue
			if self.pending_outgoing_messages[i] == False:
				self.pending_outgoing_messages[i] = True
			else:
				for j in range(0,n):
					if j == i:
						continue
					if self.own_messages[j] is None:
						break
					if j == n - 1:
						self.pending_outgoing_messages[i] = True

		return

	def send_pending_messages(self):
		n = len(self.factor_nodes)
		any_messages_sent = False
		for i in xrange(0,n):
			if self.pending_outgoing_messages[i] == True:
				''' The outgoing message is just the normalized product of messages
				    from all other neighboring factor nodes '''
				outgoing_message = [1.0, 1.0, 1.0]
				for j in xrange(0,n):
					if i == j:
						continue
					outgoing_message = [a * b for a, b in zip(outgoing_message,
						self.own_messages[j])]
				outgoing_message = [a / sum(outgoing_message) for a in outgoing_message]

				''' Send the message '''
				self.factor_nodes[i].receive_message(self, outgoing_message)
				any_messages_sent = True

				''' No message pending now '''
				self.pending_outgoing_messages[i] = False
		return any_messages_sent

class PriorMendelFactorNode:
	'Prior Factor Node (use to pass information on priors to a variable node)'

	def __init__(self, variable_node, **kwargs):
		self.variable_node = variable_node
		self.variable_node.add_factor_node(self)

		'''
		If this person was a founder, then the frequency of the recessive
		allele (f) will be an argument.
		'''
		if 'f' in kwargs:
			f = kwargs['f']
			p = [f**2, 2*f*(1-f), (1-f)**2]
		else:
			p = [1.0, 1.0, 1.0]

		''' Update based on evidence, if any exists '''
		if 'evidence' in kwargs:
			evidence = kwargs['evidence']
			if isinstance(evidence[0], str):
				'''
				Evidence has been supplied in the form of a list of possible
				genotypes (a subset of ['aa', 'Aa', 'AA'])
				'''
				for genotype, i in genotype_dict.items():
					if genotype not in evidence:
						p[i] = 0.0
			elif len(evidence) == 3 and isinstance(evidence[0], (int, float,
				long)):
				'''
				Evidence has been supplied in the form of the prior probability
				of each genotype. Just make sure the values are floats and that
				the probabilities are normalized.
				'''
				p = [float(i) for i in evidence]
			else:
				raise Exception('Did not recognize the format of evidence ' + \
					'presented.')

		''' Normalize the probabilities '''
		p = [i/sum(p) for i in p]

		''' Send a message to the variable node '''
		self.variable_node.receive_message(self, p)
		return

	'''
	Prior factor nodes are always terminal nodes. If I understand correctly,
	they do not need to do anything with messages they may receive and under
	no circumstances need to send a second message. However, we include these
	function fragments so that we can just iterate over all factor nodes w/o
	checking which type is which.
	'''
	def receive_message(self, variable_node, outgoing_message):
		return

	def send_pending_messages(self):
		''' Return value indicates that no messages were sent (used to check
			for convergence) '''
		return False

class TrioMendelFactorNode:
	'Trio Factor Node'

	def __init__(self, child_node, parent1_node, parent2_node):
		self.child_node = child_node
		self.parent1_node = parent1_node
		self.parent2_node = parent2_node
		
		self.child_node.add_factor_node(self)
		self.parent1_node.add_factor_node(self)
		self.parent2_node.add_factor_node(self)

		self.message_from_child = None
		self.message_from_parent1 = None
		self.message_from_parent2 = None

		''' To handle loopy pedigrees, we always send an initial message '''
		self.child_node.receive_message(self, [1/3.0]*3)
		self.parent1_node.receive_message(self, [1/3.0]*3)
		self.parent2_node.receive_message(self, [1/3.0]*3)

		self.pending_message_for_child = False
		self.pending_message_for_parent1 = False
		self.pending_message_for_parent2 = False

		return

	def receive_message(self, variable_node, incoming_message):
		'''
		We will store the incoming message and check whether we now need to
		send updated messages to the other two nodes. We do not need to
		send updates if (a) the message has not changed significantly or (b)
		we do not have enough information to send a messag yet.
		'''
		if variable_node is self.child_node:
			if did_message_change(self.message_from_child, incoming_message) is False:
				self.message_from_child = incoming_message
				return
			self.message_from_child = incoming_message

			if self.message_from_parent1 is not None:
				self.pending_message_for_parent2 = True
			if self.message_from_parent2 is not None:
				self.pending_message_for_parent1 = True

		elif variable_node is self.parent1_node:
			if did_message_change(self.message_from_parent1, incoming_message) is False:
				self.message_from_parent1 = incoming_message
				return
			self.message_from_parent1 = incoming_message

			if self.message_from_child is not None:
				self.pending_message_for_parent2 = True
			if self.message_from_parent2 is not None:
				self.pending_message_for_child = True	

		elif variable_node is self.parent2_node:
			if did_message_change(self.message_from_parent2, incoming_message) is False:
				self.message_from_parent2 = incoming_message
				return
			self.message_from_parent2 = incoming_message

			if self.message_from_child is not None:
				self.pending_message_for_parent1 = True
			if self.message_from_parent1 is not None:
				self.pending_message_for_child = True
		return

	def send_pending_messages(self):
		any_messages_sent = False

		if self.pending_message_for_child == True:
			parent1 = self.message_from_parent1
			parent2 = self.message_from_parent2
			temp = [0.0, 0.0 ,0.0]
			'''
			The outgoing message is calculated by summing over all possible
			genotypes of the parents and child. Probabilities for the child's
			genotype given the parents' genotypes are found in the global
			variable child_given_parents.
			'''
			for i in range(0,3):
				for j in range(0,3):
					for k in range(0,3):
						temp[k] += parent1[i] * parent2[j] * child_given_parents[k][3*i+j]
			outgoing_message = []
			for i in range(0,3):
				if temp[i] == 0:
					outgoing_message.append(0.0)
				else:
					outgoing_message.append(temp[i] / sum(temp))	

			self.child_node.receive_message(self, outgoing_message)
			any_messages_sent = True
			self.pending_message_for_child = False

		if self.pending_message_for_parent1 == True:
			'''
			Morally equivalent to the above, but we use the conditional
			probabilities in parent_given_others instead.
			'''
			child = self.message_from_child
			parent2 = self.message_from_parent2
			temp = [0, 0 ,0]
			for i in range(0,3):
				for j in range(0,3):
					for k in range(0,3):
						temp[k] += child[i] * parent2[j] * parent_given_others[k][3*i+j]
			outgoing_message = []
			for i in range(0,3):
				if temp[i] == 0:
					outgoing_message.append(0.0)
				else:
					outgoing_message.append(temp[i] / sum(temp))		

			self.parent1_node.receive_message(self, outgoing_message)
			any_messages_sent = True
			self.pending_message_for_parent1 = False

		if self.pending_message_for_parent2 == True:
			''' Same as above '''
			child = self.message_from_child
			parent1 = self.message_from_parent1
			temp = [0.0, 0.0 ,0.0]
			for i in range(0,3):
				for j in range(0,3):
					for k in range(0,3):
						temp[k] += child[i] * parent1[j] * parent_given_others[k][3*i+j]
			outgoing_message = []
			for i in range(0,3):
				if temp[i] == 0:
					outgoing_message.append(0.0)
				else:
					outgoing_message.append(temp[i] / sum(temp))		

			self.parent2_node.receive_message(self, outgoing_message)
			any_messages_sent = True
			self.pending_message_for_parent2 = False
		return any_messages_sent

def did_message_change(old_message, new_message):
	if old_message is None:
		return True
	total = 0
	for i in range(0,3):
		total += (old_message[i] - new_message[i])**2
	if total < epsilon:
		return False
	else:
		return True

''' These functions implement natural sorting (thanks Stack Overflow!) '''
def try_int(s):
    try: return int(s)
    except: return s

def natsort_key(s):
    return map(try_int, re.findall(r'(\d+|\D+)', s))

def natcmp(a, b):
    return cmp(natsort_key(a), natsort_key(b))

def natcasecmp(a, b):
    return natcmp(a.lower(), b.lower())

