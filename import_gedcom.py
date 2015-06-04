'''
import_gedcom.py

MEW 12/08/2014

Extracts name, sex, and parentage information for each individual in the
GEDCOM file. Removes individuals who are more than a defined number of
meioses from a terminal node. Retains only individuals who remain connected
to a selected person of interest. Returns:
* gedcom_dict: key: GEDCOM id; value: list
  [0] Name
  [1] Sex ('M' or 'F')
  [2] Placeholder status (True if this is a dummy person inserted where no
	    information on a parent was available)
  [3] Person of Interest status (True if name contains substring used to define
	    person of interest)
* trio_dict: key: GEDCOM id; value: two-tuple
  [0] GEDCOM id of mother
  [1] GEDCOM id of father

Python 2.7 used
PYTHONPATH=/san2/wahl/lib

regex inspired by Greg Hewgill's GEDCOM to XML parser:
http://stackoverflow.com/questions/1919593/is-there-a-gedcom-parser-written-in-python
'''
try:
	import sys, os, argparse, operator
	# For calling topological sorting program tsort
	import subprocess
	# For using a temporary directory and cleaning up when finished
	import shutil, tempfile
	# For parsing GEDCOM file entries
	import re
	# For adjacency matrix used to cull too-distant relatives
	import scipy.sparse as sps
	from scipy.sparse.csgraph import breadth_first_order, bellman_ford
except ImportError as e:
	print >> sys.stderr, 'Error importing modules:\n%s' % e
	sys.exit(1)
	
def load_gedcom(input_filename, poi_name):
	try:
		input_file = open(input_filename, 'r')
		
		mother = None
		father = None
		''' GEDCOM ids with names matching person of interest '''
		poi_ids = []
		''' Make placeholders when 1+ parents missing '''
		placeholders_used = 0
		adopted = []
		current_person = None
		
		'''
		Stores information used for pedigree visualization
		Key: GEDCOM identifier (string)
		Value: List
		[0] Name (string)
		[1] Sex (string)
		[2] Placeholder status (Boolean)
		[3] Person of interest (POI) status (Boolean)
		'''
		gedcom_dict = {}
		
		'''
		Stores information on parentage
		Key: GEDCOM identifier of child (string)
		Values: Two-tuple
		[0] Mother
		[1] Father
		'''
		trio_dict = {}
		
		for line in input_file:
			'''
			Each line should start with one or more integers ('d+'),
			representing the line "level". The second part of the line may be
			a person or family ID (e.g., '@P1@' or '@F1@'). The next word is
			a "tag" such as INDI or FAM (if the id field existed) or a
			descriptor like HUSB, WIFE, TRLR (for end of file) etc. Anything
			remaining on the line is stored  in the variable "data".
			'''
			m = re.match(r"(\d+) (@(\w+)@ )?(\w+)( (.*))?", line.strip())
			if m is None:
				continue
			id = m.group(3)
			tag = m.group(4)
			data = m.group(6)
			
			if data is None:
				if tag == 'INDI':
					'''
					This is just the first line in an entry describing this
					person. We will need to fill in the information as we
					read subsequent lines.
					'''
					gedcom_dict[id] = ["", "U", False, False]
					trio_dict[id] = (None,None)
					current_person = id
				elif tag == 'FAM':
					''' First line in an entry describing a family. '''
					mother = None
					father = None
			else:
				'''
				Check whether the data contains a person's ID: if so, this
				could be data on relationships within a family record.
				'''
				m = re.match(r"@(\w+)@", data)
				if m:
					id = m.group(1)
					'''
					If this line is describing who the parents are in the
					family, "id" contains the GEDCOM id of a parent.
					'''
					if tag == 'HUSB':
						father = id
						continue
					elif tag == 'WIFE':
						mother = id
						continue
					elif tag == 'CHIL':
						'''
						Parents always come before children in GEDCOM family
						entries, so if we reach this point, we must have a child's
						record. Skip adopted children. Create dummy parents if
						necessary. Update the child's parentage in trio_dict and
						print to the directed acyclic grant file for topo sorting.
						'''
						if id in adopted:
							continue
							
						if mother is None:
							placeholders_used += 1
							mother = 'Placeholder%d' % placeholders_used
							gedcom_dict[mother] = [mother, "F", True, False]
							trio_dict[mother] = (None ,None)
						if father is None:
							placeholders_used += 1
							father = 'Placeholder%d' % placeholders_used
							gedcom_dict[father] = [father, "M", True, False]
							trio_dict[father] = (None, None)
						
						trio_dict[id] = (mother, father)
						continue
				elif tag == 'NAME' and current_person is not None:
					current_entry = gedcom_dict[current_person]
					current_entry[0] = data
					if poi_name in data:
						current_entry[3] = True
						poi_ids.append(current_person)
					gedcom_dict[current_person] = current_entry
				elif tag == 'SEX':
					current_entry = gedcom_dict[current_person]
					current_entry[1] = data
					gedcom_dict[current_person] = current_entry
				elif tag == 'PEDI' and data == 'adopted':
					adopted.append(current_person)
		
		''' Check if we found exactly one person of interest '''
		if len(poi_ids) != 1:
			raise Exception('Expected to find one person with name match' + \
			    'ing "%s"; found %d instead' % (poi_name, len(poi_ids)))
		
	except Exception as e:
		raise Exception('Error while reading in GEDCOM file:\n%s' % e)
		return
	finally:
		input_file.close()
	return gedcom_dict, trio_dict, poi_ids[0]

def find_persons_less_than_x_meioses_from_poi(trio_dict, gedcom_dict, x, poi):
	gedcom_ids = trio_dict.keys()
	poi_index = gedcom_ids.index(poi)

	''' Build a sparse adjacency matrix '''
	n = len(gedcom_ids)
	rows =[]
	cols = []
	for i, child in enumerate(gedcom_ids):
		mother, father = trio_dict[child]
		if mother is None:
			continue
		mother = gedcom_ids.index(mother)
		father = gedcom_ids.index(father)
		rows.extend([i, i, mother, father])
		cols.extend([mother, father, i, i])
	adj_matrix = sps.coo_matrix(([1] * len(rows), (rows, cols)), shape=(n,n))
	
	''' Use Bellman-Ford to calculate distance from P.O.I. to others '''
	meioses = bellman_ford(adj_matrix, unweighted=True,
		indices=[poi_index])[0].tolist()

	''' Remove individuals who are too distant to consider. '''
	gedcom_ids_to_keep = []
	for i, gedcom_id in enumerate(gedcom_ids):
		if meioses[i] < x:
			gedcom_ids_to_keep.append(gedcom_id)

	new_gedcom_dict = {}
	new_trio_dict = {}
	for gedcom_id in gedcom_ids_to_keep:
		new_gedcom_dict[gedcom_id] = gedcom_dict[gedcom_id]
		mother, father = trio_dict[gedcom_id]
		'''
		Make sure that both parents of a child are included in the pedigree
		either parent would be included in the pedigree.
		'''
		if mother in gedcom_ids_to_keep:
			if father not in gedcom_ids_to_keep:
				new_gedcom_dict[father] = gedcom_dict[father]
				new_trio_dict[father] = (None, None)
			new_trio_dict[gedcom_id] = trio_dict[gedcom_id]
		elif father in gedcom_ids_to_keep:
			new_gedcom_dict[mother] = gedcom_dict[mother]
			new_trio_dict[mother] = (None, None)
			new_trio_dict[gedcom_id] = trio_dict[gedcom_id]
		else:
			new_trio_dict[gedcom_id] = (None, None)

	return new_trio_dict, new_gedcom_dict


def main(input_filename, x, poi):
	gedcom_dict, trio_dict, poi_id = load_gedcom(input_filename, poi)
	trio_dict, gedcom_dict = find_persons_less_than_x_meioses_from_poi(trio_dict,
		gedcom_dict, x, poi_id)
	return trio_dict, gedcom_dict, poi_id
	
if __name__ == "__main__":
	#Handle input arguments
	parser = argparse.ArgumentParser(description='Convert a GEDCOM file into' + \
	    ' a Bayesian network')
	parser.add_argument('--input-file', '-i', type=str, required=True,
		help='input GEDCOM filename', dest='input_filename', action='store')
	parser.add_argument('--max_num_generations', '-x', type=int, required=True,
		help='maximum # of generations from at least one terminal node an individual' + \
		' must be to remain included in the tree', dest='max_num_generations',
		action='store')
	parser.add_argument('--person_of_interest', '-p', type=str, required=True,
		help='substring of the name of the person of interest', dest='poi',
		action='store')
	args = parser.parse_args()

	try:
		tmp_path = tempfile.mkdtemp(prefix='gedcom_')
		trio_dict, gedcom_dict = main(args.input_filename, args.max_num_generations,
			args.poi, tmp_path)
	except Exception as e:
		print >> sys.stderr, 'Ended import_gedcom.py early:\n%s' % e
	finally:
		shutil.rmtree(tmp_path)
