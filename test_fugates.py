''' Scripts included in this repository '''
import mendelbp, import_gedcom, make_pedigree

''' For tracking runtime '''
import time

''' For visualizing the network -- not strictly necessary '''
import pygraphviz as pgv
import warnings

'''
Import the pedigree from a GEDCOM file, specifying the main person of interest.
Ben Stacy is a modern-day descendant. We receive back trio_dict (keys: personal
identifiers, values: list of parents' personal identifiers), gedcom_dict (keys:
personal identifiers, values: list of data such as name of the person), poi_id
(the personal identifier of the person of interest used in the GED file). The
dictionaries can be specified manually if no GED file is available for the
pedigree.
'''
start_time = time.time()
x = 10 # Include all persons seperated from the person of interest by <= x meioses
trio_dict, gedcom_dict, poi_id = import_gedcom.main('fugates.ged', x,
	'Ben /Stacy/')
print 'Threshold: %d meioses; %d people included' % (x, len(gedcom_dict.keys()))
print 'Time required for GEDCOM import: %0.1f seconds' % (time.time() - start_time)

'''
Here we specify any genotypes we know in evidence_dict:
key: GEDCOM id (can be found in the .ged file)
value: list of possible genotypes (subset of ['aa', 'Aa', 'AA'])
Note that this genotype specification is incomplete/inaccurate and intended for 
illustrative purposes
f is the frequency of the recessive allele 'a' in the founder population
'''
evidence_dict = {poi_id: ['aa']}
f = 0.1

'''
Now we call the believe propagation script to perform the marginalization,
and print the results to the screen.
'''
start_time = time.time()
my_network = mendelbp.MendelPedigree(f, trio_dict, evidence_dict)
rounds = my_network.bp(500)
print 'Convergence took %d rounds' % rounds
results = my_network.marginalize_all()

print 'id\tp(aa)\tp(Aa)\tp(AA)'
for person in my_network.ids:
	print '%s\t%0.4f\t%0.4f\t%0.4f' % (gedcom_dict[person][0], results[person][0],
		results[person][1], results[person][2])

print 'Time required for belief propagation network construction, conve' + \
	'rgence, and marginalization: %0.1f seconds' % (time.time() - start_time)

''' Visualize the results -- this part requires pygraphviz '''
start_time = time.time()
with warnings.catch_warnings():
	warnings.simplefilter('ignore', RuntimeWarning)
	pedigree = make_pedigree.make_pedigree_with_traits(trio_dict, gedcom_dict,
	 	evidence_dict, results, 'blue', True)
	pedigree.layout(prog='dot', args="-Goverlap=false -Gcompound=true")
	pedigree.draw('fugates_mendelian_bp.png')
print 'Time required to visualize results: %0.1f seconds' % (time.time() - start_time)

