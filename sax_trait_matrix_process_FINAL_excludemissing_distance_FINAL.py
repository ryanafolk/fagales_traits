#!/usr/bin/env python3

import csv # To process tab-delimited files
import numpy # Used for random sampling 
import sys # To process arguments
import pandas # 
import re
import random
import scipy.spatial
import itertools

infile = sys.argv[1] # Input csv
outfile = sys.argv[2] # Pre-coded output csv
outcoded = sys.argv[3] # Coded output csv
outcoded_onevalue = sys.argv[4] # Coded output csv, randomly chosen from list
codeguide = sys.argv[5] # Guide to coding done
distancematrix = sys.argv[6] # Distance matrix

traitdataset = pandas.read_csv(infile, sep='\t', header='infer', encoding = 'utf-8') # Encoding can be problematic
# If you get BOM-related errors use this:
# perl -e 's/\xef\xbb\xbf//;' -pi~ sax_traits*

# E.g., run as:
# ./sax_trait_matrix_process.py saxtraits_final.csv out.csv coded_out.csv coded_oneval_out.csv codeguide.csv distancematrix.csv

# List of terms to remove
censorlist = ["widely","broadly","narrowly","very","obliquely","subterete","semiterete","subpeltate","terete","peltate","subcylindrical","cylindrical","pale","bright","dark","deep","spots-[A-Za-z-]*","dull","dirty","stripes-[A-Za-z-]*","compact","cylindric","narrow","terminal","simple"]

# Separate term synonym list
partialmatchsynonymlist = [["circular", "orbicular"], ["subelliptic", "elliptic"], ["suboblong", "oblong"], ["round", "orbicular"], ["subflabellate", "flabellate"], ["subspatulate", "spatulate"], ["suborbicular", "orbicular"], ["rounded-ovate", "orbicular-ovate"], ["reniform-rounded", "orbicular-reniform"], ["ovate-rounded", "ovate-orbicular"], ["subovate", "ovate"], ["subreniform", "reniform"], ["subcircular", "orbicular"], ["oblong-terete", "oblong"], ["lorate", "linear"], ["ensiform", "linear"], ["flabelliorm", "flabellate"], ["orbiculate", "orbicular"], ["rounded-reniform", "orbicular-reniform"], ["suborbiculate", "orbicular"], ["subrhomboidal", "rhomboid"], ["ovate-roundish", "ovate-orbicular"], ["roundish", "orbicular"], ["ovate-circular", "ovate-orbicular"], ["rotund", "orbicular"], ["subcordate", "cordate"], ["subensate", "ensiform"], ["sublinear", "linear"], ["semicircular", "orbicular"], ["semi-circular", "orbicular"], ["ovate-suborbicular", "ovate-orbicular"], ["subcuneate", "cuneate"], ["pentagonal", "polygonal"], ["septagonal", "polygonal"], ["pentangular", "polygonal"], ["angular-orbiculate", "polygonal"], ["orbicular-pentagonal", "polygonal"], ["rectangular", "rhomboid"], ["subobovate", "obovate"], ["navicular", "cymbiform"], ["keeled", "cymbiform"], ["liguliform", "ligulate"], ["semi-orbicular", "orbicular"], ["semiorbicular", "orbicular"], ["deltate", "deltoid"], ["subdeltoid", "deltoid"], ["sublanceolate", "lanceolate"], ["suboblanceolate", "oblanceolate"]]
fullmatchsynonymlist = [["^corymbose cymes$", "corymbs"], ["^corymbose$", "corymbs"], ["^corymbiform$", "corymbs"], ["^paniculate cymes$", "panicles"], ["^spicate$", "spikes"], ["^spiciform$", "spikes"], ["^corymbose panicles$", "corymbs"], ["^spicate thyrses$", "thyrses"], ["^subcorymbose$", "corymbs"], ["^corymbose-paniculate$", "corymbs"], ["^corymbiform panicles$", "corymbs"], ["^cymose-corymbiform$", "corymbs"], ["^paniculate-corymbiform$", "corymbs"], ["^spicate-paniculate$", "panicles"], ["^racemose thyrses$", "thyrses"], ["^thyrsoid panicle$", "thyrses"], ["^paniculate$", "panicles"], ["^cymose$", "cymes"], ["^racemiform cymes$", "cymes"], ["^spicate racemes$", "spikes"], ["^corymbiform cymes$", "corymbs"], ["^thyrsoid$", "thyrses"], ["^subpaniculate$", "panicles"], ["^thyrsoid$", "thyrses"]]
colorsynonymlist = [["pinkish", "pink"], ["lilac", "purple"], ["whitish", "white"], ["yellowish", "yellow"], ["reddish", "red"], ["greenish", "green"], ["golden-yellow", "yellow"], ["crimson", "red"], ["cream", "white"], ["sulphur-yellow", "yellow"], ["rose-coloured", "pink"], ["rose", "pink"], ["brownish", "brown"], ["lemon", "yellow"], ["creamy", "cream"], ["salmon-pink", "orange-pink"], ["violet", "purple"], ["ivory-white", "white"], ["pink-violet", "pink-purple"], ["ivory", "white"], ["straw-colored", "yellow"], ["lavendar", "purple"], ["violetish", "purple"], ["purplish", "purple"], ["glaucous-pink", "pink"], ["olive-green", "green"], ["stramineous", "yellow"], ["grey-green", "green"], ["deep rose-coloured", ""], ["maroon", "red-brown"], ["gray-green", "green"], ["lavender", "purple"], ["purpleish", "purple"], ["scarlet", "red"]]

# Full synonymy list
combinedsyn = partialmatchsynonymlist + fullmatchsynonymlist + colorsynonymlist

# Process categorical text data
def cleaner(commalist): # Dataframe elements to lists, alphabetize hyphenated terms, censor, and synonymize
	y = []
	if pandas.notnull(commalist):
		commalist = commalist.split(",") # Split comma delimited text into python list
		for x in commalist: # Iterate over list element
			for censor in censorlist: # Apply censor
				x = re.sub(censor, "", x)
				x = re.sub("^ ", "", x)
				x = re.sub(" $", "", x)
			for syn in combinedsyn: # Apply synonymy
				x = re.sub(syn[0], syn[1], x)
			if x.find('-')!=-1: # Make a list for sorting given any hyphen-delimited term -- alphabetize to standardize
				x = x.split("-")
				x.sort()
				x = list(filter(None, x)) # Remove null values from censoring
				x = '-'.join(x)
				y.append(x)
			else:
				y.append(x)
	y = list(filter(None, y)) # Remove completely null values if an item was completely removed
	y = list(set(y)) # Remove any potential duplicates
	y = ",".join(y)
	return y

# Process numerical data
def numcleaner(commalist):
	y = []
	if pandas.notnull(commalist):
		commalist = commalist.split(",") # Split comma delimited text into python list
		for x in commalist: # Iterate over list element
			x = re.sub("<=", "", x)
			x = re.sub(">=", "", x)
			x = re.sub("\([0-9.-]*\)", "", x)
			x = re.sub(" ", "", x)
			if x.find('-')!=-1: # Make a list for midpoint determination
				x = x.split("-")
				x = [float(i) for i in x]
				y.append(str(numpy.mean(x)))
			else:
				y.append(x)
	y = ",".join(y)
	return y

# Further process meristic data	
def integerify(commalist):
	y = []
	if pandas.notnull(commalist):
		commalist = commalist.split(",") # Split comma delimited text into python list -- was converted back into string above
		for x in commalist: # Iterate over list element
			if x.isdigit():
				x = int(float(x))
				y.append(x)
	return y
	
print(traitdataset.columns.values)

# Format trait columns
for i in ["Leaf shape", "Sepal shape", "Petal shape", "Inflorescence type"]:
	traitdataset[i].astype(str)
	traitdataset[i] = traitdataset[i].apply(cleaner)

# Format color columns 
for i in ["Sepal color", "Petal color"]:
	traitdataset[i].astype(str)
	traitdataset[i] = traitdataset[i].apply(cleaner)

# Format flowering time 
for i in ["Flowering time"]:
	traitdataset[i].astype(str)
	traitdataset[i] = traitdataset[i].apply(cleaner)

# Format measurement columns
for i in ["Leaf length", "Plant height", "Petiole length", "Leaf width", "Sepal length", "Sepal width", "Petal length", "Petal width", "Stamen length", "Style length", "Seed length"]:
	traitdataset[i].astype(str)
	traitdataset[i] = traitdataset[i].apply(numcleaner)

# Format meristic columns
for i in ["Stamen number", "Sepal number", "Petal number"]:
	traitdataset[i].astype(str)
	traitdataset[i] = traitdataset[i].apply(numcleaner) # Take midpoint of range, if a list of ranges, will be a list of the midpoints of each
	traitdataset[i] = traitdataset[i].apply(integerify) # Make sure count is integer (average)

# Write uncoded result
traitdataset.to_csv(path_or_buf=outfile, sep=",")

codesdone = []
# Code the data in alphabetical order sequentially from zero
for w in ["Flowering time", "Leaf shape", "Sepal shape", "Petal shape", "Inflorescence type", "Sepal color", "Petal color"]:
	codes = list(traitdataset[w])
	codes = [i.split(",") for i in codes] # Split commas
	codes = [item for i in codes for item in i] # Flatten list
	codes = list(filter(None, list(set(codes))))
	codes.sort()
	
	intermediatedataset = list(traitdataset[w])
	intermediatedataset = [i.split(",") for i in intermediatedataset] # Split commas
	replacement = []
	for i in intermediatedataset:
		i = list(filter(None, list(set(i))))
		z = [re.sub("".join(["^",j,"$"]), str(codes.index(j)), j) for j in i]
		replacement.append(z)
	traitdataset[w] = replacement
	codesdone.append([w, codes])

# Write coded result
traitdataset.to_csv(path_or_buf=outcoded, sep=",")

def commaseparate(value):
	value = str(value) # Force numbers as strings for regex and split
	if value.find(',')!=-1: 
		value = value.split(",")
	return value
	
for w in ["Woody", "Perennial"]:
	traitdataset[w] = traitdataset[w].apply(commaseparate)
	#print(traitdataset[w])

for i in ["Woody", "Perennial", "Flowering time", "Leaf shape", "Sepal shape", "Petal shape", "Inflorescence type", "Sepal color", "Petal color", "Stamen number", "Sepal number", "Petal number"]:
	column = traitdataset[i]
	newcolumn = []
	for cell in column:
		if cell: # Check if non-empty
			newcolumn.append(random.choice(cell))
		else: # Else leave it empty
			newcolumn.append(None)
	traitdataset[i] = newcolumn
	
traitdataset.to_csv(path_or_buf=outcoded_onevalue, sep=",")
for i in codesdone:
	filename = "".join([i[0], "_", codeguide])
	#print(i[1])
	n = 0
	with open(filename,'w+') as file:
		wr = csv.writer(file, dialect='excel')
		wr.writerow([i[0]])
		for x in i[1]:
			wr.writerow([x, n])
			n += 1

##############
# Now we prepare for distance matrix construction

# Set taxon as row label
traitdataset.set_index("TAXON", inplace=True)
	
# Clean up dataframe
for column in traitdataset.select_dtypes([numpy.object]): # Selecting objects needed for empty columns (?)
	traitdataset[column] = traitdataset[column].str.replace(r"[a-zA-Z?]",'')
	traitdataset[column] = traitdataset[column].str.replace(r"\s",'')
traitdataset = traitdataset.drop(["Elevation","Special characters","Notes","Name source","Habitat"], axis = 1) # Axis 1 for columns

# Replace empty values as 0 for continuous data -- correct arithmetic
for a in ["Leaf length", "Plant height", "Petiole length", "Leaf width", "Sepal length", "Sepal width", "Petal length", "Petal width", "Stamen length", "Style length", "Seed length"]:
	traitdataset[a].replace('', numpy.NaN, inplace=True)

# For discrete data, 0 is a valid datum, so use unoccupied number 
for a in ["Woody", "Perennial", "Flowering time", "Leaf shape", "Sepal shape", "Petal shape", "Inflorescence type", "Sepal color", "Petal color", "Stamen number", "Sepal number", "Petal number"]:
	traitdataset[a].replace('', numpy.NaN, inplace=True)

# Jaccard distance, i.e. identical is 0, maximally different is 1
traitdataset_matrix = traitdataset.as_matrix(["Woody", "Perennial", "Flowering time", "Leaf shape", "Sepal shape", "Petal shape", "Inflorescence type", "Sepal color", "Petal color", "Stamen number", "Sepal number", "Petal number"])

# Solution from SourceForge for screening missing data values for scipy: https://stackoverflow.com/questions/24781461/compute-the-pairwise-distance-in-scipy-with-missing-values
discrete_dissimilarity_matrix = scipy.spatial.distance.squareform(numpy.array([scipy.spatial.distance.pdist(traitdataset_matrix[s][:, ~pandas.isnull(traitdataset_matrix[s]).any(axis=0)], "jaccard") for s in map(list, itertools.combinations(range(traitdataset_matrix.shape[0]), 2))]).ravel())

discrete_dissimilarity = pandas.DataFrame(discrete_dissimilarity_matrix, columns=traitdataset.index.values, index=traitdataset.index.values)	

#print(discrete_dissimilarity)

# Standardized euclidean distance
traitdataset_matrix = traitdataset.as_matrix(["Leaf length", "Plant height", "Petiole length", "Leaf width", "Sepal length", "Sepal width", "Petal length", "Petal width", "Stamen length", "Style length", "Seed length"])

# Solution from SourceForge for screening missing data values for scipy: https://stackoverflow.com/questions/24781461/compute-the-pairwise-distance-in-scipy-with-missing-values
continuous_dissimilarity_matrix = scipy.spatial.distance.squareform(numpy.array([scipy.spatial.distance.pdist(traitdataset_matrix[s][:, ~pandas.isnull(traitdataset_matrix[s]).any(axis=0)], "euclidean") for s in map(list, itertools.combinations(range(traitdataset_matrix.shape[0]), 2))]).ravel())

continuous_dissimilarity_matrix = continuous_dissimilarity_matrix/numpy.max(continuous_dissimilarity_matrix) # Normalize to max 1
continuous_dissimilarity = pandas.DataFrame(continuous_dissimilarity_matrix, columns=traitdataset.index.values, index=traitdataset.index.values)	

#print(continuous_dissimilarity)

# Sum two matrices, then divide by two to yield metric ranging from 0 to 1
combined_dataframe = (discrete_dissimilarity.add(continuous_dissimilarity))
combined_dataframe = combined_dataframe.divide(2)

# For cases where missing data filtering results in no distance
for a in combined_dataframe:
	combined_dataframe[a].replace(numpy.NaN, 0, inplace=True)
for a in combined_dataframe:
	combined_dataframe[a].replace('', 0, inplace=True)

print(type(combined_dataframe))
combined_dataframe.to_csv(path_or_buf=distancematrix, sep=",")

tree_occ_matched_taxa = ["Adromischus_roaneanus", "Adromischus_montium-klinghardtii", "Adromischus_alstonii", "Adromischus_filicaulis", "Adromischus_maximus", "Adromischus_mammillaris", "Adromischus_liebenbergii", "Adromischus_triflorus", "Adromischus_sphenophyllus", "Adromischus_fallax", "Adromischus_umbraticola", "Adromischus_schuldtianus", "Adromischus_marianiae", "Adromischus_cristatus", "Adromischus_cooperi", "Adromischus_trigynus", "Adromischus_nanus", "Adromischus_subviridis", "Adromischus_leucophyllus", "Adromischus_hemisphaericus", "Adromischus_caryophyllaceus", "Adromischus_bicolor", "Adromischus_inamoenus", "Adromischus_maculatus", "Adromischus_humilis", "Adromischus_phillipsiae", "Cotyledon_tomentosa", "Cotyledon_woodii", "Cotyledon_velutina", "Cotyledon_barbeyi", "Cotyledon_campanulata", "Cotyledon_cuneata", "Cotyledon_papillaris", "Cotyledon_eliseae", "Cotyledon_adscendens", "Tylecodon_paniculatus", "Tylecodon_hirtifolius", "Tylecodon_wallichii", "Tylecodon_sulphureus", "Tylecodon_grandiflorus", "Tylecodon_reticulatus", "Tylecodon_suffultus", "Tylecodon_cacalioides", "Tylecodon_occultans", "Tylecodon_nolteei", "Tylecodon_decipiens", "Tylecodon_schaeferianus", "Tylecodon_aurusbergensis", "Tylecodon_buchholzianus", "Tylecodon_tribblei", "Tylecodon_similis", "Tylecodon_pearsonii", "Tylecodon_rubrovenosus", "Tylecodon_bayeri", "Tylecodon_ellaphieae", "Tylecodon_pusillus", "Tylecodon_torulosus", "Tylecodon_viridiflorus", "Tylecodon_pygmaeus", "Tylecodon_tenuis", "Tylecodon_leucothrix", "Tylecodon_hallii", "Tylecodon_kritzingeri", "Tylecodon_stenocaulis", "Tylecodon_ventricosus", "Tylecodon_albiflorus", "Tylecodon_striatus", "Tylecodon_singularis", "Kalanchoe_rhombopilosa", "Kalanchoe_orgyalis", "Kalanchoe_beharensis", "Kalanchoe_millotii", "Kalanchoe_arborescens", "Kalanchoe_tomentosa", "Kalanchoe_linearifolia", "Kalanchoe_eriophylla", "Kalanchoe_integrifolia", "Kalanchoe_pumila", "Kalanchoe_gracilipes", "Kalanchoe_thyrsiflora", "Kalanchoe_peltata", "Kalanchoe_uniflora", "Kalanchoe_porphyrocalyx", "Kalanchoe_fedtschenkoi", "Kalanchoe_delagoensis", "Kalanchoe_beauverdii", "Kalanchoe_daigremontiana", "Kalanchoe_marnieriana", "Kalanchoe_rosei", "Kalanchoe_waldheimii", "Kalanchoe_synsepala", "Kalanchoe_blossfeldiana", "Kalanchoe_latisepala", "Kalanchoe_lateritia", "Kalanchoe_marmorata", "Kalanchoe_mitejea", "Kalanchoe_humilis", "Kalanchoe_elizae", "Kalanchoe_longiflora", "Kalanchoe_prittwitzii", "Kalanchoe_farinacea", "Kalanchoe_scapigera", "Kalanchoe_germanae", "Kalanchoe_nyikae", "Kalanchoe_petitiana", "Kalanchoe_densiflora", "Kalanchoe_glaucescens", "Kalanchoe_lanceolata", "Kalanchoe_crenata", "Kalanchoe_velutina", "Kalanchoe_garambiensis", "Kalanchoe_bentii", "Kalanchoe_citrina", "Kalanchoe_rotundifolia", "Kalanchoe_luciae", "Kalanchoe_brachyloba", "Kalanchoe_paniculata", "Kalanchoe_jongmansii", "Kalanchoe_campanulata", "Kalanchoe_manginii", "Kalanchoe_pubescens", "Kalanchoe_miniata", "Kalanchoe_gastonis-bonnieri", "Kalanchoe_pinnata", "Kalanchoe_prolifera", "Kalanchoe_laciniata", "Kalanchoe_streptantha", "Pseudosedum_longidentatum", "Pseudosedum_condensatum", "Rhodiola_nobilis", "Rhodiola_dumulosa", "Rhodiola_kirilowii", "Rhodiola_semenovii", "Rhodiola_fastigiata", "Rhodiola_alsia", "Rhodiola_tibetica", "Rhodiola_ishidae", "Rhodiola_tangutica", "Rhodiola_rosea", "Rhodiola_heterodonta", "Rhodiola_yunnanensis", "Rhodiola_bupleuroides", "Rhodiola_atsaensis", "Rhodiola_quadrifida", "Rhodiola_coccinea", "Rhodiola_gelida", "Rhodiola_serrata", "Rhodiola_calliantha", "Rhodiola_chrysanthemifolia", "Rhodiola_hobsonii", "Rhodiola_sinuata", "Rhodiola_prainii", "Rhodiola_stapfii", "Rhodiola_wallichiana", "Rhodiola_algida", "Rhodiola_rhodantha", "Rhodiola_integrifolia", "Rhodiola_recticaulis", "Rhodiola_crenulata", "Rhodiola_purpureoviridis", "Rhodiola_sherriffii", "Rhodiola_pamiroalaica", "Rhodiola_macrocarpa", "Rhodiola_himalensis", "Pseudosedum_lievenii", "Rhodiola_humilis", "Rhodiola_cretinii", "Rhodiola_amabilis", "Rhodiola_imbricata", "Rhodiola_smithii", "Rhodiola_atuntsuensis", "Rhodiola_discolor", "Rhodiola_primuloides", "Rhodiola_pachyclados", "Phedimus_kamtschaticus", "Phedimus_aizoon", "Phedimus_sikokianus", "Phedimus_middendorffianus", "Phedimus_hybridus", "Phedimus_spurius", "Phedimus_stoloniferus", "Phedimus_obtusifolius", "Phedimus_stellatus", "Meterostachys_sikokianus", "Orostachys_spinosa", "Orostachys_thyrsiflora", "Orostachys_fimbriata", "Orostachys_japonica", "Orostachys_schoenlandii", "Orostachys_chanetii", "Sedum_correptum", "Sinocrassula_yunnanensis", "Sinocrassula_indica", "Sinocrassula_densirosulata", "Hylotelephium_erythrostictum", "Hylotelephium_tatarinowii", "Orostachys_paradoxa", "Orostachys_gorovoii", "Orostachys_iwarenge", "Orostachys_malacophylla", "Hylotelephium_mingjinianum", "Hylotelephium_caucasicum", "Hylotelephium_sordidum", "Hylotelephium_pallescens", "Sedum_diffusum", "Hylotelephium_telephium", "Hylotelephium_angustum", "Hylotelephium_telephioides", "Hylotelephium_anacampseros", "Hylotelephium_verticillatum", "Hylotelephium_subcapitatum", "Hylotelephium_sieboldii", "Hylotelephium_cauticola", "Hylotelephium_spectabile", "Hylotelephium_viride", "Hylotelephium_ewersii", "Hylotelephium_cyaneum", "Hylotelephium_populifolium", "Hylotelephium_pluricaule", "Umbilicus_tropaeolifolius", "Umbilicus_parviflorus", "Umbilicus_erectus", "Umbilicus_oppositifolius", "Umbilicus_heylandianus", "Umbilicus_horizontalis", "Umbilicus_rupestris", "Umbilicus_botryoides", "Diamorpha_smallii", "Sedum_caeruleum", "Sedum_modestum", "Sedum_jaccardianum", "Monanthes_atlantica", "Aeonium_percarneum", "Aeonium_balsamiferum", "Aeonium_undulatum", "Aeonium_leucoblepharum", "Aeonium_arboreum", "Aeonium_korneliuslemsii", "Aeonium_simsii", "Aeonium_lancerottense", "Aeonium_canariense", "Aeonium_castello-paivae", "Aeonium_nobile", "Aeonium_haworthii", "Aeonium_pseudourbicum", "Aeonium_urbicum", "Aeonium_spathulatum", "Aeonium_ciliatum", "Aeonium_decorum", "Aeonium_davidbramwellii", "Aeonium_hierrense", "Aeonium_gomerense", "Aeonium_glandulosum", "Aeonium_tabuliforme", "Aeonium_glutinosum", "Aeonium_goochiae", "Aeonium_saundersii", "Aeonium_lindleyi", "Aeonium_aizoon", "Aeonium_diplocyclum", "Aeonium_sedifolium", "Aeonium_aureum", "Aeonium_dodrantale", "Monanthes_icterica", "Aichryson_villosum", "Aichryson_laxum", "Aichryson_bollei", "Aichryson_palmense", "Aichryson_bethencourtianum", "Aichryson_tortuosum", "Aichryson_dumosum", "Aichryson_divaricatum", "Aichryson_pachycaulon", "Aichryson_parlatorei", "Aichryson_punctatum", "Aichryson_porphyrogennetos", "Monanthes_polyphylla", "Monanthes_brachycaulos", "Monanthes_minima", "Monanthes_anagensis", "Monanthes_laxiflora", "Monanthes_muralis", "Rosularia_serrata", "Rosularia_modesta", "Rosularia_rosulata", "Rosularia_adenotricha", "Rosularia_lineata", "Sedum_pallidum", "Sedum_eriocarpum", "Sedum_cepaea", "Sedum_apoleipon", "Sedum_ursi", "Sedum_laconicum", "Sedum_acre", "Sedum_lineare", "Sedum_onychopetalum", "Sedum_mexicanum", "Sedum_multicaule", "Sedum_tosaense", "Sedum_japonicum", "Sedum_oryzifolium", "Sedum_uniflorum", "Sedum_annuum", "Sedum_meyeri-johannis", "Sedum_ruwenzoriense", "Sedum_urvillei", "Sedum_grisebachii", "Sedum_alpestre", "Sedum_baileyi", "Sedum_alfredii", "Sedum_erythrospermum", "Sedum_nokoense", "Sedum_formosanum", "Sedum_emarginatum", "Sedum_hakonense", "Sedum_morrisonense", "Sedum_makinoi", "Sedum_yabeanum", "Sedum_polytrichoides", "Sedum_bulbiferum", "Sedum_sarmentosum", "Sedum_subtile", "Sedum_platysepalum", "Sedum_triactina", "Sedum_chauveaudii", "Sedum_celatum", "Sedum_przewalskii", "Sedum_tsinghaicum", "Sedum_roborowskii", "Sedum_forrestii", "Sedum_trullipetalum", "Sedum_gagei", "Sedum_glaebosum", "Sedum_fischeri", "Sedum_henrici-robertii", "Sedum_oreades", "Sedum_obtrullatum", "Sedum_bergeri", "Sedum_sexangulare", "Sedum_anglicum", "Sedum_cymatopetalum", "Sedum_andinum", "Sedum_reniforme", "Sedum_grandyi", "Sedum_humifusum", "Sedum_oaxacanum", "Sedum_obcordatum", "Sedum_fuscum", "Sedum_vinicolor", "Sedum_pacense", "Lenophyllum_guttatum", "Lenophyllum_acutifolium", "Sedum_carinatifolium", "Sedum_bourgaei", "Sedum_chloropetalum", "Sedum_quevae", "Sedum_griseum", "Sedum_retusum", "Villadia_pringlei", "Villadia_cucullata", "Villadia_aristata", "Sedum_liebmannianum", "Sedum_wrightii", "Sedum_cockerellii", "Sedum_trichromum", "Sedum_alamosanum", "Sedum_allantoides", "Sedum_alexanderi", "Sedum_moranense", "Sedum_ulricae", "Sedum_obtusipetalum", "Sedum_jurgensenii", "Sedum_batesii", "Villadia_guatemalensis", "Villadia_misera", "Villadia_minutiflora", "Villadia_diffusa", "Villadia_imbricata", "Villadia_recurva", "Villadia_nelsonii", "Villadia_albiflora", "Villadia_incarum", "Sedum_versadense", "Sedum_stelliforme", "Sedum_jaliscanum", "Sedum_adolphii", "Sedum_commixtum", "Sedum_treleasei", "Sedum_furfuraceum", "Sedum_clavatum", "Sedum_corynephyllum", "Sedum_morganianum", "Sedum_nussbaumerianum", "Echeveria_grisea", "Pachyphytum_caesium", "Echeveria_australis", "Pachyphytum_fittkaui", "Pachyphytum_viride", "Pachyphytum_kimnachii", "Pachyphytum_glutinicaule", "Pachyphytum_hookeri", "Pachyphytum_oviferum", "Pachyphytum_compactum", "Echeveria_nuda", "Echeveria_setosa", "Echeveria_pringlei", "Echeveria_coccinea", "Echeveria_longissima", "Echeveria_purpusorum", "Echeveria_pulvinata", "Echeveria_recurvata", "Echeveria_montana", "Echeveria_megacalyx", "Echeveria_racemosa", "Echeveria_nodulosa", "Echeveria_colorata", "Echeveria_simulans", "Echeveria_chihuahuaensis", "Echeveria_rosea", "Echeveria_prolifica", "Echeveria_amoena", "Graptopetalum_filiferum", "Graptopetalum_rusbyi", "Graptopetalum_bartramii", "Graptopetalum_pachyphyllum", "Echeveria_fulgens", "Thompsonella_mixtecana", "Thompsonella_spathulata", "Thompsonella_minutiflora", "Thompsonella_nellydiegoae", "Thompsonella_platyphylla", "Thompsonella_xochipalensis", "Thompsonella_garcia-mendozae", "Graptopetalum_bellum", "Echeveria_strictiflora", "Echeveria_bifida", "Echeveria_cante", "Echeveria_gibbiflora", "Echeveria_juarezensis", "Cremnophila_nutans", "Cremnophila_linguifolia", "Graptopetalum_macdougallii", "Graptopetalum_pentandrum", "Graptopetalum_fruticosum", "Graptopetalum_grande", "Graptopetalum_amethystinum", "Graptopetalum_mendozae", "Graptopetalum_paraguayense", "Sedum_dendroideum", "Sedum_hultenii", "Sedum_caducum", "Sedum_gypsophilum", "Sedum_ternatum", "Sedum_catorce", "Sedum_palmeri", "Sedum_hintonii", "Sedum_chazaroi", "Sedum_ebracteatum", "Sedum_hemsleyanum", "Sedum_greggii", "Sedum_frutescens", "Sedum_longipes", "Sedum_oxypetalum", "Sedum_guatemalense", "Sedum_lancerottense", "Sedum_fusiforme", "Sedum_nudum", "Sedum_stefco", "Sedum_elatinoides", "Sedum_glaucophyllum", "Sedum_magellense", "Sedum_caespitosum", "Prometheum_chrysanthum", "Sedum_fragrans", "Sedum_hispanicum", "Sedum_brevifolium", "Sedum_atratum", "Sedum_gypsicola", "Sedum_album", "Sedum_dasyphyllum", "Sedum_hirsutum", "Sedum_gracile", "Sedum_subulatum", "Rosularia_platyphylla", "Sedum_obtusatum", "Sedum_oreganum", "Sedum_oregonense", "Sedum_spathulifolium", "Sedum_debile", "Sedum_divergens", "Sedum_lanceolatum", "Sedum_stenopetalum", "Sedella_pumila", "Dudleya_greenei", "Dudleya_candelabrum", "Dudleya_acuminata", "Dudleya_ingens", "Dudleya_anthonyi", "Dudleya_attenuata", "Dudleya_anomala", "Dudleya_linearis", "Dudleya_pachyphytum", "Dudleya_brittonii", "Dudleya_formosa", "Dudleya_edulis", "Dudleya_pulverulenta", "Dudleya_blochmaniae", "Dudleya_verityi", "Dudleya_variegata", "Dudleya_saxosa", "Dudleya_cymosa", "Dudleya_guadalupensis", "Dudleya_virens", "Dudleya_multicaulis", "Dudleya_farinosa", "Dudleya_densiflora", "Dudleya_abramsii", "Dudleya_stolonifera", "Dudleya_albiflora", "Sedum_microcarpum", "Sedum_tenellum", "Prometheum_sempervivoides", "Prometheum_pilosum", "Sedum_candollei", "Sedum_mucizonia", "Sedum_lagascae", "Sedum_villosum", "Pistorinia_breviflora", "Pistorinia_hispanica", "Rosularia_glabra", "Rosularia_alpestris", "Sempervivum_globiferum", "Sempervivum_heuffelii", "Sempervivum_cantabricum", "Sempervivum_pittonii", "Sempervivum_wulfenii", "Sempervivum_brevipilum", "Sempervivum_atropatanum", "Sempervivum_arachnoideum", "Sempervivum_grandiflorum", "Sempervivum_transcaucasicum", "Sempervivum_pumilum", "Sempervivum_davisii", "Sempervivum_armenum", "Sempervivum_caucasicum", "Sempervivum_minutum", "Sempervivum_calcareum", "Sempervivum_marmoreum", "Sempervivum_leucanthum", "Sempervivum_iranicum", "Sempervivum_ruthenicum", "Sempervivum_tectorum", "Sempervivum_giuseppii", "Sempervivum_montanum", "Sedum_sediforme", "Sedum_amplexicaule", "Sedum_rupestre", "Sedum_ochroleucum", "Sedum_montanum", "Sedum_pruinatum", "Sedum_forsterianum", "Sedum_churchillianum", "Sedum_drymarioides", "Sedum_stellariifolium", "Crassula_umbraticola", "Crassula_dichotoma", "Crassula_muscosa", "Crassula_manaia", "Crassula_sieberiana", "Crassula_mataikona", "Crassula_colligata", "Crassula_schimperi", "Crassula_campestris", "Crassula_thunbergiana", "Crassula_alata", "Crassula_tillaea", "Crassula_connata", "Crassula_glomerata", "Crassula_vaginata", "Crassula_setulosa", "Crassula_natalensis", "Crassula_sediflora", "Crassula_orbicularis", "Crassula_perfoliata", "Crassula_perforata", "Crassula_deceptor", "Crassula_fascicularis", "Crassula_rupestris", "Crassula_ausensis", "Crassula_cultrata", "Crassula_nudicaulis", "Crassula_clavata", "Crassula_dejecta", "Crassula_barklyi", "Crassula_globularioides", "Crassula_montana", "Crassula_capitella", "Crassula_pubescens", "Crassula_compacta", "Crassula_pseudohemisphaerica", "Crassula_brachystachya", "Crassula_flava", "Crassula_atropurpurea", "Crassula_crenulata", "Crassula_multicava", "Crassula_sarmentosa", "Crassula_inandensis", "Crassula_spathulata", "Crassula_strigosa", "Crassula_decumbens", "Crassula_vaillantii", "Crassula_papillosa", "Crassula_tenuicaulis", "Crassula_pellucida", "Crassula_alticola", "Crassula_ericoides", "Crassula_sarcocaulis", "Crassula_dependens", "Crassula_expansa", "Crassula_volkensii", "Crassula_peploides", "Crassula_gemmifera", "Crassula_multicaulis", "Crassula_kirkii", "Crassula_ruamahanga", "Crassula_hunua", "Crassula_sinclairii", "Crassula_moschata", "Crassula_inanis", "Crassula_natans", "Crassula_helmsii", "Crassula_granvikii", "Crassula_aquatica", "Crassula_peduncularis", "Crassula_venezuelensis", "Tetracarpaea_tasmanica", "Aphanopetalum_clematideum", "Aphanopetalum_resinosum", "Penthorum_sedoides", "Penthorum_chinense", "Meionectes_brownii", "Meionectes_tenuifolia", "Gonocarpus_eremophilus", "Haloragis_hamata", "Haloragis_foliosa", "Haloragis_acutangula", "Haloragis_heterophylla", "Haloragis_serra", "Haloragis_aspera", "Haloragis_erecta", "Haloragis_masatierrana", "Haloragis_stricta", "Haloragis_exalata", "Haloragis_glauca", "Haloragis_eyreana", "Haloragis_dura", "Haloragis_odontocarpa", "Haloragis_trigonocarpa", "Gonocarpus_montanus", "Gonocarpus_micranthus", "Gonocarpus_ephemerus", "Gonocarpus_acanthocarpus", "Gonocarpus_leptothecus", "Gonocarpus_chinensis", "Gonocarpus_elatus", "Gonocarpus_longifolius", "Gonocarpus_teucrioides", "Gonocarpus_mezianus", "Gonocarpus_humilis", "Gonocarpus_oreophilus", "Gonocarpus_aggregatus", "Gonocarpus_salsoloides", "Gonocarpus_tetragynus", "Gonocarpus_nodulosus", "Gonocarpus_scordioides", "Gonocarpus_urceolatus", "Gonocarpus_pusillus", "Gonocarpus_paniculatus", "Gonocarpus_benthamii", "Gonocarpus_intricatus", "Gonocarpus_cordiger", "Gonocarpus_trichostachyus", "Trihaloragis_hexandra", "Laurembergia_repens", "Myriophyllum_trifidum", "Myriophyllum_amphibium", "Myriophyllum_tillaeoides", "Myriophyllum_pedunculatum", "Myriophyllum_lapidicola", "Myriophyllum_muricatum", "Myriophyllum_dicoccum", "Myriophyllum_latifolium", "Myriophyllum_humile", "Myriophyllum_laxum", "Myriophyllum_tenellum", "Myriophyllum_farwellii", "Myriophyllum_hippuroides", "Myriophyllum_heterophyllum", "Myriophyllum_pinnatum", "Myriophyllum_trachycarpum", "Myriophyllum_coronatum", "Myriophyllum_petraeum", "Myriophyllum_propinquum", "Myriophyllum_papillosum", "Myriophyllum_crispatum", "Myriophyllum_ussuriense", "Myriophyllum_variifolium", "Myriophyllum_votschii", "Myriophyllum_lophatum", "Myriophyllum_alpinum", "Myriophyllum_simulans", "Myriophyllum_echinatum", "Myriophyllum_drummondii", "Myriophyllum_limnophilum", "Myriophyllum_filiforme", "Myriophyllum_mattogrossensis", "Myriophyllum_aquaticum", "Myriophyllum_alterniflorum", "Myriophyllum_verticillatum", "Myriophyllum_robustum", "Myriophyllum_triphyllum", "Myriophyllum_quitense", "Myriophyllum_verrucosum", "Myriophyllum_sibiricum", "Myriophyllum_spicatum", "Myriophyllum_salsugineum", "Myriophyllum_balladoniense", "Myriophyllum_caput-medusae", "Myriophyllum_decussatum", "Proserpinaca_palustris", "Proserpinaca_pectinata", "Glischrocaryon_flavescens", "Glischrocaryon_roei", "Glischrocaryon_aureum", "Glischrocaryon_behrii", "Haloragodendron_racemosum", "Haloragodendron_glandulosum", "Haloragodendron_monospermum", "Haloragodendron_baeuerlenii", "Haloragodendron_gibsonii", "Haloragodendron_lucasii", "Cynomorium_coccineum", "Cynomorium_songaricum", "Pterostemon_mexicanus", "Pterostemon_rotundifolius", "Itea_virginica", "Itea_ilicifolia", "Itea_yunnanensis", "Itea_macrophylla", "Itea_japonica", "Itea_kiukiangensis", "Itea_oldhamii", "Itea_omeiensis", "Itea_chinensis", "Itea_glutinosa", "Itea_parviflora", "Choristylis_rhamnoides", "Ribes_giraldii", "Ribes_humile", "Ribes_fasciculatum", "Ribes_orientale", "Ribes_maximowiczianum", "Ribes_laciniatum", "Ribes_komarovii", "Ribes_villosum", "Ribes_rubrisepalum", "Ribes_spicatum", "Ribes_maximowiczii", "Ribes_diacanthum", "Ribes_hirtellum", "Ribes_alpinum", "Ribes_mandshuricum", "Ribes_vilmorinii", "Ribes_rubrum", "Ribes_tenue", "Ribes_janczewskii", "Ribes_petraeum", "Ribes_multiflorum", "Ribes_himalense", "Ribes_altissimum", "Ribes_triste", "Ribes_takare", "Ribes_meyeri", "Ribes_glaciale", "Ribes_griffithii", "Ribes_moupinense", "Ribes_longeracemosum", "Ribes_velutinum", "Ribes_lobbii", "Ribes_roezlii", "Ribes_thacherianum", "Ribes_amarum", "Ribes_menziesii", "Ribes_binominatum", "Ribes_californicum", "Ribes_burejense", "Ribes_sinanense", "Ribes_stenocarpum", "Ribes_formosanum", "Ribes_aciculare", "Ribes_graveolens", "Ribes_rotundifolium", "Ribes_uva-crispa", "Ribes_alpestre", "Ribes_divaricatum", "Ribes_missouriense", "Ribes_oxyacanthoides", "Ribes_inerme", "Ribes_curvatum", "Ribes_watsonianum", "Ribes_speciosum", "Ribes_cynosbati", "Ribes_niveum", "Ribes_echinellum", "Ribes_montigenum", "Ribes_lacustre", "Ribes_andicola", "Ribes_valdivianum", "Ribes_ovalifolium", "Ribes_wolfii", "Ribes_cereum", "Ribes_affine", "Ribes_ciliatum", "Ribes_dugesii", "Ribes_horridum", "Ribes_ceriferum", "Ribes_brandegeei", "Ribes_erythrocarpum", "Ribes_malvaceum", "Ribes_indecorum", "Ribes_sanguineum", "Ribes_viscosissimum", "Ribes_aureum", "Ribes_glandulosum", "Ribes_howellii", "Ribes_laxiflorum", "Ribes_bracteosum", "Ribes_hudsonianum", "Ribes_nigrum", "Ribes_densiflorum", "Ribes_americanum", "Ribes_tortuosum", "Ribes_viburnifolium", "Peltoboykinia_tellimoides", "Chrysosplenium_griffithii", "Chrysosplenium_macrophyllum", "Chrysosplenium_lanuginosum", "Chrysosplenium_davidianum", "Chrysosplenium_carnosum", "Chrysosplenium_axillare", "Chrysosplenium_forrestii", "Chrysosplenium_nudicaule", "Chrysosplenium_rosendahlii", "Chrysosplenium_alternifolium", "Chrysosplenium_iowense", "Chrysosplenium_wrightii", "Chrysosplenium_tetrandrum", "Chrysosplenium_japonicum", "Chrysosplenium_hydrocotylifolium", "Chrysosplenium_flagelliferum", "Chrysosplenium_tosaense", "Chrysosplenium_uniflorum", "Chrysosplenium_microspermum", "Chrysosplenium_biondianum", "Chrysosplenium_valdivicum", "Chrysosplenium_oppositifolium", "Chrysosplenium_americanum", "Chrysosplenium_dubium", "Chrysosplenium_nepalense", "Chrysosplenium_sinicum", "Chrysosplenium_fauriae", "Chrysosplenium_echinus", "Chrysosplenium_kiotense", "Chrysosplenium_macrostemon", "Chrysosplenium_grayanum", "Chrysosplenium_pilosum", "Chrysosplenium_rhabdospermum", "Chrysosplenium_album", "Chrysosplenium_kamtschaticum", "Chrysosplenium_ramosum", "Micranthes_merkii", "Micranthes_tolmiei", "Micranthes_petiolaris", "Micranthes_bryophora", "Micranthes_ferruginea", "Micranthes_laciniata", "Micranthes_redofskyi", "Micranthes_foliolosa", "Micranthes_stellaris", "Micranthes_clusii", "Micranthes_japonica", "Micranthes_odontoloma", "Micranthes_calycina", "Micranthes_razshivinii", "Micranthes_spicata", "Micranthes_fusca", "Micranthes_manchuriensis", "Micranthes_nelsoniana", "Micranthes_reflexa", "Micranthes_micranthidifolia", "Micranthes_pensylvanica", "Micranthes_subapetala", "Micranthes_lyallii", "Micranthes_integrifolia", "Micranthes_eriophora", "Micranthes_rhomboidea", "Micranthes_occidentalis", "Micranthes_tempestiva", "Micranthes_rufidula", "Micranthes_aprica", "Micranthes_howellii", "Micranthes_californica", "Micranthes_apetala", "Micranthes_nidifica", "Micranthes_texana", "Micranthes_caroliniana", "Micranthes_careyana", "Micranthes_virginiensis", "Micranthes_tenuis", "Micranthes_nivalis", "Micranthes_oregana", "Micranthes_hieraciifolia", "Micranthes_davidii", "Micranthes_gageana", "Micranthes_atrata", "Micranthes_clavistaminea", "Micranthes_pallida", "Micranthes_lumpuensis", "Micranthes_melanocentra", "Micranthes_sachalinensis", "Mitella_integripetala", "Elmera_racemosa", "Mitella_caulescens", "Mitella_pentandra", "Mitella_pauciflora", "Mitella_furusei", "Mitella_yoshinagae", "Mitella_formosana", "Mitella_stylosa", "Mitella_diphylla", "Mitella_nuda", "Mitella_ovalis", "Mitella_breweri", "Mitella_japonica", "Mitella_kiusiana", "Mitella_koshiensis", "Mitella_acerina", "Mitella_diversifolia", "Mitella_stauropetala", "Mitella_trifida", "Conimitella_williamsii", "Tolmiea_menziesii", "Bensoniella_oregona", "Lithophragma_heterophyllum", "Lithophragma_cymbalaria", "Lithophragma_maximum", "Lithophragma_affine", "Lithophragma_parviflorum", "Lithophragma_glabrum", "Lithophragma_tenellum", "Lithophragma_trifoliatum", "Lithophragma_bolanderi", "Lithophragma_campanulatum", "Tellima_grandiflora", "Heuchera_grossulariifolia", "Heuchera_lakelae", "Heuchera_longipetala", "Heuchera_mexicana", "Heuchera_caespitosa", "Heuchera_hirsutissima", "Heuchera_abramsii", "Heuchera_elegans", "Heuchera_wellsiae", "Heuchera_rubescens", "Heuchera_rosendahlii", "Heuchera_brevistaminea", "Heuchera_sanguinea", "Heuchera_bracteata", "Heuchera_hallii", "Heuchera_pulchella", "Heuchera_acutifolia", "Heuchera_woodsiaphila", "Heuchera_parishii", "Heuchera_merriamii", "Heuchera_villosa", "Heuchera_missouriensis", "Heuchera_parviflora", "Heuchera_puberula", "Heuchera_cylindrica", "Heuchera_chlorantha", "Heuchera_micrantha", "Heuchera_pilosissima", "Heuchera_maxima", "Heuchera_pubescens", "Heuchera_americana", "Heuchera_alba", "Heuchera_caroliniana", "Heuchera_richardsonii", "Heuchera_longiflora", "Heuchera_wootonii", "Heuchera_parvifolia", "Heuchera_eastwoodiae", "Heuchera_glomerulata", "Heuchera_novomexicana", "Heuchera_glabra", "Tiarella_trifoliata", "Tiarella_cordifolia", "Tiarella_polyphylla", "Tanakaea_radicans", "Leptarrhena_pyrolifolia", "Saniculiphyllum_guanxiense", "Astilbe_chinensis", "Astilbe_rubra", "Astilbe_microphylla", "Astilbe_japonica", "Astilbe_shikokiana", "Astilbe_grandis", "Astilbe_philippinensis", "Astilbe_longicarpa", "Astilbe_macroflora", "Astilbe_macrocarpa", "Astilbe_simplicifolia", "Astilbe_thunbergii", "Astilbe_odontophylla", "Astilbe_biternata", "Astilbe_rivularis", "Saxifragopsis_fragarioides", "Sullivantia_sullivantii", "Sullivantia_oregana", "Sullivantia_hapemanii", "Boykinia_rotundifolia", "Boykinia_occidentalis", "Boykinia_major", "Boykinia_intermedia", "Boykinia_aconitifolia", "Boykinia_lycoctonifolia", "Suksdorfia_violacea", "Suksdorfia_ranunculifolia", "Bolandra_californica", "Bolandra_oregana", "Hieronymusia_alchemilloides", "Boykinia_richardsonii", "Telesonix_heucheriformis", "Telesonix_jamesii", "Jepsonia_malvifolia", "Jepsonia_parryi", "Rodgersia_podophylla", "Rodgersia_pinnata", "Rodgersia_aesculifolia", "Rodgersia_sambucifolia", "Astilboides_tabularis", "Darmera_peltata", "Rodgersia_nepalensis", "Oresitrophe_rupifraga", "Mukdenia_rossii", "Bergenia_crassifolia", "Bergenia_purpurascens", "Bergenia_emeiensis", "Bergenia_pacumbis", "Bergenia_ciliata", "Bergenia_stracheyi", "Saxifraga_mertensiana", "Saxifraga_fortunei", "Saxifraga_mengtzeana", "Saxifraga_sendaica", "Saxifraga_rufescens", "Saxifraga_imparilis", "Saxifraga_stolonifera", "Saxifraga_nipponica", "Saxifragella_bicuspidata", "Saxifraga_bronchialis", "Saxifraga_taylorii", "Saxifraga_tricuspidata", "Saxifraga_kruhsiana", "Saxifraga_spinulosa", "Saxifraga_firma", "Saxifraga_rebunshirensis", "Saxifraga_caulescens", "Saxifraga_nishidae", "Saxifraga_omolojensis", "Saxifraga_cherlerioides", "Saxifraga_vespertina", "Saxifraga_granulifera", "Saxifraga_sibirica", "Saxifraga_cernua", "Saxifraga_carpatica", "Saxifraga_bracteata", "Saxifraga_radiata", "Saxifraga_rivularis", "Saxifraga_hyperborea", "Saxifraga_taygetea", "Saxifraga_roylei", "Saxifraga_repanda", "Saxifraga_rotundifolia", "Saxifraga_squarrosa", "Saxifraga_caesia", "Saxifraga_bryoides", "Saxifraga_aspera", "Saxifraga_mutata", "Saxifraga_aizoides", "Saxifraga_rudolphiana", "Saxifraga_biflora", "Saxifraga_oppositifolia", "Saxifraga_retusa", "Saxifraga_duthiei", "Saxifraga_pulvinaria", "Saxifraga_likiangensis", "Saxifraga_pulchra", "Saxifraga_chionophila", "Saxifraga_poluniniana", "Saxifraga_alpigena", "Saxifraga_lilacina", "Saxifraga_stolitzkae", "Saxifraga_andersonii", "Saxifraga_decora", "Saxifraga_wendelboi", "Saxifraga_diapensioides", "Saxifraga_tombeanensis", "Saxifraga_aretioides", "Saxifraga_felineri", "Saxifraga_ferdinandi-coburgi", "Saxifraga_burseriana", "Saxifraga_subsessiliflora", "Saxifraga_georgei", "Saxifraga_afghanica", "Saxifraga_lowndesii", "Saxifraga_kotschyi", "Saxifraga_marginata", "Saxifraga_unguipetala", "Saxifraga_scardica", "Saxifraga_scleropoda", "Saxifraga_subverticillata", "Saxifraga_sancta", "Saxifraga_federici-augusti", "Saxifraga_sempervivum", "Saxifraga_porophylla", "Saxifraga_corymbosa", "Saxifraga_stribrnyi", "Saxifraga_florulenta", "Saxifraga_cochlearis", "Saxifraga_paniculata", "Saxifraga_crustata", "Saxifraga_callosa", "Saxifraga_cartilaginea", "Saxifraga_kolenatiana", "Saxifraga_hostii", "Saxifraga_longifolia", "Saxifraga_valdensis", "Saxifraga_cotyledon", "Saxifraga_cuneifolia", "Saxifraga_umbrosa", "Saxifraga_hirsuta", "Saxifraga_spathularis", "Saxifraga_tridactylites", "Saxifraga_osloensis", "Saxifraga_wahlenbergii", "Saxifraga_adscendens", "Saxifraga_conifera", "Saxifraga_aquatica", "Saxifraga_hariotii", "Saxifraga_intricata", "Saxifraga_fragilis", "Saxifraga_maderensis", "Saxifraga_babiana", "Saxifraga_trifurcata", "Saxifraga_pedemontana", "Saxifraga_cebennensis", "Saxifraga_cespitosa", "Saxifraga_granulata", "Saxifraga_sedoides", "Saxifraga_reuteriana", "Saxifraga_globulifera", "Saxifraga_corsica", "Saxifraga_erioblasta", "Saxifraga_rigoi", "Saxifraga_praetermissa", "Saxifraga_exarata", "Saxifraga_cintrana", "Saxifraga_hypnoides", "Saxifraga_magellanica", "Saxifraga_bourgaeana", "Saxifraga_gemmulosa", "Saxifraga_bulbifera", "Saxifraga_carpetana", "Saxifraga_dichotoma", "Saxifraga_latepetiolata", "Saxifraga_biternata", "Saxifraga_haenseleri", "Saxifraga_geranioides", "Saxifraga_nevadensis", "Saxifraga_nervosa", "Saxifraga_pentadactylis", "Saxifraga_losae", "Saxifraga_genesiana", "Saxifraga_pubescens", "Saxifraga_moncayensis", "Saxifraga_vayredana", "Saxifraga_rosacea", "Saxifraga_camposii", "Saxifraga_canaliculata", "Saxifraga_irrigua", "Saxifraga_seguieri", "Saxifraga_androsacea", "Saxifraga_styriaca", "Saxifraga_paradoxa", "Saxifraga_tenella", "Saxifraga_facchinii", "Saxifraga_presolanensis", "Saxifraga_aphylla", "Saxifraga_arachnoidea", "Saxifraga_cuneata", "Saxifraga_continentalis", "Saxifraga_sibthorpii", "Saxifraga_hederacea", "Saxifraga_cymbalaria", "Saxifraga_eschscholtzii", "Saxifraga_substrigosa", "Saxifraga_hispidula", "Saxifraga_strigosa", "Saxifraga_llonakhensis", "Saxifraga_balfourii", "Saxifraga_oreophila", "Saxifraga_filicaulis", "Saxifraga_gemmipara", "Saxifraga_gouldii", "Saxifraga_hypericoides", "Saxifraga_wallichiana", "Saxifraga_wardii", "Saxifraga_cinerascens", "Saxifraga_brachypoda", "Saxifraga_macrostigmatoides", "Saxifraga_versicallosa", "Saxifraga_brunonis", "Saxifraga_mucronulata", "Saxifraga_pilifera", "Saxifraga_consanguinea", "Saxifraga_angustata", "Saxifraga_microgyna", "Saxifraga_hemisphaerica", "Saxifraga_mucronulatoides", "Saxifraga_flagellaris", "Saxifraga_serpyllifolia", "Saxifraga_chrysanthoides", "Saxifraga_filifolia", "Saxifraga_aleutica", "Saxifraga_stella-aurea", "Saxifraga_signata", "Saxifraga_finitima", "Saxifraga_jacquemontiana", "Saxifraga_sanguinea", "Saxifraga_chrysantha", "Saxifraga_sediformis", "Saxifraga_candelabrum", "Saxifraga_dielsiana", "Saxifraga_heterotricha", "Saxifraga_signatella", "Saxifraga_brunneopunctata", "Saxifraga_punctulata", "Saxifraga_petraea", "Saxifraga_uninervia", "Saxifraga_glacialis", "Saxifraga_flexilis", "Saxifraga_tatsienluensis", "Saxifraga_unguiculata", "Saxifraga_aurantiaca", "Saxifraga_gemmigera", "Saxifraga_microphylla", "Saxifraga_engleriana", "Saxifraga_prattii", "Saxifraga_nanella", "Saxifraga_pasumensis", "Saxifraga_drabiformis", "Saxifraga_umbellulata", "Saxifraga_perpusilla", "Saxifraga_diapensia", "Saxifraga_isophylla", "Saxifraga_kingdonii", "Saxifraga_kingiana", "Saxifraga_lychnitis", "Saxifraga_oresbia", "Saxifraga_heleonastes", "Saxifraga_montanella", "Saxifraga_ciliatopetala", "Saxifraga_parva", "Saxifraga_tangutica", "Saxifraga_przewalskii", "Saxifraga_forrestii", "Saxifraga_cordigera", "Saxifraga_nigroglandulifera", "Saxifraga_auriculata", "Saxifraga_dianxibeiensis", "Saxifraga_gedangensis", "Saxifraga_cardiophylla", "Saxifraga_omphalodifolia", "Saxifraga_glaucophylla", "Saxifraga_insolens", "Saxifraga_implicans", "Saxifraga_subaequifoliata", "Saxifraga_yarlungzangboensis", "Saxifraga_densifoliata", "Saxifraga_litangensis", "Saxifraga_yezhiensis", "Saxifraga_peplidifolia", "Saxifraga_tsangchanensis", "Saxifraga_elliptica", "Saxifraga_caveana", "Saxifraga_culcitosa", "Saxifraga_linearifolia", "Saxifraga_giraldiana", "Saxifraga_egregia", "Saxifraga_hirculus", "Saxifraga_sphaeradena", "Saxifraga_parnassifolia", "Saxifraga_moorcroftiana", "Saxifraga_hookeri", "Saxifraga_diversifolia", "Saxifraga_sinomontana", "Saxifraga_pseudohirculus", "Saxifraga_saginoides", "Saxifraga_gonggashanensis", "Saxifraga_triaristulata", "Saxifraga_brachyphylla", "Saxifraga_hengduanensis", "Saxifraga_hirculoides", "Saxifraga_tibetica", "Saxifraga_aristulata", "Saxifraga_bergenioides", "Saxifraga_congestiflora", "Saxifraga_bulleyana", "Saxifraga_zimmermannii", "Saxifragodes_albowiana", "Cascadia_nuttallii", "Soyauxia_floribunda", "Soyauxia_ledermannii", "Soyauxia_talbotii", "Whittonia_guianensis", "Peridiscus_lucidus", "Medusandra_richardsiana", "Medusandra_mpomiana", "Paeonia_morisii", "Paeonia_obovata", "Paeonia_parnassica", "Paeonia_arietina", "Paeonia_intermedia", "Paeonia_tenuifolia", "Paeonia_lactiflora", "Paeonia_officinalis", "Paeonia_mascula", "Paeonia_cambessedesii", "Paeonia_clusii", "Paeonia_broteri", "Paeonia_coriacea", "Paeonia_mlokosewitschii", "Paeonia_peregrina", "Paeonia_veitchii", "Paeonia_emodi", "Paeonia_anomala", "Paeonia_mairei", "Paeonia_wittmanniana", "Paeonia_brownii", "Paeonia_californica", "Paeonia_decomposita", "Paeonia_ostii", "Paeonia_jishanensis", "Paeonia_szechuanica", "Paeonia_rockii", "Paeonia_delavayi", "Paeonia_ludlowii", "Paeonia_macrophylla", "Paeonia_caucasica", "Paeonia_daurica", "Daphniphyllum_laurinum", "Daphniphyllum_oldhamii", "Daphniphyllum_scortechinii", "Daphniphyllum_glaucescens", "Daphniphyllum_himalayense", "Daphniphyllum_teijsmannii", "Daphniphyllum_gracile", "Daphniphyllum_paxianum", "Daphniphyllum_phillipinensis", "Daphniphyllum_calycinum", "Daphniphyllum_borneense", "Daphniphyllum_macropodum", "Cercidiphyllum_magnificum", "Cercidiphyllum_japonicum", "Exbucklandia_longipetala", "Exbucklandia_tonkinensis", "Exbucklandia_populnea", "Rhodoleia_henryi", "Rhodoleia_championii", "Chunia_bucklandioides", "Mytilaria_laosensis", "Disanthus_cercidifolius", "Dicoryphe_stipulacea", "Trichocladus_grandiflorus", "Trichocladus_ellipticus", "Trichocladus_crinitus", "Neostrearia_fleckeri", "Ostrearia_australiana", "Noahdendron_nicholasii", "Molinadendron_guatemalense", "Molinadendron_sinaloense", "Sinowilsonia_henryi", "Fortunearia_sinensis", "Eustigma_balansae", "Eustigma_oblongifolium", "Eustigma_honbaense", "Hamamelis_japonica", "Hamamelis_mollis", "Hamamelis_vernalis", "Hamamelis_mexicana", "Hamamelis_ovalis", "Hamamelis_virginiana", "Parrotiopsis_jacquemontiana", "Distyliopsis_tutcheri", "Distylium_myricoides", "Distylium_racemosum", "Distyliopsis_dunnii", "Sycopsis_sinensis", "Parrotia_subaequalis", "Parrotia_persica", "Distyliopsis_laurifolia", "Fothergilla_gardenii", "Fothergilla_major", "Loropetalum_subcordatum", "Loropetalum_chinense", "Matudaea_trinervia", "Matudaea_colombiana", "Corylopsis_omeiensis", "Corylopsis_yunnanensis", "Corylopsis_multiflora", "Corylopsis_stenopetala", "Corylopsis_pauciflora", "Corylopsis_platypetala", "Corylopsis_sinensis", "Corylopsis_veitchiana", "Corylopsis_glandulifera", "Corylopsis_willmottiae", "Corylopsis_spicata", "Corylopsis_glabrescens", "Corylopsis_microcarpa", "Corylopsis_himalayana", "Liquidambar_orientalis", "Liquidambar_styraciflua", "Altingia_siamensis", "Altingia_excelsa", "Altingia_yunnanensis", "Altingia_poilanei", "Altingia_chinensis", "Altingia_gracilipes", "Altingia_obovata", "Liquidambar_acalycina", "Liquidambar_formosana", "Semiliquidambar_chingii", "Semiliquidambar_cathayensis"]

occurrence_not_trait = ["Orostachys_gorovoii", "Aeonium_pseudourbicum", "Sedum_yabeanum", "Sedum_pacense", "Sedum_carinatifolium", "Thompsonella_nellydiegoae", "Thompsonella_garcia-mendozae", "Sedum_chazaroi", "Aphanopetalum_clematideum", "Haloragis_masatierrana", "Myriophyllum_lapidicola", "Myriophyllum_mattogrossensis", "Haloragodendron_gibsonii", "Ribes_villosum", "Ribes_janczewskii", "Ribes_graveolens", "Ribes_andicola", "Ribes_ovalifolium", "Ribes_ceriferum", "Ribes_brandegeei", "Ribes_howellii", "Ribes_tortuosum", "Chrysosplenium_oppositifolium", "Micranthes_redofskyi", "Micranthes_clusii", "Astilbe_shikokiana", "Astilbe_philippinensis", "Saxifraga_kruhsiana", "Saxifraga_spinulosa", "Saxifraga_firma", "Saxifraga_rebunshirensis", "Saxifraga_caulescens", "Saxifraga_omolojensis", "Saxifraga_repanda", "Saxifraga_duthiei", "Saxifraga_wendelboi", "Saxifraga_felineri", "Saxifraga_scleropoda", "Saxifraga_subverticillata", "Saxifraga_sancta", "Saxifraga_cartilaginea", "Saxifraga_kolenatiana", "Saxifraga_losae", "Saxifraga_genesiana", "Saxifraga_styriaca", "Saxifraga_versicallosa", "Saxifraga_pasumensis", "Saxifraga_hengduanensis", "Medusandra_mpomiana", "Paeonia_morisii", "Paeonia_parnassica", "Paeonia_arietina", "Paeonia_mlokosewitschii", "Paeonia_veitchii", "Paeonia_wittmanniana", "Paeonia_szechuanica", "Paeonia_macrophylla", "Paeonia_caucasica", "Daphniphyllum_phillipinensis", "Molinadendron_sinaloense", "Eustigma_honbaense", "Hamamelis_mexicana", "Hamamelis_ovalis", "Matudaea_colombiana", "Corylopsis_stenopetala", "Corylopsis_himalayana", "Saxifraga_rudolphiana"]

remove_list = [x for x in tree_occ_matched_taxa if x not in occurrence_not_trait]

print(len(remove_list))
#print(list(combined_dataframe.index))

def subset_columnrow(input, subsetlist):
	sub = input[subsetlist]
	sub = sub[sub.index.isin(subsetlist)]
	if len(list(sub.index)) != len(sub.columns.values):
		print("Subset did not work.")
		sys.exit()
	return sub

combined_dataframe_subset = subset_columnrow(combined_dataframe, remove_list)
combined_dataframe_subset.to_csv(path_or_buf="_".join([distancematrix, "occ_matched.csv"]), sep=",")
