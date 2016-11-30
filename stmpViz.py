#written by Noah Friedman

#What does this script do, you may ask.
#Good question! 
#Its purpose is to take raw STMP data from the stmp TSV and properly format it for drawing by the STMP visualization program itself.
#This involves reading through the tsv, selecting values of interest, and converting them into proper "drawing vals" for the visualization.
#All that is done in here, in python.
#Python is so much easier than java, and I designed this so java does the bare minimum except drawing.e

import sys
import os
import random

#file format: each variant is listed, separated by a newline.  Within each variant, the relevant values are separated by the specififed delimiters

#file writing delimiters: (NOTE! IT IS IMPORTANT TO ENSURE THAT THESE NEVER OCCUR IN THE TSV DATA.  IF THEY DO THE SCRIPT BREAKS)
#the three sections of each line (variant info, numeric values, and string values [design as of 10.24])
lineSectionDelimiter = '\t'
#marks the limit between attributes within the same section
attributeDelimeter = '|'
#marks the limit between values within the same section
valueDelimiter = ';'

#for testing purposes 
sampleCols = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'UDN755592', 'Gene_Summary', 'Description_Summary', 'Function_Summary', 'ExonicFunction_Summary', 'Max_Allele_Freq_Summary', 'AC', 'GT', 'DP', 'AD', 'Annovar_Func_refGene', 'Annovar_Gene_refGene', 'Annovar_GeneDetail_refGene', 'Annovar_ExonicFunc_refGene', 'Annovar_AAChange_refGene', 'Annovar_Func_knownGene', 'Annovar_Gene_knownGene', 'Annovar_GeneDetail_knownGene', 'Annovar_ExonicFunc_knownGene', 'Annovar_AAChange_knownGene', 'Annovar_Func_wgEncodeGencodeBasicV19', 'Annovar_Gene_wgEncodeGencodeBasicV19', 'Annovar_GeneDetail_wgEncodeGencodeBasicV19', 'Annovar_ExonicFunc_wgEncodeGencodeBasicV19', 'Annovar_AAChange_wgEncodeGencodeBasicV19', 'snpeff_INFO', 'uk10k_AF', 'gonl_AF', 'wellderly_freq_AF', 'hg19_popfreq_all_20150413_pop_freq_max', 'hg19_popfreq_all_20150413_1000g_all', 'hg19_popfreq_all_20150413_1000g_afr', 'hg19_popfreq_all_20150413_1000g_amr', 'hg19_popfreq_all_20150413_1000g_eas', 'hg19_popfreq_all_20150413_1000g_eur', 'hg19_popfreq_all_20150413_1000g_sas', 'hg19_popfreq_all_20150413_exac_all', 'hg19_popfreq_all_20150413_exac_afr', 'hg19_popfreq_all_20150413_exac_amr', 'hg19_popfreq_all_20150413_exac_eas', 'hg19_popfreq_all_20150413_exac_fin', 'hg19_popfreq_all_20150413_exac_nfe', 'hg19_popfreq_all_20150413_exac_oth', 'hg19_popfreq_all_20150413_exac_sas', 'hg19_popfreq_all_20150413_esp6500siv2_all', 'hg19_popfreq_all_20150413_esp6500siv2_aa', 'hg19_popfreq_all_20150413_esp6500siv2_ea', 'hg19_popfreq_all_20150413_cg46', 'refseq_r_geneName', 'refseq_r_name', 'refseq_r_strand', 'refseq_r_cdsStart', 'refseq_r_cdsEnd', 'refseq_r_exonCount', 'clinvar_mut', 'clinvar_measureset_id', 'clinvar_symbol', 'clinvar_clinical_significance', 'clinvar_review_status', 'clinvar_hgvs_c', 'clinvar_hgvs_p', 'clinvar_all_submitters', 'clinvar_all_traits', 'clinvar_all_pmids', 'clinvar_pathogenic', 'clinvar_conflicted', 'clinvar_CLNSTARS', 'hg19_phastConsElements46way_r_MSA_MCE_lod', 'hg19_phastConsElements46way_r_MSA_MCE_score', 'hg19_ljb26_all_SIFT_score', 'hg19_ljb26_all_SIFT_pred', 'hg19_ljb26_all_Polyphen2_HDIV_score', 'hg19_ljb26_all_Polyphen2_HDIV_pred', 'hg19_ljb26_all_Polyphen2_HVAR_score', 'hg19_ljb26_all_Polyphen2_HVAR_pred', 'hg19_ljb26_all_LRT_score', 'hg19_ljb26_all_LRT_pred', 'hg19_ljb26_all_MutationTaster_score', 'hg19_ljb26_all_MutationTaster_pred', 'hg19_ljb26_all_MutationAssessor_score', 'hg19_ljb26_all_MutationAssessor_pred', 'hg19_ljb26_all_FATHMM_score', 'hg19_ljb26_all_FATHMM_pred', 'hg19_ljb26_all_RadialSVM_score', 'hg19_ljb26_all_RadialSVM_pred', 'hg19_ljb26_all_LR_score', 'hg19_ljb26_all_LR_pred', 'hg19_ljb26_all_VEST3_score', 'hg19_ljb26_all_CADD_raw', 'hg19_ljb26_all_CADD_phred', 'hg19_ljb26_all_GERP++_RS', 'hg19_ljb26_all_phyloP46way_placental', 'hg19_ljb26_all_phyloP100way_vertebrate', 'hg19_ljb26_all_SiPhy_29way_logOdds', 'exac03_AF', 'exac_tolerance_r_transcript', 'exac_tolerance_r_gene', 'exac_tolerance_r_n_exons', 'exac_tolerance_r_bp', 'exac_tolerance_r_syn_z', 'exac_tolerance_r_mis_z', 'exac_tolerance_r_lof_z', 'exac_tolerance_r_pLI', 'gene_info_r_synonyms', 'gene_info_r_type_of_gene', 'gene_info_r_symbol_from_nomenclature_authority', 'gene_info_r_full_name_from_nomenclature_authority', 'gene_info_r_other_designations', 'hg19_wgEncodeBroadHmmNhlfHMM_r_info', 'hg19_wgEncodeBroadHmmHmecHMM_r_info', 'hg19_wgEncodeBroadHmmH1hescHMM_r_info', 'hg19_wgEncodeBroadHmmHuvecHMM_r_info', 'hg19_wgEncodeBroadHmmNhekHMM_r_info', 'hg19_wgEncodeBroadHmmHsmmHMM_r_info', 'hg19_wgEncodeBroadHmmGm12878HMM_r_info', 'hg19_wgEncodeRegDnaseClusteredV3_r_info', 'hg19_wgEncodeRegDnaseClusteredV3_r_disease', 'hg19_geneReviews_r_bin', 'hg19_geneReviews_r_name']
variant = ['2', '83351', '', 'T', 'G', '1543.98', 'PASS', 'AC=1;AF=0.5;AN=2;BaseQRankSum=-0.106;ClippingRankSum=0.622;DP=44;FS=1.298;GQ_MEAN=173.27;GQ_STDDEV=155.57;InbreedingCoeff=-0.1111;MQ=60;MQ0=0;MQRankSum=-0.07;NCC=0;QD=10.43;ReadPosRankSum=-0.599;SOR=0.832;VQSLOD=16.59;culprit=ReadPosRankSum', 'GT:AD:DP:GQ:PL', '0/1:22,22:44:99:522,0,518', 'FAM110C,SH3YL1', '', 'intergenic', '', '', '1', '0/1', '44', '22,22', 'intergenic', 'FAM110C,SH3YL1', 'dist=36763;dist=134785', '', '', 'intergenic', 'FAM110C,SH3YL1', 'dist=36763;dist=134785', '', '', 'intergenic', 'FAM110C,AC079779.7', 'dist=36966;dist=114218', '', '', 'AC=1;AF=0.5;AN=2;BaseQRankSum=-0.106;ClippingRankSum=0.622;DP=44;FS=1.298;GQ_MEAN=173.27;GQ_STDDEV=155.57;InbreedingCoeff=-0.1111;MQ=60;MQ0=0;MQRankSum=-0.07;NCC=0;QD=10.43;ReadPosRankSum=-0.599;SOR=0.832;VQSLOD=16.59;culprit=ReadPosRankSum;ANN=G|intergenic_region|MODIFIER|FAM110C-SH3YL1|FAM110C-SH3YL1|intergenic_region|FAM110C-SH3YL1|||n.83351T>G||||||', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '', '\n']
#a dictionary of functions that maps the name of a value to a function that scales it for drawing
drawingFunctions = {
	'QUAL': lambda x: qual_score_to_drawing_val(x),
	'Max_Allele_Freq_Summary': lambda x: qual_score_to_drawing_val(x),
	'hg19_phastConsElements46way_r_MSA_MCE_lod': lambda x: const_score_to_drawing_val(x),
	'hg19_ljb26_all_CADD_raw': lambda x: cadd_score_to_drawing_val(x),
	'AD': lambda x: allele_depth_to_drawing_val(x),
	'hg19_ljb26_all_Polyphen2_HDIV_score': lambda x: polyphen2_score_to_drawing_val(x),
	'exac_tolerance_r_lof_z': lambda x: exac_tolerance_to_drawing_val(x),
	'DP': lambda x: dp_to_drawing_val(x)
}

nameFunctions = {
	'clinvar_clinical_significance': lambda x: clinvar_clinical_significance_to_val(x),
	'Gene_Summary': lambda x: gene_summary_to_val(x),
	'Function_Summary': lambda x: function_summary_to_val(x),
	'ExonicFunction_Summary': lambda x: exonic_function_summary_to_val(x)
}

#functions for translating scores to drawing values

#todo--add support for na
def dp_to_drawing_val(dp):
	return random.randrange(0, 10)

def exac_tolerance_to_drawing_val(x):
	return random.randrange(0, 10)

def polyphen2_score_to_drawing_val(x):
	return random.randrange(0, 10)

def allele_depth_to_drawing_val(ad):
	return random.randrange(0, 10)

def qual_score_to_drawing_val(qual):
	return random.randrange(0, 10)
	#return qual

def max_allele_freq_to_drawing_val(mxaf):
	return random.randrange(0, 10)
	#return mxaf

def const_score_to_drawing_val(constScore):
	return random.randrange(0, 10)
	#return constScore

def cadd_score_to_drawing_val(caddScore):
	return random.randrange(0, 10)
	#return caddScore

#----------------------------------------------------------------
#NAME type functions
#these return indicies to the colors defined in the config file
def clinvar_clinical_significance_to_val(clinvarVal):
	if clinvarVal == 'Likely pathogenic': return 1
	if clinvarVal == 'Benign': return 3
	if clinvarVal == 'Uncertain significance': return 4
	if clinvarVal == 'Pathogenic': return 0
	if clinvarVal == 'Likely Benign': return 2

def gene_summary_to_val(geneSummary):
	return random.choice([0, 1, 2, 3])

def function_summary_to_val(functionSummary):
	if functionSummary == "exonic": return 1
	if functionSummary == "intronic": return 4
	if functionSummary == "UTR3": return 3
	if functionSummary == "splicing": return 2
	if functionSummary == "ncRNA_exonic": return 0

def exonic_function_summary_to_val(exonicFunctionSummary):
	if exonicFunctionSummary == "nonsynonymous": return 0
	if exonicFunctionSummary == "synonymous": return 1
	if exonicFunctionSummary == "frameshift deletion": return 2
	if exonicFunctionSummary == "stopgain": return 3
	if exonicFunctionSummary == "frameshift insertion": return 4
	if exonicFunctionSummary == "stoploss": return 5
	
#----------------------------------------------------------------

#figures out the columns we should use as input for the file bases on what is specified in the config file
#the sister java program reads from the exact same config file
def read_config_columns(configFile):
	with open(configFile, 'r') as f:
		#skip the header line
		lines = f.readlines()
		#the first line of the config should be info, the second drawing cols, the third name cols
		#the string manipulations strip out the new lines and then split the line by commas
		#we skip the header, ergo we are 1 indexed
		print lines[1]
		print lines[2]
		infoCols = [lines[1].strip('\n').split(',')[i] for i in range(len(lines[1].strip('\n').split(',')))]
		drawingCols = [lines[2].strip('\n').split(',')[i] for i in range(len(lines[2].strip('\n').split(',')))]
		nameCols = [lines[3].strip('\n').split(',')[i] for i in range(len(lines[3].strip('\n').split(',')))]
	
	f.close()
	return infoCols, drawingCols, nameCols

def create_idx_dict(columns):
	idx_dict = {}
	cntr = 0
	for c in columns:
		idx_dict[c] = cntr
		cntr += 1
	return idx_dict


#reads out necessary fields from the annotation tsv file
def read_tsv(tsv):
	data = []
	with open(tsv) as f:
		lines = f.readlines()
		columns = lines[0].split()
		#for line in lines[1:]:
		for line in lines[1:100]:
			data.append(line.split('\t'))
	return columns, data

#iterates through the values in a 

#FREQUENCY, QUALITY, CLINVAR, CONSERVATION
def get_drawing_vals(variantLine, idx_dict, drawingCols, valToWrite):
	for col in drawingCols:
		val = variantLine[idx_dict[col]]
		#ensure that we don't append the empty string, this breaks interpretation
		#if val == '': val = 'na'
		#Values are written in this order: Name, real value, drawing value
		drawingVal = str(drawingFunctions[col](val))
		valToWrite += col
		valToWrite += valueDelimiter
		valToWrite += val
		valToWrite += valueDelimiter
		valToWrite += drawingVal
		valToWrite += attributeDelimeter

	return valToWrite

def get_name_vals(variantLine, idx_dict, nameCols, valToWrite):
	for col in nameCols:
		val = variantLine[idx_dict[col]]
		#ensure that we don't append the empty string, this breaks interpretation
		if val == '': val = 'na'

		#TEMPORARY HACK FOR DEMONSTRATION PURPOSES
		if col == 'clinvar_clinical_significance':
			r = random.uniform(0,1)
			#make pathogenic more rare
			if r < .01: val = 'Pathogenic'
			elif r < .1: val = 'Likely pathogenic'
			elif r < .5: val = 'Likely Benign'
			else: val = 'Benign' 
			#val = random.choice(['Likely pathogenic', 'Likely benign', 'Benign', 'Uncertain significance', 'Pathogenic'])
		if col == 'Function_Summary': val = random.choice(["exonic", "intronic", "UTR3", "splicing", "ncRNA_exonic"])
		if col == 'ExonicFunction_Summary': val = random.choice(["nonsynonymous", "synonymous", "frameshift deletion", "stopgain", "frameshift insertion", "stoploss"])

		nameVal = str(nameFunctions[col](val))
		valToWrite += col
		valToWrite += valueDelimiter
		valToWrite += val
		valToWrite += valueDelimiter
		valToWrite += nameVal
		valToWrite += attributeDelimeter

	return valToWrite

#gets basic variant info (i.e ref/alt etc) and writes it
def get_variant_info(variantLine, idx_dict, valToWrite):
	#CHANGE the structure of this
	for col in infoCols:
		val = variantLine[idx_dict[col]]
		#ensure that we don't append the empty string, this breaks interpretation
		if val == '': val = 'na'
		if col == 'Gene_Summary': val = random.choice(['OR2T35', "BRCA1", "AFF3", "MYO7B", "ZNF806", "NEB", "SP100", "SYN2"])

		valToWrite += val
		valToWrite += valueDelimiter

	return valToWrite

#obselete???
def prepare_header(mode, valsToWrite):
	if mode == 'single':
		valsToWrite.append(['#single', 'explanation: single mode--this data is intended to be only a single representation of data'])
	elif mode == "full":
		valsToWrite.append(['#multiple', 'explanation: mutliple mode--this data is intended to be a representation of the entire data returned by stmp'])
	else:
		print "error: improper mode specified"
		sys.exit(0) 

#sorts the data to be written by a specific value
#three column indicies
#line section columns: which tab separated column do we want?  ie do we want chrom ref info or numeric info or others?
#attribute columns: which attribute value do we want
#value columns: which value within an attribute do we want?
#you must pass 0 for a column if it doesnt exist/isnt relevant
#for example, a call sort_data_by_value(linesToWrite, 1, 0, 1) would sort by the 1st line (numeric variant attributes), the 0th attribute (quality) and the first value (raw qual score)
def sort_data_by_value(dataLines, lineSectionCol, attributeSectionCol, valueSelectionCol, sortMode):
	#the lambda function to sort the data fully splits it and extracts the value of interest
	#sort by the float value at this position

	#TODO: add error catching for NAN/empty values
	if sortMode == 'float':
		dataLines.sort(key = lambda x: 
			float(x.split(lineSectionDelimiter)[lineSectionCol]
			.split(attributeDelimeter)[attributeSectionCol]
			.split(valueDelimiter)[valueSelectionCol]
			))
	#sort by the alpha numeric value at this position
	else:
		dataLines.sort(key = lambda x: 
			x.split(lineSectionDelimiter)[lineSectionCol]
			.split(attributeDelimeter)[attributeSectionCol]
			.split(valueDelimiter)[valueSelectionCol]
			)
	print dataLines


#writes output to a file that can then be read by the graphical interface
#each variant gets its own file
def write_file(columns, linesToWrite, pos):
	savePath = "/home/noahfrie/noahfrie/devCode/stmpViz/outputFiles"
	fullName = os.path.join(savePath, pos +'viz.txt')
	f = open(fullName, 'w')
	for line in linesToWrite:
		print line
		f.write(line)
		f.write('\n')
	f.close



#--------------------MAIN CODE-------------------------------

tsv = sys.argv[1]
#columns = sampleCols
mode = 'multiple'

#columns: the names for the values that ought to be extracted from the STMP data--these columns are set by the read_config_cols function
#they are global variables 
infoCols, drawingCols, nameCols = read_config_columns(os.getcwd() + '/config.txt')
#infoCols = ['QUAL','Max_Allele_Freq_Summary','hg19_phastConsElements46way_r_MSA_MCE_lod','hg19_ljb26_all_CADD_raw','AD','hg19_ljb26_all_Polyphen2_HDIV_score','exac_tolerance_r_lof_z','DP']

columns, data = read_tsv(tsv)
valsToWrite = []

#prepare_header(mode, valsToWrite)
idx_dict = create_idx_dict(columns)

linesToWrite = []
#NOTE THIS CODE IS OBSELETE AND OUGHT TO BE DELETED
if mode == 'single':
	#n = 0
	#variant = data[n]
	get_variant_info(variant, idx_dict, valsToWrite)
	valsToWrite.append('#')
	get_drawing_vals(variant, idx_dict, drawingCols, valsToWrite)
	#append a section delimiter
	valsToWrite.append('#')
	get_name_vals(variant, idx_dict, nameCols, valsToWrite)

elif mode == 'multiple':
	cntr = 0
	for line in data:
		valToWrite = ''
		variant = line 
		valToWrite = get_variant_info(variant, idx_dict, valToWrite)
		valToWrite += lineSectionDelimiter
		valToWrite = get_drawing_vals(variant, idx_dict, drawingCols, valToWrite)
		#append a section delimiter
		valToWrite += lineSectionDelimiter
		valToWrite = get_name_vals(variant, idx_dict, nameCols, valToWrite)
		linesToWrite.append(valToWrite)

	sort_data_by_value(linesToWrite, 1, 0, 1, 'float')


write_file(columns, linesToWrite, "2")









