import numpy as np
from collections import Counter
import networkx as nx
import os
import re
import random
import itertools
from scipy import stats as stats
import sys

#########################################################################################################################
# Author: Hiren Karathia																								#
# Last time Modified: 02/08/2014-3:27PM Sunday																			#
# This script is designed and checked for Calculating Edge_Fraction values of Gene Pathways.							#
# The Edge-Fraction is calculated for the gene sets that are overlapped with both the end points of HiC interactions	#
# and the same gene sets are mapped with Pathways involved genes.														#
# Initially a file is needed as a global data set of mapped Genes and Transcripts IDs.									#
# These Genes are filtered based on how many of them are mapped in overall HiC interacting Gene sets.					#
# The filtered Genes are further reduced in number based on how many of them are mapped on each of the chromosomes.		#
# These chromosomal co-localized Genes are ultimately checked in HiC interacting Networks and find out how many   		#
# Gene-Gene interactions are there in between the Genes of the corresponding pathway.  									#
# Finally Edge-Fraction is calculated as the Number of interactions amongst genes in a pathway (Edges) divided by 		#
# number of all possible interactions amongst the Genes of the same pathway (Nodes) that had been mapped in 			#	
# Hi-C network.																											#  
# The same number of randomly chosen other Genes are selected from the same chromosome where these Nodes are localized 	#
# and perform the same Edge-fraction analysis to check probability of the pathways involved Genes Edge-Fraction.   		#
#########################################################################################################################


def getData(fl):
	data = np.genfromtxt(fl, dtype='S', delimiter='\t', usecols=(0,1,2,3,4,7,8,9,10,13,14,15,16,17,20,21,22,24,25,26,29,31,32,33))
	#data = np.genfromtxt(fl, dtype='S', delimiter='\t', comments="#", usecols=(0,1,2,3,4,7,8,9,10,13,14,15,16,17,20,21,22,24,25,26,29,31,32,33))
	data = data[1:]
	return data
	
def getFDRData(fl):
	data = np.genfromtxt(fl, dtype='S', delimiter='\t', usecols=(0,1,2,3,4,7,8,9,10,13,14,15,16,17,20,21,22,24,25,26,29,31,32,33,18))
	#data = np.genfromtxt(fl, dtype='S', delimiter='\t', comments="#", usecols=(0,1,2,3,4,7,8,9,10,13,14,15,16,17,20,21,22,24,25,26,29,31,32,33))
	INT_NO = float(data[0][24].split("based on")[1].split("total")[0].strip())
	data = data[1:]
	return data, INT_NO
	
def removeOthers(classes, tots):
	idx = np.where(classes == "Intergenic")
	tots[idx] = False
	idx = np.where(classes == "non-coding")
	tots[idx] = False
	#idx = np.where(classes == "promoter-TSS")
	#tots[idx] = False
	idx = np.where(classes == "NA")
	tots[idx] = False
	return tots
	
def getPromoter(classes, tots):
	idx = np.where(classes == "Intergenic")
	tots[idx] = False
	idx = np.where(classes == "non-coding")
	tots[idx] = False
	idx = np.where(classes == "promoter-TSS")
	tots[idx] = False
	return tots	

def getGene_Genome():
	data = np.loadtxt("/cbcbhomes/hiren/software/homer/data/genomes/hg19/hg19.annotation", dtype='S', delimiter='\t')
	terms = data[:,0]
	t = np.core.defchararray.strip(np.core.defchararray.partition(terms, "(")[:,0])
	genes = np.core.defchararray.strip(np.core.defchararray.partition(terms, "(")[:,2])
	genes = np.core.defchararray.replace(genes, ")", "")
	genes = np.core.defchararray.partition(genes, ",")[:,0]
	genes = np.core.defchararray.strip(genes)
	U_Genes = np.unique(genes[np.where(np.logical_or(t == "exon", t == "intron"))])
	return U_Genes

def getAverageDistance(K):
		KK = nx.Graph(K.edges(data=True))
		#D = np.array([K.edges(data=True)[i][2]['weight'] for i in range(len(K.edges()))], dtype='i')
		D = np.array([KK.edges(data=True)[i][2]['weight'] for i in range(len(KK.edges()))], dtype='i')
		MD = np.mean(D)
		return MD
		
def getDistZscore(R_EFs, EF):
	R_EFs = np.array(R_EFs, dtype='f')
	M = np.mean(R_EFs)
	V = np.std(R_EFs)
	if V == 0.0:
		Z = 20.0
	else:
		Z = -float(EF - M)/V
	if Z > 20.0:
		Z = 20.0
	return Z
	
def getDijkstra_path_length(a, b, G):
	try:
		l = nx.dijkstra_path_length(G, a, b)
	except:
		l = 0
	return l

def randomizeLength():
	pass
	
def randomizeDistance():
	pass
	
def randomizeBoth():
	pass	

def getLenChoice(g, STOT_C, STOT_G, STOT_L):
	
	l = GEN_LEN[g]
	
	u = l + ((l*20)/100)
	l = l - ((l*20)/100)

	#LG = TOT_G[np.where(np.logical_and(TOT_L >= l, TOT_G <= u))]	
	#LG = S_G
	#LG = S_G[np.where(np.logical_and(S_L >= l, S_L <= u))]
	
	LG = np.unique(STOT_G[np.where(np.logical_and(STOT_L >= l, STOT_L <= u))])
	LG = np.delete(LG, np.where(LG == g))
	
	if LG.size > 0:
		return LG
	else:
		LG = np.delete(STOT_G, np.where(STOT_G == g))
		return np.unique(LG)
	LG = np.delete(STOT_G, np.where(STOT_G == g))
	return np.unique(LG)
	
def getRNAgeneChoice(g, STOT_C, STOT_G, STOT_L):
	
	#LG = TOT_G[np.where(np.logical_and(TOT_L >= l, TOT_G <= u))]	
	#LG = S_G
	#LG = S_G[np.where(np.logical_and(S_L >= l, S_L <= u))]
	
	LG = np.unique(STOT_G)
	LG = np.delete(LG, np.where(LG == g))
	
	return np.unique(LG)	

def getLenAndChrChoice(g, c, STOT_C, STOT_G, STOT_L):
	
	l = GEN_LEN[g]
	
	if re.search("_", c):
		c = c.split("_")[0]
		
	u = l + ((l*20)/100)
	l = l - ((l*20)/100)
	
	
	S_G = STOT_G[np.where(STOT_C == c)]
	S_L = STOT_L[np.where(STOT_C == c)]
	#LG = TOT_G[np.where(np.logical_and(TOT_L >= l, TOT_G <= u))]	
	#LG = S_G
	LG = S_G[np.where(np.logical_and(S_L >= l, S_L <= u))]
	LG = np.delete(LG, np.where(LG == g))
	
	if LG.size > 0:
		return LG
		
	elif LG.size <= 0:
		LG = S_G[np.where(np.logical_and(S_L >= l, S_L <= u))]
		LG = np.delete(LG, np.where(LG == g))
				
		if LG.size > 0:
			return LG 
		else:
			if c == "chrY":
				LG = np.delete(S_G, np.where(S_G == g))
				if LG.size > 0:
					return LG
			else:
				if S_G.size > 0:
					LG = np.delete(S_G, np.where(S_G == g))
					return LG
	else:
		S_G = STOT_G[np.where(np.logical_and(STOT_G >= l, STOT_G <= u))]
		LG = np.delete(S_G, np.where(S_G == g))
	return S_G
	
def getChrAndLenChoice(g, STOT_C, STOT_G, STOT_L):
	
	l = GEN_LEN[g]
	
	u = l + ((l*20)/100)
	l = l - ((l*20)/100)

	#LG = TOT_G[np.where(np.logical_and(TOT_L >= l, TOT_G <= u))]	
	#LG = S_G
	#LG = S_G[np.where(np.logical_and(S_L >= l, S_L <= u))]
	
	LG = np.unique(STOT_G[np.where(np.logical_and(STOT_L >= l, STOT_L <= u))])
	LG = np.delete(LG, np.where(LG == g))
	
	if LG.size > 0:
		return LG
	else:
		LG = np.delete(STOT_G, np.where(STOT_G == g))
		return np.unique(LG)
	LG = np.delete(STOT_G, np.where(STOT_G == g))
	return np.unique(LG)	
	
def getGeneChoice(g):
	l = GEN_LEN[np.where(TOTAL_GENES == g)][0]
	u = 1000000
	l = 0
	LG = TOTAL_GENES[np.where(np.logical_and(GEN_LEN >= l, GEN_LEN <= u))]
	LG = np.delete(LG, np.where(LG == g))
	return LG
	
def getRandomInteraction(F_Nodes, Int_Chr, GD, STOT_C, STOT_L, STOT_G, siter, RNA_ANALYSIS):
	RG = []
	REF = []
	F_GENES = np.array(F_Nodes, dtype='S') 
	idx = np.where(np.logical_not(np.in1d(STOT_G, F_GENES)))
	G_STOT_G = STOT_G[idx]
	G_STOT_C = STOT_C[idx]
	G_STOT_L = STOT_L[idx]
	
	for k in range(siter):
		Nodes = []
		Edge_no = []
		Node_no = []
		Edges = []
		if chr_control == 2:
			#r_genes = np.array([random.choice(getLenChoice(F_Nodes[i], STOT_C, STOT_G, STOT_L)) for i in range(len(F_Nodes))])
			if RNA_ANALYSIS:
				#r_genes = np.array([random.choice(NC_GN) for i in range(len(F_Nodes))])
				r_genes = np.array([random.choice(getLenChoice(F_Nodes[i], G_STOT_C, G_STOT_G, G_STOT_L)) for i in range(len(F_Nodes))])
			else:
				r_genes = np.array([random.choice(getLenChoice(F_Nodes[i], G_STOT_C, G_STOT_G, G_STOT_L)) for i in range(len(F_Nodes))])
		elif chr_control == 1:
			#r_genes = np.array([random.choice(getLenAndChrChoice(F_Nodes[i], Int_Chr[i], STOT_C, STOT_G, STOT_L)) for i in range(len(F_Nodes))])
			if RNA_ANALYSIS:
				#r_genes = np.array([random.choice(NC_GN) for i in range(len(F_Nodes))])
				r_genes = np.array([random.choice(getLenAndChrChoice(F_Nodes[i], Int_Chr[i], STOT_C, STOT_G, STOT_L)) for i in range(len(F_Nodes))])
			else:
				r_genes = np.array([random.choice(getLenAndChrChoice(F_Nodes[i], Int_Chr[i], STOT_C, STOT_G, STOT_L)) for i in range(len(F_Nodes))])	
		
		(E, N, ENo, NNo, K) = getSubGraph(r_genes, GD)
		Edges.append(E)
		Nodes.append(N)
		Edge_no.append(ENo)
		Node_no.append(NNo)
		Nodes = np.array(Nodes)
		Edges = np.array(Edges)
		Edge_no = np.array(Edge_no, dtype='i')
		Node_no = np.array(Node_no, dtype='i')	
				
		if np.where(Node_no > 0)[0].size > 0:
			if RNA_ANALYSIS:
				(EF, Nodes, Int_Chr) = getEdgeExpression(Node_no, Edge_no, Nodes, Edges, Int_Chr)
			else:	
				(EF, Nodes, Int_Chr) = getEdgeFractionValue(Node_no, Edge_no, Nodes, Int_Chr)
		else:
			EF = 0.0
			
		REF.append(EF)
		
	REF = np.nan_to_num(np.array(REF))
	
	return REF, Nodes, Int_Chr

def getChrRandomInteraction(F_Nodes, Int_Chr, GD, STOT_C, STOT_L, STOT_G):
	RG = []
	REF = []
	
	F_GENES = np.array(F_Nodes, dtype='S') 
	idx = np.where(np.logical_not(np.in1d(STOT_G, F_GENES)))
	G_STOT_G = STOT_G[idx]
	G_STOT_C = STOT_C[idx]
	G_STOT_L = STOT_L[idx]	
	
	for k in range(iter):
		
		Nodes = []
		Edge_no = []
		Node_no = []
				
		#r_genes = np.array([random.choice(getLenAndChrChoice(F_Nodes, Int_Chr[j], STOT_C, STOT_G, STOT_L)) for i in range(len(F_Nodes))])
		r_genes = np.array([random.choice(getLenAndChrChoice(F_Nodes, Int_Chr[j], G_STOT_C, G_STOT_G, G_STOT_L)) for i in range(len(F_Nodes))])
			
		(E, N, ENo, NNo, K) = getSubGraph(r_genes, GD)
		Nodes.append(N)
		Edge_no.append(ENo)
		Node_no.append(NNo)
		Nodes = np.array(Nodes)
		Edge_no = np.array(Edge_no, dtype='i')
		Node_no = np.array(Node_no, dtype='i')	
				
		if np.where(Node_no > 0)[0].size > 0:
			(EF, Nodes, Int_Chr) = getEdgeFractionValue(Node_no, Edge_no, Nodes, Int_Chr)
		else:
			EF = 0.0
			
		REF.append(EF)
	REF = np.nan_to_num(np.array(REF))	
	return REF
	
def getSubGraph(genes, GD):

	K = GD.subgraph(genes)
	Edges = np.array(K.edges()).tolist()
	Nodes = K.nodes()
	Nodes_no = len(Nodes)
	Edge_no = K.size()
	
	return Edges, Nodes, Edge_no, Nodes_no, K

def getRPKMValuesOfGenes(SET):
	IDX = np.array([np.where(SET[i] == RGENES)[0].size > 0 for i in range(SET.size)])
	NSET = SET[np.where(IDX)]
	R_SET = np.array([np.max(RFPKM[np.where(RGENES == NSET[i])]) for i in range(NSET.size)])
	return R_SET	
	
def getEdgeExpression(Nodes_no, Edge_no, Nodes, Edges, Int_Chr):
	Edge_no = np.array(Edge_no, dtype='f')
	idx = np.where(Nodes_no > 0)
	Nodes = Nodes[idx][0]
	Edge_Genes = np.unique(np.ndarray.flatten(Edges))
	Node_Genes = np.unique(np.ndarray.flatten(Nodes))
	
	if rna_func_choice == 1:
		E_RPKM = getRPKMValuesOfGenes(Edge_Genes)
		N_RPKM = getRPKMValuesOfGenes(Node_Genes)
		if E_RPKM.size > 0:
			EF = float(np.sum(E_RPKM))/np.sum(N_RPKM)
			#EF = np.log2(np.mean(RPKM))
		else:
			EF = 0.0
	elif rna_func_choice == 2:
		E_RPKM = getRPKMValuesOfGenes(Edge_Genes)
		if E_RPKM.size > 0:
			EF = float(np.sum(E_RPKM))
			if EF != 0.0:
				EF = np.log2(EF)
				#EF = EF
			else:
				EF = 0.0
		else:
			EF = 0.0
	elif rna_func_choice == 3:
		N_RPKM = getRPKMValuesOfGenes(Node_Genes)
		if N_RPKM.size > 0:
			EF = float(np.sum(N_RPKM))
			EF = np.log2(EF)
		else:
			EF = 0.0
			
	else:
		print "Your Selection of the Function do not available for the analysis !! "
		print "Please choose appropriate options !! Thank you !! "
		sys.exit()
	
		
	return EF, Nodes, Int_Chr	
	
def getEdgeFractionValue(Nodes_no, Edge_no, Nodes, Int_Chr):
	Edge_no = np.array(Edge_no, dtype='f')
	#idx = np.where(Edge_no > 0.0)
	
	idx = np.where(Nodes_no > 0)
	
	Nodes = Nodes[idx][0]
	
	#Int_Chr = Int_Chr[idx]
	
	HE = np.array(Edge_no[0], dtype='f')
	HN = np.array(Nodes_no[0], dtype='f')
	MHN = HN - 1
	
	BHN = HN*MHN/2
	BHN = np.array(BHN, dtype='f')
	
	FG = np.sum(HE)
	BG = np.sum(BHN)
	
	if BG > 0:
		EF = FG/BG
	else:	
		EF = 0.0
			
	return EF, Nodes, Int_Chr		

	
def getEdgeDecision(pg_gid, p_chr, GD):
	
	u_chr = np.unique(p_chr)
	Edges = []
	Nodes = []
	Edge_no = []
	Node_no = []
	Int_Chr = []
	
	
	flag = False
	
	#for u in u_chr:
	#print u, "\t", t_genes
			
	(E, N, ENo, NNo, K) = getSubGraph(pg_gid, GD)
	Edges.append(E)
	Nodes.append(N)
	Edge_no.append(ENo)
	Node_no.append(NNo)
	
	Edges = np.array(Edges, dtype='S')
	Nodes = np.array(Nodes, dtype='S')
	idx = np.in1d(pg_gid, Nodes[0])
	Edge_no = np.array(Edge_no, dtype='i')
	Node_no = np.array(Node_no, dtype='i')
	Int_Chr = np.array(p_chr)
	Int_Chr = Int_Chr[idx]
	
	Edges_Node_No = np.ndarray.flatten(Edges).size
	Edge_Genes = np.unique(np.ndarray.flatten(Edges))
	
	#print Edge_no
	#print Node_no
	
	if np.where(Edge_no > 0)[0].size > 0:
		flag = True
	
	return Edges, Nodes, Edge_no, Node_no, Edges_Node_No, Edge_Genes, Int_Chr, flag

	
def getPathwaysInfo_With_Unique(fl):
	pflag = False
	p_chr = []
	p_genes = []
	pg_mid = []
	pg_gid = []
	
	dirs = "/fs/ftp-umiacs/pub/hiren/Hiren/HiC_Samples/Annotations/Gene_Kegg_Pathways/"
	fl = dirs + fl
	data = np.loadtxt(fl, dtype='S', delimiter='\t')
	data = data[np.unique(data[:,1], return_index=True)[1]]
	
	if data[:,0].size < 5:
		return pflag, p_chr, p_genes, pg_mid, pg_gid
		
	p_chr, pg_gid, p_genes, pg_s, pg_e = (data[:,0], data[:,1], data[:,2], data[:,4], data[:,5])
	pg_s = np.array(pg_s, dtype='i')
	pg_e = np.array(pg_e, dtype='i')

	pg_mid = np.array((pg_s + np.abs(pg_e - pg_s)/2), dtype='i')
	pflag = True
	return pflag, p_chr, p_genes, pg_mid

def getPathwaysInfo_Without_Unique(fl, dirs):
	p_chr = []
	p_genes = []
	pg_mid = []
	pg_gid = []
	
	#dirs = "/fs/ftp-umiacs/pub/hiren/Hiren/HiC_Samples/Annotations/Gene_Kegg_Pathways/"
	fl = dirs + fl
	data = np.loadtxt(fl, dtype='S', delimiter='\t')
	idx = np.unique(data[:,1], return_index=True)[1]
	
	p_chr, pg_gid, p_genes, pg_s, pg_e = (data[idx][:,0], data[idx][:,1], data[idx][:,2], data[idx][:,4], data[idx][:,5])
		
	pg_s = np.array(pg_s, dtype='i')
	pg_e = np.array(pg_e, dtype='i')
	pg_mid = np.array((pg_s + np.abs(pg_e - pg_s)/2), dtype='i')
	pg_gid = np.array(pg_gid, dtype='S')
	pflag = True
	return pflag, p_chr, p_genes, pg_mid, pg_gid

def getPvalue(EF, R_EFs, siter):
	R_EFs = np.array(R_EFs, dtype='f')
	EF = float(EF)
	
	a = R_EFs[np.where(R_EFs > EF)]
	
	pval = float(a.size)/siter
	return pval
	
def createArray_Horizontal(ARR, vals):
	
	if ARR.shape[0] == 0:
		ARR = vals
	else:
		ARR = np.hstack((ARR, vals))
	return ARR	

def chooseRandomInteractions(genes):
	#R_I = np.array([[random.choice(genes), random.choice(genes)] for i in range(1000)])
	R_I = np.array([[random.choice(genes), random.choice(genes)] for i in range(10000)])
	return R_I
	
def createArray_Vertical(ARR, vals):
	if ARR.shape[0] == 0:
		ARR = vals
	else:
		ARR = np.vstack((ARR, vals))
	return ARR	

def getZscore(R_EFs, EF):	
	
	M = np.mean(R_EFs)
	V = np.std(R_EFs)
	if V == 0.0:
		Z = 20.0
	else:
		Z = float(EF - M)/V
	if Z > 20.0:
		Z = 20.0
	return Z

def getRandomNonIntreaction(TOTAL, num):
	A = [i for i in range(TOTAL.number_of_edges())]
	if len(A) > 0:
		idx = [random.choice(A) for i in range(num)]
		NTOTAL = np.array(TOTAL.edges())[idx]
		NG = nx.Graph()
		NG.add_edges_from(NTOTAL)
	else:
		NG = TOTAL.copy()
	return NG
	
def subsetNonInteraction(E, N):
	GE = nx.Graph()
	GE.add_edges_from(E[0])
	GN = nx.Graph()
	NN = list(itertools.combinations(N[0], 2))
	GN.add_edges_from(NN)
	GN.remove_edges_from(GE.edges())
	return GN
	
def subsetNonInteractionFraction(E, N):
	GE = nx.Graph()
	GE.add_edges_from(E[0])
	GN = nx.Graph()
	NN = list(itertools.combinations(N[0], 2))
	GN.add_edges_from(NN)
	GN.remove_edges_from(GE.edges())
	if GN.number_of_edges > 10:
		GN = getRandomNonIntreaction(GN, 10)
	return GN

	
def prepareExperimentReport(OUT_PATH_GENES, PATH_EDGES, NON_PATH_EDGES):
	TGD = GD.copy()
	TOT_P = TGD.subgraph(np.unique(OUT_PATH_GENES))
	
	
	IP = nx.Graph()
	IP.add_edges_from(PATH_EDGES)
	IP_G = np.unique(np.array(IP.nodes()))
	
	NIP = nx.Graph()
	if path_choice == 6:	
		NIP.add_edges_from(NON_PATH_EDGES)
		NIP.remove_edges_from(NIP.selfloop_edges())
		#NIP = getRandomNonIntreaction(NIP, 1000)
	
	NIP_G = np.unique(np.array(NIP.nodes()))
	
	TOT_P.remove_edges_from(IP.edges())
	OP = TOT_P.copy()
	OP_G = np.setdiff1d(np.unique(np.array(OP.nodes())), np.unique(np.array(IP.nodes())))
	
	TGD.remove_edges_from(IP.edges())
	TGD.remove_edges_from(TOT_P.edges())
	
	IG = TGD.copy()
	IG_G = np.setdiff1d(np.array(TGD.nodes()), np.unique(OUT_PATH_GENES))
	
	
	TG_G = np.setdiff1d(GENOME_GENE, TOTAL_GENES)
	R_I = chooseRandomInteractions(TG_G)
	TG = nx.Graph()
	TG.add_edges_from(R_I)
	
	print "\t ===================================================================================="	
	print "\t             TOTAL NUMBER OF HiC INTERACTIONS CLASSIFIED DISJOINTED SETS "
	print "\t ===================================================================================="	
	
	print "\t\t As mentioned above the UNIQUE GENES IN GENOME             :: ", GENOME_GENE.size
	print "\t\t As mentioned above the HiC GENES IN THE TISSUE            :: ", TOTAL_GENES.size
	print "\t\t As mentioned above the HiC INTERACTIONS                   :: ", len(GD.edges())
	
	print "\t\t TOTAL SUBSET OF INTERACTIONS ONLY FOUND IN TOTAL PATHWAYS :: ", len(TOT_P.edges())
	print "\t\t THE UNIQUE GENES THAT INTERACT IN ALL THE PATHWAYS        :: ", np.unique(OUT_PATH_GENES).size
	
	print "\t\t INRA-PATHWAYS INTERACTIONS                                :: ", len(IP.edges())
	print "\t\t UNIQUE SET OF INTRA-PATHWAYS INTERACTING GENES            :: ", IP_G.size
	
	if path_choice == 6:
		print "\t\t ALL POSSIBLE INTRA PATHWAYS NON-INTERACTIONS              :: ", len(NIP.edges())
		print "\t\t UNIQUE SET OF INTRA-PATHWAYS NON INTERACTING GENES        :: ", NIP_G.size
	
	print "\t\t INTER-PATHWAYS INTERACTIONS                               :: ", len(OP.edges())
	print "\t\t UNIQUE SET OF INTER-PATHWAYS INTERACTING GENES            :: ", OP_G.size
	
	print "\t\t IN-GENOME INTERACTIONS                                    :: ", len(IG.edges())
	print "\t\t UNIQUE SET OF IN-GENOME INTERACTING GENES                 :: ", IG_G.size
	
	print "\t ===================================================================================="
	
	return IP, NIP, OP, IG, TG, IP_G, NIP_G, OP_G, IG_G, TG_G
	
	
def getPathwaysGenes(GD, OFL, STOT_C, STOT_L, STOT_G, Finished_Pathway, path_dir, RNA_ANALYSIS, t_name, siter):
		
	"""
	This function is the one that perform complete Pathways Level analysis. 
	"""
	
	print "\t 7. Actual Pathways Analysis is STARTING !!! "
	
	#dirs = "/fs/ftp-umiacs/pub/hiren/Hiren/HiC_Samples/Annotations/Gene_Kegg_Pathways/"
	#dirs = "/fs/ftp-umiacs/pub/hiren/Hiren/HiC_Samples/Annotations/COMMON_GENES_PATHWAYS/"
	if path_choice == 5:
		dirs = "/cbcb/personal-scratch/hiren/ftp_umiacs_hiren/HiC_Samples/Annotations/"+path_dir+"/"+t_name+"/"
	else:	
		dirs = "/cbcb/personal-scratch/hiren/ftp_umiacs_hiren/HiC_Samples/Annotations/"+path_dir+"/"
	
	fls = os.listdir(dirs)
	FG_PATHWAYS = []
	ABSENT = []
	
	PATH_EDGE_GENES = np.array([], dtype='S')
	PATH_EDGES = np.array([], dtype='S')
	NON_PATH_EDGES = np.array([], dtype='S')
	OUT_PATH_GENES = np.array([], dtype='S')
	INTER_PATHWAYS = {}
	INTRA_PATHWAYS_INTERACTIONS = {}
	
	for fl in fls:
		path_name = fl.split(".")[0]
		if path_name not in Finished_Pathway:
			
			#(pflag, p_chr, p_genes, pg_mid) = getPathwaysInfo_With_Unique(fl)
			(pflag, p_chr, p_genes, pg_mid, pg_gid) = getPathwaysInfo_Without_Unique(fl, dirs)
			
			if pflag:
				(Edges, Nodes, Edge_no, Node_no, Edges_Node_No, Edge_Genes, Int_Chr, flag) = getEdgeDecision(pg_gid, p_chr, GD)
				if flag:
					#if path_name != "Housekeeping_Genes":
					
					if path_choice == 6:
						NN = subsetNonInteraction(Edges, Nodes)
						#NN = subsetNonInteractionFraction(Edges, Nodes)
						NN = np.array(NN.edges(), dtype='S')
						if NN.shape[0] > 0:
							NON_PATH_EDGES = createArray_Vertical(NON_PATH_EDGES, NN)
						
					PATH_EDGE_GENES = createArray_Horizontal(PATH_EDGE_GENES, Edge_Genes)
					OUT_PATH_GENES = createArray_Horizontal(OUT_PATH_GENES, Nodes[0])
					PATH_EDGES = createArray_Vertical(PATH_EDGES, Edges[0])
					
					if path_choice == 9:
						INTER_PATHWAYS[path_name] = Nodes[0]
						INTRA_PATHWAYS_INTERACTIONS[path_name] =  Edges[0]
						
					if RNA_ANALYSIS:
						(EF, Nodes, Int_Chr) = getEdgeExpression(Node_no, Edge_no, Nodes, Edges, Int_Chr)
					else:		
						(EF, Nodes, Int_Chr) = getEdgeFractionValue(Node_no, Edge_no, Nodes, Int_Chr)
					
					FG_PATHWAYS.append(EF)
					
					ENo = np.sum(np.array(np.array(Edge_no, dtype='f')), dtype='i')
					NNo = np.sum(np.array(np.array(Node_no, dtype='f')), dtype='i')
						
					#if path_choice != 6 or path_choice != 7:
					if path_choice == 1 or path_choice == 2 or path_choice == 3 or path_choice == 4 or path_choice == 5 or path_choice == 11:

						if exp_control == 1 or exp_control == 3:
							#siter = 500
							(REF, RNODES, RCHRS) = getRandomInteraction(Nodes, Int_Chr, GD, STOT_C, STOT_L, STOT_G, siter, RNA_ANALYSIS)
						
						elif exp_control == 2:
							fiter = 1
							#siter = 100
							(REF, RNODES, RCHRS) = getRandomInteraction(Nodes, Int_Chr, GD, STOT_C, STOT_L, STOT_G, fiter, RNA_ANALYSIS)
							EF = float(REF[0])
							(REF, RNODES, RCHRS) = getRandomInteraction(RNODES, RCHRS, GD, STOT_C, STOT_L, STOT_G, siter, RNA_ANALYSIS)
						
						REFs = np.array(REF, dtype='f')
						
						if np.where(REFs > 0.0):
							pval = getPvalue(EF, REF, siter)
							Z = getZscore(REF, EF)
						else:
							pval = 0.0
							Z = -3.0
						
						e_vals = ",".join(np.array(REF, dtype='S'))
						print path_name, "\t", p_genes.size, "\t", NNo, "\t", ENo, "\t", Edges_Node_No, "\t", EF, "\t", pval, "\t", Z
						txt = path_name + "\t" + str(p_genes.size) + "\t" + str(NNo) + "\t" + str(ENo) + "\t" + str(Edges_Node_No) + "\t" + str(EF) + "\t" + str(pval) + "\t"+ str(Z) + "\t" + e_vals + "\n" 
						
						OFL.write(txt)
					
					#elif path_choice == 6 or path_choice == 7 or path_choice == 10:
					elif path_choice == 10:
						Edge_List = ",".join([str(Edges[0][i][0])+"-"+str(Edges[0][i][1]) for i in range(Edges[0].shape[0])])
						Nodes_List = ",".join(Nodes)
						txt = path_name + "\t" + str(p_genes.size) + "\t" + str(NNo) + "\t" + str(ENo) + "\t" + str(Edges_Node_No) + "\t" + str(Edge_List) + "\t" + str(Nodes_List) + "\n"	
						
						OFL.write(txt)
					
				else:
					ABSENT.append(path_name)
			else:
				ABSENT.append(path_name)
		else:
			print path_name, "\t FINISHED !!! "
	
	(IP, NIP, OP, IG, TG, IP_G, NIP_G, OP_G, IG_G, TG_G) = prepareExperimentReport(OUT_PATH_GENES, PATH_EDGES, NON_PATH_EDGES)
	
	return OFL, IP, NIP, OP, IG, TG, IP_G, NIP_G, OP_G, IG_G, TG_G, INTER_PATHWAYS, INTRA_PATHWAYS_INTERACTIONS		
		
			
	FG_PATHWAYS	= np.array(FG_PATHWAYS, dtype='f')
	ABSENT = np.array(ABSENT, dtype='S')
	
		
def getTotalInteraction_Gene(genes):
	G = np.unique(np.hstack((genes[:,0], genes[:,1])))
	return G

	
def getChrDictionary(data, TOTAL_GENES):
	CHR = {}
	LEN = {}
	MID = {}
	GEN = {}
	
	ANN_GENES = data[:,5]
	TOTAL_GENES = TOTAL_GENES[np.where(np.in1d(TOTAL_GENES, ANN_GENES))]
	data = data[np.where(np.in1d(ANN_GENES, TOTAL_GENES))]
	ANN_GENES = ANN_GENES[np.where(np.in1d(ANN_GENES, TOTAL_GENES))]
	
	Chrs = data[:,1]
	Chrs = np.core.defchararray.partition(Chrs, "_")[:,0]
	Chrs = Chrs[np.where(Chrs != "chrUn")]
	chr = np.unique(Chrs)
	
	for i in range(chr.size):
		d = data[np.where(Chrs == chr[i])][:,5]
		CHR[chr[i]] = d
	
	l = np.abs(np.array(data[:,2], dtype='i') - np.array(data[:,3], dtype='i'))
	
	m = np.array(data[:,2], dtype='i') + np.array(l/2, dtype='i')
	
	for i in range(l.size):
		LEN[ANN_GENES[i]] = l[i] 
		MID[ANN_GENES[i]] = m[i]
		GEN[ANN_GENES[i]] = data[i][1]
	
		
	return CHR, LEN, MID, GEN, TOTAL_GENES

	
def getGeneLength(TOTAL_GENES):
	
	fl = "/fs/ftp-umiacs/pub/hiren/Hiren/HiC_Samples/Annotations/TOTAL_GENES_TRANSCRIPTS.txt"
	data = np.loadtxt(fl, dtype='S', delimiter='\t')
	(CHR, LEN, MID, GEN, TOTAL_GENES) = getChrDictionary(data, TOTAL_GENES)
	print "\t 6. Collect Gene Length, Mid Point of the Gene and Coordinate Positions !!! "
	
	TOT_G = np.array(LEN.keys(), dtype='S')
	TOT_L = np.array(LEN.values(), dtype='i')
	
	return TOTAL_GENES, CHR, LEN, MID, GEN, TOT_G, TOT_L
	
def getGraph(genes1, genes2):
	genes = np.array((genes1, genes2)).T
	genes = removeNull(genes)
	GR = nx.Graph()
	GR.add_edges_from(genes)
	G.remove_edges_from(G.selfloop_edges())
	return GR
	
def getEdges(GR):
	E = np.array(GR.edges())
	return E

def removeRedundancy(TOTAL_INTERACTION):
	MG = nx.MultiDiGraph([(i[0],i[1]) for i in TOTAL_INTERACTION ])
	NEW = nx.DiGraph()
	for (u,v) in MG.edges():
		NEW.add_edge(u, v, weight=len(MG[u][v]))
	g = NEW.to_undirected()	
	MAIN_G = np.array(g.edges())
	return MAIN_G
	
def getFinishedPath(ofl):
	Finished_Pathway = np.array([], dtype='S')
	if os.path.isfile(ofl):
		if os.path.getsize(ofl) > 0:
			Finished_Pathway = np.loadtxt(ofl, dtype='S', delimiter='\t')[:,0]
		else:
			print "\t EMPTY FILE IS DELETED !!!"
			os.remove(ofl)
	return Finished_Pathway	

def checkFileMode(ofl, Finished_Pathway):
	if Finished_Pathway.shape[0] > 0:
		OFL = open(ofl, "a")
	else:
		OFL = open(ofl, "w")
	return OFL	

def filterGenes(classes):
	classes[np.where(classes == "3' UTR")] = "Gene"
	classes[np.where(classes == "5' UTR")] = "Gene"
	classes[np.where(classes == "TTS")] = "Gene"
	classes[np.where(classes == "exon")] = "Gene"
	classes[np.where(classes == "intron")] = "Gene"
	classes[np.where(classes == "promoter-TSS")] = "Gene"
	return classes

def getPromoter_Gene(data, classes1, classes2, names1, names2):
	idx = np.where(np.logical_or(np.logical_and(classes1 == "Gene", classes2 == "Gene"), np.logical_and(classes2 == "Gene", classes1 == "Gene")))
	gene1 = names1[idx]
	gene2 = names2[idx]
	pro_gene1 = classes1[idx]
	pro_gene2 = classes2[idx]
	data = data[idx]
	distance = data[:,9]
	return distance, pro_gene1, pro_gene2, gene1, gene2	 

def reorderGenes(pro_gene1, pro_gene2, gene1, gene2):
	genes = np.vstack((gene1, gene2)).T
	P_G =  np.vstack((pro_gene1, pro_gene2)).T
	return P_G, genes
	
def getGenes(data, pos1, pos2):
	classes1 = np.core.defchararray.partition(data[:,pos1], "(")[:,0]
	classes1 = np.core.defchararray.strip(classes1)
	names1 = np.core.defchararray.partition(data[:,pos1], "(")[:,2]
	names1 = np.core.defchararray.partition(names1, ",")[:,0]
	names1 = np.core.defchararray.replace(names1, ")", "")
	names1 = np.core.defchararray.replace(names1, "(", "")
	names1 = np.core.defchararray.strip(names1)
	
	
	classes2 = np.core.defchararray.partition(data[:,pos2], "(")[:,0]
	classes2 = np.core.defchararray.strip(classes2)
	names2 = np.core.defchararray.partition(data[:,pos2], "(")[:,2]
	names2 = np.core.defchararray.partition(names2, ",")[:,0]
	names2 = np.core.defchararray.replace(names2, ")", "")
	names2 = np.core.defchararray.replace(names2, "(", "")
	names2 = np.core.defchararray.strip(names2)
	
	classes1 = filterGenes(classes1)
	classes2 = filterGenes(classes2)
	
	(distance, pro_gene1, pro_gene2, gene1, gene2) = getPromoter_Gene(data, classes1, classes2, names1, names2)
	
	return gene1, gene2, distance

def removeNull(genes, distance):
	genes = np.delete(genes, np.where(np.logical_or(genes[:,0] == '', genes[:,1] == '')), axis=0)
	distance = np.delete(distance, np.where(np.logical_or(genes[:,0] == '', genes[:,1] == '')), axis=0)
	return genes, distance

def getWeights(graph):
    d={}
    for (u,v) in graph.edges():
        d[u,v]=graph[u][v]['weight']
    weights = np.array(d.values(), dtype='S')
    return weights
	
def getMainGraph(genes, distance):
	(genes, distance) = removeNull(genes, distance)
	EdgeList = [(str(genes[i][0]),str(genes[i][1]),int(distance[i])) for i in range(genes.shape[0])]
	G = nx.Graph()
	G.add_weighted_edges_from(EdgeList)
	print "\t 2. TOTAL GENE-GENE INTERACTIONS WITHOUT REMOVING SELFLOOP:: ", str(len(G.edges()))
	G.remove_edges_from(G.selfloop_edges())
	weights = getWeights(G)
	
	return G, weights

def saveGenes(G, weights, odir, t_name):
	ofl = odir + t_name + ".txt"
	T = np.array(G.edges(), dtype='S')
	F = np.array((T[:,0], T[:,1], weights), dtype='S').T
	np.savetxt(ofl, F, fmt="%s", delimiter='\t')

def filterPoints(g):
	g = np.core.defchararray.partition(g, ".")[:,0]
	g = np.core.defchararray.partition(g, "-")[:,0]
	return g

def getFilteredGenes(g1, g2):

	ag1 = np.array([G_DATA[g1[i]] for i in range(g1.size)])
	ag2 = np.array([G_DATA[g2[i]] for i in range(g2.size)])
	
	return ag1, ag2	

def filterFDRdata(data, INT_NO):
	PVAL = 10**np.array(data[:,13], dtype='f')
	IDX = np.where(data[:,0])[0] + 1
	FDRs = (PVAL*INT_NO)/IDX
	data = data[np.where(FDRs <= float(fdrval))]
	return data
	
def customerFDRfilter(alpha,data,V):
	"""
	Controls for false discovery rate using the Benjamin and Hochberg procedure.
	Given a nominal Type 1 error rate alpha and a list of nominal p-values,
	returns a corresponding list of booleans, set to True only if the null 
	hypothesis should still be rejected after controlling the false discovery rate. 
	"""
	m = float(V)
	PVAL = 10**np.array(data[:,13], dtype='f')
	PVAL = sorted(PVAL)
	reject = [False for p in PVAL]
	Q = []
	for k, p in enumerate(PVAL):
		l = k+1
		#print_(p,k*alpha/m)
		Q.append(k*alpha/m)
		if p <= k*alpha/m:
			reject[k] = True
	reject = np.array(reject)
	Q = np.array(Q)
	data = data[np.where(reject)]
	return data
	
def getTotalValues(s1, s2):
	if resolution_type == 1:
		fl1 = "/cbcb/project2-scratch/hiren/HiC_Analysis/HiC_Results/"+s1+"/HiC_Annotation/interactionAnnotation.txt"	
		fl2 = "/cbcb/project2-scratch/hiren/HiC_Analysis/HiC_Results/"+s2+"/HiC_Annotation/interactionAnnotation.txt"
		
		if fdrval == 0:
			data1_t = getData(fl1)
		elif fdrval > 0:
			data1_t, INT_NO = getFDRData(fl1)
			#data1_t = data1_t[np.where(np.array(data1_t[:,24], dtype='f') <= float(fdrval))]
			#data1_t = filterFDRdata(data1_t, INT_NO)
			data1_t = customerFDRfilter(0.1, data1_t, INT_NO)
		print "\t\t a. First of Replicate of 10K Loaded !!!! "	
		
		if fdrval == 0:
			data2_t = getData(fl2)
		elif fdrval > 0:
			data2_t, INT_NO = getFDRData(fl2)
			#data2_t = data2_t[np.where(np.array(data2_t[:,24], dtype='f') <= float(fdrval))]
			#data2_t = filterFDRdata(data2_t, INT_NO)
			data2_t = customerFDRfilter(0.1, data2_t, INT_NO) 
			
		print "\t\t b. Second of Replicate of 10K Loaded !!!! "
	elif resolution_type == 2:
		fl3 = "/cbcb/project2-scratch/hiren/HiC_Analysis/HiC_Results/"+s1+"/HiC_Annotation_100k/interactionAnnotation.txt"	
		fl4 = "/cbcb/project2-scratch/hiren/HiC_Analysis/HiC_Results/"+s2+"/HiC_Annotation_100k/interactionAnnotation.txt"
		if fdrval == 0:
			data1_t = getData(fl3)
		elif fdrval > 0:
			data1_t, INT_NO = getFDRData(fl3)
			#data1_t = data1_t[np.where(np.array(data1_t[:,24], dtype='f') <= float(fdrval))]
			#data1_t = filterFDRdata(data1_t, INT_NO)
			data1_t = customerFDRfilter(0.1, data1_t, INT_NO)
		#data1_t = getData(fl3)
		print "\t\t c. First of Replicate of 100K Loaded !!!! "
		
		if fdrval == 0:
			data2_t = getData(fl4)
		elif fdrval > 0:
			data2_t, INT_NO = getFDRData(fl4)
			#data2_t = data2_t[np.where(np.array(data2_t[:,24], dtype='f') <= float(fdrval))]
			#data2_t = filterFDRdata(data2_t, INT_NO)
			data2_t = customerFDRfilter(0.1, data2_t, INT_NO)
		#data2_t = getData(fl4)
		print "\t\t d. Second of Replicate of 100K Loaded !!!! "
	else:
		fl1 = "/cbcb/project2-scratch/hiren/HiC_Analysis/HiC_Results/"+s1+"/HiC_Annotation/interactionAnnotation.txt"	
		fl2 = "/cbcb/project2-scratch/hiren/HiC_Analysis/HiC_Results/"+s2+"/HiC_Annotation/interactionAnnotation.txt"
		fl3 = "/cbcb/project2-scratch/hiren/HiC_Analysis/HiC_Results/"+s1+"/HiC_Annotation_100k/interactionAnnotation.txt"	
		fl4 = "/cbcb/project2-scratch/hiren/HiC_Analysis/HiC_Results/"+s2+"/HiC_Annotation_100k/interactionAnnotation.txt"
		#data1 = getData(fl1)
		if fdrval == 0:
			data1 = getData(fl1)
		elif fdrval > 0:
			data1, INT_NO = getFDRData(fl1)
			#data1 = data1[np.where(np.array(data1[:,24], dtype='f') <= int(fdrval))]
			#data1 = filterFDRdata(data1, INT_NO)
			data1 = customerFDRfilter(0.1, data1, INT_NO)

		print "\t\t a. First of Replicate of 10K Loaded !!!! "
		#data2 = getData(fl2)
		if fdrval == 0:
			data2 = getData(fl2)
		elif fdrval > 0:
			data2, INT_NO = getFDRData(fl2)
			#data2 = data2[np.where(np.array(data2[:,24], dtype='f') <= int(fdrval))]
			#data2 = filterFDRdata(data2, INT_NO)
			data2 = customerFDRfilter(0.1, data2, INT_NO)
			
		print "\t\t b. Second of Replicate of 10K Loaded !!!! "
		#data3 = getData(fl3)
		if fdrval == 0:
			data3 = getData(fl3)
		elif fdrval > 0:
			data3, INT_NO = getFDRData(fl3)
			#data3 = data3[np.where(np.array(data3[:,24], dtype='f') <= int(fdrval))]
			#data3 = filterFDRdata(data3, INT_NO)
			data3 = customerFDRfilter(0.1, data3, INT_NO)
			
		print "\t\t c. First of Replicate of 100K Loaded !!!! "
		#data4 = getData(fl4)
		if fdrval == 0:
			data4 = getData(fl4)
		elif fdrval > 0:
			data4, INT_NO = getFDRData(fl4)
			#data4 = data4[np.where(np.array(data4[:,24], dtype='f') <= int(fdrval))]
			#data4 = filterFDRdata(data4, INT_NO)
			data4 = customerFDRfilter(0.1, data4, INT_NO)
			
		print "\t\t d. Second of Replicate of 100K Loaded !!!! "
	
		data1_t = np.vstack((data1, data3))
		data2_t = np.vstack((data2, data4))
	
	data = np.vstack((data1_t, data2_t))
		
	if analysis_type == 1:
		data = data[np.where(data[:,9] != "interchromosomal")]
	if analysis_type == 2:
		data = data[np.where(data[:,9] == "interchromosomal")]
	return data

def getGCTotalValues(s1, s2):
	if resolution_type == 1:
		print "We have not processed the data....Please wait untill data uploads....Thanks "
		exit()
	elif resolution_type == 2:
		fl3 = "/cbcb/project2-scratch/hiren/HiC_Analysis/HiC_Results/"+s1+"/HiC_Annotation_100k_GC/interactionAnnotation.txt"	
		fl4 = "/cbcb/project2-scratch/hiren/HiC_Analysis/HiC_Results/"+s2+"/HiC_Annotation_100k_GC/interactionAnnotation.txt"
		#data1_t = getData(fl3)
		#print "\t\t c. First of Replicate of 100K Loaded !!!! "
		#data2_t = getData(fl4)
		#print "\t\t d. Second of Replicate of 100K Loaded !!!! "
		
		if fdrval == 0:
			data1_t = getData(fl3)
		elif fdrval > 0:
			data1_t, INT_NO = getFDRData(fl3)
			data1_t = customerFDRfilter(0.1, data1_t, INT_NO)
		print "\t\t c. First of Replicate of 100K Loaded !!!! "
		
		if fdrval == 0:
			data2_t = getData(fl4)
		elif fdrval > 0:
			data2_t, INT_NO = getFDRData(fl4)
			data2_t = customerFDRfilter(0.1, data2_t, INT_NO)
		
	else:
		print "We have not processed the data....Please wait untill data uploads....Thanks "
		exit()
	data = np.vstack((data1_t, data2_t))
		
	if analysis_type == 1:
		data = data[np.where(data[:,9] != "interchromosomal")]
	if analysis_type == 2:
		data = data[np.where(data[:,9] == "interchromosomal")]
	
	return data
	
def getGenePairData(cfile):
	"""
	This function is designed to avoid heavy files of HiC to be loaded and its faster that captures Gene-Gene interactions in transformed Gene IDs!!!.
	"""	
	T = np.loadtxt(cfile, dtype='S', delimiter='\t')
	g1, g2 = T[:,0], T[:,1]
	ID = np.zeros(g1.size, dtype=np.int)
	TOTAL_GENES = np.unique(np.hstack((g1, g2)))
	print "\t 4. Number of Unique Gene-Gene Interactions :: ", T.shape[0]
	print "\t 5. Total Genes involved in Gene-Gene Interaction :: ", TOTAL_GENES.size
	GE = np.array((g1, g2)).T
	GD = nx.Graph()
	GD.add_edges_from(GE)
	GD.remove_edges_from(GD.selfloop_edges())
	return g1, g2, GD, ID, TOTAL_GENES
	
def getTotalData(s1, s2, t_name, odir, resolution):
	"""
	This function is designed if you want to get gene data from HiC annotation file. If you have created gene-gene interaction sets then you can call "getGenePairData" function in order to avoid heavy files to be loaded and its faster!!!.
	"""
	#cdir = "/cbcb/project2-scratch/hiren/HiC_Analysis/HiC_Results/COMPARATIVE/10_KB_RESULTS/FINAL_"+resolution+"KB/ALL_INTERACTIONS/"
	#cdir = "/cbcb/project2-scratch/hiren/HiC_Analysis/HiC_Results/COMPARATIVE/10_KB_RESULTS/REVISIONS/100KB_REVISION/FINAL_INTERACTIONS/"
	if fdrval == 0:
		cdir = odir+"FINAL_INTERACTIONS/"
		cfile = cdir+t_name+"_"+resolution+"KB_Gene_Gene.txt"
	elif fdrval > 0:
		cdir = odir+"FINAL_INTERACTIONS/"
		cfile = cdir+t_name+"_"+resolution+"_CG_FDR_"+str(fdrval)+"_KB_Gene_Gene.txt"
	print "IS FILE THERE:: ", cfile
	if os.path.isdir(cdir) and os.path.isfile(cfile) and os.stat(cfile).st_size > 0:
		(g1, g2, GD, ID, TOTAL_GENES) = getGenePairData(cfile)
	else:
		if not os.path.isdir(cdir):
			os.mkdir(cdir)
		
		if gcval == 1:
			data = getTotalValues(s1, s2)
		elif gcval == 2:
			data = getGCTotalValues(s1, s2)
		
		(g1, g2, ID) = getGenes(data, 16, 20)

		print "\t 2. Unfiltered Transcript-Transcript Interactions :: ", g1.size
			
		GD = nx.Graph()
		GD.add_edges_from(np.array((g1, g2)).T)
		GD.remove_edges_from(GD.selfloop_edges())
		
		print "\t 3. Number of Unique Transcript-Transcript Interactions :: ", len(GD.edges())
		
		(g1, g2) = np.array(GD.edges())[:,0], np.array(GD.edges())[:,1]
			
		Tg1 = filterPoints(g1)
		Tg2 = filterPoints(g2)
		
		(g1, g2) = getFilteredGenes(Tg1, Tg2)
		
		GD = nx.Graph()
		GD.add_edges_from(np.array((g1, g2)).T)
		GD.remove_edges_from(GD.selfloop_edges())
		
		T = np.array(GD.edges())
		(g1, g2) = (T[:,0], T[:,1]) 
		
		print "\t 4. Number of Unique Gene-Gene Interactions :: ", T.shape[0]
		
		TOTAL_GENES = np.unique(np.hstack((T[:,0], T[:,1])))
				
		print "\t 5. Unique Genes involved in the above Interactions :: ", TOTAL_GENES.size
		
		
		GE = np.array((g1, g2)).T
		GD = nx.Graph()
		GD.add_edges_from(GE)
		GD.remove_edges_from(GD.selfloop_edges())
		np.savetxt(cfile, np.array(GD.edges()), fmt="%s", delimiter='\t')
		
	return g1, g2, GD, ID, TOTAL_GENES

def getExpressionData(rdir, rs):
	#fl = rdir+rs+"/AVERAGE_REPLICATE_EXPRESSION.txt"
	#fl = rdir+rs+"_Expression/AVERAGE_REPLICATE_EXPRESSION_WITH_NOVEL_GENES.txt"
	fl = rdir+rs+"/"+rs+"_FINAL_GENES_FPKM.txt"
	#fl = rdir+rs+"/"+rs+"_REVISED_FINAL_GENES_FPKM_ONLY_CUFFLINK_ISO.txt"
	data = np.loadtxt(fl, dtype='S', delimiter='\t')
	rgene = np.array(data[:,0], dtype='S')
	rfpkm = np.array(data[:,1], dtype='f')
	return rgene, rfpkm
	

def getval(gdata, M, i, TOT_C, TOT_L):
	G_DATA[gdata[i][0]] = gdata[i][5]
	T_DATA[gdata[i][5]] = gdata[i][0]
	GEN_LEN[gdata[i][5]] = TOT_L[i]
	CHR_GEN[TOT_C[i]] = gdata[i][5]
	GEN_CHR[gdata[i][5]] = TOT_C[i]
	GEN_MID[gdata[i][5]] = M[i]
	
def getGeneData():
	#fl = "/fs/ftp-umiacs/pub/hiren/Hiren/HiC_Samples/Annotations/TOTAL_GENES.txt"
	#fl = "/fs/ftp-umiacs/pub/hiren/Hiren/HiC_Samples/Annotations/TOTAL_GENES_TRANSCRIPTS.txt"
	fl = "/cbcb/personal-scratch/hiren/ftp_umiacs_hiren/HiC_Samples/Annotations/TOTAL_GENES_TRANSCRIPTS.txt"
	#fl = "/cbcb/personal-scratch/hiren/ftp_umiacs_hiren/HiC_Samples/Annotations/TOTAL_GENES_TRANSCRIPTS_OLD.txt"
	gdata = np.loadtxt(fl, dtype='S', delimiter='\t')
	TOT_L = np.abs(np.array(gdata[:,2], dtype='i') - np.array(gdata[:,3], dtype='i'))
	
	M = np.array(gdata[:,2], dtype='i') + TOT_L/2
	TOT_C = gdata[:,1]
	TOT_C = np.core.defchararray.partition(TOT_C, "_")[:,0]
	TOT_G = gdata[:,5]
	[getval(gdata, M, i, TOT_C, TOT_L) for i in range(gdata.shape[0])]
	idx = np.where(np.core.defchararray.find(np.array(G_DATA.keys()), "NM_") == 0)
	NC_TR = np.unique(np.array(G_DATA.keys())[idx])
	NC_GN = np.unique(np.array(G_DATA.values())[idx])
	
	return G_DATA, T_DATA, GEN_LEN, CHR_GEN, GEN_MID, GEN_CHR, TOT_C, TOT_G, TOT_L, NC_TR, NC_GN

def getProteinData():
	if protein_db == 1:
		fl = "/cbcb/personal-scratch/hiren/ftp_umiacs_hiren/HiC_Samples/Annotations/Protein_Interactions/PROTEIN_PROTEIN_INTERACTION_ALL.txt"
	elif protein_db == 2:
		fl = "/cbcb/personal-scratch/hiren/ftp_umiacs_hiren/HiC_Samples/Annotations/Protein_Interactions/HPRD_GENE_INTERACTION.txt"
	elif protein_db == 3:
		fl = "/cbcb/personal-scratch/hiren/ftp_umiacs_hiren/HiC_Samples/Annotations/Protein_Interactions/HPRD_AND_STRING_GENE_INTERACTION.txt"
	else:
		fl = "/cbcb/personal-scratch/hiren/ftp_umiacs_hiren/HiC_Samples/Annotations/Protein_Interactions/HPRD_GENE_INTERACTION.txt"
	data = np.loadtxt(fl, dtype='S', delimiter='\t')
	PG = nx.Graph()
	PG.add_edges_from(data)
	return PG	

def getSubProteinNetwork(GLarge, GSmall):
	#S = GLarge.subgraph(GSmall)
	S = GSmall.subgraph(GLarge)
	return S

def getSubProteinNetworkOther(PPE, OG):
	#S = GLarge.subgraph(GSmall)
	#S = PPE.copy()
	#S.remove_edges_from(OG.edges())
	S1 = PPE.copy()
	S2 = PPE.copy()
	
	S1.remove_edges_from(OG.edges())
	S2.remove_edges_from(S1.edges())
	
	return S2	



def getCommonSubnetwork(PPE, IP, NIP, OP, IG, TG, OFL):
	
	
	IP_P = getSubProteinNetworkOther(PPE, IP)
	
	if path_choice == 6:
		NIP_P = getSubProteinNetworkOther(PPE, NIP)
		D2 = len(NIP.edges()) - len(NIP_P.edges())
		X2 = float(NIP_P.number_of_edges())
		Y2 = NIP.number_of_edges()
		FR2 = (X2/Y2)*1000
		AP, AO = stats.fisher_exact([[len(IP_P.edges()), len(NIP_P.edges())], [D1, D2]])
		EP, EO = stats.fisher_exact([[len(NIP_P.edges()), len(OP_P.edges())], [D2, D3]])
		FP, FO = stats.fisher_exact([[len(NIP_P.edges()), len(IG_P.edges())], [D2, D4]])
		GP, GO = stats.fisher_exact([[len(NIP_P.edges()), len(TG_P.edges())], [D2, D5]])
	else:
		NIP_P = getSubProteinNetworkOther(PPE, NIP)
		D2 = 0; X2 = 0; Y2 = 0; FR2 = 0; AP = 0; AO = 0; EP = 0; EO = 0;FP = 0; FO = 0; GP = 0; GO = 0;
	
	OP_P = getSubProteinNetworkOther(PPE, OP)
	IG_P = getSubProteinNetworkOther(PPE, IG)
	TG_P = getSubProteinNetworkOther(PPE, TG)
	
	D1 = len(IP.edges()) - len(IP_P.edges())
	D3 = len(OP.edges()) - len(OP_P.edges())
	D4 = len(IG.edges()) - len(IG_P.edges())
	D5 = len(TG.edges()) - len(TG_P.edges())
	
	X1 = float(IP_P.number_of_edges())
	X3 = float(OP_P.number_of_edges())
	X4 = float(IG_P.number_of_edges())
	X5 = float(TG_P.number_of_edges())
	
	Y1 = IP.number_of_edges()
	Y3 = OP.number_of_edges()
	Y4 = IG.number_of_edges()
	Y5 = TG.number_of_edges()
	
	FR1 = (X1/Y1)*1000
	FR3 = (X3/Y3)*1000
	FR4 = (X4/Y4)*1000
	FR5 = (X5/Y5)*1000	
	
	
	BP, BO = stats.fisher_exact([[len(IP_P.edges()), len(OP_P.edges())], [D1, D3]])
	CP, CO = stats.fisher_exact([[len(IP_P.edges()), len(IG_P.edges())], [D1, D4]])
	DP, DO = stats.fisher_exact([[len(IP_P.edges()), len(TG_P.edges())], [D1, D5]])
	HP, HO = stats.fisher_exact([[len(OP_P.edges()), len(IG_P.edges())], [D3, D4]])
	JP, JO = stats.fisher_exact([[len(OP_P.edges()), len(TG_P.edges())], [D3, D5]])
	KP, KO = stats.fisher_exact([[len(IG_P.edges()), len(TG_P.edges())], [D4, D5]])
	
	
	return IP_P, NIP_P, OP_P, IG_P, TG_P, D1, D2, D3, D4, D5, X1, X2, X3, X4, X5, Y1, Y2, Y3, Y4, Y5, FR1, FR2, FR3, FR4, FR5, AP, AO, BP, BO, CP, CO, DP, DO, EP, EO, FP, FO, GP, GO, HP, HO, JP, JO, KP, KO
	
	
def makeCorrelationReport(GCE, IP, NIP, OP, IG, TG, OFL, corr_term):
		
	(IP_P, NIP_P, OP_P, IG_P, TG_P, D1, D2, D3, D4, D5, X1, X2, X3, X4, X5, Y1, Y2, Y3, Y4, Y5, FR1, FR2, FR3, FR4, FR5, AP, AO, BP, BO, CP, CO, DP, DO, EP, EO, FP, FO, GP, GO, HP, HO, JP, JO, KP, KO) = getCommonSubnetwork(GCE, IP, NIP, OP, IG, TG, OFL)
	
	print IP_P.edges()[0:10]
	
	print "\t\t In-Pathways Total Interactions :: ", len(IP.edges()), "\t", corr_term, " Co-expression Found :: ", len(IP_P.edges()), "\t Co-expression Not found :: ", D1, "\t", "FRACTION :: ", FR1
	
	print "\t\t Out-Pathways Total Interactions :: ", len(OP.edges()), "\t", corr_term, " Co-expression Found :: ", len(OP_P.edges()), "\t Co-expression Not found :: ", D3, "\t", "FRACTION :: ", FR3
	print "\t\t In-Genome Total Interactions :: ", len(IG.edges()), "\t", corr_term, " Co-expression Found :: ", len(IG_P.edges()), "\t Co-expression Not found :: ", D4, "\t", "FRACTION :: ", FR4
	print "\t\t Total-Genome Non HiC Interactions :: ", len(TG.edges()), "\t", corr_term, " Co-expression Found :: ", len(TG_P.edges()), "\t Co-expression Not found :: ", D5, "\t", "FRACTION :: ", FR5
	
	print "in Pathways Interaction Vs. in Pathways Non-Interaction P-value & ODD Ratio :: ", AP, "\t", AO
	print "in Pathways Interaction Vs. out Pathways Interaction P-value & ODD Ratio :: ", BP, "\t", BO
	print "in Pathways Interaction Vs. in-Genome Interaction P-value & ODD Ratio :: ", CP, "\t", CO
	print "in Pathways Interaction Vs. tot-Genome Interaction P-value & ODD Ratio :: ", DP, "\t", DO
	
	print "out Pathways Interaction Vs. in Genome Interaction P-value & ODD Ratio :: ", HP, "\t", HO
	print "out Pathways Interaction Vs. tot Genome Interaction P-value & ODD Ratio :: ", JP, "\t", JO
	
	print "in Genome Interaction Vs. tot Genome Interaction P-value & ODD Ratio :: ", KP, "\t", KO
	
	
	txt = str(FR1) + "\t" + str(AP) + "\t" + str(BP) + "\t" + str(CP) + "\t" + str(DP) + "\t" + str(AO) + "\t" + str(BO) + "\t" + str(CO) + "\t" + str(DO) + "\n"
	OFL.write(txt)
	txt = str(FR3) + "\t" + str(HP) + "\t" + str(JP) + "\t" + str(HO) + "\t" + str(JO) + "\n"
	OFL.write(txt)
	txt = str(FR4) + "\t" + str(KP) + "\t" + str(KO) + "\n"
	OFL.write(txt)
	txt = str(FR5) + "\n"
	OFL.write(txt)
	
	IP_Coex = np.array([w['weight'] for n, v, w in IP_P.edges(data=True)])
	txt = "\t".join(np.array(IP_Coex, dtype='S')) + "\n"
	OFL.write(txt)
	OP_Coex = np.array([w['weight'] for n, v, w in OP_P.edges(data=True)])
	txt = "\t".join(np.array(OP_Coex, dtype='S')) + "\n"
	OFL.write(txt)
	IG_Coex = np.array([w['weight'] for n, v, w in IG_P.edges(data=True)])
	txt = "\t".join(np.array(IG_Coex, dtype='S')) + "\n"
	OFL.write(txt)
	TG_Coex = np.array([w['weight'] for n, v, w in TG_P.edges(data=True)])
	txt = "\t".join(np.array(TG_Coex, dtype='S')) + "\n"
	OFL.write(txt)
		
	return OFL
	
def makeProtProtReport(PIE, IP, NIP, OP, IG, TG, OFL):
	
	IP_P = getSubProteinNetworkOther(PIE, IP)
	NIP_P = getSubProteinNetworkOther(PIE, NIP)
	OP_P = getSubProteinNetworkOther(PIE, OP)
	IG_P = getSubProteinNetworkOther(PIE, IG)
	TG_P = getSubProteinNetworkOther(PIE, TG)
	
	
	D1 = len(IP.edges()) - len(IP_P.edges())
	D2 = len(NIP.edges()) - len(NIP_P.edges())
	D3 = len(OP.edges()) - len(OP_P.edges())
	D4 = len(IG.edges()) - len(IG_P.edges())
	D5 = len(TG.edges()) - len(TG_P.edges())
	
	X1 = float(IP_P.number_of_edges())
	X2 = float(NIP_P.number_of_edges())
	X3 = float(OP_P.number_of_edges())
	X4 = float(IG_P.number_of_edges())
	X5 = float(TG_P.number_of_edges())
	
	Y1 = IP.number_of_edges()
	Y2 = NIP.number_of_edges()
	Y3 = OP.number_of_edges()
	Y4 = IG.number_of_edges()
	Y5 = TG.number_of_edges()
	
	FR1 = (X1/Y1)*1000
	FR2 = (X2/Y2)*1000
	FR3 = (X3/Y3)*1000
	FR4 = (X4/Y4)*1000
	FR5 = (X5/Y5)*1000	
	
	print "In-Pathways :: ", IP_P.number_of_edges()
	print "in-N Pathways :: ", NIP_P.number_of_edges()
	print "out Pathways :: ", OP_P.number_of_edges()
	print "in-Genome :: ", IG_P.number_of_edges()
	print "Total Genome :: ", TG_P.number_of_edges()
	
	print "\t\t In-Pathways Total Interactions :: ", len(IP.edges()), "\t Protein-Protein Found :: ", len(IP_P.edges()), "\t Protein-Protein Not found :: ", D1, "\t", "FRACTION :: ", FR1
	print "\t\t In-Pathways Total Non Interactions :: ", len(NIP.edges()), "\t Protein-Protein Found :: ", len(NIP_P.edges()), "\t Protein-Protein Not found :: ", D2, "\t", "FRACTION :: ", FR2
	print "\t\t Out-Pathways Total Interactions :: ", len(OP.edges()), "\t Protein-Protein Found :: ", len(OP_P.edges()), "\t Protein-Protein Not found :: ", D3, "\t", "FRACTION :: ", FR3
	print "\t\t In-Genome Total Interactions :: ", len(IG.edges()), "\t Protein-Protein Found :: ", len(IG_P.edges()), "\t Protein-Protein Not found :: ", D4, "\t", "FRACTION :: ", FR4
	print "\t\t Total-Genome Non HiC Interactions :: ", len(TG.edges()), "\t Protein-Protein Found :: ", len(TG_P.edges()), "\t Protein-Protein Not found :: ", D5, "\t", "FRACTION :: ", FR5
	"""
	AP, AO = stats.fisher_exact([[len(IP_P.edges()), len(NIP_P.edges())], [D1, D2]])
	BP, BO = stats.fisher_exact([[len(IP_P.edges()), len(OP_P.edges())], [D1, D3]])
	CP, CO = stats.fisher_exact([[len(IP_P.edges()), len(IG_P.edges())], [D1, D4]])
	DP, DO = stats.fisher_exact([[len(IP_P.edges()), len(TG_P.edges())], [D1, D5]])
	
	
	EP, EO = stats.fisher_exact([[len(NIP_P.edges()), len(OP_P.edges())], [D2, D3]])
	FP, FO = stats.fisher_exact([[len(NIP_P.edges()), len(IG_P.edges())], [D2, D4]])
	GP, GO = stats.fisher_exact([[len(NIP_P.edges()), len(TG_P.edges())], [D2, D5]])
	
	HP, HO = stats.fisher_exact([[len(OP_P.edges()), len(IG_P.edges())], [D3, D4]])
	IP, IO = stats.fisher_exact([[len(OP_P.edges()), len(TG_P.edges())], [D3, D5]])
	
	JP, JO = stats.fisher_exact([[len(IG_P.edges()), len(TG_P.edges())], [D4, D5]])
	"""
	
	"""	
	print "in Pathways Interaction Vs. in Pathways Non-Interaction P-value & ODD Ratio :: ", AP, "\t", AO
	print "in Pathways Interaction Vs. out Pathways Interaction P-value & ODD Ratio :: ", BP, "\t", BO
	print "in Pathways Interaction Vs. in-Genome Interaction P-value & ODD Ratio :: ", CP, "\t", CO
	print "in Pathways Interaction Vs. tot-Genome Interaction P-value & ODD Ratio :: ", DP, "\t", DO
	
	print "in Pathways Non-Interaction Vs. out Pathways Interaction P-value & ODD Ratio :: ", EP, "\t", EO
	print "in Pathways Non-Interaction Vs. in Genome Interaction P-value & ODD Ratio :: ", FP, "\t", FO
	print "in Pathways Non-Interaction Vs. tot-Genome Interaction P-value & ODD Ratio :: ", GP, "\t", GO
	
	print "out Pathways Interaction Vs. in Genome Interaction P-value & ODD Ratio :: ", HP, "\t", HO
	print "out Pathways Interaction Vs. tot Genome Interaction P-value & ODD Ratio :: ", IP, "\t", IO
	
	print "in Genome Interaction Vs. tot Genome Interaction P-value & ODD Ratio :: ", JP, "\t", JO
	"""
	
	AP = stats.chi2_contingency([[len(IP_P.edges()), len(NIP_P.edges())], [D1, D2]])[1]
	BP = stats.chi2_contingency([[len(IP_P.edges()), len(OP_P.edges())], [D1, D3]])[1]
	CP = stats.chi2_contingency([[len(IP_P.edges()), len(IG_P.edges())], [D1, D4]])[1]
	DP = stats.chi2_contingency([[len(IP_P.edges()), len(TG_P.edges())], [D1, D5]])[1]
	
	
	EP = stats.chi2_contingency([[len(NIP_P.edges()), len(OP_P.edges())], [D2, D3]])[1]
	FP = stats.chi2_contingency([[len(NIP_P.edges()), len(IG_P.edges())], [D2, D4]])[1]
	GP = stats.chi2_contingency([[len(NIP_P.edges()), len(TG_P.edges())], [D2, D5]])[1]
	
	HP = stats.chi2_contingency([[len(OP_P.edges()), len(IG_P.edges())], [D3, D4]])[1]
	IP = stats.chi2_contingency([[len(OP_P.edges()), len(TG_P.edges())], [D3, D5]])[1]
	
	JP = stats.chi2_contingency([[len(IG_P.edges()), len(TG_P.edges())], [D4, D5]])[1]
	
	print "in Pathways Interaction Vs. in Pathways Non-Interaction P-value :: ", AP
	print "in Pathways Interaction Vs. out Pathways Interaction P-value :: ", BP
	print "in Pathways Interaction Vs. in-Genome Interaction P-value    :: ", CP
	print "in Pathways Interaction Vs. tot-Genome Interaction P-value :: ", DP
	
	print "in Pathways Non-Interaction Vs. out Pathways Interaction P-value :: ", EP
	print "in Pathways Non-Interaction Vs. in Genome Interaction P-value :: ", FP
	print "in Pathways Non-Interaction Vs. tot-Genome Interaction P-value :: ", GP
	
	print "out Pathways Interaction Vs. in Genome Interaction P-value :: ", HP
	print "out Pathways Interaction Vs. tot Genome Interaction P-value :: ", IP
	
	print "in Genome Interaction Vs. tot Genome Interaction P-value :: ", JP
	
	txt = str(FR1) + "\t" + str(AP) + "\t" + str(BP) + "\t" + str(CP) + "\t" + str(DP) + "\n"
	OFL.write(txt)
	txt = str(FR2) + "\t" + str(EP) + "\t" + str(FP) + "\t" + str(GP) + "\n"
	OFL.write(txt)
	txt = str(FR3) + "\t" + str(HP) + "\t" + str(IP) + "\n"
	OFL.write(txt)
	txt = str(FR4) + "\t" + str(JP) + "\n"
	OFL.write(txt)
	txt = str(FR5) + "\n"
	OFL.write(txt)
	return OFL

	
def commonEdgeCollector(p_genes, chrs, GD):
	(Edges, Nodes, Edge_no, Node_no, Edges_Node_No, Edge_Genes, Int_Chr, flag) = getEdgeDecision(p_genes, chrs, GD)
	return Edges[0], Nodes[0]

def interPathwaysRoutine(oEdges, iEdges):
	iG = nx.Graph()
	oG = nx.Graph()
	iG.add_edges_from(iEdges)
	oG.add_edges_from(oEdges)
	oG.remove_edges_from(iG.edges())
	
	Node_no = np.array(len(oG.nodes()), dtype='i')
	Edge_no = np.array(len(oG.edges()), dtype='i')
	Nodes = np.array(oG.nodes(), dtype='S')
		
	Int_Chr = [GEN_CHR[i] for i in Nodes]
	Edges = np.array([oG.edges()], dtype='S')
	return Node_no, Edge_no, Nodes, Edges, Int_Chr
	
def interRandomPathwaysRoutine(F_Nodes1, F_Nodes2, Int_Chr1, Int_Chr2, GD, STOT_C, STOT_L, STOT_G, siter, RNA_ANALYSIS, ta, tb):
	REF = []
	F_Nodes1 = np.array(F_Nodes1, dtype='S')
	F_Nodes2 = np.array(F_Nodes2, dtype='S')
	Int_Chr1 = np.array(Int_Chr1, dtype='S')
	Int_Chr2 = np.array(Int_Chr2, dtype='S')
	
	F_Nodes1 = F_Nodes1[np.where(np.in1d(F_Nodes1, TOT_G))]
	F_Nodes2 = F_Nodes2[np.where(np.in1d(F_Nodes2, TOT_G))]
	Int_Chr1 = Int_Chr1[np.where(np.in1d(F_Nodes1, TOT_G))]
	Int_Chr2 = Int_Chr2[np.where(np.in1d(F_Nodes2, TOT_G))]
	
	Int_Chr1 = Int_Chr1[np.where(Int_Chr1 != "chrMT")]
	F_Nodes1 = F_Nodes1[np.where(Int_Chr1 != "chrMT")]
	
	Int_Chr2 = Int_Chr2[np.where(Int_Chr2 != "chrMT")]
	F_Nodes2 = F_Nodes2[np.where(Int_Chr2 != "chrMT")]
	
	Int_Chr1 = Int_Chr1[np.where(Int_Chr1 != "chrY")]
	F_Nodes1 = F_Nodes1[np.where(Int_Chr1 != "chrY")]
	Int_Chr2 = Int_Chr2[np.where(Int_Chr2 != "chrY")]
	F_Nodes2 = F_Nodes2[np.where(Int_Chr2 != "chrY")]
	
	F_Nodes1 = F_Nodes1[np.in1d(F_Nodes1, STOT_G)]
	F_Nodes2 = F_Nodes2[np.in1d(F_Nodes2, STOT_G)]
	Int_Chr1 = Int_Chr1[np.in1d(F_Nodes1, STOT_G)]
	Int_Chr2 = Int_Chr2[np.in1d(F_Nodes2, STOT_G)]
	
	f_finish = False
	r_temp = np.array([], dtype='S')
	
	F_GENES_1 = np.array(F_Nodes1, dtype='S') 
	idx = np.where(np.logical_not(np.in1d(STOT_G, F_GENES_1)))
	G_STOT_G_1 = STOT_G[idx]
	G_STOT_C_1 = STOT_C[idx]
	G_STOT_L_1 = STOT_L[idx]
	
	F_GENES_2 = np.array(F_Nodes2, dtype='S') 
	idx = np.where(np.logical_not(np.in1d(STOT_G, F_GENES_2)))
	G_STOT_G_2 = STOT_G[idx]
	G_STOT_C_2 = STOT_C[idx]
	G_STOT_L_2 = STOT_L[idx]

	for k in range(siter):
		Nodes1 = []
		Edge_no1 = []
		Node_no1 = []
		Edges1 = []
		
		Nodes2 = []
		Edge_no2 = []
		Node_no2 = []
		Edges2 = []
		
		
		if chr_control == 2:
			#r_genes1 = np.array([random.choice(getLenChoice(F_Nodes1[i], STOT_C, STOT_G, STOT_L)) for i in range(len(F_Nodes1))])
			#r_genes2 = np.array([random.choice(getLenChoice(F_Nodes2[i], STOT_C, STOT_G, STOT_L)) for i in range(len(F_Nodes2))])
			r_genes1 = np.array([random.choice(getLenChoice(F_Nodes1[i], G_STOT_C_1, G_STOT_G_1, G_STOT_L_1)) for i in range(len(F_Nodes1))])
			r_genes2 = np.array([random.choice(getLenChoice(F_Nodes2[i], G_STOT_C_2, G_STOT_G_2, G_STOT_L_2)) for i in range(len(F_Nodes2))])
		elif chr_control == 1:
			if ta == "Housekeeping_Genes" or tb == "Housekeeping_Genes":
				if ta == "Housekeeping_Genes":
					if f_finish:
						r_genes1 = r_temp
					else:	
						r_temp = np.array([random.choice(getLenAndChrChoice(F_Nodes1[i], Int_Chr1[i], G_STOT_C_1, G_STOT_G_1, G_STOT_L_1)) for i in range(len(F_Nodes1))])
						f_finish = True
						r_genes1 = r_temp
					r_genes2 = np.array([random.choice(getLenAndChrChoice(F_Nodes2[i], Int_Chr2[i], G_STOT_C_2, G_STOT_G_2, G_STOT_L_2)) for i in range(len(F_Nodes2))])
				
				else:
					r_genes1 = np.array([random.choice(getLenAndChrChoice(F_Nodes1[i], Int_Chr1[i], G_STOT_C_1, G_STOT_G_1, G_STOT_L_1)) for i in range(len(F_Nodes1))])	
					if f_finish:
						r_genes2 = r_temp
					else:
						r_temp = np.array([random.choice(getLenAndChrChoice(F_Nodes2[i], Int_Chr2[i], G_STOT_C_2, G_STOT_G_2, G_STOT_L_2)) for i in range(len(F_Nodes2))])
						f_finish = True
						r_genes2 = r_temp
			else:
				r_genes1 = np.array([random.choice(getLenAndChrChoice(F_Nodes1[i], Int_Chr1[i], G_STOT_C_1, G_STOT_G_1, G_STOT_L_1)) for i in range(len(F_Nodes1))])
				r_genes2 = np.array([random.choice(getLenAndChrChoice(F_Nodes2[i], Int_Chr2[i], G_STOT_C_2, G_STOT_G_2, G_STOT_L_2)) for i in range(len(F_Nodes2))])
		
		r_genes1 = r_genes1[np.in1d(r_genes1, STOT_G)]
		r_genes2 = r_genes2[np.in1d(r_genes2, STOT_G)]
			
		p_chr1 = Int_Chr1
		p_chr2 = Int_Chr2
				
		(Edges1, Nodes1, Edge_no1, Node_no1, Edges_Node_No1, Edge_Genes1, Int_Chr1, flag1) = getEdgeDecision(r_genes1, p_chr1, GD)
		(Edges2, Nodes2, Edge_no2, Node_no2, Edges_Node_No2, Edge_Genes2, Int_Chr2, flag2) = getEdgeDecision(r_genes2, p_chr2, GD)
		
		if np.array(Edges1[0]).shape[0] > 0 and np.array(Edges2[0]).shape[0] > 0:
			iEDGES = np.vstack((np.array(Edges1[0], dtype='S'), np.array(Edges2[0], dtype='S')))
		else:
			if np.array(Edges1[0]).shape[0] > 0:
				iEDGES = np.array(Edges1[0], dtype='S')
			else:
				iEDGES = np.array(Edges2[0], dtype='S')
		
		#dg1 = np.setdiff1d(Nodes1, Nodes2)
		#dg2 = np.setdiff1d(Nodes2, Nodes1)
		dg1 = np.unique(Nodes1)
		dg2 = np.unique(Nodes2)
		
		p_genes = np.unique(np.hstack((dg1, dg2)))
		p_chr = [GEN_CHR[i] for i in p_genes]
		(Edges, Nodes, Edge_no, Node_no, Edges_Node_No, Edge_Genes, Int_Chr, flag) = getEdgeDecision(p_genes, p_chr, GD)
		oEDGES = np.array(Edges[0], dtype='S')
		
		if oEDGES.shape[0] > 0:
			(Node_no, Edge_no, Nodes, Edges, Int_Chr) = interPathwaysRoutine(oEDGES, iEDGES)
			EF = Edge_no
		else:
			EF = 0
			
		REF.append(EF)
	REF = np.array(REF, dtype='i')
	
	return REF
	

def getInterPathwaysData(t_name, INTER_PATHWAYS, INTRA_PATHWAYS_INTERACTIONS, GD, OFL, STOT_C, STOT_L, STOT_G, RNA_ANALYSIS, path_dir, Finished_Pathway):
	
	Path_names = []
	dirs = "/cbcb/personal-scratch/hiren/ftp_umiacs_hiren/HiC_Samples/Annotations/"+path_dir+"/"
	fls = os.listdir(dirs)
	for fl in fls:
		path_name = fl.split(".")[0]
		Path_names.append(path_name)
	Path_names = np.array(Path_names, dtype='S')
	#Path_names = np.delete(Path_names, np.where(Path_names == "Housekeeping_Genes"))
	TOTAL = list(itertools.combinations(Path_names, 2))
	
	for a, b in TOTAL:
		fl1 = a + ".txt"
		fl2 = b + ".txt"
		
		path_name = a + "@" + b
		if path_name not in Finished_Pathway:
			
			(pflag1, p_chr1, p_genes1, pg_mid1, pg_gid1) = getPathwaysInfo_Without_Unique(fl1, dirs)
			(pflag2, p_chr2, p_genes2, pg_mid2, pg_gid2) = getPathwaysInfo_Without_Unique(fl2, dirs)
			(Edges1, Nodes1, Edge_no1, Node_no1, Edges_Node_No1, Edge_Genes1, Int_Chr1, flag1) = getEdgeDecision(pg_gid1, p_chr1, GD)
			(Edges2, Nodes2, Edge_no2, Node_no2, Edges_Node_No2, Edge_Genes2, Int_Chr2, flag2) = getEdgeDecision(pg_gid2, p_chr2, GD)
			
			if np.array(Edges1[0]).shape[0] > 0 and np.array(Edges2[0]).shape[0] > 0:
				iEDGES = np.vstack((np.array(Edges1[0], dtype='S'), np.array(Edges2[0], dtype='S')))
			else:
				if np.array(Edges1[0]).shape[0] > 0:
					iEDGES = np.array(Edges1[0], dtype='S')
				else:
					iEDGES = np.array(Edges2[0], dtype='S')
					
			#dg1 = np.setdiff1d(Nodes1, Nodes2)
			#dg2 = np.setdiff1d(Nodes2, Nodes1)
			dg1 = np.unique(Nodes1)
			dg2 = np.unique(Nodes2)
			
			p_genes = np.unique(np.hstack((dg1, dg2)))
			p_genes = p_genes[np.where(np.in1d(p_genes, TOT_G))]
			p_chr = [GEN_CHR[i] for i in p_genes]
			
			(Edges, Nodes, Edge_no, Node_no, Edges_Node_No, Edge_Genes, Int_Chr, flag) = getEdgeDecision(p_genes, p_chr, GD)
			oEDGES = np.array(Edges[0], dtype='S')
			
			
			(oNode_no, oEdge_no, oNodes, oEDGES, Int_Chr) = interPathwaysRoutine(oEDGES, iEDGES)
			
			if oEDGES.size > 0:
				#(Node_no, Edge_no, Nodes, Edges, Int_Chr) = interPathwaysRoutine(oEDGES, iEDGES)
				
							
				EF = np.sum(np.array(np.array(oEdge_no, dtype='f')), dtype='i')
				NNo = np.sum(np.array(np.array(oNode_no, dtype='f')), dtype='i')
				
				siter = 100
				REF = interRandomPathwaysRoutine(pg_gid1, pg_gid2, p_chr1, p_chr2, GD, STOT_C, STOT_L, STOT_G, siter, RNA_ANALYSIS, a, b)
				
				if np.where(REF > 0.0)[0].shape[0] > 0:
					pval = getPvalue(EF, REF, siter)
					Z = getZscore(REF, EF)
				else:
					pval = 0.0
					Z = 20.0
				
					
				e_vals = ",".join(np.array(REF, dtype='S'))
				
				print path_name, "\t", str(pg_gid1.size), "\t", str(pg_gid2.size), "\t", dg1.size, "\t", dg2.size, "\t", NNo, "\t", EF, "\t", pval, "\t", Z
				#txt = path_name + "\t" + str(pg_gid1.size) + "\t" + str(pg_gid2.size) + "\t" + str(dg1.size) + "\t" + str(dg2.size) + "\t" + str(NNo) + "\t" + str(EF) + "\t" + str(pval) + "\t"+ str(Z) + "\t" + e_vals + "\n" 
				txt = path_name + "\t" + str(pg_gid1.size) + "\t" + str(pg_gid2.size) + "\t" + str(dg1.size) + "\t" + str(dg2.size) + "\t" + str(NNo) + "\t" + str(EF) + "\t" + str(pval) + "\t"+ str(Z) + "\n" 
				
			else:
				print path_name, "\t", str(pg_gid1.size), "\t", str(pg_gid2.size), "\t", dg1.size, "\t", dg2.size, "\t PATHWAYS ABSENT "
				#txt = path_name + "\t" + str(pg_gid1.size) + "\t" + str(pg_gid2.size) + "\t" + str(dg1.size) + "\t" + str(dg2.size) + "\t" + str(0) + "\t" + str(0) + "\t" + str(1.0) + "\t"+ str(-3.0) + "\t" + "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0" + "\n" 
				txt = path_name + "\t" + str(pg_gid1.size) + "\t" + str(pg_gid2.size) + "\t" + str(dg1.size) + "\t" + str(dg2.size) + "\t" + str(0) + "\t" + str(0) + "\t" + str(1.0) + "\t"+ str(-3.0) + "\n" 
			OFL.write(txt)
		else:
			print path_name, "\t IS FINISHED !!! "
	return OFL	
			
			
def createInteractingGeneSets(IP, IP_G, OP, OP_G, OFL):
	IP1 = np.array(IP.edges(), dtype='S')
	IP2 = np.array(OP.edges(), dtype='S')
	
	A = ",".join(IP_G)
	B = ["-".join(IP1[i]) for i in range(IP1.shape[0])]
	B = ",".join(B)
	C = ",".join(OP_G)
	D = ["-".join(IP2[i]) for i in range(IP2.shape[0])]
	D = ",".join(D)
	
	txt1 = "IN-PATHWAYS GENES\t" + A + "\n"
	txt2 = "IN-PATHWAYS INTERACTIONS\t" + B + "\n"
	txt3 = "OUT-PATHWAYS GENES\t" + C + "\n"
	txt4 = "OUT-PATHWAYS INTERACTIONS\t" + D + "\n"
	
	OFL.write(txt1)
	OFL.write(txt2)
	OFL.write(txt3)
	OFL.write(txt4)
	
	return OFL
	
	
def getWeightVector(data):	
	A = [(str(data[i][0]), str(data[i][1]), float(data[i][2])) for i in range(data.shape[0])]
	G = nx.Graph()
	G.add_weighted_edges_from(A)
	return G

def getOneCoExpressionData(t):
	dirs = "/fs/ftp-umiacs/pub/hiren/Hiren/HiC_Samples/Annotations/Coexpression/"
	fl = dirs + t + "_HiC_Coexpressed.txt"
	data = np.loadtxt(fl, dtype='S', delimiter='\t')
	G = getWeightVector(data)
	return G
	
def getCoExpressionData(t):
	
	dirs = "/fs/ftp-umiacs/pub/hiren/Hiren/HiC_Samples/Annotations/Coexpression/"
	fl1 = dirs+t+"_POS_GID.coexp.txt"
	fl2 = dirs+t+"_NEG_GID.coexp.txt"
	fl3 = dirs+t+"_ZERO_GID.coexp.txt"
		
	if corr_type == 1:
		corr_term = "Positive correlated"
		data = np.loadtxt(fl1, dtype='S', delimiter='\t')[0:10000]
		Corr_G = getWeightVector(data)
		print "\t TOTAL NUMBER OF POSITIVE CORRELATED GENES OF THE TISSUE 			:: ", Corr_G.number_of_edges()
		print "\t TOTAL NUMBER OF POSITIVE CORRELATED UNIQUE GENES OF THE TISSUE 	:: ", np.unique(np.array(Corr_G.nodes())).size
	elif corr_type == 2:
		corr_term = "Negative correlated"
		data = np.loadtxt(fl2, dtype='S', delimiter='\t')
		Corr_G = getWeightVector(data)
		print "\t TOTAL NUMBER OF NEGATIVE CORRELATED GENES OF THE TISSUE 			:: ", Corr_G.number_of_edges()
		print "\t TOTAL NUMBER OF POSITIVE CORRELATED UNIQUE GENES OF THE TISSUE 	:: ", np.unique(np.array(Corr_G.nodes())).size
	elif corr_type == 3:
		corr_term = "Zero correlated"
		data = np.loadtxt(fl3, dtype='S', delimiter='\t')
		Corr_G = getWeightVector(data)
		print "\t TOTAL NUMBER OF NEGATIVE CORRELATED GENES OF THE TISSUE 			:: ", Corr_G.number_of_edges()
		print "\t TOTAL NUMBER OF POSITIVE CORRELATED UNIQUE GENES OF THE TISSUE 	:: ", np.unique(np.array(Corr_G.nodes())).size
	elif corr_type == 4:
		corr_term = "All type correlated"
		data1 = np.loadtxt(fl1, dtype='S', delimiter='\t')
		data2 = np.loadtxt(fl2, dtype='S', delimiter='\t')
		data1 = np.loadtxt(fl3, dtype='S', delimiter='\t')
		data = np.vstack((data1, data2, data3))
		Corr_G = getWeightVector(data)
	else:
		corr_term = "!!ERROR !!"
		print "\t YOU HAVE NOT CHOSEN AN APPROPRIATE CORRELATION OPTION :-( THANK YOU !! BYE "
		sys.exit()		
	
	return Corr_G	

	
def getFilteredData(TOTAL_GENES, TOT_C, TOT_G, TOT_L):
	
	STOT_C = TOT_C[np.in1d(TOT_G, TOTAL_GENES)]
	STOT_L = TOT_L[np.in1d(TOT_G, TOTAL_GENES)]
	STOT_G = TOT_G[np.in1d(TOT_G, TOTAL_GENES)]
	return STOT_C, STOT_L, STOT_G
	
def getGlobalArrays():
	G_DATA = {}
	T_DATA = {}
	GEN_LEN = {}
	CHR_GEN = {}
	GEN_MID = {}
	GEN_CHR = {}
	TOT_C = np.array([], dtype='S')
	TOT_G = np.array([], dtype='S')
	TOT_L = np.array([], dtype='i')
	RGENES = np.array([], dtype='i')
	RFPKM = np.array([], dtype='i')
	NC_TR = np.array([], dtype='S') 
	NC_GN = np.array([], dtype='S')
	return G_DATA, T_DATA, GEN_LEN, CHR_GEN, GEN_MID, GEN_CHR, TOT_C, TOT_G, TOT_L, RGENES, RFPKM, NC_TR, NC_GN

def choosePathwayType():
	if path_choice == 1:
		path_dir = "Gene_Kegg_Pathways"
		if chr_control == 1:
			oterm = "KEGG"
		elif chr_control == 2:
			oterm = "KEGG_Chr_Controlled"
			
	if path_choice == 2:
		path_dir = "Receptor_Induced_Pathways"
		if chr_control == 1:
			oterm = "Receptor_Induced"
		elif chr_control == 2:
			oterm = "Receptor_Induced_Chr_Controlled"
				
	if path_choice == 3:
		if sub_path_choice == 1:
			if corr_type == 2:
				path_dir = "Directred_Gene_Kegg_Pathways"
			else:
				path_dir = "Gene_Kegg_Pathways"
		elif sub_path_choice == 2:
			path_dir = "Receptor_Induced_Pathways"
			
		if chr_control == 1:
			oterm = "KEGG"
		elif chr_control == 2:
			oterm = "KEGG_Chr_Controlled"	
			
	if path_choice == 4:
		path_dir = "COMMON_GENES_PATHWAYS"
		if chr_control == 1:
			oterm = "KEGG"
		elif chr_control == 2:
			oterm = "KEGG_Chr_Controlled"
			
	if path_choice == 5:
		path_dir = "Joined_Pathways"
		if chr_control == 1:
			oterm = "KEGG"
		elif chr_control == 2:
			oterm = "KEGG_Chr_Controlled"
			
	if path_choice == 6:
		if sub_path_choice == 1:
			path_dir = "Gene_Kegg_Pathways"
		elif sub_path_choice == 2:
			path_dir = "Receptor_Induced_Pathways"
		if chr_control == 1:
			oterm = "KEGG"
		elif chr_control == 2:
			oterm = "KEGG_Chr_Controlled"
	
	if path_choice == 7:
		if sub_path_choice == 1:
			path_dir = "Gene_Kegg_Pathways"
		elif sub_path_choice == 2:
			path_dir = "Receptor_Induced_Pathways"
		
		if chr_control == 1:
			oterm = "KEGG"
		elif chr_control == 2:
			oterm = "KEGG_Chr_Controlled"
			
	if path_choice == 8 or path_choice == 9 or path_choice == 10:
		if sub_path_choice == 1:
			if path_choice == 10 and corr_type == 2:
				path_dir = "Directred_Gene_Kegg_Pathways"
			else:	
				path_dir = "Gene_Kegg_Pathways"
		elif sub_path_choice == 2:
			path_dir = "Receptor_Induced_Pathways"
		

		if chr_control == 1:
			oterm = "KEGG"
		elif chr_control == 2:
			oterm = "KEGG_Chr_Controlled"

	if path_choice == 11:
		path_dir = "WIKI_PATHWAYS"
		oterm = "NETPATH"		
			
	return oterm, path_dir

def getParameters():
	os.system('clear')
	resolution_type = raw_input("\t What Resolution you want to analyse the data ?? [10KB/100KB/Both] --> Type 1 OR 2 OR 3 for your selection !! :: ")
	analysis_type = raw_input("\t What type of Analysis do you want to do ?? [Intra-Chromosomal/Inter-Chromosomal/Both] --> Type 1 OR 2 OR 3 for your selection !! :: ")
	iterval = raw_input("\t What number of iteration you want to run for random simulation ?? [100/500/1000] --> Type 1 OR 2 OR 3 for your selection !! :: ")
	gcval = int(raw_input("\t Do you want to analyse the data with GC check OR without GC check ?? [without/with] --> Type 1 OR 2 for your selection !! :: "))
	fdrcont = int(raw_input("\t Do you want to Control FDR ?? [Yes/No] --> Type 1 OR 2 for your selection !! :: "))
	if fdrcont == 1:
		tempval = int(raw_input("\t\t\t At what cut-off FDR ?? [0.05/0.1/0.3] --> Type 1 OR 2 OR 3 for your selection !! :: "))
		if tempval == 1:
			fdrval = 0.05
		elif tempval == 2:
			fdrval = 0.1
		elif tempval == 3:
			fdrval = 0.3
	else:
		fdrval = 10.0
		
	print "\t What Pahtways Genes you want to work on ?? See the options below :: "
	print "\t\t 1. KEGG Pathways "
	print "\t\t 2. Receptor Induced Pathways "
	print "\t\t 3. RNA Gene Expression "
	print "\t\t 4. Grouped Pathways "
	print "\t\t 5. Paired-wise Pathways Proximity "
	print "\t\t 6. Protein-Protein Interaction Correlation "
	print "\t\t 7. Co-expression Analysis "
	print "\t\t 8. Make HiC interactions Gene Sets "
	print "\t\t 9. Inter-Pathways Proximity "
	print "\t\t 10. Pathways Nodes and Edges Genes Capturing "
	print "\t\t 11. NETPATH Edge Fraction "
	print "\t\t 12. Calculate GC Content of Pathways "
		
	corr_type = 0
	corr_term = ""
	sub_path_choice = 0
	path_choice = int(raw_input("\t Choose Only One Number from above :: "))
	
	if path_choice == 1:
		print "\t\t What type of Graph data you want to use in the KEGG analysis :: "
		print "\t\t\t 1. Undirected "
		print "\t\t\t 2. Directed "
		sub_path_choice = int(raw_input("\t\t\t Choose Only One Number from above :: "))
		
	if path_choice == 3:
		print "\t\t What factor of Expression function you want to analyse :: "
		print "\t\t\t 1. Fractional Expression [i.e., Sum(Max(RPKM(in-Pathways)))/Sum(Max(RPKM(All-HiC Genes)))] "
		print "\t\t\t 2. Only HiC interacting Genes only [i.e., Sum(Max(RPKM(in-Pathways)))] "
		print "\t\t\t 3. All interacting and non-interacting Genes "
		rna_func_choice = int(raw_input("\t\t\t Choose Only One Number from above :: "))
		
		print "\t\t What Genes of Pathways you want to find Level of Expressions :: "
		print "\t\t\t 1. KEGG "
		print "\t\t\t 2. Receptor Induced Pathways "
		sub_path_choice = int(raw_input("\t\t\t Choose Only One Number from above :: "))
		
		if sub_path_choice == 1:
			print "\t\t What type of Graph data you want to use in the KEGG analysis :: "
			print "\t\t\t 1. Undirected "
			print "\t\t\t 2. Directed "
			corr_type = int(raw_input("\t\t\t Choose Only One Number from above :: "))

	else:
		rna_func_choice = int(0)
	
	if path_choice == 6:
		print "\t\t\t What Protein Protein Interaction you want to select for the Interaction Analysis ?? See the options below :: "
		print "\t\t\t\t 1. STRING All Interaction "
		print "\t\t\t\t 2. HPRD Only Physical & Validated Interactions "
		print "\t\t\t\t 3. HPRD and STRING Combined Physical & Validated Interactions "
		protein_db = int(raw_input("\t\t\t\t Choose Only One Number from above :: "))
		
		print "\t\t What Genes of Pathways you want to find PPI :: "
		print "\t\t\t 1. KEGG "
		print "\t\t\t 2. Receptor Induced Pathways "
		sub_path_choice = int(raw_input("\t\t\t Choose Only One Number from above :: "))
	else:
		protein_db = int(0)
		
		
	if path_choice == 7:
		print "\t\t\t What type of Correlated Genes you want to analyse for the Co-expression :: "
		print "\t\t\t\t 1. Positive Correlated "
		print "\t\t\t\t 2. Negative Correlated "
		print "\t\t\t\t 3. Zero Correlated "
		print "\t\t\t\t 4. All types "
		corr_type = int(raw_input("\t\t\t\t Choose Only One Number from above :: "))
		corr_term = "Positive Co-expression"
		
		print "\t\t What Genes of Pathways you want to find PPI :: "
		print "\t\t\t 1. KEGG "
		print "\t\t\t 2. Receptor Induced Pathways "
		sub_path_choice = int(raw_input("\t\t\t Choose Only One Number from above :: "))
	
	if path_choice == 8:
		print "\t\t What Genes of Pathways you want to Create the Gene Sets :: "
		print "\t\t\t 1. KEGG "
		print "\t\t\t 2. Receptor Induced Pathways "
		sub_path_choice = int(raw_input("\t\t\t Choose Only One Number from above :: "))
	else:
		#corr_type = int(0)
		corr_term = "ERROR"
	
	if path_choice == 9:
		print "\t\t What Genes of Pathways you want to Create the Gene Sets :: "
		print "\t\t\t 1. KEGG "
		print "\t\t\t 2. Receptor Induced Pathways "
		print "\t\t\t 2. RNA Expression Inter Pathways "
		sub_path_choice = int(raw_input("\t\t\t Choose Only One Number from above :: "))
	else:
		#corr_type = int(0)
		corr_term = "ERROR"	
		
	if path_choice == 10:
		print "\t\t What Genes of Pathways you want to Create the Gene Sets :: "
		print "\t\t\t 1. KEGG "
		print "\t\t\t 2. Receptor Induced Pathways "
		sub_path_choice = int(raw_input("\t\t\t Choose Only One Number from above :: "))
		
		if sub_path_choice == 1:
			print "\t\t What type of Graph data you want to use in the KEGG analysis :: "
			print "\t\t\t 1. Undirected "
			print "\t\t\t 2. Directed "
			corr_type = int(raw_input("\t\t\t Choose Only One Number from above :: "))
		
	else:
		#corr_type = int(0)
		corr_term = "ERROR"		

	
	print "\t Do you want to select the random Genes from the same Chromosome of the Pathways Genes ?? :: "
	print "\t\t 1. YES "
	print "\t\t 2. NO "
	chr_control = raw_input("\t Choose Only One Number from above :: ")

	print "\t Do you want to do this analysis for Positive Control OR Negative Control  OR Grouped Pathways ?? :: "
	print "\t\t 1. Positive Control "
	print "\t\t 2. Negative Control "
	exp_control = raw_input("\t Choose Only One Number from above :: ")

	analysis_type = int(analysis_type)
	resolution_type = int(resolution_type)
	path_choice = int(path_choice)
	chr_control = int(chr_control)
	exp_control = int(exp_control)
	iterval = int(iterval)
	
	if iterval == 1:
		siter = 100
	elif iterval == 2:
		siter = 500
	elif iterval == 3:
		siter = 1000
	else:
		print "Please select correct iteration number. We perform only three iteration analysis !!! "
		exit()
		
	
	if resolution_type == 1:
		resolution = "10"
		#odir = "/cbcb/project2-scratch/hiren/HiC_Analysis/HiC_Results/COMPARATIVE/10_KB_RESULTS/FINAL_10KB/"	
		odir = "/cbcb/project2-scratch/hiren/HiC_Analysis/HiC_Results/COMPARATIVE/10_KB_RESULTS/REVISIONS/10KB_REVISION/"
	elif resolution_type == 2:
		resolution = "100"
		if gcval == 1:
			if path_choice == 1 and sub_path_choice == 2:
				odir = "/cbcb/project2-scratch/hiren/HiC_Analysis/HiC_Results/COMPARATIVE/10_KB_RESULTS/REVISIONS/100KB_REVISION/"
			elif path_choice == 3 and sub_path_choice == 1 and corr_type == 2:
				odir = "/cbcb/project2-scratch/hiren/HiC_Analysis/HiC_Results/COMPARATIVE/10_KB_RESULTS/REVISIONS/100KB_REVISION/"	
			elif path_choice == 10 and sub_path_choice == 1 and corr_type == 2:
				odir = "/cbcb/project2-scratch/hiren/HiC_Analysis/HiC_Results/COMPARATIVE/10_KB_RESULTS/REVISIONS/100KB_REVISION/"
			else:
				odir = "/cbcb/project2-scratch/hiren/HiC_Analysis/HiC_Results/COMPARATIVE/10_KB_RESULTS/REVISIONS/100KB_REVISION/"
		elif gcval == 2:
			odir = "/cbcb/project2-scratch/hiren/HiC_Analysis/HiC_Results/COMPARATIVE/10_KB_RESULTS/REVISIONS/GC_COMPARISON/"
		
	elif resolution_type == 3:
		resolution = "10_100"
		odir = "/cbcb/project2-scratch/hiren/HiC_Analysis/HiC_Results/COMPARATIVE/10_KB_RESULTS/REVISIONS/10KB_100KB_REVISION/"
	else:
		print "Your selection of resolution is not in our analysed data..... Please use either 1 OR 2 OR 3 options."
		print "========== Run again the script and choose appropriate option .... THANK YOU !!!! ================="
		sys.exit()
	
	if path_choice == 1:	
		pn = "KEGG"
	elif path_choice == 2:
		pn = "RI"
	elif path_choice == 3:
		if rna_func_choice == 1:
			pn = "RNA"
		if rna_func_choice == 2:
			pn = "RNA" + "_PEDGES"
		if rna_func_choice == 3:
			pn = "RNA" + "_PNODES"
	elif path_choice == 4:
		pn = "GROUPED"
	elif path_choice == 5:
		pn = "PAIRED_COMMON"
	elif path_choice == 6:
		pn = "PROT_PROT"
	elif path_choice == 7:
		pn = "CO-EXPRESSED"
	elif path_choice == 8:
		pn = "INTERACTING_GENES"
		odir = odir + "ALL_INTERACTIONS/"
	elif path_choice == 9:
		pn = "INTER-PATHWAYS"
	elif path_choice == 10:
		pn = "GENES-ENTRY"
	elif path_choice == 11:
		pn = "NETPATH"	
		
	if chr_control == 1:
		ch = "CHR_LEN"
	elif chr_control == 2:
		ch = "LEN"
	
	if exp_control == 1:
		ex = "POSITIVE"
		iter = 500
	elif exp_control == 2:
		ex = "NEGATIVE"
		iter = 100
			
	if path_choice == 8:
		ofl_back = "_" + pn + "_" + ch + "_" + ex + "_Controlled.txt"
		if fdrval > 0:
			ofl_back = "_" + str(fdrval)+ "_" + pn + "_" + ch + "_" + ex + "_Controlled.txt"
	else:	
		ofl_back = "_" + pn + "_EDGE_Fraction_With_Zscore_" + ch + "_" + ex + "_Controlled" + "_" + str(siter) + "_temp.txt"
		if fdrval > 0:
			ofl_back = "_UPDATED_" + str(fdrval)+ "_" + pn + "_EDGE_Fraction_With_Zscore_" + ch + "_" + ex + "_Controlled" + "_" + str(siter) + "_temp.txt"
	
	
	"""
		if chr_control == 1:
			ofl_back = "_KEGG_EDGE_Fraction_With_Zscore_CHR_LEN_Controlled.txt"
		elif chr_control == 2:
			ofl_back = "_KEGG_EDGE_Fraction_With_Zscore_LEN_Controlled.txt"
	elif path_choice == 2:
		if chr_control == 1:
			ofl_back = "_RI_EDGE_Fraction_With_Zscore_CHR_LEN_Controlled.txt"
		elif chr_control == 2:
			ofl_back = "_RI_EDGE_Fraction_With_Zscore_LEN_Controlled.txt"			
	"""
	return analysis_type, resolution_type, path_choice, sub_path_choice, corr_type, corr_term, rna_func_choice, chr_control, exp_control, resolution, odir, ofl_back, iter, protein_db, siter, gcval, fdrval
	
global iter, Total_Genes, TOTAL_GENES, TOTAL_TRANSCRIPT, GEN_LEN, CHR_GEN, GEN_MID, GEN_CHR, TOT_G, TOT_L, TOT_C, G_DATA, T_DATA, RGENES, RFPKM, analysis_type, resolution_type, path_choice, sub_path_choice, corr_type, corr_term, chr_control, exp_control, rna_func_choice, protein_db, GENOME_GENE, siter, gcval, NC_TR, NC_GN, fdrval

#TISSUES = ["GSM1340638_GSM1340637", "GSM1081530_GSM1081531", "GSM1055800_GSM1055801", "GSM862723_GSM892306", "GSM1340639_GSM1340640", "GSM927076_GSM927076"]
#T_NAMES = ["BT483", "HEK293", "IMR90", "hESC", "GM06990", "RWPE1"]


TISSUES = ["GSM1340638_GSM1340637", "GSM1081530_GSM1081531", "GSM1055800_GSM1055801", "GSM862723_GSM892306", "GSM1340639_GSM1340640", "GSM927076_GSM927076"]
T_NAMES = ["BT483", "HEK293", "IMR90", "hESC", "GM06990", "RWPE1"]
#R_NAMES = ["HEK293", "IMR90", "RWPE1"]
#R_TISSUES = ["GSM1081530_GSM1081531", "GSM1055800_GSM1055801", "GSM927076_GSM927076"]
R_TISSUES = ["GSM1340638_GSM1340637", "GSM1081530_GSM1081531", "GSM1055800_GSM1055801", "GSM862723_GSM892306", "GSM1340639_GSM1340640", "GSM927076_GSM927076"]
R_NAMES = ["BT483", "HEK293", "IMR90", "hESC", "GM06990", "RWPE1"]

"""
TISSUES = ["GSM927076_GSM927076"]
T_NAMES = ["RWPE1"]
"""

(analysis_type, resolution_type, path_choice, sub_path_choice, corr_type, corr_term, rna_func_choice, chr_control, exp_control, resolution, odir, ofl_back, iter, protein_db, siter, gcval, fdrval) = getParameters()

if path_choice == 3:
	#rdir = "/cbcb/project2-scratch/hiren/HiC_Analysis/HiC_Results/RNA_Results/"
	rdir = "/cbcb/project2-scratch/hiren/Justin/Chromatin_Project/RNA_Seq/"
	T_NAMES = R_NAMES
	TISSUES = R_TISSUES

if gcval == 2:
	TISSUES = ["GSM1081530_GSM1081531"]
	T_NAMES = ["HEK293"]
	R_NAMES = ["HEK293"]	
	
(G_DATA, T_DATA, GEN_LEN, CHR_GEN, GEN_MID, GEN_CHR, TOT_C, TOT_G, TOT_L, RGENES, RFPKM, NC_TR, NC_GN) = getGlobalArrays()
#Total_Genes = getGene_Genome()
(G_DATA, T_DATA, GEN_LEN, CHR_GEN, GEN_MID, GEN_CHR, TOT_C, TOT_G, TOT_L, NC_TR, NC_GN) = getGeneData()

GENOME_GENE = np.unique(TOT_G)
print "\t TOTAL UNIQUE GENOMIC GENES ( With Transcripts ) :: ", GENOME_GENE.size

if path_choice == 6:
	PIE = getProteinData()
	if protein_db == 1:
		print "\t TOTAL NUMBER OF INTERACTIONS IN THE STRING DATABASE 				:: ", len(PIE.edges())
		print "\t TOTAL NUMBER OF UNIQUE GENES IN THE STRING DATABASE 				:: ", np.unique(np.array(PIE.nodes())).size
	if protein_db == 2:
		print "\t TOTAL NUMBER OF INTERACTIONS IN THE HPRD DATABASE 				:: ", len(PIE.edges())
		print "\t TOTAL NUMBER OF UNIQUE GENES IN THE HPRD DATABASE 				:: ", np.unique(np.array(PIE.nodes())).size
	if protein_db == 3:
		print "\t TOTAL NUMBER OF INTERACTIONS IN BOTH HPRD & STRING DATABASE 				:: ", len(PIE.edges())
		print "\t TOTAL NUMBER OF UNIQUE GENES IN BOTH HPRD & STRING DATABASE 				:: ", np.unique(np.array(PIE.nodes())).size
		
	"""	
	else:
		print "===>", protein_db
		print "\t YOU HAVE NOT CHOSEN AN APPROPRIATE DATABASE OPTION :-( THANK YOU !! BYE "
		sys.exit()
	"""	

for i in range(len(TISSUES)):
	(s1, s2) = TISSUES[i].split("_")
	RNA_ANALYSIS = False
	t_name = T_NAMES[i]
	
	if path_choice == 7:
		if t_name not in R_NAMES:
			continue
	
	
	(oterm, path_dir) = choosePathwayType()
	
	#ofl = odir+t_name+"_"+oterm+"_EDGE_Fraction_With_Zscore.txt"
	#ofl = odir+t_name+"_"+oterm+ofl_back
	
	
	ofl = odir + t_name + ofl_back
	#print ofl
	
	if path_choice == 1 or path_choice == 2 or path_choice == 3 or path_choice == 4 or path_choice == 5 or path_choice == 9:
		Finished_Pathway = getFinishedPath(ofl)
	else:
		Finished_Pathway = np.array([], dtype='S')
	
	OFL = checkFileMode(ofl, Finished_Pathway)
	
	print t_name, "\t analysis is Starting !!! "
	(g1, g2, GD, ID, TOTAL_GENES) = getTotalData(s1, s2, t_name, odir, resolution)
	
	#(TOTAL_GENES, CHR_GEN, GEN_LEN, GEN_MID, GEN_CHR, TOT_G, TOT_L) = getGeneLength(TOTAL_GENES)
	print "\t 3. UNIFIED SET OF GENES ARE FILTERED WITH TOTAL ANNOTATED GENES !!! "
	STOT_C, STOT_L, STOT_G = getFilteredData(TOTAL_GENES, TOT_C, TOT_G, TOT_L)
	
	if path_choice == 3 or sub_path_choice == 3:
		if t_name in R_NAMES:
			(RGENES, RFPKM) = getExpressionData(rdir, t_name)
			print "\t 4. GENES IN THE SAMPLE THAT HAS RNA SEQ DATA :: ", RGENES.size
			RNA_ANALYSIS = True
			
	
	(OFL, IP, NIP, OP, IG, TG, IP_G, NIP_G, OP_G, IG_G, TG_G, INTER_PATHWAYS, INTRA_PATHWAYS_INTERACTIONS) = getPathwaysGenes(GD, OFL, STOT_C, STOT_L, STOT_G, Finished_Pathway, path_dir, RNA_ANALYSIS, t_name, siter)
	
	
	if path_choice == 8:
		OFL = createInteractingGeneSets(IP, IP_G, OP, OP_G, OFL)
	
	if path_choice == 6:
		OFL = makeProtProtReport(PIE, IP, NIP, OP, IG, TG, OFL)
	
	if path_choice == 7:
		GCE = getOneCoExpressionData(t_name)
		OFL = makeCorrelationReport(GCE, IP, NIP, OP, IG, TG, OFL, corr_term)
	
	if path_choice == 9:
		OFL = getInterPathwaysData(t_name, INTER_PATHWAYS, INTRA_PATHWAYS_INTERACTIONS, GD, OFL, STOT_C, STOT_L, STOT_G, RNA_ANALYSIS, path_dir, Finished_Pathway)
	
	OFL.close()
	
	print "Congratulations !!! ", s1, "\t", s2, "\t is Completed !!! "
	print "See your output in the file :: ", ofl
	
	
	