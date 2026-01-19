
# Script to identify event types for DEAD-dependent differentially expressed genes
# Created by Maximilian Sack, Bioinformatics Group, Department of Computer Science and Interdisciplinary Centre for Bioinformatics, Leipzig University
# Part of the Study:
# A structured RNA balances DEAD-box RNA helicase function in plant alternative splicing control (2026)
# Rica Burgardt, Julia Bauer, Maren Reinhardt, Natalie Rupp, Christoph Engel, Lukas Hellmann, Maximilian Sack, Zasha Weinberg, and Andreas Wachter
# Cite that paper if you use this script

#Example command line
#python3 findSpliceEvents.py 
#	AtRTDv2_QUASI_19April2016.gtf 
#	DTU_transcripts_cpm4_sample1_3h_6h_all_transcripts_3h.tab 
#	--output_all SpliceEvents_all.tab 
#	--output_regulated SpliceEvents_regulated.tab 
#	--output_reciprocal SpliceEvents_reciprocal.tab

#Differential gene expression table:
#AT1G01260_P1,X03h.EST-X03h.mock,0.00395097538878148,0.171207076302256,0.369400232338061,up-regulated
#AT1G01260_P3,X03h.EST-X03h.mock,0.00395097538878148,-0.171207076302256,-1.19466698377057,down-regulated

#Transcript table from AtRTD:
#Chr1    Araport11       exon    108946  111699  .       +       .       transcript_id "AT1G01260_P1"; gene_id "AT1G01260"; Note "basic";
#Chr1    Araport11       exon    199527  199763  .       +       .       transcript_id "AT1G01550_P1"; gene_id "AT1G01550"; Note "BPS1-like";
#Chr1    Araport11       exon    199890  199959  .       +       .       transcript_id "AT1G01550_P1"; gene_id "AT1G01550"; Note "BPS1-like";
#Chr1    Araport11       exon    200511  201775  .       +       .       transcript_id "AT1G01550_P1"; gene_id "AT1G01550"; Note "BPS1-like";

#Chr1    Araport11       exon    108946  109330  .       +       .       transcript_id "AT1G01260_P3"; gene_id "AT1G01260"; Note "basic";
#Chr1    Araport11       exon    109406  111699  .       +       .       transcript_id "AT1G01260_P3"; gene_id "AT1G01260"; Note "basic";

#Basic method:
#read csv, get pairs and up-/down-regulated	#Gene -> ((ID,up),(ID,down))
#read gtf, get lists of exons (start,end) and strand 
#for each pair: compare the exons
#	discard identical
#	compare the ones that are not
#	only one remaining for one transcript: CE (retention/skipping)
#	two exons and one exon, with same boundaries: IR (skipping/retention)
#	one exon each, different start/end: alt5p/alt3p, also dependent on strand


import argparse

upString = "up-regulated"
downString = "down-regulated"
regCon = (upString,downString)
fStrandS = "+"
rStrandS = "-"

ess = "exon_skip"
irs = "intron_ret"
a5s = "alt_5prime"
a3s = "alt_3prime"
mxes = "MXE"
ses = "SE"
eventIDs = {ess:0,irs:1,a5s:2,a3s:3,mxes:4,ses:5}

def getStrand(strandString):
	if strandString == fStrandS: return True
	elif strandString == rStrandS: return False
	else: return None

def loadTranscripts(fileName):
	geneDict = dict()
	with open(fileName,"r") as gtf_reader:
		for i,line in enumerate(gtf_reader):
			lsp = line.split("\t")
			if lsp [2] != "exon":
				print("Not an exon: "+lsp[2]+" at line "+str(i))
				continue
			
			eStart = int(lsp[3])
			eEnd = int(lsp[4])
			
			fStrand = getStrand(lsp[6])
			if fStrand is None:
				print("Error, unknown strand-string: "+strandString+" at line "+str(i))
				continue
			
			infos = lsp[8].split(";")	#transcript_id "AT1G01550_P1"; gene_id "AT1G01550"; Note "BPS1-like";
			transcriptID = None
			geneID = None
			for info in infos:
				if info.strip().startswith("transcript_id"):
					transcriptID = info.split("\"")[1]
				elif info.strip().startswith("gene_id"):
					geneID = info.split("\"")[1]
				else:
					pass
			
			if transcriptID is None or geneID is None:
				print("Error, cant find gene and transcript IDs: "+infos+" at line "+str(i))
				continue
			
			
			if geneID not in geneDict:geneDict[geneID] = [fStrand,list()]
			
			transcriptPos = -1
			for i in range(len(geneDict[geneID][1])):
				if geneDict[geneID][1][i][0] == transcriptID:
					transcriptPos = i
					break
			if transcriptPos == -1:
				transcriptPos = len(geneDict[geneID][1])
				geneDict[geneID][1].append([transcriptID,list()])
			
			geneDict[geneID][1][transcriptPos][1].append((eStart,eEnd))
	return geneDict

def loadRegulatedTranscripts(regulatedTable):
	transcriptRegDict = dict()	#dict of transcriptID -> reg-type
	with open (regulatedTable,"r") as table_reader:
		for i,line in enumerate(table_reader):
			lsp = line.strip().split("\t")
			if lsp[0]=="target":continue	#header
			transcriptID = lsp[0]
			upReg = None
			if lsp[-1] == upString:
				upReg = 0
			elif lsp[-1] == downString:
				upReg = 1
			else:
				print("Error, unknown reg-string: "+lsp[-1]+" at line "+str(i))
				continue
			transcriptRegDict[transcriptID] = upReg	#TODO
	return transcriptRegDict

def writeOutput():
	with open(sys.argv[3],"w") as out_writer:
		unknownEvents=0
		writtenEvents=0
		missingTranscripts=0
		#Gene	eventType	up-reg-ID	up-reg-type	down-reg-ID	down-reg-type
		out_writer.write("gene\teventType\tup-reg-ID\tup-reg-eventType\tdown-reg-ID\tdown-reg-eventType")
		for gene,(uReg,dReg,fStrand) in geneDict.items():
			if uReg is None or dReg is None:
				print("Error, transcript is missing for "+gene)
				missingTranscripts+=1
				continue
			
			eventType = None
			uRegType = None
			dRegType = None
			#only one remaining for one transcript: CE (retention/skipping)
			#two exons and one exon, with same boundaries: IR (skipping/retention)
			#one exon each, different start/end: alt5p/alt3p, also dependent on strand
			uRegExonSet = set(uReg[1])
			dRegExonSet = set(dReg[1])
			for exon in uReg[1]:dRegExonSet.discard(exon)
			for exon in dReg[1]:uRegExonSet.discard(exon)
			
			if len(uRegExonSet) + len(dRegExonSet) == 1:		#ES
				eventType = "exon_skip"
				if len(uRegExonSet) == 0 and len(dRegExonSet) == 1:
					uRegType = "exon_skip"
					dRegType = "exon_retention"
				elif len(uRegExonSet) == 1 and len(dRegExonSet) == 0:
					uRegType = "exon_retention"
					dRegType = "exon_skip"
				else:
					print("Error, exon skip event type could not be identified: "+gene+"\n"+str(uRegExonSet)+" "+str(dRegExonSet))
					eventType = "Unknown"
					uRegType = "Unknown"
					dRegType = "Unknown"
					unknownEvents+=1
					
			elif len(uRegExonSet) + len(dRegExonSet) == 3:		#IR
				eventType = "intron_ret"
				if len(uRegExonSet) == 1 and len(dRegExonSet) == 2:
					irEx = list(uRegExonSet)[0]
					isEx1,isEx2 = sorted(list(dRegExonSet))
					if irEx[0]!= isEx1[0] or irEx[1] != isEx2[1]:	#Compare coordinates to validate event
						print("Error, intron retention event could not be validated: "+gene+"\n"+str(uRegExonSet)+" "+str(dRegExonSet))
						eventType = "Unknown"
						uRegType = "Unknown"
						dRegType = "Unknown"
						unknownEvents+=1
					uRegType = "intron_retention"
					dRegType = "intron_skip"
				elif len(uRegExonSet) == 2 and len(dRegExonSet) == 1:
					irEx = list(dRegExonSet)[0]
					isEx1,isEx2 = sorted(list(uRegExonSet))
					if irEx[0]!= isEx1[0] or irEx[1] != isEx2[1]:
						print("Error, intron retention event could not be validated: "+gene+"\n"+str(uRegExonSet)+" "+str(dRegExonSet))
						eventType = "Unknown"
						uRegType = "Unknown"
						dRegType = "Unknown"
						unknownEvents+=1
					uRegType = "intron_skip"
					dRegType = "intron_retention"
				else:
					print("Error, intron retention event type could not be identified: "+gene+"\n"+str(uRegExonSet)+" "+str(dRegExonSet))
					eventType = "Unknown"
					uRegType = "Unknown"
					dRegType = "Unknown"
					unknownEvents+=1
				
			elif len(uRegExonSet) == 1 and len(dRegExonSet) == 1:	#alt
				uEx = list(uRegExonSet)[0]
				dEx = list(dRegExonSet)[0]
				
				#strand=-, [1]==[1] => alt5p
				#	[0]<[0] => 0 is 3p, 1 is 5p
				#strand=+, [0]==[0] => alt5p
				#	[1]<[1] => 0 is 5p, 1 is 3p
				
				if fStrand:	#is on the forward strand
					if uEx[0] == dEx[0]:
						eventType = "alt_5prime"
						if uEx[1] < dEx[1]:
							uRegType = "alt_5prime_5"
							dRegType = "alt_5prime_3"
						else:
							uRegType = "alt_5prime_3"
							dRegType = "alt_5prime_5"
					elif uEx[1] == dEx[1]:
						eventType = "alt_3prime"
						if uEx[0] < dEx[0]:
							uRegType = "alt_3prime_5"
							dRegType = "alt_3prime_3"
						else:
							uRegType = "alt_3prime_3"
							dRegType = "alt_3prime_5"
					else:
						print("Error, alt prime event type could not be identified: "+gene+"\n"+str(uRegExonSet)+" "+str(dRegExonSet))
						eventType = "Unknown"
						uRegType = "Unknown"
						dRegType = "Unknown"
						unknownEvents+=1
						#continue
				else:
					if uEx[0] == dEx[0]:
						eventType = "alt_3prime"
						if uEx[1] < dEx[1]:
							uRegType = "alt_3prime_3"
							dRegType = "alt_3prime_5"
						else:
							uRegType = "alt_3prime_5"
							dRegType = "alt_3prime_3"
					elif uEx[1] == dEx[1]:
						eventType = "alt_5prime"
						if uEx[0] < dEx[0]:
							uRegType = "alt_5prime_3"
							dRegType = "alt_5prime_5"
						else:
							uRegType = "alt_5prime_5"
							dRegType = "alt_5prime_3"
					else:
						print("Error, alt prime event type could not be identified: "+gene+"\n"+str(uRegExonSet)+" "+str(dRegExonSet))
						eventType = "Unknown"
						uRegType = "Unknown"
						dRegType = "Unknown"
						unknownEvents+=1
						#continue
			else:
				print("Error, splice event type could not be identified: "+gene+" "+str(len(uRegExonSet))+" "+str(len(dRegExonSet)))
				eventType = "Unknown"
				uRegType = "Unknown"
				dRegType = "Unknown"
				unknownEvents+=1
				#continue
			
			out_writer.write("\n"+gene+"\t"+eventType+"\t"+uReg[0]+"\t"+uRegType+"\t"+dReg[0]+"\t"+dRegType)
			writtenEvents+=1

def getEvent(i_exonList,j_exonList,fStrand):
	#only one remaining for one transcript: CE (retention/skipping)
	#two exons and one exon, with same boundaries: IR (skipping/retention)
	#one exon each, different start/end: alt5p/alt3p, also dependent on strand
	iRegExonSet = set(i_exonList)
	jRegExonSet = set(j_exonList)
	for exon in i_exonList:jRegExonSet.discard(exon)
	for exon in j_exonList:iRegExonSet.discard(exon)
	
	if len(iRegExonSet) + len(jRegExonSet) == 1:		#ES
		foundEventType = "exon_skip"
		if len(iRegExonSet) == 0 and len(jRegExonSet) == 1:
			iType = "exon_skip"
			jType = "exon_retention"
			foundEvent = list(jRegExonSet)[0]
		elif len(iRegExonSet) == 1 and len(jRegExonSet) == 0:
			iType = "exon_retention"
			jType = "exon_skip"
			foundEvent = list(iRegExonSet)[0]
		else:
			return None,None,None,None
		return foundEventType,foundEvent,iType,jType
			
	elif len(iRegExonSet) + len(jRegExonSet) == 3:		#IR
		foundEventType = "intron_ret"
		if len(iRegExonSet) == 1 and len(jRegExonSet) == 2:
			irEx = list(iRegExonSet)[0]
			isEx1,isEx2 = sorted(list(jRegExonSet))
			if irEx[0]!= isEx1[0] or irEx[1] != isEx2[1]:
				return None,None,None,None
			iType = "intron_retention"
			jType = "intron_skip"
			foundEvent = (irEx,(isEx1,isEx2))
		elif len(iRegExonSet) == 2 and len(jRegExonSet) == 1:
			irEx = list(jRegExonSet)[0]
			isEx1,isEx2 = sorted(list(iRegExonSet))
			if irEx[0]!= isEx1[0] or irEx[1] != isEx2[1]:
				return None,None,None,None
			iType = "intron_skip"
			jType = "intron_retention"
			foundEvent = (irEx,(isEx1,isEx2))
		else:
			return None,None,None,None
		return foundEventType,foundEvent,iType,jType
		
	elif len(iRegExonSet) == 1 and len(jRegExonSet) == 1:	#alt	#TODO or MXE or shifted exon
		iEx = list(iRegExonSet)[0]
		jEx = list(jRegExonSet)[0]
		foundEvent = (iEx,jEx)
		
		if iEx[1] < jEx[0] or iEx[0] > jEx[1]:	#MXE
			return "MXE",foundEvent,"MXE","MXE"
		elif iEx[0]==jEx[0] or iEx[1]==jEx[1]:	#alt
			#strand=-, [1]==[1] => alt5p
			#	[0]<[0] => 0 is 3p, 1 is 5p
			#strand=+, [0]==[0] => alt5p
			#	[1]<[1] => 0 is 5p, 1 is 3p
			
			if fStrand:	#is on the forward strand
				if iEx[0] == jEx[0]:
					foundEventType = "alt_5prime"
					if iEx[1] < jEx[1]:
						iType = "alt_5prime_5"
						jType = "alt_5prime_3"
					else:
						iType = "alt_5prime_3"
						jType = "alt_5prime_5"
				elif iEx[1] == jEx[1]:
					foundEventType = "alt_3prime"
					if iEx[0] < jEx[0]:
						iType = "alt_3prime_5"
						jType = "alt_3prime_3"
					else:
						iType = "alt_3prime_3"
						jType = "alt_3prime_5"
				else:
					return None,None,None,None
				return foundEventType,foundEvent,iType,jType
			else:
				if iEx[0] == jEx[0]:
					foundEventType = "alt_3prime"
					if iEx[1] < jEx[1]:
						iType = "alt_3prime_3"
						jType = "alt_3prime_5"
					else:
						iType = "alt_3prime_5"
						jType = "alt_3prime_3"
				elif iEx[1] == jEx[1]:
					foundEventType = "alt_5prime"
					if iEx[0] < jEx[0]:
						iType = "alt_5prime_3"
						jType = "alt_5prime_5"
					else:
						iType = "alt_5prime_5"
						jType = "alt_5prime_3"
				else:
					return None,None,None,None
				return foundEventType,foundEvent,iType,jType
				
		else:	#shifted
			return ses,foundEvent,ses,ses
	else:
		return None,None,None,None	#not clearly identifyable!

def evaluateGene(gene,fStrand,transcripts,trancriptReg):
	geneLine = list()	#GeneID	nTranscripts	ES-fraction	IR-fraction	alt5-fraction	alt3-fraction	MXE-fraction	shiftedExon-fraction	Events
	geneLine.append(gene)
	geneLine.append(str(len(transcripts)))
	eventCounts = [0]*len(eventIDs)
	foundEventTypes = list()	#list of all clear types
	foundEvents = set()		#set of all events, to prevent duplicates
	regSet = set()
	regIDs = set()
	reciprocalList = list()
	reciprocalIDs = set()
	
	for i in range(len(transcripts)-1):
		i_transcriptID,i_exonList = transcripts[i]
		for j in range(i+1,len(transcripts)):
			#print(i,j)
			j_transcriptID,j_exonList = transcripts[j]
			
			foundEventType,foundEvent,iType,jType = getEvent(i_exonList,j_exonList,fStrand)
			
			if i_transcriptID in trancriptReg:
				regSet.add((gene,i_transcriptID,j_transcriptID,regCon[trancriptReg[i_transcriptID]],str(foundEventType),str(iType)))
				regIDs.add(i_transcriptID)
			if j_transcriptID in trancriptReg:
				regSet.add((gene,j_transcriptID,i_transcriptID,regCon[trancriptReg[j_transcriptID]],str(foundEventType),str(jType)))
				regIDs.add(j_transcriptID)
			
			if i_transcriptID in trancriptReg and j_transcriptID in trancriptReg:
				if trancriptReg[i_transcriptID] != trancriptReg[j_transcriptID]:
					#Gene	eventType	up-reg-ID	up-reg-type	down-reg-ID	down-reg-type
					if trancriptReg[i_transcriptID] == 0 and trancriptReg[j_transcriptID] == 1:
						upRegID     = i_transcriptID
						upRegType   = iType
						downRegID   = j_transcriptID
						downRegType = jType
					elif trancriptReg[i_transcriptID] == 1 and trancriptReg[j_transcriptID] == 0:
						upRegID     = j_transcriptID
						upRegType   = jType
						downRegID   = i_transcriptID
						downRegType = iType
					
					reciprocalList.append([None,gene,str(foundEventType),upRegID,str(upRegType),downRegID,str(downRegType)])
					reciprocalIDs.add(i_transcriptID)
					reciprocalIDs.add(j_transcriptID)
			
			if foundEventType is None:	#either none found (impossible) or not crealy identifyable, ignore in that case 
				continue
			if foundEvent in foundEvents:	#event has already been recorded
				continue
			#print(foundEventType)
			foundEvents.add(foundEvent)
			foundEventTypes.append(foundEventType)
			eventCounts[eventIDs[foundEventType]]+=1
	
	#print(eventCounts)
	geneLine.extend([str(round(e/len(foundEvents) if len(foundEvents)>0 else 0,4)) for e in eventCounts])
	geneLine.append(",".join(foundEventTypes))
	
	regList = list()
	for l in regSet:
		tmp = [str(len(regIDs))]
		tmp.extend(list(l))
		regList.append(tmp)
	
	for l in reciprocalList:l[0] = str(len(reciprocalIDs))
	return geneLine,regList,reciprocalList

def evaluateGenes(geneDict,trancriptReg	):
	fullEvaluation = list()
	regulatedEvaluation = list()
	reciprocalEvaluation = list()
	
	for gene,(fStrand,transcripts) in geneDict.items():
		
		geneLine,regList,reciprocalList = evaluateGene(gene,fStrand,transcripts,trancriptReg)
		fullEvaluation.append(geneLine)
		regulatedEvaluation.extend(regList)
		reciprocalEvaluation.extend(reciprocalList)
	
	return fullEvaluation,regulatedEvaluation,reciprocalEvaluation

def writeOutputAll(output_all,fullEvaluation):
	with open(output_all,"w") as out_writer:
		out_writer.write("GeneID\tnTranscripts\tES-fraction\tIR-fraction\talt5-fraction\talt3-fraction\tMXE-fraction\tshiftedExon-fraction\tEvents")
		for line in fullEvaluation:
			out_writer.write("\n"+"\t".join(line))

def writeOutputRegulated(output_regulated,regulatedEvaluation):
	with open(output_regulated,"w") as out_writer:
		out_writer.write("nRegTranscForGene\tgene\ttranscriptID\tpartner-transcriptID\tregulation\teventType\ttranscriptEventType")
		for line in regulatedEvaluation:
			out_writer.write("\n"+"\t".join(line))

def writeOutputReciprocal(output_reciprocal,reciprocalEvaluation):
	with open(output_reciprocal,"w") as out_writer:
		out_writer.write("nPairsForGene\tgene\teventType\tup-reg-ID\tup-reg-eventType\tdown-reg-ID\tdown-reg-eventType")
		for line in reciprocalEvaluation:
			out_writer.write("\n"+"\t".join(line))

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('transcriptTable', type=str, help="Table of transcript exons in .gtf format")
	parser.add_argument('regulatedTable', type=str, help="Table of transcriptIDs that are affected")
	parser.add_argument('--output_all', type=str, help="Output")
	parser.add_argument('--output_regulated', type=str, help="Output")
	parser.add_argument('--output_regulated_Breakdown', type=str, help="Output")
	parser.add_argument('--output_reciprocal', type=str, help="Output")
	args = parser.parse_args()
	
	geneDict = loadTranscripts(args.transcriptTable)
	
	trancriptReg = loadRegulatedTranscripts(args.regulatedTable)
	fullEvaluation,regulatedEvaluation,reciprocalEvaluation = evaluateGenes(geneDict,trancriptReg)
	
	if args.output_all is not None:
		writeOutputAll(args.output_all,fullEvaluation)
	if args.output_regulated is not None:
		writeOutputRegulated(args.output_regulated,regulatedEvaluation)
	if args.output_regulated_Breakdown is not None:
		writeOutputRegulatedBreakdown(args.output_regulated_Breakdown,regulatedEvaluation)
	if args.output_reciprocal is not None:
		writeOutputReciprocal(args.output_reciprocal,reciprocalEvaluation)










