#!/usr/bin/env python
################################################################
# Orthology mapping-
#
# (c) Hasiba Asma July 2020

################################################################


import os
import sys
import argparse
import pprint
import re
import csv


# This dictionary is created to save id mapping file of Spec_X. The format of this dictionary is something like this: [OFAS00001] = 'OFAS2:0001'
def idMap_dict(nameOfDict, file):
	with open(file, 'r') as fi:
		rows = (line.split('\t') for line in fi)
		nameOfDict = {row[0]: row[1].strip('\n\r') for row in rows}
	return (nameOfDict)


def mapping(fileName):
	specList = []
	with open(fileName, 'r') as file1:
		for line in file1:
			cols = line.split(' ')
			scaf = cols[1]
			specList.append(scaf)

	return specList


# write a fucnction that would take feature table as an input and create a dictionary outta it
def idMap_dict_FT(nameOfDict, file):
	with open(file, 'r') as fi:
		rows = (line.split('\t') for line in fi)
		nameOfDict = {row[10]: row[14] for row in rows}
	return (nameOfDict)


#write a function that would go through each item of species X and see if its the ortholog so like AA23189: hh

def getOrthologs(dict_spX,dict_spD,diction2):
	#gbidlocal=''
	#gbidXlocal=''
	finalDict={}
	for item in dict_spX: # for each gene-loc
		gbidlocal=''
		gbidXlocal=''
		gbidlocal = dict_spX[item]  # gbid=Xps
		c = 0
		# check if this XP in above dictionary, if so write others here
		for key in diction2:
			if gbidlocal in diction2[key]:
				c = 1
				gbidXlocal = diction2[key][gbidlocal]  # [0][XP]=NP's
				if not gbidXlocal:
					gbidXlocal = ''

				#### these would be one or more DMEL's NP ids ..time to break them and convert them to their gene symbols
				# check if one or more
				# if msultiple orthologs
				if len(gbidXlocal) == 1:
					for item2 in gbidXlocal:
						if item2 in dict_spD:  # check if this NP in feature table > if so > write symbol else no symbol found thus orig NP id
							gbidXlocal = dict_spD[item2]
						else:
							gbidXlocal = item2
				# print('FOUND item yay ',dict_spDid[item])
				elif len(gbidXlocal) > 1:
					oneI = []
					commas = 0
					for item3 in gbidXlocal:
						if item3 in dict_spD:
							if commas == 0:

								oneI.append(dict_spD[item3])
							else:
								if dict_spD[item3] not in oneI:
									oneI.append(dict_spD[item3])
									#oneI = oneI + ',' + dict_spD[item3]
							if commas == len(gbidXlocal) - 1:
								gbidXlocal = oneI
							commas += 1
						# if no symbol
						else:
							if commas == 0:
								oneI.append(dict_spD[item3])
								#oneI = oneI + item3
							else:
								if dict_spD[item3] not in oneI:
								#	oneI = oneI + ',' + item3
									oneI.append(dict_spD[item3])
							if commas == len(gbidXlocal) - 1:
								gbidXlocal = oneI
							commas += 1

				# print('Found it0')
				break
			elif c == 1:
				print(gbidlocal, 'not in dictionary--no drosophile ortholog')

		#print(item,gbidlocal,gbidXlocal)
		if gbidXlocal!='':
			if item not in finalDict:
				finalDict[item]=str(gbidXlocal).strip("[]").replace('"','').replace("'","").replace(" ","")

	return finalDict


# -----------------------------MAIN FUNCTION-------------------
def main():
	# global d1
	parser = argparse.ArgumentParser()
	parser.add_argument('-sp1id', '--sp1idmap', help='species 1 id map text file', required=True)
	parser.add_argument('-ft', '--featureTable', help='species DMEL id map feature table', required=True)
	parser.add_argument('-spec', '--speciesName', help='species name')

	parser.add_argument('-mD', '--mappingDMEL', help='mapping file of DMEL')
	parser.add_argument('-mX', '--mappingX', help='mapping file of spec X')
	parser.add_argument('-og', '--ogfile', help='result file from orthologer')
	args = parser.parse_args()
	# absolute path
	#	namesp1 = args.namesp1
	sp1idmap = os.path.abspath(args.sp1idmap)
	featureTable = os.path.abspath(args.featureTable)
	speciesName = (args.speciesName)
#	scrmshawOutputPath = os.path.abspath(scrmshawOutput)

	mappingDMEL = os.path.abspath(args.mappingDMEL)
	mappingX = os.path.abspath(args.mappingX)
	ogfile = os.path.abspath(args.ogfile)

	# mapping files
	listDMEL = mapping(mappingDMEL)
	listX = mapping(mappingX)
	diction = {}

	###creating dictionary fron mapping ids
	# result file
	with open(ogfile, 'r') as infile:
		for line in infile:
			if not line.startswith('#'):  # == 'FALSE':
				# print(line)
				cols = line.split(' ')
				srNo = cols[0]
				chrID = cols[1]
				id = int(cols[2])

				if srNo not in diction:
					diction[srNo] = chrID
				else:
					diction[srNo] = diction[srNo] + ',' + chrID

	#	pprint.pprint(diction)
	# for i in range(0,int(srNo)):
	#	quit()

	diction2 = {}
	for key in diction:

		dmelIDs = []
		notDmelIDs = []
		# print(key)
		# print(diction[key])
		listOfItems = diction[key].split(',')
		numOfItems = len(listOfItems)
		# print('number of items in this list are', numOfItems)
		count = 0
		for item in listOfItems:

			if item in listDMEL:
				dmelIDs.append(item)
			# print('dmel item ',item)
			else:
				notDmelIDs.append(item)
			# print('not dmel ',item)

			# when it hits the last item in the list, add all the dmel ids items to non dmel ids items
			if count == (numOfItems - 1):
				# print(dmelIDs)
				# print(notDmelIDs)
				for eachNonDmelID in notDmelIDs:
					# print(eachNonDmelID)

					if key not in diction2:
						diction2[key] = {}
					if eachNonDmelID not in diction2[key]:
						# print(diction2)
						diction2[key][eachNonDmelID] = dmelIDs
					else:
						print('already exist, adding up')
			# 						diction2[eachNonDmelID]=diction2[eachNonDmelID]+','+dmelIDs
			count += 1

	#pprint.pprint(diction2)
	# creating a dictionary of SpecX using idmap file (format: AGAP : XP )
	dict_sp1id = idMap_dict('sp1id', sp1idmap)

	# creating a dictionary2 of SpecX using feature table file (format: NP: symb )
	dict_spDid = idMap_dict_FT('spDid', featureTable)

	dictOrthologs= getOrthologs(dict_sp1id,dict_spDid,diction2)

	#pprint.pp(dictOrthologs)
	#saving final dict to a text file
	file = open(speciesName+"_ortholog_dictfile.txt", "w")

	for key in dictOrthologs.keys():
		file.write(str(key) + "\t" + str(dictOrthologs[key]))
		file.write("\n")

	file.close()

	#pprint.pprint(dict_spDid)
	#pprint.pprint(dict_sp1id)
	# squit()
	#orthologOutput = 'SO_' + scrmshawOutput

	# # need to fix it at some point, exclude anything useless-- op stuff (June 22)
	# with open(scrmshawOutputPath, 'r') as so, open(orthologOutput, 'w') as fo:
	# 	for line in so:
	# 		print(line)
	# 		cols = line.split('\t')
	#
	# 		# Check col 5 and col 10
	# 		# To reduce calculation
	#
	# 		######			#if both right and left are similar
	# 		if cols[5] == cols[10]:
	# 			######				# if one gene------------
	# 			if cols[5].find(',') == -1:
	# 				gbid = ''
	# 				gbidX = ''
	# 				# ccheck to see if mapping of X species exist-- AGAP > XP
	# 				if cols[5] in dict_sp1id:# if gene-loc in gene-loc: XPs
	# 					gbid = dict_sp1id[cols[5]]# gbid=Xps
	# 					c = 0
	# 					# check if this XP in above dictionary, if so write others here
	# 					for key in diction2:
	# 						if gbid in diction2[key]:
	# 							c = 1
	# 							gbidX = diction2[key][gbid]  # [0][XP]=NP's
	# 							if not gbidX:
	# 								gbidX = ''
	#
	# 							#### these would be one or more DMEL's NP ids ..time to break them and convert them to their gene symbols
	# 							# check if one or more
	# 							# if multiple orthologs
	# 							if len(gbidX) == 1:
	# 								for item in gbidX:
	# 									if item in dict_spDid:  # check if this NP in feature table > if so > write symbol else no symbol found thus orig NP id
	# 										gbidX = dict_spDid[item]
	# 									else:
	# 										gbidX = item
	# 							# print('FOUND item yay ',dict_spDid[item])
	# 							elif len(gbidX) > 1:
	# 								oneI = ''
	# 								commas = 0
	# 								for item in gbidX:
	# 									if item in dict_spDid:
	# 										if commas == 0:
	# 											oneI = oneI + dict_spDid[item]
	# 										else:
	# 											oneI = oneI + ',' + dict_spDid[item]
	# 										if commas == len(gbidX) - 1:
	# 											gbidX = oneI
	# 										commas += 1
	# 									# if no symbol
	# 									else:
	# 										if commas == 0:
	# 											oneI = oneI + item
	# 										else:
	# 											oneI = oneI + ',' + item
	# 										if commas == len(gbidX) - 1:
	# 											gbidX = oneI
	# 										commas += 1
	#
	# 							# print('Found it0')
	# 							break
	# 						elif c == 1:
	# 							print(gbid, 'not in dictionary')
	#
	# 				# if no mapping found
	# 				else:
	# 					gbid = 'No_ID_mapped'
	# 					gbidX = 'No_ID_mapped'
	# 				# fo.write(cols[0] + '\t' + cols[1] + '\t' + cols[2] + '\t' + cols[3] + '\t' + cols[4] + '\t' + gbid + '\t' + cols[6] + '\t' + cols[7] + '\t' + cols[8] + '\t' + cols[9] + '\t' + cols[
	# 				#			 10] + '\t' + gbid + '\t' + cols[12] + '\t' + cols[13] + '\t' + cols[14] + '\t' +
	# 				#		 cols[15] + '\t' + cols[16] + '\t' + cols[17])
	#
	# 				if gbidX != '':
	# 					fo.write(
	# 						cols[0] + '\t' + cols[1] + '\t' + cols[2] + '\t' + cols[3] + '\t' + cols[4] + '\t' + cols[
	# 							6] + '\t' + str(gbidX) + '\t' + cols[7] + '\t' + cols[8] + '\t' + cols[9] + '\t' + cols[
	# 							10] + '\t' + str(gbidX) + '\t' + cols[12] + '\t' + cols[13] + '\t' + cols[14] + '\t' +
	# 						cols[15] + '\t' + cols[16] + '\t' + cols[17])
	# 				elif gbidX == '' or gbidX == '[]':
	# 					fo.write(
	# 						cols[0] + '\t' + cols[1] + '\t' + cols[2] + '\t' + cols[3] + '\t' + cols[4] + '\t' + cols[
	# 							6] + '\t' + str("No_OrthoPara") + '\t' + cols[7] + '\t' + cols[8] + '\t' + cols[
	# 							9] + '\t' + cols[
	# 							10] + '\t' + str("No_OrthoPara") + '\t' + cols[12] + '\t' + cols[13] + '\t' + cols[
	# 							14] + '\t' +
	# 						cols[15] + '\t' + cols[16] + '\t' + cols[17])
	#
	# 			# if multiple genes and both are similar
	#
	# 			if cols[5].find(',') != -1:
	# 				# ortho_para1 = ''
	# 				gbid = ''
	# 				gbidX = ''
	# 				genes1 = cols[5].split(',')
	# 				for i in range(0, len(genes1)):
	# 					# to put comma between
	# 					if i < len(genes1) - 1:
	# 						if genes1[i] in dict_sp1id:
	# 							gbid = gbid + dict_sp1id[genes1[i]] + ','
	# 							c = 0
	# 							# check if this XP in above dictionary, if so write others here
	# 							for key in diction2:
	# 								if dict_sp1id[genes1[i]] in diction2[key]:
	# 									c = 1
	# 									# print('Found it1')
	# 									gbidX = gbidX + str(diction2[key][dict_sp1id[genes1[i]]]) + ','
	# 									# this gbidx[1 agap] has a list of NP's now ..try to convert these to symbols
	# 									# may be should do it at the end not here
	#
	# 									break
	# 								elif c == 1:
	# 									print(dict_sp1id[genes1[i]], 'not in dictionary')
	# 						else:
	# 							gbid = 'No_ID_mapped,'
	# 							gbidX = 'No_ID_mapped,'
	#
	# 					elif i == len(genes1) - 1:
	# 						if genes1[i] in dict_sp1id:
	# 							gbid = gbid + dict_sp1id[genes1[i]]
	# 							c = 0
	# 							for key in diction2:
	# 								if dict_sp1id[genes1[i]] in diction2[key]:
	# 									c = 1
	# 									# print('Found it1a')
	# 									gbidX = gbidX + str(diction2[key][dict_sp1id[genes1[i]]])
	# 									# this is the last one, here gbidx is at its final round
	# 									# work on symbols
	# 									gbidXSymbList = gbidX.split(',')
	# 									if len(gbidXSymbList) == 1:
	# 										for item in gbidXSymbList:
	# 											if item in dict_spDid:  # check if this NP in feature table > if so > write symbol else no symbol found thus orig NP id
	# 												gbidX = dict_spDid[item]
	# 											else:
	# 												gbidX = item
	# 									# print('FOUND item yay ',dict_spDid[item])
	# 									elif len(gbidXSymbList) > 1:
	# 										oneI = ''
	# 										oneIList = []
	# 										commas = 0
	# 										for item in gbidXSymbList:
	# 											item = item.replace('"', '').replace('\'', '').replace('[', '').replace(
	# 												']', '').replace('\\', '').replace(' ', '')
	# 											if item in dict_spDid:
	# 												if dict_spDid[item] not in oneIList:
	# 													oneIList.append(dict_spDid[item])
	# 												# if commas == 0:
	# 												#	oneI = oneI + dict_spDid[item]
	#
	# 												#	else:
	# 												#		oneI = oneI + ',' + dict_spDid[item]
	#
	# 												if commas == len(gbidXSymbList) - 1:
	# 													#	gbidX = oneI
	# 													gbidX = ",".join(map(str, oneIList))
	# 												commas += 1
	# 											# if no symbol
	# 											else:
	# 												if item not in oneIList:
	# 													oneIList.append(item)
	# 												#	if commas == 0:
	# 												#		oneI = oneI + item
	# 												#	else:
	# 												#		oneI = oneI + ',' + item
	# 												if commas == len(gbidXSymbList) - 1:
	# 													#	gbidX = oneI
	# 													gbidX = ",".join(map(str, oneIList))
	# 												commas += 1
	#
	# 									# uptil here symbol work
	#
	# 									break
	# 								elif c == 1:
	# 									print(dict_sp1id[genes1[i]], 'not in dictionary')
	#
	# 						else:
	# 							# if last gene is not mapped then have to do symbol work for prev ones
	# 							# work on symbols
	# 							gbidXSymbList = gbidX.split(',')
	# 							if len(gbidXSymbList) == 1:
	# 								for item in gbidXSymbList:
	# 									if item in dict_spDid:  # check if this NP in feature table > if so > write symbol else no symbol found thus orig NP id
	# 										gbidX = dict_spDid[item]
	# 									else:
	# 										gbidX = item
	# 							# print('FOUND item yay ',dict_spDid[item])
	# 							elif len(gbidXSymbList) > 1:
	# 								oneI = ''
	# 								oneIList = []
	# 								commas = 0
	# 								for item in gbidXSymbList:
	# 									item = item.replace('"', '').replace('\'', '').replace('[', '').replace(']',
	# 																											'').replace(
	# 										'\\', '').replace(' ', '')
	# 									if item in dict_spDid:
	# 										if dict_spDid[item] not in oneIList:
	# 											oneIList.append(dict_spDid[item])
	# 										# if commas == 0:
	# 										#	oneI = oneI + dict_spDid[item]
	#
	# 										#	else:
	# 										#		oneI = oneI + ',' + dict_spDid[item]
	#
	# 										if commas == len(gbidXSymbList) - 1:
	# 											#	gbidX = oneI
	# 											gbidX = ",".join(map(str, oneIList))
	# 										commas += 1
	# 									# if no symbol
	# 									else:
	# 										if item not in oneIList:
	# 											oneIList.append(item)
	# 										#	if commas == 0:
	# 										#		oneI = oneI + item
	# 										#	else:
	# 										#		oneI = oneI + ',' + item
	# 										if commas == len(gbidXSymbList) - 1:
	# 											#	gbidX = oneI
	# 											gbidX = ",".join(map(str, oneIList))
	# 										commas += 1
	#
	# 							# uptil here symbol work
	# 							gbid = gbid + ',No_ID_mapped'
	# 							gbidX = gbidX + ',No_ID_mapped'
	# 				print(gbid, gbidX)
	# 				if gbidX != '' or gbidX != '[]':
	# 					fo.write(
	# 						cols[0] + '\t' + cols[1] + '\t' + cols[2] + '\t' + cols[3] + '\t' + cols[4] + '\t' + cols[
	# 							6] + '\t' + str(gbidX) + '\t' + cols[7] + '\t' + cols[8] + '\t' + cols[9] + '\t' + cols[
	# 							10] + '\t' + str(gbidX) + '\t' + cols[12] + '\t' + cols[13] + '\t' + cols[14] + '\t' +
	# 						cols[15] + '\t' + cols[16] + '\t' + cols[17])
	# 				elif gbidX == '' or gbidX == '[]':
	# 					fo.write(
	# 						cols[0] + '\t' + cols[1] + '\t' + cols[2] + '\t' + cols[3] + '\t' + cols[4] + '\t' + cols[
	# 							6] + '\t' + str("No_OrthoPara") + '\t' + cols[7] + '\t' + cols[8] + '\t' + cols[
	# 							9] + '\t' + cols[
	# 							10] + '\t' + str("No_OrthoPara") + '\t' + cols[12] + '\t' + cols[13] + '\t' + cols[
	# 							14] + '\t' +
	# 						cols[15] + '\t' + cols[16] + '\t' + cols[17])
	#
	#
	# 		# if both are different
	# 		elif cols[5] != cols[10]:
	# 			gbid = ''
	# 			gbidX = ''
	# 			gbid2 = ''
	# 			gbidX2 = ''
	# 			# CASE A
	# 			# Neither col 5 nor col 10 has more than one gene on it---
	# 			if cols[5].find(',') == -1 and cols[10].find(',') == -1:
	# 				if cols[5] in dict_sp1id:
	# 					gbid = dict_sp1id[cols[5]]
	# 					c = 0
	# 					# check if this XP in above dictionary, if so write others here
	# 					for key in diction2:
	# 						if gbid in diction2[key]:
	# 							c = 1
	# 							gbidX = diction2[key][gbid]
	# 							if not gbidX:
	# 								gbidX = ''
	#
	# 							# print('Found it3')
	# 							#### these would be one or more DMEL's NP ids ..time to break them and convert them to their gene symbols
	# 							# check if one or more
	# 							# if one
	# 							# if len(gbidX)==1:
	# 							# print(str(gbidX))
	# 							if len(gbidX) == 1:
	# 								for item in gbidX:
	# 									# if it has symbol
	# 									if item in dict_spDid:
	# 										gbidX = dict_spDid[item]
	# 									# if no symbol
	# 									else:
	# 										gbidX = item
	#
	#
	# 							elif len(gbidX) > 1:
	# 								symbList = []
	# 								commas = 0
	# 								# for each item of gbidx that is Np1, Np2
	# 								for item in gbidX:
	# 									# check if it has symbol from feature table
	# 									if item in dict_spDid:
	# 										if dict_spDid[item] not in symbList:
	# 											symbList.append(dict_spDid[item])
	#
	# 										if commas == len(gbidX) - 1:
	# 											# gbidX=symbList
	# 											# gbidX=*symbList,sep =','
	# 											gbidX = ",".join(map(str, symbList))
	# 										commas += 1
	# 									# if no symbol found
	# 									else:
	# 										if item not in symbList:
	# 											symbList.append(item)
	#
	# 										if commas == len(gbidX) - 1:
	# 											gbidX = ",".join(map(str, symbList))
	# 										commas += 1
	# 							# oneI=''
	#
	# 							break
	# 						elif c == 1:
	# 							print(gbid, 'not in dictionary')
	# 				else:
	# 					gbid = 'No_ID_mapped'
	# 					gbidX = 'No_ID_mapped'
	#
	# 				# fo.write(cols[0]+'\t'+cols[1]+'\t'+cols[2]+'\t'+cols[3]+'\t'+cols[4]+'\t'+cols[5]+'\t'+ortho_para+'\t'+cols[7]+'\t'+cols[8]+'\t'+cols[9]+'\t'+cols[10]+'\t'+cols[11]+'\t'+cols[12]+'\t'+cols[13]+'\t'+cols[14]+'\t'+cols[15]+'\t'+cols[16]+'\t'+cols[17])
	# 				if cols[10] in dict_sp1id:
	# 					gbid2 = dict_sp1id[cols[10]]
	# 					c = 0
	# 					# check if this XP in above dictionary, if so write others here
	# 					for key in diction2:
	# 						if gbid2 in diction2[key]:
	# 							c = 1
	# 							gbidX2 = diction2[key][gbid2]
	# 							if not gbidX2:
	# 								gbidX2 = ''
	# 							symbList = []
	# 							commas = 0
	# 							# for each item of gbidx that is Np1, Np2
	# 							for item in gbidX2:
	# 								# check if it has symbol from feature table
	# 								if item in dict_spDid:
	# 									if dict_spDid[item] not in symbList:
	# 										symbList.append(dict_spDid[item])
	#
	# 									if commas == len(gbidX2) - 1:
	# 										# gbidX2=symbList
	# 										# gbidX2=*symbList,sep =','
	# 										gbidX2 = ",".join(map(str, symbList))
	# 									commas += 1
	# 							# oneI=''
	# 							break
	# 						elif c == 1:
	# 							print(gbid2, 'not in dictionary')
	#
	# 				else:
	# 					gbid2 = 'No_ID_mapped'
	# 					gbidX2 = 'No_ID_mapped'
	#
	# 				# print(gbidX,gbidX2)
	# 				if (gbidX != '' and gbidX != '[]' and gbidX2 != '' and gbidX2 != '[]'):
	# 					fo.write(
	# 						cols[0] + '\t' + cols[1] + '\t' + cols[2] + '\t' + cols[3] + '\t' + cols[4] + '\t' + cols[
	# 							5] + '\t' +
	# 						str(gbidX) + '\t' + cols[7] + '\t' + cols[8] + '\t' + cols[9] + '\t' + cols[10] + '\t' +
	# 						str(gbidX2) + '\t' + cols[12] + '\t' + cols[13] + '\t' + cols[14] + '\t' + cols[
	# 							15] + '\t' + cols[16] + '\t' +
	# 						cols[17])
	# 				elif (gbidX != '' and gbidX != '[]') and (gbidX2 == '' or gbidX2 == '[]'):
	# 					fo.write(
	# 						cols[0] + '\t' + cols[1] + '\t' + cols[2] + '\t' + cols[3] + '\t' + cols[4] + '\t' + cols[
	# 							5] + '\t' +
	# 						str(gbidX) + '\t' + cols[7] + '\t' + cols[8] + '\t' + cols[9] + '\t' + cols[10] + '\t' +
	# 						str("No_OrthoPara") + '\t' + cols[12] + '\t' + cols[13] + '\t' + cols[14] + '\t' + cols[
	# 							15] + '\t' + cols[16] + '\t' +
	# 						cols[17])
	# 				elif (gbidX == '' or gbidX == '[]') and (gbidX2 != '' and gbidX2 != '[]'):
	# 					fo.write(
	# 						cols[0] + '\t' + cols[1] + '\t' + cols[2] + '\t' + cols[3] + '\t' + cols[4] + '\t' + cols[
	# 							5] + '\t' +
	# 						str("No_OrthoPara") + '\t' + cols[7] + '\t' + cols[8] + '\t' + cols[9] + '\t' + cols[
	# 							10] + '\t' +
	# 						str(gbidX2) + '\t' + cols[12] + '\t' + cols[13] + '\t' + cols[14] + '\t' + cols[
	# 							15] + '\t' + cols[16] + '\t' +
	# 						cols[17])
	# 				elif (gbidX == '' or gbidX == '[]') and (gbidX2 == '' or gbidX2 == '[]'):
	# 					fo.write(
	# 						cols[0] + '\t' + cols[1] + '\t' + cols[2] + '\t' + cols[3] + '\t' + cols[4] + '\t' + cols[
	# 							5] + '\t' +
	# 						str("No_OrthoPara") + '\t' + cols[7] + '\t' + cols[8] + '\t' + cols[9] + '\t' + cols[
	# 							10] + '\t' +
	# 						str("No_OrthoPara") + '\t' + cols[12] + '\t' + cols[13] + '\t' + cols[14] + '\t' + cols[
	# 							15] + '\t' + cols[16] + '\t' +
	# 						cols[17])
	#
	# 			# CASE B
	# 			# Column 5 has more than 1 gene but col 10 has 1 gene only
	# 			elif cols[5].find(',') != -1 and cols[10].find(',') == -1:
	# 				gbid = ''
	# 				gbidX = ''
	# 				gbid2 = ''
	# 				gbidX2 = ''
	# 				genes1 = cols[5].split(',')
	# 				for i in range(0, len(genes1)):
	# 					if i < len(genes1) - 1:
	# 						if genes1[i] in dict_sp1id:
	# 							gbid = gbid + dict_sp1id[genes1[i]] + ','
	# 							c = 0
	# 							for key in diction2:
	# 								if dict_sp1id[genes1[i]] in diction2[key]:
	# 									c = 1
	# 									# print('Found it4')
	# 									gbidX = gbidX + str(diction2[key][dict_sp1id[genes1[i]]]) + ','
	#
	# 									break
	# 								elif c == 1:
	# 									print(dict_sp1id[genes1[i]], 'not in dictionary')
	#
	# 						else:
	# 							gbid = 'No_ID_mapped'
	# 							gbidX = 'No_ID_mapped'
	# 					elif i == len(genes1) - 1:
	# 						if genes1[i] in dict_sp1id:
	# 							gbid = gbid + dict_sp1id[genes1[i]]
	# 							c = 0
	# 							for key in diction2:
	# 								if dict_sp1id[genes1[i]] in diction2[key]:
	# 									c = 1
	# 									# print('Found it4a')
	# 									gbidX = gbidX + str(diction2[key][dict_sp1id[genes1[i]]])
	# 									# this is the last one, here gbidx is at its final round
	# 									# work on symbols
	# 									gbidXSymbList = gbidX.split(',')
	# 									if len(gbidXSymbList) == 1:
	# 										for item in gbidXSymbList:
	# 											if item in dict_spDid:  # check if this NP in feature table > if so > write symbol else no symbol found thus orig NP id
	# 												gbidX = dict_spDid[item]
	# 											else:
	# 												gbidX = item
	# 									# print('FOUND item yay ',dict_spDid[item])
	# 									elif len(gbidXSymbList) > 1:
	# 										oneI = ''
	# 										oneIList = []
	# 										commas = 0
	# 										for item in gbidXSymbList:
	# 											item = item.replace('"', '').replace('\'', '').replace('[', '').replace(
	# 												']', '').replace('\\', '').replace(' ', '')
	# 											if item in dict_spDid:
	# 												if dict_spDid[item] not in oneIList:
	# 													oneIList.append(dict_spDid[item])
	# 												# if commas == 0:
	# 												#	oneI = oneI + dict_spDid[item]
	#
	# 												#	else:
	# 												#		oneI = oneI + ',' + dict_spDid[item]
	#
	# 												if commas == len(gbidXSymbList) - 1:
	# 													#	gbidX = oneI
	# 													gbidX = ",".join(map(str, oneIList))
	# 												commas += 1
	# 											# if no symbol
	# 											else:
	# 												if item not in oneIList:
	# 													oneIList.append(item)
	# 												#	if commas == 0:
	# 												#		oneI = oneI + item
	# 												#	else:
	# 												#		oneI = oneI + ',' + item
	# 												if commas == len(gbidXSymbList) - 1:
	# 													#	gbidX = oneI
	# 													gbidX = ",".join(map(str, oneIList))
	# 												commas += 1
	#
	# 									# uptil here symbol work
	# 									break
	# 								elif c == 1:
	# 									print(dict_sp1id[genes1[i]], 'not in dictionary')
	# 						else:
	# 							# this is the last one, here gbidx is at its final round
	# 							# work on symbols
	# 							gbidXSymbList = gbidX.split(',')
	# 							if len(gbidXSymbList) == 1:
	# 								for item in gbidXSymbList:
	# 									if item in dict_spDid:  # check if this NP in feature table > if so > write symbol else no symbol found thus orig NP id
	# 										gbidX = dict_spDid[item]
	# 									else:
	# 										gbidX = item
	# 							# print('FOUND item yay ',dict_spDid[item])
	# 							elif len(gbidXSymbList) > 1:
	# 								oneI = ''
	# 								oneIList = []
	# 								commas = 0
	# 								for item in gbidXSymbList:
	# 									item = item.replace('"', '').replace('\'', '').replace('[', '').replace(']',
	# 																											'').replace(
	# 										'\\', '').replace(' ', '')
	# 									if item in dict_spDid:
	# 										if dict_spDid[item] not in oneIList:
	# 											oneIList.append(dict_spDid[item])
	# 										# if commas == 0:
	# 										#	oneI = oneI + dict_spDid[item]
	#
	# 										#	else:
	# 										#		oneI = oneI + ',' + dict_spDid[item]
	#
	# 										if commas == len(gbidXSymbList) - 1:
	# 											#	gbidX = oneI
	# 											gbidX = ",".join(map(str, oneIList))
	# 										commas += 1
	# 									# if no symbol
	# 									else:
	# 										if item not in oneIList:
	# 											oneIList.append(item)
	# 										#	if commas == 0:
	# 										#		oneI = oneI + item
	# 										#	else:
	# 										#		oneI = oneI + ',' + item
	# 										if commas == len(gbidXSymbList) - 1:
	# 											#	gbidX = oneI
	# 											gbidX = ",".join(map(str, oneIList))
	# 										commas += 1
	#
	# 							# uptil here symbol work
	# 							gbid = gbid + ',No_ID_mapped'
	# 							gbidX = gbidX + ',No_ID_mapped'
	# 				if cols[10] in dict_sp1id:
	# 					gbid2 = dict_sp1id[cols[10]]
	# 					c = 0
	# 					# check if this XP in above dictionary, if so write others here
	# 					for key in diction2:
	# 						if gbid2 in diction2[key]:
	# 							c = 1
	# 							gbidX2 = diction2[key][gbid2]
	# 							if not gbidX2:
	# 								gbidX2 = ''
	#
	# 							# print('Found it4b')
	# 							if len(gbidX2) == 1:
	# 								for item in gbidX2:
	# 									if item in dict_spDid:
	# 										gbidX2 = dict_spDid[item]
	# 							# print('FOUND item yay ',dict_spDid[item])
	# 							elif len(gbidX2) > 1:
	# 								symbList = []
	# 								commas = 0
	# 								# for each item of gbidx that is Np1, Np2
	# 								for item in gbidX2:
	# 									# check if it has symbol from feature table
	# 									if item in dict_spDid:
	# 										if dict_spDid[item] not in symbList:
	# 											symbList.append(dict_spDid[item])
	#
	# 										if commas == len(gbidX2) - 1:
	# 											# gbidX2=symbList
	# 											# gbidX2=*symbList,sep =','
	# 											gbidX2 = ",".join(map(str, symbList))
	#
	# 										commas += 1
	#
	# 							break
	# 						elif c == 1:
	# 							print(gbid2, 'not in dictionary')
	#
	# 				else:
	# 					gbid2 = 'No_ID_mapped'
	# 					gbidX2 = 'No_ID_mapped'
	#
	# 				# check if col 10 has ortholog or paralog or both
	# 				# if (cols[10] in dict_orthologs) and (cols[10] not in dict_paralogs):
	# 				print(gbidX, gbidX2)
	# 				if (gbidX != '' and gbidX != '[]' and gbidX2 != '' and gbidX2 != '[]'):
	# 					fo.write(
	# 						cols[0] + '\t' + cols[1] + '\t' + cols[2] + '\t' + cols[3] + '\t' + cols[4] + '\t' + cols[
	# 							5] + '\t' +
	# 						str(gbidX) + '\t' + cols[7] + '\t' + cols[8] + '\t' + cols[9] + '\t' + cols[10] + '\t' +
	# 						str(gbidX2) + '\t' + cols[12] + '\t' + cols[13] + '\t' + cols[14] + '\t' + cols[15] + '\t' +
	# 						cols[16] + '\t' +
	# 						cols[17])
	# 				elif (gbidX != '' and gbidX != '[]') and (gbidX2 == '' or gbidX2 == '[]'):
	# 					fo.write(
	# 						cols[0] + '\t' + cols[1] + '\t' + cols[2] + '\t' + cols[3] + '\t' + cols[4] + '\t' + cols[
	# 							5] + '\t' +
	# 						str(gbidX) + '\t' + cols[7] + '\t' + cols[8] + '\t' + cols[9] + '\t' + cols[10] + '\t' +
	# 						str("No_OrthoPara") + '\t' + cols[12] + '\t' + cols[13] + '\t' + cols[14] + '\t' + cols[
	# 							15] + '\t' + cols[16] + '\t' +
	# 						cols[17])
	# 				elif (gbidX == '' or gbidX == '[]') and (gbidX2 != '' and gbidX2 != '[]'):
	# 					fo.write(
	# 						cols[0] + '\t' + cols[1] + '\t' + cols[2] + '\t' + cols[3] + '\t' + cols[4] + '\t' + cols[
	# 							5] + '\t' +
	# 						str("No_OrthoPara") + '\t' + cols[7] + '\t' + cols[8] + '\t' + cols[9] + '\t' + cols[
	# 							10] + '\t' +
	# 						str(gbidX2) + '\t' + cols[12] + '\t' + cols[13] + '\t' + cols[14] + '\t' + cols[15] + '\t' +
	# 						cols[16] + '\t' +
	# 						cols[17])
	# 				elif (gbidX == '' or gbidX == '[]') and (gbidX2 == '' or gbidX2 == '[]'):
	# 					fo.write(
	# 						cols[0] + '\t' + cols[1] + '\t' + cols[2] + '\t' + cols[3] + '\t' + cols[4] + '\t' + cols[
	# 							5] + '\t' +
	# 						str("No_OrthoPara") + '\t' + cols[7] + '\t' + cols[8] + '\t' + cols[9] + '\t' + cols[
	# 							10] + '\t' +
	# 						str("No_OrthoPara") + '\t' + cols[12] + '\t' + cols[13] + '\t' + cols[14] + '\t' + cols[
	# 							15] + '\t' + cols[16] + '\t' +
	# 						cols[17])
	# 			# if gbidX!='' and gbidX2!='':
	# 			# 						fo.write(cols[0] + '\t' + cols[1] + '\t' + cols[2] + '\t' + cols[3] + '\t' + cols[4] + '\t' + cols[5]+ '\t' +
	# 			# 							str(gbidX)  + '\t' + cols[7] + '\t' + cols[8] + '\t' + cols[9] +'\t' + cols[10] + '\t' +
	# 			# 							str(gbidX2) + '\t' + cols[12] +  '\t' + cols[13] + '\t' + cols[14] + '\t' +cols[15] + '\t' + cols[16] + '\t' +
	# 			# 							cols[17])
	# 			# 					elif gbidX!='' and gbidX2=='':
	# 			# 						fo.write(cols[0] + '\t' + cols[1] + '\t' + cols[2] + '\t' + cols[3] + '\t' + cols[4] + '\t' + cols[5]+ '\t' +
	# 			# 							str(gbidX)  + '\t' + cols[7] + '\t' + cols[8] + '\t' + cols[9] +'\t' + cols[10] + '\t' +
	# 			# 							str("No_OrthoPara") + '\t' + cols[12] +  '\t' + cols[13] + '\t' + cols[14] + '\t' +cols[15] + '\t' + cols[16] + '\t' +
	# 			# 							cols[17])
	# 			# 					elif gbidX=='' and gbidX2!='':
	# 			# 						fo.write(cols[0] + '\t' + cols[1] + '\t' + cols[2] + '\t' + cols[3] + '\t' + cols[4] + '\t' + cols[5]+ '\t' +
	# 			# 							str("No_OrthoPara")  + '\t' + cols[7] + '\t' + cols[8] + '\t' + cols[9] +'\t' + cols[10] + '\t' +
	# 			# 							str(gbidX2) + '\t' + cols[12] +  '\t' + cols[13] + '\t' + cols[14] + '\t' +cols[15] + '\t' + cols[16] + '\t' +
	# 			# 							cols[17])
	# 			# 					elif gbidX=='' and gbidX2=='':
	# 			# 						fo.write(cols[0] + '\t' + cols[1] + '\t' + cols[2] + '\t' + cols[3] + '\t' + cols[4] + '\t' + cols[5]+ '\t' +
	# 			# 							str("No_OrthoPara")  + '\t' + cols[7] + '\t' + cols[8] + '\t' + cols[9] +'\t' + cols[10] + '\t' +
	# 			# 							str("No_OrthoPara") + '\t' + cols[12] +  '\t' + cols[13] + '\t' + cols[14] + '\t' +cols[15] + '\t' + cols[16] + '\t' +
	# 			# 							cols[17])
	#
	# 			# CASE C
	# 			# 			#Column 5 has 1 gene but col 10 has more than 1 gene
	# 			#
	# 			elif cols[5].find(',') == -1 and cols[10].find(',') != -1:
	# 				gbid = ''
	# 				gbidX = ''
	# 				gbid2 = ''
	# 				gbidX2 = ''
	# 				if cols[5] in dict_sp1id:
	# 					gbid = dict_sp1id[cols[5]]
	# 					c = 0
	# 					# check if this XP in above dictionary, if so write others here
	# 					for key in diction2:
	# 						if gbid in diction2[key]:
	# 							c = 1
	# 							gbidX = diction2[key][gbid]
	# 							if not gbidX:
	# 								gbidX = ''
	#
	# 							#
	# 							### these would be one or more DMEL's NP ids ..time to break them and convert them to their gene symbols
	# 							# check if one or more
	# 							# if one
	# 							if len(gbidX) == 1:
	# 								for item in gbidX:
	# 									if item in dict_spDid:
	# 										gbidX = dict_spDid[item]
	# 							# print('FOUND item yay ',dict_spDid[item])
	# 							elif len(gbidX) > 1:
	# 								symbList = []
	# 								commas = 0
	# 								# for each item of gbidx that is Np1, Np2
	# 								for item in gbidX:
	# 									# check if it has symbol from feature table
	# 									if item in dict_spDid:
	# 										if dict_spDid[item] not in symbList:
	# 											symbList.append(dict_spDid[item])
	#
	# 										if commas == len(gbidX) - 1:
	# 											# gbidX=symbList
	# 											# gbidX=*symbList,sep =','
	# 											gbidX = ",".join(map(str, symbList))
	# 										commas += 1
	#
	# 							break
	# 						elif c == 1:
	# 							print(gbid, 'not in dictionary')
	# 				else:
	# 					gbid = 'No_ID_mapped'
	# 					gbidX = 'No_ID_mapped'
	# 				genes2 = cols[10].split(',')
	# 				for i in range(0, len(genes2)):
	# 					# for first few genes
	# 					if i < len(genes2) - 1:
	# 						if genes2[i] in dict_sp1id:
	# 							gbid2 = gbid2 + dict_sp1id[genes2[i]] + ','
	# 							c = 0
	# 							for key in diction2:
	# 								if dict_sp1id[genes2[i]] in diction2[key]:
	# 									c = 1
	# 									# print('Found it5a')
	# 									# gbidX2.append(diction2[key][dict_sp1id[genes2[i]]])
	# 									gbidX2 = gbidX2 + ','.join(diction2[key][dict_sp1id[genes2[i]]])
	#
	#
	#
	# 									break
	# 								elif c == 1:
	# 									print(dict_sp1id[genes2[i]], 'not in dictionary')
	#
	# 						else:
	# 							gbid2 = 'No_ID_mapped'
	# 							gbidX2 = 'No_ID_mapped'
	# 					# for last gene
	# 					elif i == len(genes2) - 1:
	# 						if genes2[i] in dict_sp1id:
	# 							gbid2 = gbid2 + dict_sp1id[genes2[i]]
	# 							c = 0
	# 							for key in diction2:
	# 								if dict_sp1id[genes2[i]] in diction2[key]:
	# 									c = 1
	# 									# print('Found it5b')
	# 									gbidX2 = gbidX2 +',' +str(diction2[key][dict_sp1id[genes2[i]]])
	#
	# 									# check it out
	# 									# this is the last one, here gbidx is at its final round
	# 									# work on symbols
	# 									gbidX2SymbList = gbidX2.split(',')
	# 									if len(gbidX2SymbList) == 1:
	# 										for item in gbidX2SymbList:
	# 											if item in dict_spDid:  # check if this NP in feature table > if so > write symbol else no symbol found thus orig NP id
	# 												gbidX2 = dict_spDid[item]
	# 											else:
	# 												gbidX2 = item
	# 									# print('FOUND item yay ',dict_spDid[item])
	# 									elif len(gbidX2SymbList) > 1:
	# 										oneI = ''
	# 										oneIList = []
	# 										commas = 0
	# 										for item in gbidX2SymbList:
	# 											item = item.replace('"', '').replace('\'', '').replace('[', '').replace(
	# 												']', '').replace('\\', '').replace(' ', '')
	# 											if item in dict_spDid:
	# 												if dict_spDid[item] not in oneIList:
	# 													oneIList.append(dict_spDid[item])
	# 												# if commas == 0:
	# 												#	oneI = oneI + dict_spDid[item]
	#
	# 												#	else:
	# 												#		oneI = oneI + ',' + dict_spDid[item]
	#
	# 												if commas == len(gbidX2SymbList) - 1:
	# 													#	gbidX = oneI
	# 													gbidX2 = ",".join(map(str, oneIList))
	# 												commas += 1
	# 											# if no symbol
	# 											else:
	# 												if item not in oneIList:
	# 													oneIList.append(item)
	# 												#	if commas == 0:
	# 												#		oneI = oneI + item
	# 												#	else:
	# 												#		oneI = oneI + ',' + item
	# 												if commas == len(gbidX2SymbList) - 1:
	# 													#	gbidX = oneI
	# 													gbidX2 = ",".join(map(str, oneIList))
	# 												commas += 1
	#
	# 									# uptil here symbol work
	# 									break
	# 							if c == 0:
	# 								print(dict_sp1id[genes2[i]], 'not in dictionary')
	# 								gbidX2 = gbidX2 + str(',No_OrthoPara')
	# 								# since second gene did not have any orthopara mapped you need to do symbol thing for first genes NP ids
	# 								# this is the last one, here gbidx is at its final round
	# 								# work on symbols
	# 								gbidX2SymbList = gbidX2.split(',')
	# 								if len(gbidX2SymbList) == 1:
	# 									for item in gbidX2SymbList:
	# 										if item in dict_spDid:  # check if this NP in feature table > if so > write symbol else no symbol found thus orig NP id
	# 											gbidX2 = dict_spDid[item]
	# 										else:
	# 											gbidX2 = item
	# 								# print('FOUND item yay ',dict_spDid[item])
	# 								elif len(gbidX2SymbList) > 1:
	# 									oneI = ''
	# 									oneIList = []
	# 									commas = 0
	# 									for item in gbidX2SymbList:
	# 										item = item.replace('"', '').replace('\'', '').replace('[', '').replace(
	# 											']', '').replace('\\', '').replace(' ', '')
	# 										if item in dict_spDid:
	# 											if dict_spDid[item] not in oneIList:
	# 												oneIList.append(dict_spDid[item])
	# 											# if commas == 0:
	# 											#	oneI = oneI + dict_spDid[item]
	#
	# 											#	else:
	# 											#		oneI = oneI + ',' + dict_spDid[item]
	#
	# 											if commas == len(gbidX2SymbList) - 1:
	# 												#	gbidX = oneI
	# 												gbidX2 = ",".join(map(str, oneIList))
	# 											commas += 1
	# 										# if no symbol
	# 										else:
	# 											if item not in oneIList:
	# 												oneIList.append(item)
	# 											#	if commas == 0:
	# 											#		oneI = oneI + item
	# 											#	else:
	# 											#		oneI = oneI + ',' + item
	# 											if commas == len(gbidX2SymbList) - 1:
	# 												#	gbidX = oneI
	# 												gbidX2 = ",".join(map(str, oneIList))
	# 											commas += 1
	#
	# 						# uptil here symbol work
	#
	# 						else:
	# 							# check it out
	# 							# this is the last one, here gbidx is at its final round
	# 							# work on symbols
	# 							gbidX2SymbList = gbidX2.split(',')
	# 							if len(gbidX2SymbList) == 1:
	# 								for item in gbidX2SymbList:
	# 									if item in dict_spDid:  # check if this NP in feature table > if so > write symbol else no symbol found thus orig NP id
	# 										gbidX2 = dict_spDid[item]
	# 									else:
	# 										gbidX2 = item
	# 							# print('FOUND item yay ',dict_spDid[item])
	# 							elif len(gbidX2SymbList) > 1:
	# 								oneI = ''
	# 								oneIList = []
	# 								commas = 0
	# 								for item in gbidX2SymbList:
	# 									item = item.replace('"', '').replace('\'', '').replace('[', '').replace(']',
	# 																											'').replace(
	# 										'\\', '').replace(' ', '')
	# 									if item in dict_spDid:
	# 										if dict_spDid[item] not in oneIList:
	# 											oneIList.append(dict_spDid[item])
	# 										# if commas == 0:
	# 										#	oneI = oneI + dict_spDid[item]
	#
	# 										#	else:
	# 										#		oneI = oneI + ',' + dict_spDid[item]
	#
	# 										if commas == len(gbidX2SymbList) - 1:
	# 											#	gbidX = oneI
	# 											gbidX2 = ",".join(map(str, oneIList))
	# 										commas += 1
	# 									# if no symbol
	# 									else:
	# 										if item not in oneIList:
	# 											oneIList.append(item)
	# 										#	if commas == 0:
	# 										#		oneI = oneI + item
	# 										#	else:
	# 										#		oneI = oneI + ',' + item
	# 										if commas == len(gbidX2SymbList) - 1:
	# 											#	gbidX = oneI
	# 											gbidX2 = ",".join(map(str, oneIList))
	# 										commas += 1
	#
	# 							# uptil here symbol work
	# 							gbid2 = gbid2 + ',No_ID_mapped'
	# 							gbidX2 = gbidX2 + ',No_ID_mapped'
	# 				print(gbidX, gbidX2)
	# 				if (gbidX != '' and gbidX != '[]' and gbidX2 != '' and gbidX2 != '[]'):
	# 					fo.write(cols[0] + '\t' + cols[1] + '\t' + cols[2] + '\t' + cols[3] + '\t' + cols[4] + '\t' + cols[
	# 						5] + '\t' +
	# 							 str(gbidX) + '\t' + cols[7] + '\t' + cols[8] + '\t' + cols[9] + '\t' + cols[10] + '\t' +
	# 							 str(gbidX2) + '\t' + cols[12] + '\t' + cols[13] + '\t' + cols[14] + '\t' + cols[15] + '\t' +
	# 							 cols[16] + '\t' +
	# 							 cols[17])
	# 				elif (gbidX != '' and gbidX != '[]') and (gbidX2 == '' or gbidX2 == '[]'):
	# 					fo.write(cols[0] + '\t' + cols[1] + '\t' + cols[2] + '\t' + cols[3] + '\t' + cols[4] + '\t' + cols[
	# 						5] + '\t' +
	# 							 str(gbidX) + '\t' + cols[7] + '\t' + cols[8] + '\t' + cols[9] + '\t' + cols[10] + '\t' +
	# 							 str("No_OrthoPara") + '\t' + cols[12] + '\t' + cols[13] + '\t' + cols[14] + '\t' + cols[
	# 								 15] + '\t' + cols[16] + '\t' +
	# 							 cols[17])
	# 				elif (gbidX == '' or gbidX == '[]') and (gbidX2 != '' and gbidX2 != '[]'):
	# 					fo.write(cols[0] + '\t' + cols[1] + '\t' + cols[2] + '\t' + cols[3] + '\t' + cols[4] + '\t' + cols[
	# 						5] + '\t' +
	# 							 str("No_OrthoPara") + '\t' + cols[7] + '\t' + cols[8] + '\t' + cols[9] + '\t' + cols[
	# 								 10] + '\t' +
	# 							 str(gbidX2) + '\t' + cols[12] + '\t' + cols[13] + '\t' + cols[14] + '\t' + cols[15] + '\t' +
	# 							 cols[16] + '\t' +
	# 							 cols[17])
	# 				elif (gbidX == '' or gbidX == '[]') and (gbidX2 == '' or gbidX2 == '[]'):
	# 					fo.write(cols[0] + '\t' + cols[1] + '\t' + cols[2] + '\t' + cols[3] + '\t' + cols[4] + '\t' + cols[
	# 						5] + '\t' +
	# 							 str("No_OrthoPara") + '\t' + cols[7] + '\t' + cols[8] + '\t' + cols[9] + '\t' + cols[
	# 								 10] + '\t' +
	# 							 str("No_OrthoPara") + '\t' + cols[12] + '\t' + cols[13] + '\t' + cols[14] + '\t' + cols[
	# 								 15] + '\t' + cols[16] + '\t' +
	# 							 cols[17])
	# 	# if gbidX!='' and gbidX2!='':
	# 	# 						fo.write(cols[0] + '\t' + cols[1] + '\t' + cols[2] + '\t' + cols[3] + '\t' + cols[4] + '\t' + cols[5]+ '\t' +
	# 	# 							str(gbidX)  + '\t' + cols[7] + '\t' + cols[8] + '\t' + cols[9] +'\t' + cols[10] + '\t' +
	# 	# 							str(gbidX2) + '\t' + cols[12] +  '\t' + cols[13] + '\t' + cols[14] + '\t' +cols[15] + '\t' + cols[16] + '\t' +
	# 	# 							cols[17])
	# 	# 					elif gbidX!='' and gbidX2=='':
	# 	# 						fo.write(cols[0] + '\t' + cols[1] + '\t' + cols[2] + '\t' + cols[3] + '\t' + cols[4] + '\t' + cols[5]+ '\t' +
	# 	# 							str(gbidX)  + '\t' + cols[7] + '\t' + cols[8] + '\t' + cols[9] +'\t' + cols[10] + '\t' +
	# 	# 							str("No_OrthoPara") + '\t' + cols[12] +  '\t' + cols[13] + '\t' + cols[14] + '\t' +cols[15] + '\t' + cols[16] + '\t' +
	# 	# 							cols[17])
	# 	# 					elif gbidX=='' and gbidX2!='':
	# 	# 						fo.write(cols[0] + '\t' + cols[1] + '\t' + cols[2] + '\t' + cols[3] + '\t' + cols[4] + '\t' + cols[5]+ '\t' +
	# 	# 							str("No_OrthoPara")  + '\t' + cols[7] + '\t' + cols[8] + '\t' + cols[9] +'\t' + cols[10] + '\t' +
	# 	# 							str(gbidX2) + '\t' + cols[12] +  '\t' + cols[13] + '\t' + cols[14] + '\t' +cols[15] + '\t' + cols[16] + '\t' +
	# 	# 							cols[17])
	# 	# 					elif gbidX=='' and gbidX2=='':
	# 	# 						fo.write(cols[0] + '\t' + cols[1] + '\t' + cols[2] + '\t' + cols[3] + '\t' + cols[4] + '\t' + cols[5]+ '\t' +
	# 	# 							str("No_OrthoPara")  + '\t' + cols[7] + '\t' + cols[8] + '\t' + cols[9] +'\t' + cols[10] + '\t' +
	# 	# 							str("No_OrthoPara") + '\t' + cols[12] +  '\t' + cols[13] + '\t' + cols[14] + '\t' +cols[15] + '\t' + cols[16] + '\t' +
	# 	# 							cols[17])
	#
	# 	# CASE D
	# 	# Column 5 and col 10 both have more than 1 gene  but different
	# 			elif cols[5].find(',') != -1 and cols[10].find(',') != -1:
	# 				gbid = ''
	# 				gbidX = ''
	# 				gbid2 = ''
	# 				gbidX2 = ''
	# 				genes1 = cols[5].split(',')
	# 				for i in range(0, len(genes1)):
	# 					if i < len(genes1) - 1:
	# 						if genes1[i] in dict_sp1id:
	# 							gbid = gbid + dict_sp1id[genes1[i]] + ','
	# 							c = 0
	# 							for key in diction2:
	# 								if dict_sp1id[genes1[i]] in diction2[key]:
	# 									c = 1
	# 									# print('Found it6')
	# 									gbidX = gbidX + str(diction2[key][dict_sp1id[genes1[i]]]) + ','
	# 									break
	# 								elif c == 1:
	# 									print(dict_sp1id[genes1[i]], 'not in dictionary')
	# 						else:
	# 							gbid = 'No_ID_mapped'
	# 							gbidX2 = 'No_ID_mapped'
	# 					elif i == len(genes1) - 1:
	# 						if genes1[i] in dict_sp1id:
	# 							gbid = gbid + dict_sp1id[genes1[i]]
	# 							c = 0
	# 							for key in diction2:
	# 								if dict_sp1id[genes1[i]] in diction2[key]:
	# 									c = 1
	# 									# print('Found it6a')
	# 									gbidX = gbidX + str(diction2[key][dict_sp1id[genes1[i]]])
	# 									# add here
	# 									# this is the last one, here gbidx is at its final round
	# 									# work on symbols
	# 									gbidXSymbList = gbidX.split(',')
	# 									if len(gbidXSymbList) == 1:
	# 										for item in gbidXSymbList:
	# 											if item in dict_spDid:  # check if this NP in feature table > if so > write symbol else no symbol found thus orig NP id
	# 												gbidX = dict_spDid[item]
	# 											else:
	# 												gbidX = item
	# 									# print('FOUND item yay ',dict_spDid[item])
	# 									elif len(gbidXSymbList) > 1:
	# 										oneI = ''
	# 										oneIList = []
	# 										commas = 0
	# 										for item in gbidXSymbList:
	# 											item = item.replace('"', '').replace('\'', '').replace('[', '').replace(']',
	# 																													'').replace(
	# 												'\\', '').replace(' ', '')
	# 											if item in dict_spDid:
	# 												if dict_spDid[item] not in oneIList:
	# 													oneIList.append(dict_spDid[item])
	# 												# if commas == 0:
	# 												#	oneI = oneI + dict_spDid[item]
	#
	# 												#	else:
	# 												#		oneI = oneI + ',' + dict_spDid[item]
	#
	# 												if commas == len(gbidXSymbList) - 1:
	# 													#	gbidX = oneI
	# 													gbidX = ",".join(map(str, oneIList))
	# 												commas += 1
	# 											# if no symbol
	# 											else:
	# 												if item not in oneIList:
	# 													oneIList.append(item)
	# 												#	if commas == 0:
	# 												#		oneI = oneI + item
	# 												#	else:
	# 												#		oneI = oneI + ',' + item
	# 												if commas == len(gbidXSymbList) - 1:
	# 													#	gbidX = oneI
	# 													gbidX = ",".join(map(str, oneIList))
	# 												commas += 1
	#
	# 									# uptil here symbol work
	# 									break
	# 								elif c == 1:
	# 									print(dict_sp1id[genes1[i]], 'not in dictionary')
	# 						else:
	# 							# work on symbols
	# 							gbidXSymbList = gbidX.split(',')
	# 							if len(gbidXSymbList) == 1:
	# 								for item in gbidXSymbList:
	# 									if item in dict_spDid:  # check if this NP in feature table > if so > write symbol else no symbol found thus orig NP id
	# 										gbidX = dict_spDid[item]
	# 									else:
	# 										gbidX = item
	# 							# print('FOUND item yay ',dict_spDid[item])
	# 							elif len(gbidXSymbList) > 1:
	# 								oneI = ''
	# 								oneIList = []
	# 								commas = 0
	# 								for item in gbidXSymbList:
	# 									item = item.replace('"', '').replace('\'', '').replace('[', '').replace(']',
	# 																											'').replace(
	# 										'\\', '').replace(' ', '')
	# 									if item in dict_spDid:
	# 										if dict_spDid[item] not in oneIList:
	# 											oneIList.append(dict_spDid[item])
	# 										# if commas == 0:
	# 										#	oneI = oneI + dict_spDid[item]
	#
	# 										#	else:
	# 										#		oneI = oneI + ',' + dict_spDid[item]
	#
	# 										if commas == len(gbidXSymbList) - 1:
	# 											#	gbidX = oneI
	# 											gbidX = ",".join(map(str, oneIList))
	# 										commas += 1
	# 									# if no symbol
	# 									else:
	# 										if item not in oneIList:
	# 											oneIList.append(item)
	# 										#	if commas == 0:
	# 										#		oneI = oneI + item
	# 										#	else:
	# 										#		oneI = oneI + ',' + item
	# 										if commas == len(gbidXSymbList) - 1:
	# 											#	gbidX = oneI
	# 											gbidX = ",".join(map(str, oneIList))
	# 										commas += 1
	#
	# 							# uptil here symbol work
	# 							gbid = gbid + ',No_ID_mapped'
	# 							gbidX = gbidX + ',No_ID_mapped'
	# 				genes2 = cols[10].split(',')
	# 				for i in range(0, len(genes2)):
	# 					if i < len(genes2) - 1:
	# 						if genes2[i] in dict_sp1id:
	# 							gbid2 = gbid2 + dict_sp1id[genes2[i]] + ','
	# 							c = 0
	# 							for key in diction2:
	# 								if dict_sp1id[genes2[i]] in diction2[key]:
	# 									c = 1
	# 									# print('Found it6b')
	# 									gbidX2 = gbidX2 + str(diction2[key][dict_sp1id[genes2[i]]])
	# 									break
	# 								elif c == 1:
	# 									print(dict_sp1id[genes2[i]], 'not in dictionary')
	#
	# 						else:
	# 							gbid2 = 'No_ID_mapped'
	# 							gbidX2 = 'No_ID_mapped'
	# 					elif i == len(genes2) - 1:
	# 						if genes2[i] in dict_sp1id:
	# 							gbid2 = gbid2 + dict_sp1id[genes2[i]]
	# 							c = 0
	# 							for key in diction2:
	# 								if dict_sp1id[genes2[i]] in diction2[key]:
	# 									c = 1
	# 									# print('Found it6c')
	# 									gbidX2 = gbidX2 + str(diction2[key][dict_sp1id[genes2[i]]])
	# 									# add here
	# 									# this is the last one, here gbidx is at its final round
	# 									# work on symbols
	# 									gbidX2SymbList = gbidX2.split(',')
	# 									if len(gbidX2SymbList) == 1:
	# 										for item in gbidX2SymbList:
	# 											if item in dict_spDid:  # check if this NP in feature table > if so > write symbol else no symbol found thus orig NP id
	# 												gbidX2 = dict_spDid[item]
	# 											else:
	# 												gbidX2 = item
	# 									# print('FOUND item yay ',dict_spDid[item])
	# 									elif len(gbidX2SymbList) > 1:
	# 										oneI = ''
	# 										oneIList = []
	# 										commas = 0
	# 										for item in gbidX2SymbList:
	# 											item = item.replace('"', '').replace('\'', '').replace('[', '').replace(']',
	# 																													'').replace(
	# 												'\\', '').replace(' ', '')
	# 											if item in dict_spDid:
	# 												if dict_spDid[item] not in oneIList:
	# 													oneIList.append(dict_spDid[item])
	# 												# if commas == 0:
	# 												#	oneI = oneI + dict_spDid[item]
	#
	# 												#	else:
	# 												#		oneI = oneI + ',' + dict_spDid[item]
	#
	# 												if commas == len(gbidX2SymbList) - 1:
	# 													#	gbidX = oneI
	# 													gbidX2 = ",".join(map(str, oneIList))
	# 												commas += 1
	# 											# if no symbol
	# 											else:
	# 												if item not in oneIList:
	# 													oneIList.append(item)
	# 												#	if commas == 0:
	# 												#		oneI = oneI + item
	# 												#	else:
	# 												#		oneI = oneI + ',' + item
	# 												if commas == len(gbidX2SymbList) - 1:
	# 													#	gbidX = oneI
	# 													gbidX2 = ",".join(map(str, oneIList))
	# 												commas += 1
	#
	# 									# uptil here symbol work
	# 									break
	# 							if c == 0:
	# 								print(dict_sp1id[genes2[i]], 'not in dictionary')
	# 						 #	   print(dict_sp1id[genes2[i]], 'not in dictionary')
	# 								gbidX2 = gbidX2 + str(',No_OrthoPara')
	# 								# since second gene did not have any orthopara mapped you need to do symbol thing for first genes NP ids
	# 								# this is the last one, here gbidx is at its final round
	# 								# work on symbols
	# 								gbidX2SymbList = gbidX2.split(',')
	# 								if len(gbidX2SymbList) == 1:
	# 									for item in gbidX2SymbList:
	# 										if item in dict_spDid:  # check if this NP in feature table > if so > write symbol else no symbol found thus orig NP id
	# 											gbidX2 = dict_spDid[item]
	# 										else:
	# 											gbidX2 = item
	# 								# print('FOUND item yay ',dict_spDid[item])
	# 								elif len(gbidX2SymbList) > 1:
	# 									oneI = ''
	# 									oneIList = []
	# 									commas = 0
	# 									for item in gbidX2SymbList:
	# 										item = item.replace('"', '').replace('\'', '').replace('[', '').replace(
	# 											']', '').replace('\\', '').replace(' ', '')
	# 										if item in dict_spDid:
	# 											if dict_spDid[item] not in oneIList:
	# 												oneIList.append(dict_spDid[item])
	# 											# if commas == 0:
	# 											#	oneI = oneI + dict_spDid[item]
	#
	# 											#	else:
	# 											#		oneI = oneI + ',' + dict_spDid[item]
	#
	# 											if commas == len(gbidX2SymbList) - 1:
	# 												#	gbidX = oneI
	# 												gbidX2 = ",".join(map(str, oneIList))
	# 											commas += 1
	# 										# if no symbol
	# 										else:
	# 											if item not in oneIList:
	# 												oneIList.append(item)
	# 											#	if commas == 0:
	# 											#		oneI = oneI + item
	# 											#	else:
	# 											#		oneI = oneI + ',' + item
	# 											if commas == len(gbidX2SymbList) - 1:
	# 												#	gbidX = oneI
	# 												gbidX2 = ",".join(map(str, oneIList))
	# 											commas += 1
	#
	# 						# uptil here symbol work
	# 						else:
	# 							# work on symbols
	# 							gbidX2SymbList = gbidX2.split(',')
	# 							if len(gbidX2SymbList) == 1:
	# 								for item in gbidX2SymbList:
	# 									if item in dict_spDid:  # check if this NP in feature table > if so > write symbol else no symbol found thus orig NP id
	# 										gbidX2 = dict_spDid[item]
	# 									else:
	# 										gbidX2 = item
	# 							# print('FOUND item yay ',dict_spDid[item])
	# 							elif len(gbidX2SymbList) > 1:
	# 								oneI = ''
	# 								oneIList = []
	# 								commas = 0
	# 								for item in gbidX2SymbList:
	# 									item = item.replace('"', '').replace('\'', '').replace('[', '').replace(']',
	# 																											'').replace(
	# 										'\\', '').replace(' ', '')
	# 									if item in dict_spDid:
	# 										if dict_spDid[item] not in oneIList:
	# 											oneIList.append(dict_spDid[item])
	# 										# if commas == 0:
	# 										#	oneI = oneI + dict_spDid[item]
	#
	# 										#	else:
	# 										#		oneI = oneI + ',' + dict_spDid[item]
	#
	# 										if commas == len(gbidX2SymbList) - 1:
	# 											#	gbidX = oneI
	# 											gbidX2 = ",".join(map(str, oneIList))
	# 										commas += 1
	# 									# if no symbol
	# 									else:
	# 										if item not in oneIList:
	# 											oneIList.append(item)
	# 										#	if commas == 0:
	# 										#		oneI = oneI + item
	# 										#	else:
	# 										#		oneI = oneI + ',' + item
	# 										if commas == len(gbidX2SymbList) - 1:
	# 											#	gbidX = oneI
	# 											gbidX2 = ",".join(map(str, oneIList))
	# 										commas += 1
	#
	# 							# uptil here symbol work
	# 							gbid2 = gbid2 + ',No_ID_mapped'
	# 							gbidX2 = gbidX2 + ',No_ID_mapped'
	#
	# 				if (gbidX != '' and gbidX != '[]' and gbidX2 != '' and gbidX2 != '[]'):
	# 					fo.write(
	# 						cols[0] + '\t' + cols[1] + '\t' + cols[2] + '\t' + cols[3] + '\t' + cols[4] + '\t' + cols[5] + '\t' +
	# 						str(gbidX) + '\t' + cols[7] + '\t' + cols[8] + '\t' + cols[9] + '\t' + cols[10] + '\t' +
	# 						str(gbidX2) + '\t' + cols[12] + '\t' + cols[13] + '\t' + cols[14] + '\t' + cols[15] + '\t' + cols[
	# 							16] + '\t' +
	# 						cols[17])
	# 				elif (gbidX != '' and gbidX != '[]') and (gbidX2 == '' or gbidX2 == '[]'):
	# 					fo.write(
	# 						cols[0] + '\t' + cols[1] + '\t' + cols[2] + '\t' + cols[3] + '\t' + cols[4] + '\t' + cols[5] + '\t' +
	# 						str(gbidX) + '\t' + cols[7] + '\t' + cols[8] + '\t' + cols[9] + '\t' + cols[10] + '\t' +
	# 						str("No_OrthoPara") + '\t' + cols[12] + '\t' + cols[13] + '\t' + cols[14] + '\t' + cols[15] + '\t' +
	# 						cols[16] + '\t' +
	# 						cols[17])
	# 				elif (gbidX == '' or gbidX == '[]') and (gbidX2 != '' and gbidX2 != '[]'):
	# 					fo.write(
	# 						cols[0] + '\t' + cols[1] + '\t' + cols[2] + '\t' + cols[3] + '\t' + cols[4] + '\t' + cols[5] + '\t' +
	# 						str("No_OrthoPara") + '\t' + cols[7] + '\t' + cols[8] + '\t' + cols[9] + '\t' + cols[10] + '\t' +
	# 						str(gbidX2) + '\t' + cols[12] + '\t' + cols[13] + '\t' + cols[14] + '\t' + cols[15] + '\t' + cols[
	# 							16] + '\t' +
	# 						cols[17])
	# 				elif (gbidX == '' or gbidX == '[]') and (gbidX2 == '' or gbidX2 == '[]'):
	# 					fo.write(
	# 						cols[0] + '\t' + cols[1] + '\t' + cols[2] + '\t' + cols[3] + '\t' + cols[4] + '\t' + cols[5] + '\t' +
	# 						str("No_OrthoPara") + '\t' + cols[7] + '\t' + cols[8] + '\t' + cols[9] + '\t' + cols[10] + '\t' +
	# 						str("No_OrthoPara") + '\t' + cols[12] + '\t' + cols[13] + '\t' + cols[14] + '\t' + cols[15] + '\t' +
	# 						cols[16] + '\t' +
	# 						cols[17])
	# # if gbidX!='' and gbidX2!='':


# 						fo.write(cols[0] + '\t' + cols[1] + '\t' + cols[2] + '\t' + cols[3] + '\t' + cols[4] + '\t' + cols[5]+ '\t' +
# 							str(gbidX)  + '\t' + cols[7] + '\t' + cols[8] + '\t' + cols[9] +'\t' + cols[10] + '\t' + 
# 							str(gbidX2) + '\t' + cols[12] +  '\t' + cols[13] + '\t' + cols[14] + '\t' +cols[15] + '\t' + cols[16] + '\t' + 
# 							cols[17])
# 					elif gbidX!='' and gbidX2=='':
# 						fo.write(cols[0] + '\t' + cols[1] + '\t' + cols[2] + '\t' + cols[3] + '\t' + cols[4] + '\t' + cols[5]+ '\t' +
# 							str(gbidX)  + '\t' + cols[7] + '\t' + cols[8] + '\t' + cols[9] +'\t' + cols[10] + '\t' + 
# 							str("No_OrthoPara") + '\t' + cols[12] +  '\t' + cols[13] + '\t' + cols[14] + '\t' +cols[15] + '\t' + cols[16] + '\t' + 
# 							cols[17])
# 					elif gbidX=='' and gbidX2!='':
# 						fo.write(cols[0] + '\t' + cols[1] + '\t' + cols[2] + '\t' + cols[3] + '\t' + cols[4] + '\t' + cols[5]+ '\t' +
# 							str("No_OrthoPara")  + '\t' + cols[7] + '\t' + cols[8] + '\t' + cols[9] +'\t' + cols[10] + '\t' + 
# 							str(gbidX2) + '\t' + cols[12] +  '\t' + cols[13] + '\t' + cols[14] + '\t' +cols[15] + '\t' + cols[16] + '\t' + 
# 							cols[17])
# 					elif gbidX=='' and gbidX2=='':								
# 						fo.write(cols[0] + '\t' + cols[1] + '\t' + cols[2] + '\t' + cols[3] + '\t' + cols[4] + '\t' + cols[5]+ '\t' +
# 							str("No_OrthoPara")  + '\t' + cols[7] + '\t' + cols[8] + '\t' + cols[9] +'\t' + cols[10] + '\t' + 
# 							str("No_OrthoPara") + '\t' + cols[12] +  '\t' + cols[13] + '\t' + cols[14] + '\t' +cols[15] + '\t' + cols[16] + '\t' + 
# 							cols[17])


# ./orthologyMapping.py -np1 tcas -sp1id TCAST.idmap.txt -sp1pc TCAST.97pc.txt -brh TCAST_DMELA.brh -sp2pc DMELA.97pc.txt -sp2id DMELA.idmap.txt -geneSet OSG2geneSet.txt -symb fb_synonym_fb_2020_02.tsv -conv true -setConv OGS3toOGS2_conversion.txt -so scrmshawOutput_peaksCalled_antennal_lobe_imm_1388_peaks.bed -sep ':' -flip TRUE
# ../orthologyMapping.py -np1 agamb -sp1id AGAMB.idmap.txt -sp1pc AGAMB.97pc.txt -brh AGAMB_DMELA.brh -sp2pc DMELA.97pc.txt -sp2id DMELA.idmap.txt -geneSet AGAM_4.9_genes -symb fb_synonym_fb_2020_02.tsv -conv false -so scrmshawOutput_peaksCalled_adult_circulatory_imm_MedianPointAmplitudeCurve_633_peaks.bed
# ./orthologyMapping.py -np1 apis -sp1id AMELL.idmap.txt -sp1pc AMELL.97pc.txt -brh AMELL_DMELA.brh -sp2pc DMELA.97pc.txt -sp2id DMELA.idmap.txt -geneSet OSG2geneSet.txt -symb fb_synonym_fb_2020_02.tsv -conv true -setConv OSG3toOSG2_conversion.txt -so scrmshawOutput_peaksCalled_adult_circulatory_imm_MedianPointAmplitudeCurve_575_peaks.bed -sep '|' -flip false

main()
