#Import the packages we need
import allel;
import time;
import sys;
import numpy;

#@getWindowsByVariants Documentation
#$pos : An ordered list of variant positions.
#$winSize : The physical size of each window.
#$stpSize : The step size to take between each window.
#$genSize : The total size of the genome.
#Gets the start index and length of window sizes in terms of variant loci
#<-Returns:
#$winStarts : A list of start indices for each window of size $winSize for the list of variants, $pos.
#$winVarSizes : A list of the number of variants in each window.
def getWindowsByVariants(pos, winSize, stpSize, genSize):
	winStarts = [];
	winVarSizes = [];
	curLoc = 1;
	curEnd = curLoc + winSize;
	while (curLoc < genSize):
		curEnd = curLoc + winSize;
		ctr = 0;
		varStart = 0;
		foundFirst = False;
		for i in range(len(pos)):
			if(pos[i]>=curLoc and pos[i] < curEnd):
				ctr +=1;
				if(not foundFirst):
					varStart = i;
					foundFirst = True;
		winStarts.append(varStart);
		winVarSizes.append(ctr);
		if((curLoc+winSize) > genSize):
			curLoc = genSize;
		else:
			curLoc += stpSize;
	return [winStarts, winVarSizes];

#@printTimeLapsed Documentation
#$start : The clock time when the script began running
#Prints the elapsed time since a given start time.
def printTimeLapsed(start):
	print("Time Lapsed:" + str((time.perf_counter()-start)));

#@getThetaH Documentation
#$positions : An ordered list of variant positions
#$ac : An ordered list of the allele counts of each variant in the population
#$start : The start position for the window being examined
#$winSize : The size (in base pairs) of the window.
#$sampSize : The number of genomes in the population
#Calcultates thetaH for a window of size $winSize beginning at $start
#<-Returns:
#$thetaSum/winSize : The calculated value of ThetaH for the window.
def getThetaH(positions, ac, start, winSize, sampSize):

	thetaSum = 0;
	denom = sampSize * (sampSize-1);
	#Because subsamples could have variants that are "fixed" but were not fixed 
	#in the populaiton; the Si array (counter for number of variants with n copies in the sample)
	#needs to include an index for all samples having that allele. This value will not
	#be used in calculating thetaH
	Si = [0]*(sampSize);
	for i in range(len(positions)):
		if(positions[i] >= start and positions[i] < (start+winSize)):
			Si[(ac[i][1]-1)] += 1;
	for i in range(sampSize-1):
		thetaSum += ((2 * Si[i] * ((i+1)**2)) / denom);
	return thetaSum/winSize;

#@getThetaPi Documentation
#$positions : An ordered list of variant positions
#$ac : An ordered list of the allele counts of each variant in the population
#$start : The start position for the window being examined
#$winSize : The size (in base pairs) of the window.
#$sampSize : The number of genomes in the population
#Calcultates thetaPi for a window of size $winSize beginning at $start; unused in script, but useful for error checking getThetaH, etc
#<-Returns:
#$thetaSum/winSize : The calculated value of ThetaPi for the window.
def getThetaPi(positions, ac, start, winSize, sampSize):

	thetaSum = 0;
	denom = sampSize * (sampSize-1);
	#Because subsamples could have variants that are "fixed" but were not fixed 
	#in the populaiton; the Si array (counter for number of variants with n copies in the sample)
	#needs to include an index for all samples having that allele. This value will not
	#be used in calculating thetaPi
	Si = [0]*(sampSize);
	for i in range(len(positions)):
		if(positions[i] >= start and positions[i] < (start+winSize)):
			Si[(ac[i][1]-1)] += 1;
	for i in range(sampSize-1):
		thetaSum += ((2 * Si[i] * (i+1) * (sampSize - i - 1)) / denom);
	return thetaSum/winSize;

#@getOutputLineBegin Documentation
#$modParams : A list of numbers which detail the model parameters used: model type, required parameters, recombination parameters, progeny skew parameters, and DFE parameters.
#$sampScheme : A string describing the type of sampling done (full or partial).
#$sampSize : Number of genomes in the sample
#Generates the output that begins every line of the .csv output if parameter information is to be included in the output. 
#This method has low reusability as a large portion is specifically encoded with parameters from this project.
#<-Returns:
#$outBeg : The parameter information as a comma-seperated string.
def getOutputLineBegin(modParams, sampScheme, sampSize):
	#Add the model type and replicate number first
	outBeg = str(modParams[0]);
	outBeg += f",{modParams[len(modParams)-1]}";

	#Parameter Distributions
	muDist = ["2.135e-6", "2.135e-5", "2.135e-4"];
	reproDist = ["1"];
	kDist = ["5e3", "5e4", "1e5"];
	initDist = ["1", "30", "1000"];
	runtimeDist = ["336", "672", "1008"];
	rDist = ["1e-5", "5.5e-5", "10e-5"];
	xiDist = ["0.0001", "0.003", "0.1"];
	burstNDist = ["20", "100", "200"];
	dfeDist = ["4:1", "1:1", "1:4"];

	#Required Parameters
	reqdChew = int(modParams[1]);
	maxReqd = ((3**4)-1);
	reqdParam = [];
	while maxReqd > 0:
		reqdParam.append((reqdChew % 3));
		reqdChew = reqdChew // 3;
		maxReqd = maxReqd // 3;
	outBeg += f",{muDist[reqdParam[0]]},{reproDist[0]},{kDist[reqdParam[1]]},{initDist[reqdParam[2]]},{runtimeDist[reqdParam[3]]}";

	#recomb params
	recombParam = int(modParams[2]);
	if recombParam < 0:
		outBeg += ",na";
	else:
		outBeg += f",{rDist[recombParam]}";
	
	#progeny skew params
	progSkewParam = int(modParams[3]);
	if progSkewParam < 0:
		outBeg += ",na,na";
	else:
		outBeg += f",{xiDist[(progSkewParam % 3)]},{burstNDist[((progSkewParam//3)%3)]}";
	#dfe params
	dfeParam = int(modParams[4]);
	if dfeParam<0:
		outBeg += ",na";
	else:
		outBeg += f",{dfeDist[dfeParam]}";
	outBeg += f",{sampScheme},{str(sampSize)}";

	return outBeg;

#@countFixed Documentation
#$fixedMut : A 2 dimensional list of fixed variants. Information on the position, generation tick, and fixation tick is provided for each variant. 
#$burntime : The time the burnin was run for.
#$start : The first position (in base pairs) to count number of fixed alleles.
#$winSize : The size of the window within which fixed alleles are being counted
#Calculates number of fixed burnin mutations, fixed new mutations, and fixed in sample mutations
#<-Returns:
#$fixCount : A list of size 3 with the number of variants from the burnin fixed, the number of variants that arose during the simulation and fixed, and the number of variants that appear fixed in the sample.
def countFixed(fixedMut, burntime, start, winSize):
	fixCount = [0, 0, 0];
	if(len(fixedMut) > 0):
		for position, mutBirth in zip(fixedMut[0], fixedMut[1]):
			if (position >= start and position < (start+winSize)):
				if (mutBirth == None):
					fixCount[2] += 1;
				elif(mutBirth > burntime):
					fixCount[1] += 1;
				else:
					fixCount[0] +=1;
	else:
		if debug: print("No fixed alleles in population or sample.");
	return fixCount;

#Main method, runs the script.
def main():
	print(debug);

	#Get command line arguments; if none use defaults.
	if len(sys.argv) <= 1:
		fileStem = "sc2-rep1_15.39.0.8.0.1_partial_1000";
		inputMS = True;
		genomeSize = 30000;
	elif len(sys.argv) == 4:
		fileStem = sys.argv[1];
		inputMS = (sys.argv[2]=='m');
		genomeSize = int(sys.argv[3]);
	else:
		sys.exit("Expecting zero (default) or three arguments: the stem of the input/output files, the main input file type, and the genome size. Cancelling script.");
	
	haplotypeInfo = True;

	#Defining and Initiating Variables
	startTime = time.perf_counter();

	#Read main input file
	if (inputMS):
		#Preparing variables for reading .ms input file
		inputFile = fileStem + ".output.ms";
		positions = [];
		genomes = [];
		if debug: print(inputFile);

		#Reading the input .ms file.
		with open(inputFile, 'r+') as msFile:
			if debug: print(msFile);
			counter = 0;
			for curLine in msFile.readlines():
				curLine = curLine.rstrip();
				if(counter >= 3):
					curGenome = list(curLine);
					genomes.append(curGenome);
				elif(counter == 2):
					positions = curLine.split(' ');
					positions.pop(0);
					positions = list(map(float, positions));
					positions = [round((i*(genomeSize-1))+1) for i in positions];
				elif(counter == 1):
					tmp = curLine.split(' ');
					nvariants = int(tmp[1]);
					del tmp;
				elif(curLine != "//"):
					print("First line of .ms file is not valid.");
					break;
				counter += 1;
			else:
				if debug: print(".ms file read successfully.");
				nSamp = counter - 3;
				del counter, curLine;
	
		if debug: printTimeLapsed(startTime);
		#Genomes in the .ms file are read as a 2D list with the first dimension corresponding to sample and the second dimension corresponding to variant.
		#This needs to be transposed to be used by the functions in scikit-allel.
		genomes = list(zip(*genomes));
		haplArray = allel.HaplotypeArray(genomes, dtype='i1');
		alleleCount = haplArray.count_alleles();
	else:
		inputFile = fileStem + ".output.vcf";
		callset = allel.read_vcf(f"{inputFile}", numbers={'GT':1});
		positions = callset['variants/POS'];
		haplArray = allel.HaplotypeArray(callset['calldata/GT']);
		nSamp= len(haplArray[0]);
		if(nSamp <=1):
			nSamp = 100;
			with open(inputFile, 'r+') as vcfFile:
				if debug: print("vcf had only one sample; attempting to get allele counts directly from vcf file.");
				haplotypeInfo = False;
				aC = [];
				refCtTg = "RO";
				altCtTg = "AO";
				for curLine in vcfFile.readlines():
					curLine = curLine.rstrip();
					if(curLine[0] != '#'):
						infoAry = curLine.split("\t")[7].split(";");
						curRefCt = 0;
						curAltCt = 0;
						for info in infoAry:
							if info.split("=")[0] == refCtTg: curRefCt = info.split("=")[1];
							if info.split("=")[0] == altCtTg: curAltCt = info.split("=")[1];
						aC.append([curRefCt, curAltCt]);
			alleleCount = allel.AlleleCountsArray(aC, dtype='u1')
		else:
			alleleCount = haplArray.count_alleles();
	
	#Prepare output related strings and files
	#Get output file name
	outputCSV = f"{fileStem}.output.csv";
	with open(outputCSV, 'a', encoding="utf-8") as outFile:
		if (outputParameterColumns):
			#Get model parameters and the associated output that begins each line
			modelParams = (fileStem.split("_")[1]).split(".");
			sampType = (fileStem.split("_")[2]);
			outputBegin = getOutputLineBegin(modelParams, sampType, nSamp);
			#Write Header Line
			outFile.write("modelNum,rep,mu,repro,K,init,runtime,r,xi,burstN,dfeDist,sampScheme,nSamp,win/full,filter,start,size,stat,value\n");
		else:
			outputBegin = f"{fileStem}";
			#Write Header Line
			outFile.write("fileName,win/full,filter,start,size,stat,value\n");
	outputBegin_Full = f"{outputBegin},full,0,1,{genomeSize}";
	if debug: print(outputBegin_Full);

	#Do stats work requiring output file from SLIM detailing # of fixed mutations
	if (doNFixed):
		inputFile = "_".join(fileStem.split("_")[0:2]) + ".output.fix";
		runtime=int(outputBegin.split(',')[10]) if outputParameterColumns else 1;
		if debug: print(inputFile);
		fixedAlleles = [];
		with open(inputFile, 'r+') as fixFile:
			if debug: print(fixFile);
			counter = 0;
			for curLine in fixFile.readlines():
				curLine = curLine.rstrip();
				if(counter >=2):
					fixedAlleles.append([int(curLine.split(' ')[i]) for i in [3, 7, 8]]);
				elif(counter == 1 and curLine != "Mutations:"):
					print("Second line of .fix file is not valid.");
					break;
				elif(counter == 0 and curLine.split(' ')[0] != "#OUT:"):
					print("First line of .fix file is not valid");
					break;
				counter += 1;
			else:
				if debug: print(".fix file read successfully.");
				del counter, curLine;

		for position, ac in zip(positions, list(zip(*alleleCount))[1]):
			if ac == nSamp:
				fixedAlleles.append([position, None, burnintime+runtime]);
		fixedAlleles = list(zip(*sorted(fixedAlleles, key=lambda x:x[0])));
		fixCount = countFixed(fixedAlleles, burnintime, 1, genomeSize);
		with open(outputCSV, 'a', encoding="utf-8") as outFile:
			#fixBurnMuts = number of variants that arose during the burnin period and became fixed
			outFile.write(f"{outputBegin_Full},fixBurnMuts,{fixCount[0]}\n");
			#fixNewMuts = number of variants that arose and became fixed during the simulation
			outFile.write(f"{outputBegin_Full},fixNewMuts,{fixCount[1]}\n");
			#fixSampMuts = number of variants that were not fixed in the population, but appear fixed in the sample
			outFile.write(f"{outputBegin_Full},fixSampMuts,{fixCount[2]}\n");
		for windowSize, stepSize in zip(windowSizes, stepSizes):
			if debug: print(f"{windowSize}, {stepSize}");
			curStart = 1;
			while(curStart < genomeSize):
				windowedFixCount = countFixed(fixedAlleles, burnintime, curStart, windowSize);
				with open(outputCSV, 'a', encoding="utf-8") as outFile:
					outFile.write(f"{outputBegin},windowed,{curStart},{windowSize},fixBurnMuts,{windowedFixCount[0]}\n");
					outFile.write(f"{outputBegin},windowed,{curStart},{windowSize},fixNewMuts,{windowedFixCount[1]}\n");
					outFile.write(f"{outputBegin},windowed,{curStart},{windowSize},fixSampMuts,{windowedFixCount[2]}\n");
				if((curStart+windowSize) > genomeSize):
					curStart = genomeSize;
				else:
					curStart += stepSize;
		else:
			del windowSize, stepSize, curStart, fixCount, windowedFixCount, fixedAlleles, runtime;
			if debug: print("Variables used in calculating number of fixed alleles deleted.");
		if debug: print("Fixed allele counts computed");
	else:
		if debug: print("Skipping fixed allele counts.");

	#No more input files to be read so this variable can be deleted
	del inputFile;
	if debug: printTimeLapsed(startTime);
	
	#Get counts of total variants; then, for each filtering cutoff, get the number of variants higher than the cutoff,
	#the number of variants lower than the cutoff, the number of variants that pass filtering critera, and the number
	##of unique haplotypes (if haplotype information is available). 
	if(doNFiltVar):
		if debug: print(f"Total variants = {len(alleleCount)}");
		with open(outputCSV, 'a', encoding="utf-8") as outFile:
			outFile.write(f"{outputBegin},full,0,1,{genomeSize},totalVars,{len(alleleCount)}\n");
		for cutoff in filterFreqs:
			if haplotypeInfo:
				#Because of the structure of numpy arrays, it is much quicker if we construct an empty array
				#of maximum size and then slice out the unneeded rows later.
				filteredHaplArray = numpy.empty(shape = [len(alleleCount),1000]);
				curFiltVar  = 0;
			minCt = (nSamp * cutoff)//1;
			if (minCt>0):
				maxCt = nSamp - minCt;
				if debug: print(f"Cutoff Frequency {cutoff}: Minimum Count = {minCt}, Maximum Count = {maxCt}");
				cnt = 0;
				pvtCnt = 0;
				for ind in range(len(alleleCount)-1, -1, -1):
					refCt = alleleCount[ind][0];
					if(refCt <= minCt or refCt >= maxCt):
						cnt += 1;
						if (refCt >= maxCt):
							pvtCnt +=1;
					elif haplotypeInfo:
						filteredHaplArray[curFiltVar] = haplArray[ind];
						curFiltVar += 1;
				if haplotypeInfo:
					filteredHaplArray = filteredHaplArray[0:curFiltVar-1];
					filteredHaplArray = allel.HaplotypeArray(filteredHaplArray, dtype='i1');
					filteredUniqueHaplotypes = len(filteredHaplArray.distinct_counts());

				with open(outputCSV, 'a', encoding="utf-8") as outFile:
					#fixFiltMuts = number of variant sites that become fixed, new SNPs due to the MAF being above the max threshold (e.g., above 98% for a 2% threshold)
					outFile.write(f"{outputBegin},full,{cutoff},1,{genomeSize},fixFiltMuts,{cnt-pvtCnt}\n");
					#errFiltMuts = number of variant sites that erroneously become "fixed" as the reference allele due to the MAF being under the threshold
					outFile.write(f"{outputBegin},full,{cutoff},1,{genomeSize},errFiltMuts,{pvtCnt}\n");
					#filtSNPs = number of segregating variants clearing the filtering criteria
					outFile.write(f"{outputBegin},full,{cutoff},1,{genomeSize},filtSNPs,{len(alleleCount)-cnt}\n");
					#filtUniqHapl = number of unique haplotypes after filtering
					if haplotypeInfo: outFile.write(f"{outputBegin},full,{cutoff},1,{genomeSize},filtUniqHapl,{filteredUniqueHaplotypes}\n");
					if debug: print(f"Total to remove = {cnt}; Total private variants = {pvtCnt}");
			elif (debug):
				print(f"Filter frequency of {cutoff} is too small; there would be no filtering based on sample size of {nSamp}.");
	elif debug: print("Skipping filtered variant stats");

	#Do stats work other than LD stats requiring O(n^2) time algorithms.
	if (doNonLD):
		thetaPi = allel.sequence_diversity(positions, alleleCount);
		thetaW = allel.watterson_theta(positions, alleleCount);
		thetaH = getThetaH(positions, alleleCount, 1, genomeSize, nSamp);
		fayWuH = thetaPi - thetaH;
		tajimasD = allel.tajima_d(alleleCount, pos=positions);
		if haplotypeInfo: 
			nUniqueHaplotypes = len(haplArray.distinct_counts());
			haplotypeDiversity = allel.haplotype_diversity(haplArray);
			garudsH = allel.garud_h(haplArray);
		
		#Write full genome stats to output file
		with open(outputCSV, 'a', encoding="utf-8") as outFile:
			outFile.write(f"{outputBegin_Full},thetaPi,{thetaPi}\n");
			outFile.write(f"{outputBegin_Full},thetaW,{thetaW}\n");
			outFile.write(f"{outputBegin_Full},thetaH,{thetaH}\n");
			outFile.write(f"{outputBegin_Full},fayWuH,{fayWuH}\n");
			outFile.write(f"{outputBegin_Full},tajimasD,{tajimasD}\n");
			if haplotypeInfo: 
				outFile.write(f"{outputBegin_Full},nUniqHapl,{nUniqueHaplotypes}\n");
				outFile.write(f"{outputBegin_Full},haplDiv,{haplotypeDiversity}\n");
				outFile.write(f"{outputBegin_Full},gH1,{garudsH[0]}\n");
				outFile.write(f"{outputBegin_Full},gH12,{garudsH[1]}\n");
				outFile.write(f"{outputBegin_Full},gH123,{garudsH[2]}\n");
				outFile.write(f"{outputBegin_Full},gH2/H1,{garudsH[3]}\n");

		#Delete variables now that they are written to output file.
		del thetaPi, thetaW, thetaH, fayWuH, tajimasD;
		if haplotypeInfo: del nUniqueHaplotypes, haplotypeDiversity, garudsH;
		if debug: print("Variables used in calculating full genome non-LD stats deleted.");

		#Calculate windowed statistics for all window size-step size pairs.
		for windowSize, stepSize in zip(windowSizes, stepSizes):
			if debug: print(f"{windowSize}, {stepSize}");
			windowedThetaPi = allel.windowed_diversity(positions, alleleCount, size=windowSize, start=1, stop=genomeSize, step=stepSize, fill=-1);
			if debug: print("Windowed ThetaPi Calculated");
			windowedThetaW = allel.windowed_watterson_theta(positions, alleleCount, size=windowSize, start=1, stop=genomeSize, step=stepSize, fill=-1)
			if debug: print("Windowed ThetaW Calculated");
			windowedThetaH = [];
			curStart = 1;
			while(curStart < genomeSize):
				windowedThetaH.append(getThetaH(positions, alleleCount, curStart, windowSize, nSamp));
				if((curStart+windowSize) > genomeSize):
					curStart = genomeSize;
				else:
					curStart += stepSize;
			del curStart;
			if debug: print("Windowed ThetaH Calculated");
			windowedFayWuH = [];
			for i in range(len(windowedThetaPi[0])):
				windowedFayWuH.append(windowedThetaPi[0][i]-windowedThetaH[i]);
			if debug: print("Windowed FayWeH Calculated");
			windowedTajimasD = allel.windowed_tajima_d(positions, alleleCount, size=windowSize, start=1, stop=genomeSize, step=stepSize);
			if debug: print("Windowed Tajima's D Calculated");
			if haplotypeInfo:
				varWindows = getWindowsByVariants(positions, windowSize, stepSize, genomeSize);
				if debug: print("Variant Windows Found");
				windowedHaplotypeDiversity = [];
				windowedGarudsH = [[],[],[],[]];
				for i in range(len(varWindows[0])):
					if debug: print(f"size:{varWindows[1][i]}, start:{varWindows[0][i]}, stop:{(varWindows[0][i]+varWindows[1][i])}");
					if(varWindows[1][i]<=0):
						if debug: print("Current window has no variation; assigning stat to 'na'");
						curGarudsH = [["na"],["na"],["na"],["na"]];
						windowedHaplotypeDiversity.append("na");
					else:
						curGarudsH = allel.moving_garud_h(haplArray, size=varWindows[1][i], start=varWindows[0][i], stop=(varWindows[0][i] + varWindows[1][i]));
						windowedHaplotypeDiversity.append(allel.moving_haplotype_diversity(haplArray, size=varWindows[1][i], start=varWindows[0][i], stop=(varWindows[0][i] + varWindows[1][i]))[0]);
					for j in range(len(windowedGarudsH)):
						windowedGarudsH[j].append(curGarudsH[j][0]);
				if debug: print("Windowed Garud's H and Haplotype Diversity Calculated");
				#Output windowed statistics to output .csv file
				for thPi, thW, thH, fwH, tajD, hapD, gH1, gH12, gH123, gH2H1, winStart in zip(windowedThetaPi[0], windowedThetaW[0], windowedThetaH, windowedFayWuH, windowedTajimasD[0], windowedHaplotypeDiversity, windowedGarudsH[0], windowedGarudsH[1], windowedGarudsH[2], windowedGarudsH[3], list(zip(*windowedTajimasD[1]))[0]):
					with open(outputCSV, 'a', encoding="utf-8") as outFile:
						outFile.write(f"{outputBegin},windowed,0,{winStart},{windowSize},thetaPi,{thPi}\n");
						outFile.write(f"{outputBegin},windowed,0,{winStart},{windowSize},thetaW,{thW}\n");
						outFile.write(f"{outputBegin},windowed,0,{winStart},{windowSize},thetaH,{thH}\n");
						outFile.write(f"{outputBegin},windowed,0,{winStart},{windowSize},fayWuH,{fwH}\n");
						outFile.write(f"{outputBegin},windowed,0,{winStart},{windowSize},tajimasD,{tajD}\n");
						outFile.write(f"{outputBegin},windowed,0,{winStart},{windowSize},haplDiv,{hapD}\n");
						outFile.write(f"{outputBegin},windowed,0,{winStart},{windowSize},gH1,{gH1}\n");
						outFile.write(f"{outputBegin},windowed,0,{winStart},{windowSize},gH12,{gH12}\n");
						outFile.write(f"{outputBegin},windowed,0,{winStart},{windowSize},gH123,{gH123}\n");
						outFile.write(f"{outputBegin},windowed,0,{winStart},{windowSize},gH2/H1,{gH2H1}\n");
				else:
					del thPi, thW, thH, fwH, tajD, hapD, gH1, gH12, gH123, gH2H1, winStart;
					print("Current Window-Step Size stats saved.");
			else:
				for thPi, thW, thH, fwH, tajD, winStart in zip(windowedThetaPi[0], windowedThetaW[0], windowedThetaH, windowedFayWuH, windowedTajimasD[0], list(zip(*windowedTajimasD[1]))[0]):
					with open(outputCSV, 'a', encoding="utf-8") as outFile:
						outFile.write(f"{outputBegin},windowed,0,{winStart},{windowSize},thetaPi,{thPi}\n");
						outFile.write(f"{outputBegin},windowed,0,{winStart},{windowSize},thetaW,{thW}\n");
						outFile.write(f"{outputBegin},windowed,0,{winStart},{windowSize},thetaH,{thH}\n");
						outFile.write(f"{outputBegin},windowed,0,{winStart},{windowSize},fayWuH,{fwH}\n");
						outFile.write(f"{outputBegin},windowed,0,{winStart},{windowSize},tajimasD,{tajD}\n");
				else:
					del thPi, thW, thH, fwH, tajD, winStart;
					print("Current Window-Step Size stats saved.");

		else:
			if debug: print("Windowed stats done with no issues.");
			del windowedThetaPi, windowedThetaW, windowedThetaH, windowedFayWuH, windowedTajimasD
			if haplotypeInfo: del windowedHaplotypeDiversity, windowedGarudsH, varWindows;
			if debug: print("Variables used in calculating windowed non-LD stats deleted.");
		if debug: print("Non LD stats computed");
	elif debug: print("Skipping non LD stats.");
	
	if debug: printTimeLapsed(startTime);

	#Do stats work for LD that requires O(n^2) time algorithms.
	if(doLD and haplotypeInfo):
		if debug: print("Beginning computations of LD stats; this may take some time.");
		genomeAlts = haplArray.to_genotypes(ploidy=1).to_n_alt(fill=-1);
		
		#Loci that are fixed in the sample cannot be used to calculate LD and must be removed.
		#The position value at the same index must also be removed.
		sampleFixedIndices = [];
		for i in range(len(genomeAlts)):
			if(sum(genomeAlts[i]) == nSamp):
				sampleFixedIndices.insert(0, i);
		positions_LD = positions;
		for index in sampleFixedIndices:
			genomeAlts = numpy.delete(genomeAlts, index, axis=0);
			del positions_LD[index];

		#Get median r^2 for entire genome.
		fullRSquared = allel.windowed_r_squared(positions_LD, genomeAlts, start=1, stop=genomeSize, size=genomeSize, step=genomeSize, fill=-1);
		
		#Write median r^2 for entire genome to output .csv file.
		with open(outputCSV, 'a', encoding="utf-8") as outFile:
			outFile.write(f"{outputBegin_Full},r2,{fullRSquared[0][0]}\n");
		
		del fullRSquared;
		if debug: print("Variables used in calculating full genome LD stats deleted.");
		
		if debug: printTimeLapsed(startTime);

		#Calculate windowed LD statistics for all window size-step size pairs.
		for windowSize, stepSize in zip(windowSizes, stepSizes):
			if debug: print(f"{windowSize}, {stepSize}");
			windowedRSquared = allel.windowed_r_squared(positions_LD, genomeAlts, size=windowSize, start=1, step=stepSize, fill=-1);
			for rSquared, winStart in zip(windowedRSquared[0], list(zip(*windowedRSquared[1]))[0]):
				with open(outputCSV, 'a', encoding="utf-8") as outFile:
					outFile.write(f"{outputBegin},windowed,0,{winStart},{windowSize},r2,{rSquared}\n");
			else:
				del rSquared, winStart;
			if debug: printTimeLapsed(startTime);
		else:
			if debug: print("Windowed LD stats done with no issues.");
			del windowSize, stepSize, windowedRSquared, genomeAlts, positions_LD;
			if debug: print("Variables used in calculating windowed LD stats deleted.");
		if debug: print("LD stats computed.");
	elif debug: print("Skipping LD stats.");

	if debug: printTimeLapsed(startTime);

#If this script is being run (and not just sourced for the methods), then variables need to be initialized and the @main() function called.
if __name__ == "__main__":
	#Initialize hard-coded variables; changes to these need to be done in this source code.
	#debug : Boolean value for whether or not debug output should be written.
	debug = True;
	#doNFiltVar : Boolean value for whether the number of filtered variants should be counted (variants that occur above an allele frequency in $filterFreqs).
	doNFiltVar = True;
	#doNfixed : Boolean value for wheter the number of fixed mutations should be counted.
	doNFixed = True;
	#doNonLD : Boolean value for wheter non linkage disequilibrium statistics should be calculated.
	doNonLD = True;
	#doLD : Boolean value for whether linkage disequilibrium statistics should be calculated. These calculations can be time intensive.
	doLD = False;
	#outputParameterColumns : Boolean value for if output file should include columns with the paramter values used in the analyzed simulation. Requires filename to be formatted in a specific way and can greatly increase the size of the output file.
	outputParameterColumns = False;
	#filterFreqs : A list of MAF (minor allele frequencies) that a variant must exceed to be counted if $doNFiltVar is being used.
	filterFreqs = [0.005, 0.01, 0.02, 0.05];
	#windowSizes : List of window sizes to use for calculating sliding windowed statistics.
	windowSizes = [100, 1000, 2000];
	#stepSizes : List of step sizes to use alongside $windowSizes in calculating sliding windowed statistics.
	stepSizes = [50, 500, 1000];
	#burnintime : The length (in ticks) that the burnin was run.
	burnintime = 25000;
	main();
