

import sys

inFile = open(sys.argv[1])
binLength = int(sys.argv[2])
outFile = open(sys.argv[3], "w")

MIN_COV = 30
MIN_VAF = 0.2


# create the data structure for bin depth, which will be dict(Chrom) -> dict(bin) -> sumCoverage; that is a dict of dicts for which the second layer of dicts points to the sum of coverage values
chromosomes = [str(a) for a in range(1, 23)]
chromosomes.append("X")

chromDict = {}
binIds = []

for thisChrom in chromosomes: # this loop will set up the data structure storing bins (dict of dicts)
        chromDict[thisChrom] = {}

        binStart = 0

        while binStart < 300000000: # no chromosome is longer than 300,000,000 bp, so we don't need bins beyond that
                chromDict[thisChrom][binStart] = 0

                if thisChrom == "1":
                        binIds.append(binStart)

                binStart += binLength


# get header 
header = inFile.readline().rstrip("\n").split("\t")

for line in inFile:
	lineList = line.rstrip("\n").split("\t")

	chrom = lineList[header.index("Chr")].replace("chr", "")
	pos = int(lineList[header.index("Pos")])

	refNormal = float(lineList[header.index("reference_normal_count")])
	refTumor = float(lineList[header.index("reference_tumor_count")])
	altNormal = float(lineList[header.index("alternative_normal_count")])
	altTumor = float(lineList[header.index("alternative_tumor_count")])

	MinN = altNormal
	TinN = altNormal + refNormal

	MinD = altTumor
	TinD = altTumor + refTumor

	if TinN > 0 and TinD > 0:
		nVaf = MinN / TinN
		tVaf = MinD / TinD
		
		if TinD >= MIN_COV and TinN >= MIN_COV:
			if nVaf >= 0.2 and tVaf >= 0.2:
				whichBin = ""

				if chrom in chromDict: # if this is a valid chromosome, find out which chromosome bin this position is in
					# find out which bin this position is in
					whichBin = binIds[pos // binLength]

					#if counter % 1000000 == 0:
					#	print("On chromosome " + str(chrom) + " position " + str(pos))
					#counter += 1
				else:
					continue # move on to next line if not a valid chromosome

				chromDict[chrom][whichBin] = str(nVaf) + ";" + str(tVaf)

# make output matrix
outMatrix = []
outHeader = ["Bin"]
outHeader.extend(chromosomes)
outMatrix.append(outHeader)

print "Regarding bins..."
print len(binIds)

for thisBin in binIds:
        outRow = []
        outRow.append("Bin" + str(thisBin))

        for thisChrom in chromosomes:
                outRow.append(chromDict[thisChrom][thisBin])

        outMatrix.append(outRow)

print("Length of outMatrix is " + str(len(outMatrix)))

# output to file
for line in outMatrix:
        outFile.write("\t".join([str(a) for a in line]) + "\n")

outFile.close()


inFile.close()
outFile.close()

