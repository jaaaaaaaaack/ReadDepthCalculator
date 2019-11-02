'''
ReadDepthCalculator.py
A module for calculating the read depth at provided loci of interest.

Calculating coverage of genomic sequencing data is a crucial step in validating that the sequencing has produced a signal of high enough confidence for the needs of a given application. In this script we demonstrate a basic implementation of an algorithm to calculate the genomic read depth across a set of next-generation sequencing read data.

by Jack Beal
October 2019

'''

import csv
from statistics import mean
from operator import itemgetter

class ReadDepthCalculator:
    '''A class that enables the calculation of genomic read depths at various points of interest along the genome.
    '''

    def __init__(self, reads=[], loci=[]):
        # Build homes for all of our numbers
        # This is primarily so that users may add data from several files
        self.reads = reads
        self.loci = loci
        self.depths = []


    def addReads(self, readslist):
        '''Add an existing list of reads to the calculator. Expects a list of tuples like [(pos_i, len_i)..(pos_n, len_n)]'''

        self.reads += readslist


    def addLoci(self, locilist):
        '''Add an existing list of loci of interest to the calculator. Expects a list of integers.'''

        self.loci += locilist


    def getReads(self):
        return self.reads


    def getLoci(self):
        return self.loci


    def getDepths(self):
        return self.depths


    def addReadsFromCSV(self, filename):
        '''Add CSV data of position-length pairs to the calculator. Assumes presence of a header row.'''

        # Open CSV at provided location and read in the data
        with open(filename, newline='') as csvFile:
            reader = csv.reader(csvFile, delimiter=',')
            # Lose the header and int-ify the numbers from every other row
            rows = [(int(row[0]), int(row[1])) for row in list(reader)[1:]]
        # Accumulate it into our little database
        self.addReads(rows)


    def addLociFromCSV(self, filename):
        '''Add CSV data of loci at which to calculate read depth. Assumes presence of a header row.'''

        with open(filename, newline='') as csvFile:
            # Open up our CSV
            reader = csv.reader(csvFile, delimiter=',')
            # Lose the header and int-ify the numbers from every other row
            rows = [int(row[0]) for row in list(reader)[1:]]
        self.addLoci(rows)


    def outputAllDepthsToCSV(self, filename):
        '''Output the entire genome's coverage to a CSV file. The file will be structured as integers in two columns: Position, Coverage.'''

        with open(filename, 'w', newline='') as csvFile:
            writer = csv.writer(csvFile, quoting=csv.QUOTE_MINIMAL)
            for i in range(len(self.depths)):
                writer.writerow([i, self.depths[i]])


    def populateLociCSVDepthField(self, inputFilename, outputFilename):
        '''Read loci positions from the first column of a CSV file, and output each locus with its read depth in the second column, optionally in a different file than the input.'''

        with open(inputFilename, newline='') as infile:
            reader = csv.reader(infile, delimiter=',')
            rows = [(int(row[0]), 0) for row in list(reader)[1:]]

        with open(outputFilename, 'w', newline='') as outfile:
            writer = csv.writer(outfile, quoting=csv.QUOTE_MINIMAL)
            writer.writerow(['position','coverage'])
            for row in rows:
                populated = [row[0], self.depths[row[0]]]
                writer.writerow(populated)


    def preprocessReadData(self):
        '''Preprocess read data so that it's nice and pretty for the calculator. Use sparingly!'''

        self.reads.sort(key=itemgetter(0,1))
        # Any other desirable preprocessing steps can happen in here, too.


    def __getMaxPosition(self, reads):
        '''Determines the maximum read position in a list of read data. Expects a list of tuples like [(pos_i, len_i)..(pos_n, len_n)].'''

        maximum = 0
        for pair in reads:
            if pair[0]+pair[1] > maximum:
                maximum = pair[0]+pair[1]
        return maximum


    def calculateDepth(self):
        '''Calculate the read depth across the genomic data.'''

        # Preprocess data once here, rather than many times elsewhere
        self.preprocessReadData()

        # Allocate some space in which to do our dyamic programming
        tally = [0] * (self.__getMaxPosition(self.reads)+1)
        sums = tally

        # For every position in the list, track read start/stop points
        for pair in self.reads:
            tally[pair[0]] += 1
            tally[(pair[0])+(pair[1])] -= 1

        # Calculate cumulative sum for each position in the tally
        prior = 0
        for i in range(len(tally)):
            sums[i] = tally[i] + prior
            prior = sums[i]

        # Save it all for later
        self.depths = sums


def main():
    '''A little ignition test.'''

    reads = [(10,20),(20,40),(15,15)]
    loci = [0,50,15,30]

    rdc = ReadDepthCalculator()

    rdc.addReads(reads)
    rdc.addLoci(loci)

    print("Calculating read depths...")
    rdc.calculateDepth()

    avgCoverage = mean(rdc.getDepths())
    print("Average coverage across entire genome:", round(avgCoverage, 3))

    print("Done.")

if __name__=="__main__":
    main()
