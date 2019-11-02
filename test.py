'''
A module for testing ReadDepthCalculator.

by Jack Beal
October 2019

'''

import unittest
import os
from filecmp import cmp
from ReadDepthCalculator import *

class TestReadDepthCalculator(unittest.TestCase):
    '''Test the ReadDepthCalculator's functionality.'''

    def test_add_reads(self):

        '''Test that read data can be put into and pulled out of the RDC.'''
        r1 = (10,30)
        r2 = (20,40)
        reads = [r1,r2]

        rdc = ReadDepthCalculator([],[])
        rdc.addReads(reads)

        self.assertFalse(rdc.getReads() == [])
        self.assertEqual(reads, rdc.getReads())


    def test_add_loci(self):
        '''Test that locus data can be put into and pulled out of the RDC.'''

        l1 = 5
        l2 = 15
        loci = [l1, l2]

        rdc = ReadDepthCalculator([],[])
        rdc.addLoci(loci)

        self.assertFalse(rdc.getLoci == [])
        self.assertEqual(loci, rdc.getLoci())


    def test_reads_CSV(self):
        '''Test that CSV data makes it into and out of the RDC correctly.'''

        testReads = [(10,30),(20,40)]

        # Write our test data out to a CSV that we'll then import
        with open("temp-in.csv", 'w', newline='') as outfile:
            writer = csv.writer(outfile, quoting=csv.QUOTE_MINIMAL)
            writer.writerow(["start","length"])
            for row in testReads:
                writer.writerow(row)

        # Import the data and verify that it arrived unharmed
        rdc = ReadDepthCalculator([],[])
        rdc.addReadsFromCSV("temp-in.csv")
        self.assertEqual(testReads, rdc.getReads())

        # Export the data and make sure that it matches the input data
        with open("temp-out.csv", 'w', newline='') as outfile:
            writer = csv.writer(outfile, quoting=csv.QUOTE_MINIMAL)
            writer.writerow(["start","length"])
            for row in rdc.getReads():
                writer.writerow(row)

        self.assertTrue(cmp("temp-in.csv","temp-out.csv"))

        # Clean up
        os.remove("temp-in.csv")
        os.remove("temp-out.csv")


    def test_loci_CSV(self):
        '''Test that CSV data makes it into and out of the RDC correctly.'''

        testLoci = [5,15,30]

        # Write our test data out to a CSV that we'll then import
        with open("temp-in.csv", 'w', newline='') as outfile:
            writer = csv.writer(outfile, quoting=csv.QUOTE_MINIMAL)
            writer.writerow(["position","depth"])
            for item in testLoci:
                writer.writerow([item])

        # Import the data and verify that it arrived unharmed
        rdc = ReadDepthCalculator([],[])
        rdc.addLociFromCSV("temp-in.csv")
        self.assertEqual(testLoci, rdc.getLoci())

        # Export the data and make sure that it matches the input data
        with open("temp-out.csv", 'w', newline='') as outfile:
            writer = csv.writer(outfile, quoting=csv.QUOTE_MINIMAL)
            writer.writerow(["position","depth"])
            for row in rdc.getLoci():
                writer.writerow([row])

        self.assertTrue(cmp("temp-in.csv","temp-out.csv"))

        # Clean up
        os.remove("temp-in.csv")
        os.remove("temp-out.csv")


    def test_data_preprocessing(self):
        '''Test that data gets sorted as expected.'''

        testReads = [(20,40),(10,30),(20,10)]
        sortedReads = [(10,30),(20,10),(20,40)]

        # Add the unsorted reads to the RDC and cue the preprocessing step
        rdc = ReadDepthCalculator([],[])
        rdc.addReads(testReads)
        rdc.preprocessReadData()

        # Check to see that the reads match the hand-sorted ones
        self.assertEqual(sortedReads, rdc.getReads())


    def test_depth_calculations(self):
        '''Test that the RDC is correctly calculating depths.'''

        testReads = [(10,30),(20,40)]
        testLoci = [5,15,30]
        expectedDepths = [0,1,2]

        rdc = ReadDepthCalculator(testReads, testLoci)
        rdc.calculateDepth()
        calculatedDepthsAtLoci = [rdc.getDepths()[l] for l in testLoci]
        self.assertEqual(expectedDepths, calculatedDepthsAtLoci)


if __name__ == '__main__':
    unittest.main()
