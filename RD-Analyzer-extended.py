## Copyright (C) 2016 Xia Eryu (xiaeryu@u.nus.edu).
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 3
## of the License, or (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program; if not, see
## http://www.opensource.org/licenses/gpl-3.0.html

## --------------------------------
## Please report bugs to:
## xiaeryu@u.nus.edu


import os
import sys
import gzip
import subprocess
from optparse import OptionParser


######################################################
## Options and arguments
######################################################
usage = "usage: %prog [options] REF.FASTA FASTQ_1 FASTQ_2(optional)"
parser = OptionParser(usage=usage, version="%prog 1.0")

# Output arguments
parser.add_option("-O","--outdir",action="store",type="string",dest="outdir",default=".",help="output directory [Default: running directory]")
parser.add_option("-o","--output",action="store",type="string",dest="output",default="RD-Analyzer",help="basename of output files [Default: RD-Analyzer]")

# Debug mode
parser.add_option("-d","--debug",action="store_true",dest="debug",help="enable debug mode, keeping all intermediate files")

# Parse options
(options, args) = parser.parse_args()

outdir = options.outdir
output = options.output
debug = options.debug


#########################################################
## Input check
#########################################################

## Check input fastq files
narg = len(args)
if narg < 2:
    print "usage: python " + os.path.realpath(__file__) + " [options] REF.FASTA FASTQ_1 FASTQ_2(optional)"
    sys.exit()
    
reference = args[0] # Input reference sequences
if not os.path.isfile(reference):
    print "Invalid Input reference file!"
    sys.exit()

input1 = args[1]    # Input fastq file 1
if not os.path.isfile(input1):
    print "Invalid FASTQ 1 file!"
    sys.exit()

if narg == 3:
    input2 = args[2]    # Input fastq file 2
    if not os.path.isfile(input2):
        print "Invalid FASTQ_2 file!"
        sys.exit()


## Check output options
outprefix = outdir+"/"+output
if os.path.isfile(outprefix+'.sam') or os.path.isfile(outprefix+'.bam') or os.path.isfile(outprefix+'.sort.bam') or os.path.isfile(outprefix+'.depth') or os.path.isfile(outprefix+'.result'):
    print "Files exist in the output directory with the output prefix specified, refuse to overwrite..."
    sys.exit()

if not os.path.isdir(outdir):
    subprocess.call(["mkdir", outdir])


############################################################
## Main class
############################################################
class Main:
    def isFloat(self, num):
        try:
            float(num)
        except ValueError:
            return False
        return True

    def readRef(self, ref):
        '''
        This dunction reads in the reference.fasta file and returns the information.
        
        @param ref: the input reference.fasta file.
        @return {name:[depthCut, covCut, info]}
        '''
        info = {}
        trace = ""

        inH = open(ref)
        for line in inH:
            line = line.strip('\n')
            if line.startswith('>'):
                line = line.lstrip('>')
                tmp = line.split('-')
                trace = tmp[0]
                info[trace] = [0.09, 0.5, tmp[3]]
                if tmp[1].lower() != 'default' and self.isFloat(tmp[1]):
                    info[trace][0] = float(tmp[1])
                if tmp[2].lower() != 'default' and self.isFloat(tmp[2]) and float(tmp[2]) <= 1:
                    info[trace][1] = float(tmp[2])

        return info

    def calcThroughput(self, infile):
        '''
        This function calculates the number of bases of a fastq file
        '''
        if infile.endswith(".gz"):
            in_handle = gzip.open(infile, 'rb')
        else:
            in_handle = open(infile)

        total = 0
        count = 0
        for line in in_handle:
            if count % 4 == 1:
                line = line.strip('\n')
                total += len(line)
            count = (count + 1) % 4
        in_handle.close()
        return total

    def recLength(self, inbam):
        '''
        This function extracts the length of the reference sequences from bam file
        '''
        output = {}
        process = subprocess.Popen(["samtools", "view", "-H", inbam], stdout=subprocess.PIPE)
        for line in iter(process.stdout.readline, ''):
            if line.startswith('@SQ'):
                line = line.strip('\n')
                tmp = line.split('\t')
                # @SQ     SN:RD711_2      LN:888
                header = (tmp[1].split(':')[1]).split('-')[0]
                output[header] = int(tmp[2].split(':')[1])
        process.terminate()
        return output

    def median(self, inarr, totalNum):
        '''
        This function calculates the median of a list
        '''
        tmp = inarr[:]
        if len(tmp) < totalNum:
            for i in range(len(tmp),totalNum):
                tmp.append(0)

        sorts = sorted(tmp)
        l = len(sorts)
        if l % 2 == 0:
            return (int(sorts[l/2]) + int(sorts[l/2-1])) / 2
        return sorts[l / 2]

    def calcStats(self, inarr, cutoff, totalNum):
        '''
        This function returns the [max, min, median, #lager_than_cutoff]
        '''
        if len(inarr) == 0:
            return [0, 0, 0, 0]

        least = int(inarr[0])
        most = int(inarr[0])
        pf = 0

        if len(inarr) < totalNum:
            least = 0

        for item in inarr:
            item = int(item)
            if item < least:
                least = item
            if item > most:
                most = item
            if item > cutoff:
                pf += 1
        return [most, least, self.median(inarr, totalNum), pf]


############################################################
## Code starts here
############################################################
if __name__ == "__main__":
    t = Main()
    # Index reference file
    try:
        subprocess.call(["bwa", "index", reference])
    except OSError:
        print "Please check the installation of BWA."
        raise

    try:
        subprocess.call(["samtools", "faidx", reference])
    except OSError:
        print "Please check the installation of samtools."
        raise

    # BWA alignment
    if narg == 2:
        os.system("bwa mem -R '@RG\\tID:RD\\tSM:RD\\tLB:RD\\tPL:Illumina' %s %s > %s" % ( reference, input1, outprefix+".sam"))

    if narg == 3:
        os.system("bwa mem -R '@RG\\tID:RD\\tSM:RD\\tLB:RD\\tPL:Illumina' %s %s %s > %s" % ( reference, input1, input2, outprefix+".sam"))

    # Sam to Bam
    os.system("samtools view -bS %s -o %s" % (outprefix+".sam", outprefix+".bam"))

    # Sort Bam file
    os.system("samtools sort %s %s" % (outprefix+".bam", outprefix+".sort"))
    os.system("samtools index %s" % outprefix+".sort.bam")

    # Calculate depth
    os.system("samtools depth %s > %s" % (outprefix+".sort.bam", outprefix+".depth"))

    # Throughput calculation:
    throughput = t.calcThroughput(input1)
    if narg == 2:
        throughput += t.calcThroughput(input2)
    depth = throughput *1.0 / 4500000

    # RDs included in the analysis
    info = t.readRef(reference)
    RDs = info.keys()

    # Record length of the reference sequences
    RD_length = t.recLength(outprefix+".sort.bam")

    # Summary of the depth file
    storage = {}
    for RD in RDs:
        storage[RD] = []

    depthH = open(outprefix+".depth")
    for line in depthH:
        line = line.strip('\n')
        tmp = line.split('\t')
        header = tmp[0].split('-')[0]
        if header in storage:
            storage[header].append(tmp[2])
    depthH.close()

    # Result matrix
    resultMat = []
    for RD in RDs:
        outStats = t.calcStats(storage[RD], info[RD][0]*depth, RD_length[RD])
        
        percentage = outStats[3]/1.0/RD_length[RD]
        prediction = 'A'
        if percentage > info[RD][1]:
            prediction = 'P'
            resultMat.append([RD, outStats[0], outStats[1], outStats[2], outStats[3], RD_length[RD], "%.2f" % percentage, prediction])
        else:
            resultMat.append([RD, outStats[0], outStats[1], outStats[2], outStats[3], RD_length[RD], "%.2f" % percentage, prediction, info[RD][2]])

    # Format output file
    outH = open(outprefix+".result", 'w')
    outH.write("# Input: " + input1 + '\n')
    if narg == 3:
        outH.write("# Input: " + input2 + '\n')
    outH.write("# Number of bases in the input file(s): %d\n" % throughput)
    outH.write("# Estimated read depth (#Bases/(4.5Mbp)): %.2f\n\n" % depth)
    outH.write("# RD_name\tMaximum_depth\tMinimum_depth\tMedian_depth\t#pos_pass_cutoff\t#total_pos\t%pos_pass_cutoff\tPrediction\tInformation\n")

    for i in range(len(resultMat)):
        outH.write(str(resultMat[i][0]))
        for j in range(1, len(resultMat[i])):
             outH.write("\t%s" % str(resultMat[i][j]))
        outH.write('\n')

    # Cleaning up
    if not debug:
        post = ['.sam','.bam','.sort.bam','.sort.bam.bai']
        for i in post:
            os.remove("%s%s" % (outprefix,i))
