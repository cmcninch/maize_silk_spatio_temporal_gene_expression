#B73 & Mo17 Slice and Dice RNA-Seq Protocol
---
####Colton McNinch (November 3, 2018)
* The following protocol highlights the individual steps taken to analyze all RNA-Seq libraries in this study, which includes those novel to this study and the Maize B73 Stepflug and Walley expression atlases.
	
* Version four of the Maize reference genome was used [(link to paper)](https://doi.org/10.1038/nature22971). This [genome fasta file](ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/plant/Zea_mays/latest_assembly_versions/GCA_000005005.6_B73_RefGen_v4/GCA_000005005.6_B73_RefGen_v4_genomic.fna.gz) and this [gene model gff file](ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/plant/Zea_mays/latest_assembly_versions/GCA_000005005.6_B73_RefGen_v4/GCA_000005005.6_B73_RefGen_v4_genomic.gff.gz) were used. 
* The CAU Mo17 reference genome was used [(link to paper)](https://doi.org/10.1038/s41588-018-0182-0). This [genome fasta file](https://ftp.maizegdb.org/MaizeGDB/FTP/Mo17-CAU/Zm-Mo17-REFERENCE-CAU-1.0.fsa.gz) and this [gene model gff file](https://ftp.maizegdb.org/MaizeGDB/FTP/Mo17-CAU/Zm-Mo17-REFERENCE-CAU-1.0_mgdb.gff3.gz) were used. 
* All genome files used were downloaded at the above links on 11/3/2018. 

* The protocol is subdivided into four main sections:
	1. **Read Quality Control** with `FastQC` [(Documentation)] (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and `MultiQC` [(Documentation)](http://multiqc.info/)
	2. **Read Alignment** with `HISAT2` [(Documentation)] (https://ccb.jhu.edu/software/hisat2/manual.shtml)
	3. **Transcript Assembly and Counting** with `StringTie` [(Documentation)] (http://ccb.jhu.edu/software/stringtie/index.shtml)
	4. **Count Matrix Generation** with a custom perl script named `prepDE.py` 

---
#Read Quality Control
* In this section we will determine the quality of individual samples that comprise the data set. To accomplish this we will use `FastQC v0.11.7`, [(download link)](https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc) and `MultiQC v1.1.dev`, [(download link)](http://multiqc.info/). These steps should be done in the directory holding the `.fastqc` illumina sequencing files.

* We will run a single command using `FastQC` to output a `FASTQC` file for all the libraries in the data set with the following code:

```
$ fastqc *.fastq.gz
```

* Next, organize the current directory by placing all the `fastqc.html` and `fastqc.zip` files into their own directories with the following command:

```
$ mkdir fastqc_html_files fastqc_zip_files

$ mv *fastqc.html fastqc_html_files

$ mv *fastqc.zip fastqc_zip_files
```
* Next, run `MultiQC` to compile the fastqc files into one report.

```
$ multiqc fastqc_zip_files -o ./
```
---
#Read Alignment
* For the read alignment portions of this protocol we will be using `HISAT2` [(Documentation)] (https://ccb.jhu.edu/software/hisat2/manual.shtml).

* This is an extremely fast alignment tool that works by aligning to indexes of the reference genome. So the original reference genome fasta files will be used to generate these indexes and then the alignment process will work by using the resulting indexes. This should be done on all genome fasta files. To begin, make your way into the directory holding the fasta file of the genome about to be aligned to. The B73 Version 4 genome will be used as an example.

####Step 1. Build indexes of the reference genome using `hisat2-build`.

```
$ hisat2-build -f -p8 GCA_000005005.6_B73_RefGen_v4_genomic.fna B73_V4

$ mkdir genome_indexes
$ mv *.ht2 genome_indexes
```

####Step 2. Create a bash script, using a text editor, of the `hisat2` commands soon to be used and save as `read_alignment.sh`. The input following the `-x` argument should be the file path leading to the directory where the previously made genomic indexes are held in relation to your current working directory. A `/` and the basename of the genomic indexes, or `B73_V4` in this case will then follow that path. The `-U` argument will change if paired-end reads are used (see hisat2's documentation listed above for more details).
```
#!/usr/bin/bash
hisat2 --dta -p 8 -x genome_indexes/B73_V4 -U $1 -S "`basename $1 .fastq.gz`.sam" --summary-file "`basename $1 .fastq.gz`.txt"

```

####Step 3. Align the RNA-Seq reads for all libraries.
```
$ ls -1 *.fastq.gz | xargs -n 1 nohup bash read_alignment.sh
```

####Step 4. Retrieve the alignment summaries for all the libraries by extracting the information from the `.txt` files generated in the last step and creating a `<library_name>.csv` file for all the libraries. Then gather all the alignment information for the libraries into one file named `compiled_alignment_summary.csv`. This can be done by: making a `compile_alignment_summaries.sh` bash script using a text editor of the following code:
```
#!/usr/bin/bash
for File in `ls *.txt`
 	do
 		library="$(echo $File | sed 's/_.*//')"
 		total_reads="$(awk 'NR==1' $File | sed 's/r.*//')"
 		unaligned_reads="$(awk 'NR==3' $File | sed 's/(.*//')"
 		single_aligned_reads="$(awk 'NR==4' $File | sed 's/(.*//')"
 		multiple_aligned_reads="$(awk 'NR==5' $File | sed 's/(.*//')"
 		
 		echo "Library_Name,Total_Reads,Unaligned_Reads,Uniquely_Aligned_Reads,Multiple_Aligned_Reads" >> "$library.csv"
 		echo "$library" "," "$total_reads" "," "$unaligned_reads" "," "$single_aligned_reads" "," "$multiple_aligned_reads" >> "$library.csv"
 	done
 	
 echo "Library_Name,Total_Reads,Unaligned_Reads,Uniquely_Aligned_Reads,Multiple_Aligned_Reads" > alignment_summary.csv
for File in `ls *.csv`
 	do
 		tail n +2 $File >> alignment_summary.csv
 	done
 	
awk '!seen[$0]++' alignment_summary.csv > compiled_alignment_summary.csv
rm alignment_summary.csv
```

```
$ bash compile_alignment_summaries.sh
```

####Step 5. Convert `.sam` files to `.bam` files by first creating a bash script of the following code, using a text editor, and save as `sam_to_bam_convert.sh`
```
#!/usr/bin/bash
samtools sort -@ 12 -o "`basename $1 .sam`.bam" $1
```
####Step 6. Covert all `.sam` files by running this code.
```
$ ls -1 *.sam | xargs -n 1 nohup bash sam_to_bam_convert.sh
```

---
#Transcript Assembly and Counting
* In this section we will use `Stringtie` to assemble the aligned transcripts and quantify their abundance. 

####Step 1. Assemble and count transcript abundances by first creating a bash script, using a text editor, and save as `transcript_assemble_and_count.sh`. The input following the `-G` argument should be the path to the gff file of the genome being used.

```
#!/usr/bin/bash
stringtie -p 4 -M 0.5 -e -v -G GCA_000005005.6_B73_RefGen_v4_genomic.gff -o "`basename $1 .bam`.gtf" $1 -A "`basename $1 .bam`.tab"
```
####Step 2. Next execute the following command to assemble and count transcripts for all the libraries. 
```
$ ls *.bam | xargs -n 1 nohup bash transcript_assemble_and_count.sh
```

####Step 3. Generate FPKM and TPM tables for all the libraries by creating a bash script, using a text editor, and save as `FPKM_TPM_Compiler.sh`.
```
#!/usr/bin/bash
for File in `ls -1 *.tab`
    do
    Library_Name="$(echo $File | awk -F "_" '{print $1}')"
	
	tail -n +2 $File | cut -f2,8,9 |sed 's/$/\t'$Library_Name'/' | tr "\\t" "," >> FPKM_and_TPM.csv
	
	done
```

####Step 4. Execute the created bash script and add a header to the file.
```
$ bash FPKM_TPM_Compiler.sh
$ echo -e "GeneID,RPKM,TPM,Library_Name" | cat - FPKM_and_TPM.csv > FPKM_and_TPM_final.csv
```
---
# Count Matrix Generation
* In this section we will create a count matrix that contains interger count estimates, (normalized by library size), for all the libraries. This is an important step before differential expression analysis with programs such as DESeq, which require count matrices. This will be accomplished by general bash scripts and commands.

####Step 1. Create a directory for each of the libraries and then move the transcript count `.gtf` file for that particular library into the created directory (This is needed to make the script in the next step to work). 
```
$ for File in `ls *.gtf`;
 	do
		mkdir "$(echo $File | sed 's/.gtf//')"
		mv $File "$(echo $File | sed 's/.gtf//')";
	done
```

####Step 2. Create count a count matrix for all RNA-Seq libraries by first creating a bash script of the following code, using a text editor, and save as `prepDE.py`.
```
#!/usr/bin/env python2
import re, csv, sys, os, glob, warnings, itertools
from math import ceil
from optparse import OptionParser
from operator import itemgetter
#note that the gtf files in the sample folders have same # of lines, just different order(?)

parser=OptionParser(description='Generates two CSV files containing the count matrices for genes and transcripts, using the coverage values found in the output of `stringtie -e`')
parser.add_option('-i', '--input', '--in', default='ballgown', help="the parent directory of the sample sub-directories or a textfile listing the paths to GTF files [default: %default]")
parser.add_option('-g', default='gene_count_matrix.csv', help="where to output the gene count matrix [default: %default")
parser.add_option('-t', default='transcript_count_matrix.csv', help="where to output the transcript count matrix [default: %default]")
parser.add_option('-l', '--length', default=75, type='int', help="the average read length [default: %default]")
parser.add_option('-p', '--pattern', default=".", help="a regular expression that selects the sample subdirectories")
parser.add_option('-c', '--cluster', action="store_true", help="whether to cluster genes that overlap with different gene IDs, ignoring ones with geneID pattern (see below)")
parser.add_option('-s', '--string', default="MSTRG", help="if a different prefix is used for geneIDs assigned by StringTie [default: %default]")
parser.add_option('-k', '--key', default="prepG", help="if clustering, what prefix to use for geneIDs assigned by this script [default: %default]")
parser.add_option('--legend', default="legend.csv", help="if clustering, where to output the legend file mapping transcripts to assigned geneIDs [default: %default]")
(opts, args)=parser.parse_args()

samples = [] # List of tuples. If sample list, (first column, path). Else, (subdirectory name, path to gtf file in subdirectory)
if (os.path.isfile(opts.input)):
    # gtfList = True
    try:
        fin = open(opts.input, 'r')
        for line in fin:
            if line[0] != '#':
                lineLst = tuple(line.strip().split())
                if (len(lineLst) != 2):
                    print "Error: Text file with sample ID and path invalid (%s)" % (line.strip())
                    exit(1)
                if lineLst[0] in samples:
                    print "Error: Sample ID duplicated (%s)" % (lineLst[0])
                    exit(1)
                if not os.path.isfile(lineLst[1]):
                    print "Error: GTF file not found (%s)" % (lineLst[1])
                    exit(1)
                samples.append(lineLst)
    except IOError:
        print "Error: List of .gtf files, %s, doesn't exist" % (opts.input)
        exit(1)
else:
    # gtfList = False
    ## Check that opts.input directory exists
    if not os.path.isdir(opts.input):
      parser.print_help()
      print " "
      print "Error: sub-directory '%s' not found!" % (opts.input)
      sys.exit(1)

    #####
    ## Collect all samples file paths and if empty print help message and quit
    #####
    samples = [(i,glob.iglob(os.path.join(opts.input,i,"*.gtf")).next()) for i in next(os.walk(opts.input))[1] if re.search(opts.pattern,i)]

if len(samples) == 0:
  parser.print_help()
  print " "
  print "Error: no GTF files found under ./%s !" % (opts.input)
  sys.exit(1)

RE_GENE_ID=re.compile('gene_id "([^"]+)"')
RE_GENE_NAME=re.compile('gene_name "([^"]+)"')
RE_TRANSCRIPT_ID=re.compile('transcript_id "([^"]+)"')
RE_COVERAGE=re.compile('cov "([\-\+\d\.]+)"')
RE_STRING=re.compile(re.escape(opts.string))

#####
## Sort the sample names by the sample ID
#####
samples.sort()

#####
## Checks whether a given row is a transcript 
## other options: ex. exon, transcript, mRNA, 5'UTR
#####
def is_transcript(x):
  return len(x)>2 and x[2]=="transcript"

def getGeneID(s, ctg, tid):
  r=RE_GENE_ID.search(s)
  if r: return r.group(1)
  r=RE_GENE_NAME.search(s)
  if r: return ctg+'|'+r.group(1)
  return tid

def getCov(s):
  r=RE_COVERAGE.search(s)
  if r:
    v=float(r.group(1))
    if v<0.0: v=0.0
    return v
  return 0.0

def is_overlap(x,y): #NEEDS TO BE INTS!
  return x[0]<=y[1] and y[0]<=x[1]


def t_overlap(t1, t2): #from badGenes: chromosome, strand, cluster, start, end, (e1start, e1end)...
    if t1[0] != t2[0] or t1[1] != t2[1] or t1[5]<t2[4]: return False
    for i in range(6, len(t1)):
        for j in range(6, len(t2)):
            if is_overlap(t1[i], t2[j]): return True
    return False

## Average Readlength
read_len=opts.length

## Variables/Matrices to store t/g_counts
t_count_matrix, g_count_matrix=[],[]

##Get ready for clustering, stuff is once for all samples##
geneIDs={} #key=transcript, value=cluster/gene_id


## For each of the sorted sample paths
for s in samples:
    badGenes=[] #list of bad genes (just ones that aren't MSTRG)
    try:
        ## opts.input = parent directory of sample subdirectories
        ## s = sample currently iterating through
        ## os.path.join(opts.input,s,"*.gtf") path to current sample's GTF
        ## split = list of lists: [[chromosome, ...],...]

        #with open(glob.iglob(os.path.join(opts.input,s,"*.gtf")).next()) as f:
        #    split=[l.split('\t') for l in f.readlines()]
#        if not gtfList:
#            f = open(glob.iglob(os.path.join(opts.input,s[1],"*.gtf")).next())
#        else:
#            f = open(s[1])
        with open(s[1]) as f:
            split=[l.split('\t') for l in f.readlines()]

        ## i = numLine; v = corresponding i-th GTF row
        for i,v in enumerate(split):
            if is_transcript(v):
                t_id=RE_TRANSCRIPT_ID.search(v[len(v)-1]).group(1)
                try:
                  g_id=getGeneID(v[len(v)-1], v[0], t_id)
                except:
                  print "Problem at line:\n:%s\n" % (v)
                  print "i='%s', len(v)=%s" % (i, len(v));
                  sys.exit(1)
                geneIDs.setdefault(t_id, g_id)
                if not RE_STRING.match(g_id):
                    badGenes.append([v[0],v[6], t_id, g_id, min(int(v[3]),int(v[4])), max(int(v[3]),int(v[4]))]) #chromosome, strand, cluster/transcript id, start, end
                    j=i+1
                    while j<len(split) and split[j][2]=="exon":
                        badGenes[len(badGenes)-1].append((min(int(split[j][3]), int(split[j][4])), max(int(split[j][3]), int(split[j][4]))))
                        j+=1

    except StopIteration:
        warnings.warn("Didn't get a GTF in that directory. Looking in another...")

    else: #we found the "bad" genes!
        break

##THE CLUSTERING BEGINS!##
if opts.cluster and len(badGenes)>0:
    clusters=[] #lists of lists (could be sets) or something of transcripts
    badGenes.sort(key=itemgetter(3)) #sort by start coord...?
    i=0
    while i<len(badGenes): #rather un-pythonic
        temp_cluster=[badGenes[i]]

        k=0
        while k<len(temp_cluster):
            j=i+1
            while j<len(badGenes):
                if t_overlap(temp_cluster[k], badGenes[j]):
                    temp_cluster.append(badGenes[j])
                    del badGenes[j]
                else:
                    j+=1
            k+=1
        if len(temp_cluster)>1:
            clusters.append([t[2] for t in temp_cluster])
        i+=1

    print len(clusters)

    for c in clusters:
        c.sort()

    clusters.sort(key=itemgetter(0))
    legend=[]
    for u,c in enumerate(clusters):
        my_ID=opts.key+str((u+1))
        legend.append(list(itertools.chain.from_iterable([[my_ID],c]))) #my_ID, clustered transcript IDs
        for t in c:
            geneIDs[t]=my_ID
##            geneIDs[t]="|".join(c) #duct-tape transcript IDs together, disregarding ref_gene_names and things like that

    with open(opts.legend, 'w') as l_file:
        my_writer=csv.writer(l_file)
        my_writer.writerows(legend)

geneDict={} #key=gene/cluster, value=dictionary with key=sample, value=summed counts
t_dict={}
for q, s in enumerate(samples):
    print q, s[0]

    try:
        #with open(glob.iglob(os.path.join(opts.input,s,"*.gtf")).next()) as f: #grabs first .gtf file it finds inside the sample subdirectory
#        if not gtfList:
#            f = open(glob.iglob(os.path.join(opts.input,s[1],"*.gtf")).next())
#        else:
        f = open(s[1])
        
            # s = s.split('/')[len(s.split('/')) - 1].split('.gtf')[0].split('_')[0]
            # s = sample_IDs[q]

##        split=[t[:len(t)-1]+t[len(t)-1].split(";") for t in split]
##        split=[t[:len(t)-1] for t in split] #eliminate '\n' at end
##        split=[[e.lstrip() for e in t] for t in split]
        #should consider making stuff into dictionaries, maybe each split line

##            transcriptList=[]
        transcript_len=0
        for l in f:
            if l.startswith("#"):
                continue
            v=l.split('\t')
            if v[2]=="transcript":
                if transcript_len>0:
##                        transcriptList.append((g_id, t_id, int(ceil(coverage*transcript_len/read_len))))
                    t_dict.setdefault(t_id, {})
                    t_dict[t_id].setdefault(s[0], int(ceil(coverage*transcript_len/read_len)))
                t_id=RE_TRANSCRIPT_ID.search(v[len(v)-1]).group(1)
                #g_id=RE_GENE_ID.search(v[len(v)-1]).group(1)
                g_id=getGeneID(v[len(v)-1], v[0], t_id)
                #coverage=float(RE_COVERAGE.search(v[len(v)-1]).group(1))
                coverage=getCov(v[len(v)-1])
                transcript_len=0
            if v[2]=="exon":
                transcript_len+=int(v[4])-int(v[3])+1 #because end coordinates are inclusive in GTF

##            transcriptList.append((g_id, t_id, int(ceil(coverage*transcript_len/read_len))))
        t_dict.setdefault(t_id, {})
        t_dict[t_id].setdefault(s[0], int(ceil(coverage*transcript_len/read_len)))

    except StopIteration:
#        if not gtfList:
#            warnings.warn("No GTF file found in " + os.path.join(opts.input,s[1]))
#        else:
        warnings.warn("No GTF file found in " + s[1])


##        transcriptList.sort(key=lambda bla: bla[1]) #gene_id

    for i,v in t_dict.iteritems():
##        print i,v
        geneDict.setdefault(geneIDs[i],{}) #gene_id
        geneDict[geneIDs[i]].setdefault(s[0],0)
        geneDict[geneIDs[i]][s[0]]+=v[s[0]]


with open(opts.t, 'w') as csvfile:
    my_writer = csv.DictWriter(csvfile, fieldnames = ["transcript_id"] + [x for x,y in samples])
    my_writer.writerow(dict((fn,fn) for fn in my_writer.fieldnames))
    for i in t_dict:
        t_dict[i]["transcript_id"] = i
        my_writer.writerow(t_dict[i])

with open(opts.g, 'w') as csvfile:
    my_writer = csv.DictWriter(csvfile, fieldnames = ["gene_id"] + [x for x,y in samples])
##    my_writer.writerow([""]+samples)
##    my_writer.writerows(geneDict)
    my_writer.writerow(dict((fn,fn) for fn in my_writer.fieldnames))
    for i in geneDict:
        geneDict[i]["gene_id"] = i #add gene_id to row
        my_writer.writerow(geneDict[i])
```
####Step 3. Generate a count matrix by running the following line of code:
#####`StringTie_gtf_output` is the directory containing the `.gtf` files that were generated by the StringTie program. 
#####The `-l` option was set to 49 to represent the mean length of the reads.
```
$ prepDE.py -i StringTie_gtf_output -l 49  
```
####Step 4 (varies by reference genome gff files)
* In some instances the resulting expression files will not have ZmXXXXXaXXXXXX gene identifiers associated with them. They should be in the original gff file but were not transferred because the above code was unable to parse the exact identifier out of the name in the gff. The below code will obtain those identifiers from the gff file and will create a file that can be used for later use to match gene name and gene identifiers. 

`
$ grep -o '^[^#]*' GCA_000005005.6_B73_RefGen_v4_genomic.gff | cut -f9 | sed 's/\;/\t/g' | awk '$1~/ID=gene/' | cut -f1,2 | sed 's/ID\=//g' | sed 's/Name\=ZEAMMB73\_//g' >> B73_V4_Gene_Name_Table.txt
`






















 			