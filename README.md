# kwmtools
----

This tool is a set of command line tools to manipulate variant call format file (VCF).



## Install & Requirement
----



For the tools to run properly, you need to have Python 3.x installed.  
And you must install python modules: sys, re, math, time, scipy, mne, argparse.  
If you carry out the above, you can use this tools as it is.  



## Description
----


### Overall

To test that you can run this tools, run the following command in your terminal application.  
```
python kwmtools.py -h  
```
You should see a complete list of all the functions in this toolkit.  
The tools, which are all listed further below, are invoked as follows:  
```
python kwmtools.py -T ToolName -OPTION1 [value1] -OPTION2 [value2]...  
```
### By functions

#### VcfIndexChangeByFai

Change vcf file index based on fasta index file.  
Usage example:
```
python kwmtools.py -T VcfIndexChangeByFai -I [input.vcf] -O [output prefix] -fai [input.fai]
```
#### VcfClassifyUPR

Classify variant records of vcf file as unplaced or unlocalized region variant or other(normal) region variant.  
Usage example:
```
python kwmtools.py -T VcfClassifyUPR -I [input.vcf] -O [output prefix]
```
#### VcfClassifyMultiAlt

Classify variant records of vcf file as multi allele or single allele.  
Usage example:
```
python kwmtools.py -T VcfClassifyMultiAlt -I [input.vcf] -O [output prefix]
```
#### VcfClassifyInv

Classify variant records of vcf file as inversion or other indel based on reference fasta file.  
Usage example:
```
python kwmtools.py -T VcfClassifyInv -I [input.vcf] -O [output prefix] -R [reference.fa]
```
#### VcfClassifyDup

Classify variant records of vcf file as duplication or other indel based on reference fasta file.  
Usage example:
```
python kwmtools.py -T VcfClassifyDup -I [input.vcf] -O [output prefix] -R [reference.fa]
```
#### VcfClassifyRef

Classify variant records of vcf file as alternative or reference.  
Usage example:
```
python kwmtools.py -T VcfClassifyRef -I [input.vcf] -O [output prefix]
```
#### playBH

If you have detected variants using VarScan and then combine two vcf files using GATK CombineVariants, you can run this funcion. The function compare two seeds statistically using Benjamini-Hochberg test.  
Usage example:
```
python kwmtools.py -T playBH -I [input.vcf] -O [output  prefix] -pv [p-value]
```
I,O are required options. pv sets default value(0.05).  

#### filterVarscanVariants

If you have detected variants using VarScan, you can run this function.  
Select variant records which range from selected frequency.  
Usage example:
```
python kwmtools.py -T filterVarscanVariants -I [input.vcf] -O [output prefix] --minVarFreq [minimum Frequency percentage] --minCov [minimum coverage number] --maxCov [maximum coverage number]
```
I,O are required options. Other options set default value; minVarFreq: 0.0, minCov: 14, maxCov: 100.  

#### filterSV

If you have detected variants using Manta, you can run this function.  
Select variant records which range from selected frequency.  
Usage example:
```
python kwmtools.py -T fiterSV -I [input.vcf] -O [output prefix] --minVarFreq [minimum Frequency percentage] --minCov [minimum coverage number] --maxCov [maximum coverage number]
```
Option status is same with filterVarscanVariants.

#### VcfIntersection
  
Get intersection between vcf files.  
Usage example:
```
python kwmtools.py -T VcfIntersection -I [input.vcf] -O [output prefix] -V [compared input vcf file]
```
#### VcfExtractFieldData
  
Extract selected field data of vcf file.  
If you have annotated using snpEff, you can set as field items ANN.hoge, LOF.hoge, and NMD.hoge on field option.  
If you have detected using VarScan or Manta, you can set as field items VAR.hoge or MNT.hoge on field option.  
Usage example:
```
python kwmtools.py -T VcfExtractFieldData -I [input.vcf] -O [output.txt] --field [comma separated field name]
```
#### GetFastaSeq
  
Get sequence data from fasta file using bed format file.  
Usage example:
```
python kwmtools.py -T GetFastaSeq -I [input.bed] -O [output.fa] -R [reference.fa]
```
---
## Note  
### Made by  
Masaki Kawamoto  

### Notification  
Even if any disadvantage occurs due to this work, we do not take responsibility.    
