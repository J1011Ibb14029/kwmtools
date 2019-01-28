import sys
import re
import math
import time
import scipy.stats as stats
from mne.stats import fdr_correction
from argparse import ArgumentParser
from joblib import Parallel, delayed
import multiprocessing as MP
from multiprocessing import Pool

#vcfファイルのchromosomeのindexを希望のindexへ修正する機能.	
def VcfIndexChangeByFai(wf,rf,faif):
	faiList = []
	with open(wf+'.Changed.vcf','w') as cf:
		with open(wf+'.Skipped.vcf','w') as sf:
			with open(rf,'r') as f2:
				for line in f2:
					if line[0] == '#':#コメント部分の処理
						if line[0:6] == '#CHROM':
							for contigs in open(faif,"r"):
								data=contigs.split()
								cf.write("##contig=<ID="+str(data[0])+",length="+str(data[1])+",assembly=null>\n")
								faiList.append(data[0])
							cf.write(line)
							sf.write(line)
						elif line[1] != "#":
							cf.write(line)
							sf.write(line)
						elif line[1] == "#" and line[0:8] != "##contig":
							cf.write(line)
							sf.write(line)
						elif line[0:8] == "##contig":	sf.write(line)
					else:#vcfレコード部分の処理
						comp=line.split()
						if str(comp[0]) in faiList:	cf.write(line)
						else:
							if line[0] == 'c':#UCSC format
								if re.match("chr[1-9XY][1-9]?",str(comp[0])) and len(str(comp[0])) < 6:
									newChr=comp[0][3:]
									if newChr in faiList:	cf.write(str(newChr)+"\t"+"\t".join(map(str,comp[1:]))+"\n")
									else:	sf.write(line)
								elif comp[0] == "chrM":
									newChr="MT"
									if newChr in faiList:	cf.write(str(newChr)+"\t"+"\t".join(map(str,comp[1:]))+"\n")
									else:	sf.write(line)
								else:#Unplaced Region
									elem = comp[0].split('_')[1].replace('v','.')
									newChr='.'.join(map(str,elem))
									if newChr in faiList:	cf.write(str(newChr)+"\t"+"\t".join(map(str,comp[1:]))+"\n")
									else:	sf.write(line)
							elif re.match("[1-9XYM][0-9T]?",str(comp[0])) and len(comp[0]) < 3:#NCBI Normal
								if comp[0] == 'MT':	newChr = 'chrM'
								else:	newChr = 'chr'+str(comp[0])
								if newChr in faiList: cf.write(str(newChr)+"\t"+"\t".join(map(str,comp[1:]))+"\n")
								else:	sf.write(line)
							else:#NCBI UPR
								newChr = [elem for elem in faiList if elem.find(comp[0].replace('.','v')) != -1]
								if len(newChr) == 1:	cf.write(str(newChr[0])+"\t"+"\t".join(map(str,comp[1:]))+"\n")
								else:	sf.write(line)

#unplaced regionと通常の染色体で分類する機能.
def VcfClassifyUPR(wf,rf):
	with open(wf+'NormalChr.vcf','w') as ncf:	#normal chromosome record vcf file
		with open(wf+'UPR.vcf','w') as uprf:	#unplaced region record vcf file
			with open(rf,'r') as f1:	#input vcf file
				for row in f1:
					if row[0] == "#":
						ncf.write(row)
						uprf.write(row)
					else:
						rowData = row.split()
						if row[0] == 'c':#UCSC Format
							if re.match("[1-9XYM][0-9]?",str(rowData[0][3:])) and len(rowData[0][3:]) < 3:	ncf.write(row)
							else: uprf.write(row)
						else:
							if re.match("[1-9XYM][0-9T]?",str(rowData[0])) and len(rowData[0]) < 3:	ncf.write(row)
							else: uprf.write(row)
						
#VcfからMulti allele recordを取得する機能
def VcfClassifyMultiAlt(wf,rf):
	with open(wf+'MultiAllele.vcf','w') as maf:	#multi allele vcf file
		with open(wf+'SingleAllele.vcf','w') as saf:	#single allele vcf file
			with open(rf,'r') as f2:
				for row in f2:
					if row[0] != '#':
						rowData = row.split()
						if rowData[4].find(',') == -1: saf.write(row)
						else:	maf.write(row)
					else:
						maf.write(row)
						saf.write(row)
						
#共通変異に対してBH法を用いて検定する機能
def playBH(wf,rf,alpha):
	infoArray = []
	stayArray = []
	pValue = []
	pToDict = []
	VcfColDict = {}
	with open(wf+'Exclude.vcf','w') as otherf:
		with open(rf,"r") as f1:
			for row in f1:
				rowData = row.split()
				if rowData[0] == '#CHROM':
					posF = rowData.index('FORMAT')
					Snum = len(rowData)-(posF+1)
				if row[0] == '#':
					otherf.write(row)
					infoArray.append(row)
				else:
					StayFlag = 0
					if rowData[0] in VcfColDict:	VcfColDict[rowData[0]][rowData[1]] = row
					else:
						VcfColDict[rowData[0]] = {}
						VcfColDict[rowData[0]][rowData[1]] = row
					#データのあるテーブルの位置を取得
					Dpos = []
					for i in range(posF+1,len(rowData)):
						if rowData[i].find('./.') == -1: Dpos.append(i)
					#各SampleのADとRDを取得[AD,RD]．ただしAltが複数種類あることによりADが表記されない場合はADF＋ADRをADとして代用
					if len(Dpos) == 2:	#データのあるテーブルが2個の時のみχ二乗値計算およびBHを適用する．
						if rowData[posF].find(':AD:') != -1:
							Sample1AD = int(rowData[int(Dpos[0])].split(':')[2])
							Sample2AD = int(rowData[int(Dpos[1])].split(':')[2])
						else:
							Sample1AD = int(rowData[int(Dpos[0])].split(':')[2]) + int(rowData[int(Dpos[0])].split(':')[3])
							Sample2AD = int(rowData[int(Dpos[1])].split(':')[2]) + int(rowData[int(Dpos[1])].split(':')[3])
						Sample1Data = [Sample1AD,int(rowData[int(Dpos[0])].split(':')[-4])]
						Sample2Data = [Sample2AD,int(rowData[int(Dpos[1])].split(':')[-4])]
						crossM = [Sample1Data,Sample2Data]
						#print(row)
						Sample1VariantRatio = Sample1AD/(Sample1AD+int(rowData[int(Dpos[0])].split(':')[-4]))
						Sample2VariantRatio = Sample2AD/(Sample2AD+int(rowData[int(Dpos[1])].split(':')[-4]))
						if float(Sample1VariantRatio) < float(Sample2VariantRatio):
							Condition = "Up"
						elif float(Sample1VariantRatio) > float(Sample2VariantRatio):
							Condition = "Down"
						else:
							#Condition is "Stay"
							StayFlag = 1
							stayArray.append(row)
							del VcfColDict[rowData[0]][rowData[1]]
						if StayFlag == 0:
							x2, p, dof, expected = stats.chi2_contingency(crossM)
							pValue.append(p)
							SArray = [rowData[0],rowData[1],Condition,p]
							pToDict.append(SArray)
					else:	#除外データは都度保存
						otherf.write(row)
	#Benjamini-Hochberg法での検定（Benjamini-Yekutieli法の場合は'indep'を'negcorr'とする）
	ResultFDR,PValFDR=fdr_correction(pValue,float(alpha),'indep')
	"""
	wf1:有意かつHuh7.5.1-8で増加，wf2:有意かつHuh7.5.1-8で減少
	wf3:非有意かつHuh7.5.1-8で増加，wf4:非有意かつHuh7.5.1-8で減少
	wf5：HuH-7とHuh7.5.1-8で変化なし
	"""
	with open(wf+'CommonBHSig.vcf','w') as asf:
		with open(wf+'CommonBHNotSig.vcf','w') as ansf:
			with open(wf+'CommonBHSigInc.vcf',"w") as wf1:
				with open(wf+'CommonBHSigDec.vcf',"w") as wf2:
					with open(wf+'CommonBHNotSigInc.vcf',"w") as wf3:
						with open(wf+'CommonBHNotSigDec.vcf',"w") as wf4:
							with open(wf+'CommonBHNotSigStay.vcf',"w") as wf5:
								for row in infoArray:
									asf.write(row)
									ansf.write(row)
									wf1.write(row)
									wf2.write(row)
									wf3.write(row)
									wf4.write(row)
									wf5.write(row)
								for i in range(len(ResultFDR)):
									if ResultFDR[i] == True and str(pToDict[i][2]) == "Up":
										wf1.write(VcfColDict[pToDict[i][0]][pToDict[i][1]])
										asf.write(VcfColDict[pToDict[i][0]][pToDict[i][1]])
									elif ResultFDR[i] == True and str(pToDict[i][2]) == "Down":
										wf2.write(VcfColDict[pToDict[i][0]][pToDict[i][1]])
										asf.write(VcfColDict[pToDict[i][0]][pToDict[i][1]])
									elif ResultFDR[i] == False and str(pToDict[i][2]) == "Up":
										wf3.write(VcfColDict[pToDict[i][0]][pToDict[i][1]])
										ansf.write(VcfColDict[pToDict[i][0]][pToDict[i][1]])
									elif ResultFDR[i] == False and str(pToDict[i][2]) == "Down":
										wf4.write(VcfColDict[pToDict[i][0]][pToDict[i][1]])
										ansf.write(VcfColDict[pToDict[i][0]][pToDict[i][1]])
								for row in stayArray:
									wf5.write(row)
									ansf.write(row)
														
#MantaVCFをスクリーニングする機能．頻度とカバレッジでのスクリーニングが可能．
def filterSV(wf,rf,minVarFreq,minCov,maxCov):
	minFreq = float(minVarFreq/100)
	includeFile = open(wf+"Pass.vcf","w")
	excludeFile = open(wf+"Skipped.vcf","w")
	excludeFlag = 0
	with open(rf,"r") as f1:
		for row in f1:
			data = row.split()
			if row[0] == "#":
				if row[1] == "#":
					includeFile.writelines(str(row))
					excludeFile.writelines(str(row))
				#独自のINFOの説明をメタタグに加筆
				elif data[0] == "#CHROM":
					elementNum = len(data)-9
					for i in range(elementNum):
						d = -1-i
						excludeFile.writelines("##INFO=<ID=element"+str(i+1)+"filtered,Number=0,Type=Flag,Description=\""+str(data[d])+" is filtered\">\n")
						includeFile.writelines("##INFO=<ID=PRAF"+str(i+1)+",Number=1,Type=Float,Description=\"Allele Frequency of "+str(data[d])+" from Paired end mapping\">\n")
						includeFile.writelines("##INFO=<ID=SRAF"+str(i+1)+",Number=1,Type=Float,Description=\"Allele Frequency of "+str(data[d])+" from Split read mapping\">\n")
					excludeFile.writelines("##INFO=<ID=PRZD,Number=0,Type=Flag,Description=\"the number of paired end read is zero\">\n")
					includeFile.writelines("##INFO=<ID=PRZD,Number=0,Type=Flag,Description=\"the number of paired end read is zero\">\n")
					excludeFile.writelines("##INFO=<ID=SRZD,Number=0,Type=Flag,Description=\"the number of split end read is zero\">\n")
					includeFile.writelines(str(row))
					excludeFile.writelines(str(row))
			elif data[6] == "PASS":
				elementNum = len(data)-9
				excludeFlag = 0
				pr_qualFlag = 0
				for i in range(elementNum):
					d = -1-i
					readData = data[d]
					readDetail = readData.split(':')
					if len(readDetail) == 1:#Paired End Read Dataのみの場合の処理
						ratioElement = readDetail[0].split(',')
						if int(ratioElement[0]+ratioElement[1]) >= int(minCov) and int(ratioElement[0]+ratioElement[1]) <= int(maxCov):
							prReadRatio = int(ratioElement[1])/int(ratioElement[0]+ratioElement[1])
							if prReadRatio < float(minFreq):
								excludeFlag = 1
								data[7] = data[7]+";element"+str(i+1)+"filtered"
							else:
								data[7] = data[7]+";PRAF"+str(i+1)+"="+str(prReadRatio)
						else:
							data[7] = data[7]+";element"+str(i+1)+"filtered"
							excludeFlag = 1
					else:#Paired End, Single End　両方のデータがある場合の処理
						ratioElement = readDetail[0].split(',')
						if int(ratioElement[0]+ratioElement[1]) >= 1 and int(ratioElement[1]) > 0:
							prReadRatio = int(ratioElement[1])/int(ratioElement[0]+ratioElement[1])
							if prReadRatio >= float(minFreq) and int(ratioElement[0]+ratioElement[1]) >= int(minCov) and int(ratioElement[0]+ratioElement[1]) <= int(maxCov):
								pr_qualFlag == 1
								data[7] = data[7]+";PRAF"+str(i+1)+"="+str(prReadRatio)
							elif int(ratioElement[0]+ratioElement[1]) < int(minCov) or int(ratioElement[0]+ratioElement[1]) > int(maxCov):
								data[7] = data[7]+";PRAF"+str(i+1)+"="+str(prReadRatio)
							else:
								excludeFlag = 1
								data[7] = data[7]+";element"+str(i)+"filtered"
						else:
							excludeFlag = 1
							data[7] = data[7] +";element"+str(i+1)+"filtered"
						if excludeFlag == 0:
							ratioElement = readDetail[1].split(',')
							if int(ratioElement[0]+ratioElement[1]) >= int(minCov) and int(ratioElement[0]+ratioElement[1]) <= int(maxCov):
								srReadRatio = int(ratioElement[1])/int(ratioElement[0]+ratioElement[1])
								if srReadRatio < float(minFreq):
									excludeFlag = 1
									data[7] = data[7] +";element"+str(i+1)+"filtered"
								else:
									data[7] = data[7]+";SRAF"+str(i+1)+"="+str(srReadRatio)
							elif pr_qualFlag == 1:
								srReadRatio = int(ratioElement[1])/int(ratioElement[0]+ratioElement[1])
								data[7] = data[7]+";SRAF"+str(i+1)+"="+str(srReadRatio)
							else:
								data[7] = data[7]+";element"+str(i+1)+"filtered"
								excludeFlag = 1
								
				newRow = "\t".join(data)+"\n"
				if excludeFlag == 1:
					excludeFile.writelines(str(newRow))
				else:
					includeFile.writelines(str(newRow))
					print(str(readData))
				

	includeFile.close()
	excludeFile.close()

#VarScanで抽出したファイルから特定の割合以上の変異のみを抽出してvcfファイルに保存する機能	
def filterVarscanVariants(wf,rf,minFreq,minCov,maxCov):
	with open(wf+'Pass.vcf','w') as f1:
		with open(wf+'Skipped.vcf','w') as f3:
			with open(rf,'r') as f2:
				for row in f2:
					rowData = row.split()
					if rowData[0] == '#CHROM':
						posF = rowData.index('FORMAT')
						Snum = len(rowData)-(posF+1)
					if row[0] == '#':
						f1.write(row)
						f3.write(row)
					else:
						SkipFlag = 0
						Fpos = []
						formatElem = rowData[posF].split(':')
						for i in range(len(formatElem)):
							if formatElem[i].find('FREQ')!=-1: Fpos.append(i)
							elif formatElem[i].find('SDP')!=-1: Fpos.append(i)
						#データのあるテーブルの位置を取得
						Dpos = []
						for i in range(posF+1,len(rowData)):
							if rowData[i].find('./.') == -1: Dpos.append(i)
						for elem in Dpos:
							elemData = rowData[elem].split(':')
							if minFreq != -1 and float(elemData[Fpos[0]][:-1]) < minFreq: SkipFlag = 1
							if minCov != -1 and float(elemData[Fpos[1]]) < minCov: SkipFlag = 1
							if maxCov != -1 and int(elemData[Fpos[1]]) > maxCov: SkipFlag = 1
						if SkipFlag == 1:	f3.write(row)
						else:	f1.write(row)

#usage=1:指定区間のゲノムデータをFasta fileから取得する機能
#usage=0:Fasta fileを読み込み，Dict Fileに格納し，値を返す機能
def GetFastaSeq(wf,rf1,rf2,usage):
	FastaDict = {}	#Fasta file saving dict
	with open(rf1,"r") as fasta:	#fasta reading
		IsInit = 1
		for row in fasta:
			#メタ情報行においてこれまでの染色体情報をdictに格納し，新たなdictを作成
			if row[0] == '>':
				#１行目のメタ情報であればスキップし，それ以外であればdictに格納
				if IsInit == 0:		FastaDict[chrFlag] = nSeq
				elif IsInit == 1:	IsInit = 0
				#新たな子dictを作成
				FastaDict[row.split()[0][1:]] = {}
				chrFlag = row.split()[0][1:]
				#print(chrFlag)
				isDataExist = 0
			#ゲノム情報行では当該染色体のゲノム情報をnSeqに随時追記
			else:
				#当該染色体においてnSeqを宣言していれば追記，していなければ新たな値を代入．
				if isDataExist == 0:
					nSeq = row[:-1].upper()
					isDataExist = 1
				else:
					nSeq = nSeq + row[:-1].upper()
		FastaDict[chrFlag] = nSeq	#最後の染色体のdict格納用処理
	if usage == 0:	return FastaDict
	elif usage == 1:
		with open(wf,"w") as outFasta:
			with open(rf2,"r") as bedFile:
				for row in bedFile:
					bedData = row.split()
					outFasta.write('>'+str(bedData[0])+' deletion:Pos'+str(bedData[3])+':'+str(bedData[1])+':'+str(bedData[2])+"\n")
					outSeq = FastaDict[bedData[0]][int(bedData[1])-1:int(bedData[2])]
					for i in range(int(len(outSeq)/60)+1):
						if (i+1)*60 > len(outSeq):	outFasta.write(str(outSeq[i*60:])+"\n")
						else:	outFasta.write(str(outSeq[i*60:(i+1)*60])+"\n")
						
def VcfClassifyInv(wf,rf1,rf2):
	complement = {'A':'T','C':'G','G':'C','T':'A'}
	FastaDict = {}	#Fasta file 格納用辞書
	print("Start!")
	FastaDict = GetFastaSeq(wf,rf1,rf2,0)	#Fastaを読込
	#結果の書き込みと判定したいvcfの読込
	with open(wf+'Inv.vcf',"w") as Invf:	#inversion判定vcf
		with open(wf+'NoInv.vcf','w') as NoInvf:	#非inv判定vcf
			with open(rf2,"r") as vcff:	#判定したいvcf
				for row in vcff:
					#メタ情報を両ファイルに書き込み
					if row[0] == '#':
						Invf.write(row)
						NoInvf.write(row)
					#変異行の処理
					else:
						rowData = row.split()
						InvFlag = 0
						#唯一アーティファクトになり得るInsertionのみを判定．
						if len(rowData[3])<len(rowData[4]):
							refLen = len(rowData[3])
							#altが複数あるかの処理
							if rowData[4].find(',') != -1:
								altData = []
								altElems = rowData[4].split(',')
								#純粋ins部分のみを取得する
								for elem in altElems:	altData.append(elem[refLen:])
							else:
								#純粋ins部分のみを取得する
								altData = []
								altData.append(rowData[4][refLen:])
							for elem in altData:
								dist = len(elem)
								#純粋ins部の長さが2以上のものについて判定．
								if dist >= 2:
									start = int(rowData[1])-1+refLen
									end = int(rowData[1])+dist
									refData = FastaDict[rowData[0]][start:end]
									#Inversionかどうかの判定．refの配列を反転させた相補鎖と同じ並びになればinvと判定．
									elemInv = "".join([complement[i] for i in elem])
									if str(elemInv) == refData[::-1]:	InvFlag = 1
						#ファイル書き込み部分
						if InvFlag == 1:	Invf.write(row)
						else:	NoInvf.write(row)
	
def VcfClassifyDup(wf,rf1,rf2):
	FastaDict = {}	#Fasta file 格納用辞書
	print("Start!")
	FastaDict = GetFastaSeq(wf,rf1,rf2,0)	#Fastaを読込
	#結果の書き込みと判定したいvcfの読込
	with open(wf+'Dup.vcf',"w") as Dupf:	#inversion判定vcf
		with open(wf+'NoDup.vcf','w') as NoDupf:	#非inv判定vcf
			with open(rf2,"r") as vcff:	#判定したいvcf
				for row in vcff:
					#メタ情報を両ファイルに書き込み
					if row[0] == '#':
						Dupf.write(row)
						NoDupf.write(row)
					#変異行の処理
					else:
						rowData = row.split()
						DupFlag = 0
						#唯一アーティファクトになり得るDuplicationのみを判定．
						if len(rowData[3])<len(rowData[4]):
							refLen = len(rowData[3])
							#altが複数あるかの処理
							if rowData[4].find(',') != -1:
								altData = []
								altElems = rowData[4].split(',')
								#純粋ins部分のみを取得する
								for elem in altElems:	altData.append(elem[refLen:])
							else:
								#純粋ins部分のみを取得する
								altData = []
								altData.append(rowData[4][refLen:])
							for elem in altData:
								dist = len(elem)
								start = int(rowData[1])-1+refLen
								end = int(rowData[1])+dist
								refData = FastaDict[rowData[0]][start:end]
								#Duplicationかどうかの判定．refの配列と同じになればdupと判定．
								if elem == refData:	DupFlag = 1
						#ファイル書き込み部分
						if DupFlag == 1:	Dupf.write(row)
						else:	NoDupf.write(row)
	
def VcfClassifyRef(wf,rf):
	with open(wf+'.Mut.vcf','w') as mutf:
		with open(wf+'.WT.vcf','w') as wtf:
			with open(rf,'r') as f1:
				for row in f1:
					if row[0] == '#':
						mutf.write(row)
						wtf.write(row)
					else:
						rowData = row.split()
						if rowData[4] == '.':
							wtf.write(row)
						else:
							mutf.write(row)
				
def VcfIntersection(wf,rf1,rf2):
	remainVcfDict = {}
	with open(rf1,'r') as remainf:
		for row in remainf:
			if row[0] != '#':
				rowData = row.split()
				if rowData[0] in remainVcfDict:	remainVcfDict[rowData[0]][rowData[1]] = row
				else:
					remainVcfDict[rowData[0]] = {}
					remainVcfDict[rowData[0]][rowData[1]] = row
	with open(wf,'w') as f1:
		with open(rf2,'r') as leftf:
			for row in leftf:
				if row[0] != '#':
					rowData = row.split()
					if rowData[0] in remainVcfDict and rowData[1] in remainVcfDict[rowData[0]]: f1.write(remainVcfDict[rowData[0]][rowData[1]])
				else:	f1.write(row)					
			
#vcfレコード中の欲しいフィールドの情報のみ取得する機能
def VcfExtractFieldData(wf,rf,fields):
	#各フィールドのdictionary
	mainColumnDict = {'CHROM':0,'POS':1,'ID':2,'REF':3,'ALT':4,'QUAL':5,'FILTER':6,'INFO':7,'FORMAT':8,'VARIANT':9}
	AnnotationDict = {'ANN.ALLELE':0,'ANN.EFFECT':1,'ANN.IMPACT':2,'ANN.GENE':3,'ANN.GENEID':4,'ANN.FEATURE':5,'ANN.FEATUREID':6,'ANN.BIOTYPE':7,'ANN.RANK':8,'ANN.HGVS.c':9,'ANN.HGVS.p':10,'ANN.cDNA.pos/cDNA.length':11,'ANN.CDS.pos/CDS.length':12,'ANN.AA.pos/AA.length':13,'ANN.Distance':14,'ANN.ERRORS/WARNINGS/INFO':15,'ANN.LOF':16,'ANN.NMD':17}
	LOFDict = {'LOF.GENE':0, 'LOF.GENEID':1, 'LOF.NUMTR':2,'LOF.PERC':3}
	NMDDict = {'NMD.GENE':0, 'NMD.GENEID':1, 'NMD.NUMTR':2,'NMD.PERC':3}
	VarScanDict = {'VAR.GT':0,'VAR.ABQ':1,'VAR.AD':2,'VAR.ADF':3,'VAR.ADR':4,'VAR.DP':5,'VAR.FREQ':6,'VAR.GQ':7,'VAR.PVAL':8,'VAR:RBQ':9,'VAR:RD':10,'VAR:RDF':11,'VAR:RDR':12,'VAR:SDP':13}
	MantaDict={'MNT.READ':0}
	resultArray = []	#結果格納行列
	#上記dictから実際に抽出するデータのためのdict
	getDataDict = {'Main':-1,'INFO':-1,'ANN':-1,'LOF':-1,'NMD':-1,'VAR':-1,'MNT':-1}
	#実際に抽出するデータを引数から探索
	for field in fields:
	#フィールド引数は属性.カラム名となっているため，属性別で探索．
		fieldtype = field.split('.')
		if fieldtype[0] == 'ANN':
			if getDataDict['ANN'] == -1:	getDataDict['ANN'] = str(AnnotationDict[field])
			else:	getDataDict['ANN'] = getDataDict['ANN'] + ',' + str(AnnotationDict[field])
		elif fieldtype[0] == 'LOF':
			if getDataDict['LOF'] == -1:	getDataDict['LOF'] = str(LOFDict[field])
			else:	getDataDict['LOF'] = getDataDict['LOF'] + ',' + str(LOFDict[field])
		elif fieldtype[0] == 'NMD':
			if getDataDict['NMD'] == -1:	getDataDict['NMD'] = str(NMDDict[field])
			else:	getDataDict['NMD'] = getDataDict['NMD'] + ',' + str(NMDDict[field])
		elif fieldtype[0] == 'VAR':
			if getDataDict['VAR'] == -1:	getDataDict['VAR'] = str(VarScanDict[field])
			else:	getDataDict['VAR'] = getDataDict['VAR'] + ',' + str(VarScanDict[field])
		elif fieldtype[0] == 'MNT':
			if getDataDict['MNT'] == -1:	getDataDict['MNT'] = str(MantaDict[field])
			else:	getDataDict['MNT'] = getDataDict['MNT'] + ',' + str(MantaDict[field])
		elif len(fieldtype) == 2 and fieldtype[0] == 'INFO':
			if getDataDict['INFO'] == -1:	getDataDict['INFO'] = str(fieldtype[1])
			else:	getDataDict['INFO'] = getDataDict['INFO'] + ',' + str(fieldtype[1])
		else:
			if getDataDict['Main'] == -1:	getDataDict['Main'] = str(mainColumnDict[field])
			else:	getDataDict['Main'] = getDataDict['Main'] + ',' + str(mainColumnDict[field])
	with open(rf,"r") as f1:#ファイル読込
		for row in f1:
			if row[0] != '#':	#コメント部分は無視
				vcfDataArray = []
				mainField = row.split()
				#初期化処理
				sampleDetailField = ['NA']*999
				MantaField = ['NA']*999
				VarScanField = ['NA']*999
				LOField = ['NA']*4
				NMDField = ['NA']*4
				if len(mainField) >= 8:
					sampleDetailField = mainField[9:]#サンプルデータ取得
					AnnotationField = []
					INFOField = mainField[7].split(';')#INFOデータ取得
					#snpEffアノテーションデータ取得
					for infoD in INFOField:
						if infoD.startswith('ANN='):
							AnnotationRecord = infoD[0].split(',')
							for agr in AnnotationRecord:	AnnotationField.append(agr.split('|'))
					for infoD in INFOField:#LOF,NMDデータ取得
						if infoD.startswith('LOF='):	LOField = infoD[5:-1].split('|')
						if infoD.startswith('NMD='):	NMDField = infoD[5:-1].split('|')
				if row.find('Manta') != -1:#Manta Variant Data
					MantaField = sampleDetailField
				else:	#VarScan Variant Data
					#空データ以外を取得
					VarScanField = [sdf.split(':') for sdf in sampleDetailField if sdf != './.']
				annDataArray = []
				annDataPos = []
				#フィールド引数のデータをレコードから取得
				for k,v in getDataDict.items():
					if v != -1:	#要求のある属性のデータのみを取得
						getFieldNum = v.split(',')
						for n in getFieldNum:#属性別の取得処理
							if k == 'Main':	vcfDataArray.append(mainField[int(n)])
							if k == 'INFO':	vcfDataArray.append([infoD[infoD.find('=')+1:] for infoD in INFOField if infoD.startswith(n) is True][0])
							#if k == 'INFO':	
							if k == 'ANN':
								#複数個のアノテーションデータから取得
								annResult = [annD[int(n)] for annD in AnnotationField]
								for i in range(len(annResult)):#空のセルにはNAを代入
									if len(annResult[i]) == 0: annResult[i] = 'NA'
								#後の処理のためにデータと位置を保留し，ダミーデータをvcfDataArrayに代入
								annDataArray.append(annResult)
								annDataPos.append(len(vcfDataArray))
								vcfDataArray.append(len(vcfDataArray))
							if k == 'LOF':	vcfDataArray.append(LOField[int(n)])
							if k == 'NMD':	vcfDataArray.append(NMDField[int(n)])
							if k == 'MNT':	vcfDataArray.append(MantaField[int(n)])
							if k == 'VAR':
								#ADが無い場合ことを想定した処理
								varArray = []
								nowKey = [k2 for k2,v2 in VarScanDict.items() if v2 == int(n)][0]
								nowFormat = mainField[8].split(':')
								for i in range(len(nowFormat)):
									if nowFormat[i].find(nowKey.split('.')[1]) != -1 and len(nowFormat[i]) == len(nowKey.split('.')[1]):
										for vsD in VarScanField:	varArray.append(vsD[i])
								vcfDataArray += varArray
				if len(annDataArray) >= 1 and len(annDataPos) >= 1:
					annLength = max([len(annd) for annd in annDataArray])
					for i in range(annLength):#アノテーションデータが複数個の場合にそれぞれ行を分けて出力
						for j in range(len(annDataPos)):
							vcfDataArray[annDataPos[j]] = annDataArray[j][i]
						newArray = vcfDataArray[:]#別のidとして識別させるための対策
						resultArray.append(newArray)
				else:	resultArray.append(vcfDataArray)
	outputTitle = []#1行目に対応するフィールド名を出力するための処理
	for gKey,gValue in getDataDict.items():
		if gValue != -1: values = gValue.split(',')
		if gKey == 'ANN' and gValue != -1:
			outputTitle += [tk for tk,tv in AnnotationDict.items() if str(tv) in values]
		elif gKey == 'LOF' and gValue != -1:
			outputTitle += [tk for tk,tv in LOFDict.items() if str(tv) in values]
		elif gKey == 'NMD' and gValue != -1:
			outputTitle += [tk for tk,tv in NMDDict.items() if str(tv) in values]
		elif gKey == 'VAR' and gValue != -1:
			for i in range(len(VarScanField)):#データのあるサンプル数分繰り返す．
				outputTitle += [tk+".S"+str(i+1) for tk,tv in VarScanDict.items() if str(tv) in values]
		elif gKey == 'MNT' and gValue != -1:
			outputTitle += [tk for tk,tv in MantaDict.items() if str(tv) in values]
		elif gKey == 'INFO' and gValue != -1:
			for tv in values:
				tmpTitle = str(gKey)+'.'+str(tv)
				outputTitle.append(tmpTitle)
		elif gKey == 'Main' and gValue != -1:
			outputTitle += [tk for tk,tv in mainColumnDict.items() if str(tv) in values]
	with open(wf,'w') as f2:#ファイル書込
		f2.write("\t".join(outputTitle)+"\n")
		for result in resultArray:
			f2.write("\t".join(result)+"\n")
			
def useTool(args):
	I = args.input
	O = args.output
	PV = args.pv
	R = args.R
	Fai = args.fai
	minFreq = args.minVarFreq
	minCov = args.minCov
	maxCov = args.maxCov
	V = args.variant
	fields = str(args.field).split(',')
	toolDict = {'VcfClassifyUPR':1,'VcfClassifyMultiAlt':2,'playBH':3,'VcfIndexChangeByFai':4,'filterVarscanVariants':5,'filterSV':6,'VcfClassifyInv':7,'VcfClassifyDup':8, 'VcfClassifyRef':9, 'VcfIntersection':10,'VcfExtractFieldData':11,'GetFastaSeq':12}
	if toolDict[args.T]==1:
		VcfClassifyUPR(O,I)
	elif toolDict[args.T]==2:
		VcfClassifyMultiAlt(O,I)
	elif toolDict[args.T]==3:
		playBH(O,I,PV)
	elif toolDict[args.T]==4:
		VcfIndexChangeByFai(O,I,Fai)
	elif toolDict[args.T]==5:
		filterVarscanVariants(O,I,minFreq,minCov,maxCov)
	elif toolDict[args.T]==6:
		filterSV(O,I,minFreq,minCov,maxCov)
	elif toolDict[args.T]==7:
		VcfClassifyInv(O,R,I)
	elif toolDict[args.T]==8:
		VcfClassifyDup(O,R,I)
	elif toolDict[args.T]==9:
		VcfClassifyRef(O,I)
	elif toolDict[args.T]==10:
		VcfIntersection(O,I,V)
	elif toolDict[args.T]==11:
		VcfExtractFieldData(O,I,fields)
	elif toolDict[args.T]==12:
		GetFastaSeq(O,R,I,1)
	
def getArgs():
	parser = ArgumentParser()
	parser.add_argument('-T', help="Tool name to use.   Tool List(option): VcfIndexChangeByFai(O,I,fai), VcfClassifyUPR(O,I), VcfClassifyMultiAlt(O,I), VcfClassifyInv(O,I,R), VcfClassifyDup(O,I,R), VcfClassifyRef(O,I) playBH(O,I,pv), filterVarscanVariants(O,I,minFreq,minCov,maxCov), filterSV(O,I,minFreq,minCov,maxCov), VcfIntersection(O,I,V), VcfExtractFieldData(O,I,fields), GetFastaSeq(O,R,I)",metavar='ToolName', required=True)
	'''
	VcfExtractFieldData:extract selected field data of vcf file.
	VcfIndexChangeByFai:change input vcf file index based on input fasta index file.
	VcfIntersection:get intersection between vcf files.
	VcfClassifyUPR:classify variant records of vcf file as unplaced region variant or normal region variant.
	VcfClassifyMultiAlt:classify variant records of vcf file as multi allele or single allele.
	VcfClassifyInv:classify variant records of vcf file as inversion or normal indel.
	VcfClassifyDup:classify variant records of vcf file as duplication or normal indel.
	VcfClassifyRef:classify records of vcf file as variant or reference.
	playBH:play Benjamini-Hochberg test.
	GetFastaSeq:get sequence data from fasta file using bed format file:
	filterVarscanVariants: select varscan variant records which range from selected frequency.
	filterSV:select manta variant records which range from selected frequency.
	'''
	parser.add_argument('-I', '--input', help='main input variant file path', required=True, metavar='[input.vcf]')
	parser.add_argument('-O', '--output', help='output file path', required=True, metavar='[output.vcf]')
	parser.add_argument('-V', '--variant', help='sub input variant vcf file path', metavar='[variant.vcf]')
	parser.add_argument('-R', help='Fasta file path',metavar='[Ref.fa]')
	parser.add_argument('-fai', help='Fasta index file path',metavar='[Ref.fa.fai]')
	parser.add_argument('-pv', help='P value', type=float, default=0.05, metavar='P-value [0-1]')
	parser.add_argument('--minVarFreq', help='minimum variant frequency. default=0.0', type=float, default=0.0, metavar='Frequency [0-100]')
	parser.add_argument('--minCov', help='minimum coverage. default=14', type=int, default=14, metavar='Number [0-]')
	parser.add_argument('--maxCov', help='maximum coverage. default=100', type=int, default=100, metavar='Number [0-]')
	parser.add_argument('--field', help='Extract fields name (comma determined)',metavar='Field1,Field2,...')
	
	
	return parser.parse_args()

if __name__ == "__main__":
	t1 = time.time()
	args = getArgs()
	useTool(args)
	elapsed = time.time()-t1
	print(elapsed)
	print('End.')