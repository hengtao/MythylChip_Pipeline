library("ChAMP")
library("getopt")
library("doParallel")
args = commandArgs(T)
command = matrix(c("import","F",1,"character","数据类型是否为β矩阵",
                   "pdfile","PD",1,"character","表型数据文件",
                   "matrixfile","M",1,"character","β矩阵文件，如果从TCGA或者GEO上下载数据，很可能只有β矩阵文件",
                   "arraytype","T",1,'character','数据类型,450k或者EPIC(850k)',
                   "directory","D",1,"character","文件目录,idat文件所在的目录",
                   "loadmethod","LM",1,"character","数据导入方法,ChAMP或者minfi,minfi属于旧方法",
                   "probecutoff","PC",1,"double","探针在样本中NA的比例，高于该比例的探针将会被删除",
                   "samplecutoff","SC",1,"double","样本中失败的探针的比例，高于该阈值的样本将会被删除",
                   "detPcut","DP",1,"double","探针被认为检测到的P值",
                   "population","P",1,"character","人群信息，默认为T",
                   "filterXY","FXY",1,"logical","是否过滤X和Y染色体探针，T or F",
                   "detPfile", "DF",1,"character","如果没有idat文件，除了提供β矩阵文件外，还应该有P值文件",
                   "analysis", "A",1,"character","分析名称",
                   "normmethod","NM",1,"character","标准化方法，如果导入方法是minfi，可选择SWAN，PBC，BMIQ，FunctionalNormaliazation;否则只能选BMIQ和PBC",
                   "cores","C",2,"integer","核数",
                   "diffexprprobeP","DEPP",1,"double","差显探针的Padj值，默认为0.05",
                   "dmrmethod","DMRM",1,"character","DMR检测方法，可选Bumphunter,DMRcate,ProbeLasso",
                   "minProbe","MP",1,"integer","DMR检测时一个region最小的探针数目，默认为7",
                   "diffexprregionP","DERP",1,"double","差显区域的Padj值，默认为0.05",
                   "dmrcores","DC",1,"integer","DMR计算核数，只对Bumphunter和DMRcate参数有效",
                   "maxGap","MG",1,"integer","DMR最大长度，大于该长度的region将会被忽略",
                   "adjPprobe","APP",1,"double","探针是否被包含在DMR中的最小P值",
                   "bpScan","BS",1,"integer","DMB的最大范围，超出此值得DMB将被忽略",
                   "compareGroup","CG",2,'character',"比较组",
                   "gseaMethod","GM",1,'character','gometh,fisher,ebayes',
                   "gseaPvalue",'GP',1,'double','显著富集的Padj阈值',
                   "controlGroup","CGP",1,'character','CNV变异检测的比较组',
                   "freqThreshold","FT",1,"double","CNV检测的阈值",
                   "sampleType","ST",1,'character',"样本类型",
                   "help","h",0,"logical", "Help information"),
                   byrow=T, ncol=5)
args = getopt(command)
if (!is.null(args$help)) {
  cat(paste(getopt(command, usage = T), "\n"))
  q()
}

arraytype <- args$arraytype

## 处理表型数据比较组信息
compareGroup <- args$compareGroup
#### 比较组之间以空格分隔
compGroups <- strsplit(compareGroup," ")
#### 比较组数目
compNums <- length(compGroups[[1]])

## 加载idat格式的数据
import = args$import
print(import)

if(!is.null(args$loadmethod)){
  method = args$loadmethod
}else{
  method = "NA"
}
probecutoff = args$probecutoff
samplecutoff = args$samplecutoff
population = as.logical(args$population)
filterXY = as.logical(args$filterXY)
directory = args$directory
detpcut = args$detPcut
if(import=="false"){
  myLoad <- champ.load(directory = directory, method = method, methValue = "B",
             autoimpute = T, filterDetP = T, ProbeCutoff = probecutoff,
             SampleCutoff = samplecutoff, population = population,
             detPcut = detpcut, filterXY = filterXY, arraytype = arraytype)
}else {
  matrixfile = args$matrixfile
  pdfile = args$pdfile
  file <- paste(directory, "/", matrixfile, sep = "")
  detPfile <- paste(directory, "/", args$detPfile, sep = "")
  matrix <- read.table(file, header = T,row.names = 1, sep = "\t", check.names = F)
  myLoad <- champ.filter(beta = matrix, pd = pdfile, depP = detPfile,
                         SampleCutoff = samplecutoff, ProbeCutoff = probecutoff,
                         detPcut = detpcut, population = population, filterXY = filterXY,
                         arraytype = arraytype)
}

options("shiny.host" = "127.0.0.1")
options("shiny.port" = 5686)
cpggui <- CpG.GUI(CpG = rownames(myLoad$beta), arraytype = arraytype)


## 质控步骤
analysis = args$analysis
analysisdir = paste(directory,"/raw_",analysis,"/",sep = "")
champ.QC(beta = myLoad$beta, pheno = myLoad$pd$Sample_Group, 
         resultsDir = analysisdir)
#### 除了使用champ.QC外，还可以使用QC.GUI来进行交互式可视化，它会展示5个图
QC.GUI(beta = myLoad$beta, pheno = myLoad$pd$Sample_Group, arraytype = arraytype)

## 对II型探针产生的β值进行标准化。标准化和数据导入方法有关系
normmethod = args$normmethod
resultdir = paste(directory, "/CHAMP_Normalization/", sep = "")
if (method == "minfi"){
  if (normmethod == "BMIQ"){
    cores = args$cores
    myNorm <- champ.norm(beta = myLoad$beta, arraytype = arraytype, cores = cores,
                         resultdir = resultdir, plotBMIQ = T)
  }else if (normmethod == "PBC" || normmethod == "SWAN" || normmethod == "FunctionalNormalize"){
    myNorm <- champ.norm(beta = myLoad$beta, arraytype = arraytype, resultsDir = resultsDir,
                         method = normmethod)
  }else{
    print("WRONG NORMALIZE METHOD,PLEASE CHECK YOUR CONFIG FILE")
    stop()
  }
}else if(method == "ShAMP" || method == "NA"){
  if(method == "BMIQ"){
    cores = args$cores
    myNorm <- champ.norm(beta = myLoad$beta, arraytype = arraytype, cores = cores,
                         resultsDir = resultdir, plotBMIQ = T)
  }else if(method == "PBC"){
    myNorm <- champ.norm(beta = myLoad$beta, arraytype = arraytype, resultsDir = resultsDir,
                         method = normmethod)
  }else{
    print("WRONG NORMALIZE METHOD,CASE YOU'RE USING SHAMP LOADING METHOD, ONLY 'PBC' or 'BMIQ' METHOD IS SUPPORTED. PLEASE CHECK YOUR CONFIG FILE")
    stop()
  }
}
write.table(myNorm, paste(directory, "/", "MyNorm.txt", sep = ""),sep="\t", quote=F)
## SVD校正,画出图形
resultdir = paste(directory, "/CHAMP_SVDimages/", sep = "")
if (method == "minfi"){
  champ.SVD(beta = myNorm, rgSet = myLoad$rgSet, pd = myLoad$pd,
            RGEffect = T, resultsDir = resultir)
}else if(method == "ChAMP" || method == "NA"){
  champ.SVD(beta = myNorm, pd = myLoad$pd, resultsDir = resultdir)
}else{
  print("SVD plotting ERROR!!!")
  print("PLEASE CHECK YOUR DATA LOADING METHOD!!!")
  stop()
}

## 校正批效应
#### logitTrans设置为T，但如果是M矩阵，该值应该设置为F
#### pd文件的Slide列必不可少
myCombat <- champ.runCombat(beta = myNorm, pd = myLoad$pd, variablename = "Sample_Group",
                            batchname = c("Slide"), logitTrans = T)
## 差显探针鉴定
diffexprprobeP = args$diffexprprobeP
myDMP <- champ.DMP(beta = myNorm, pheno = myLoad$pd$Sample_Group, adjPVal = diffexprprobeP,
                   arraytype = arraytype)
DMR.GUI(DMR = myDMP[[1]], beta = myNorm, pheno = myLoad$pd$Sample_Group, cutgroupnumber = 4)
#### 保存每组的差显结果
myDMPnames = names(myDMP)
for(i in 1:length(myDMP)){
  filename = paste(directory, "/", myDMPnames[i], ".DMPresult.txt", sep = "")
  write.table(myDMP[[1]],filename, sep ='\t', quote=F)
}

## 差显甲基化区域鉴定
dmrmethod = args$dmrmethod
minprobe = args$minProbe
diffexprregionP = args$diffexprregionP
adjPprobe = args$adjPprobe
maxGap = args$maxGap
dmrresultdir = paste(directory, "/CHAMP_ProbeLasso/", sep = "")
if (dmrmethod == "ProbeLasso"){
  phenonum <- length(levels(factor(myLoad$pd$Sample_Group)))
  if(phenonum != 2){
    print("YOU ARE USING 'ProbeLasso' method to calculate DMR...")
    print("ONLY 2 phenotypes are supported by this method...")
    stop()
  }
  myDMR <- champ.DMR(beta = myNorm, pheno = myLoad$pd$Sample_Group, arraytype = arraytype,
                   method = dmrmethod, minProbes = minprobe, adjPvalDmr = diffexprregionP,
                   maxGap = maxGap, adjPvalProbe = adjPprobe, resultDir = dmrresultdir)
}else if(dmrmethod == "Bumphunter" || dmrmethod == "DMRcate"){
  cores = args$dmrcores
  myDMR <- champ.DMR(beta = myNorm, pheno = myLoad$pd$Sample_Group, arraytype = arraytype,
                     method = dmrmethod, minProbes = minprobe, adjPvalDmr = diffexprregionP,
                     cores = dmrcores, maxGap = maxGap, adjPvalProbe = adjPprobe,
                     resultDir = dmrresultdir)
}else{
  print("MDR Calculating ERROR!!!")
  print("PLEASE CHECK YOUR DATA LOADING METHOD!!!")
  stop()
}

#### DMR可视化
if(compNums == 1){
  DMR.GUI(DMR = myDMR,beta = myNorm, pheno = myLoad$pd$Sample_Name,
          rumDMP = T, arraytype = arraytype)
}else{
  for(i in 1:compNums){
    comp = strsplit(compGroups[[1]][i],"_vs_")
    print(paste(comp, "group was plotted by shiny..."), sep = " ")
    group1 = comp[[1]][1]
    group2 = comp[[1]][2]
    DMR.GUI(DMR = myDMR,beta = myNorm, pheno = myLoad$pd$Sample_Name,
            rumDMP = T, arraytype = arraytype, compare.group = c(group1, group2))
  }
}
#### 保存DMR结果
myDMRnames = names(myDMR)
for(i in 1:length(myDMR)){
  filename = paste(directory, "/", myDMRnames[i], ".DMRresult.txt", sep = "")
  write.table(myDMR[[1]],filename, sep ='\t', quote=F)
}

## 差显block
bpSpan = args$bpSpan
cores = args$dmrcores
myBlock <- champ.Block(beta = mNorm, pheno = myLoad$pd$Sample_Group, arraytype = arraytype,
                       bpSpan = bpSpan, cores = cores)
#### 差显block的可视化
if(compNums == 1){
  Block.GUI(Block = myBlock, beta = myNorm, pheno = myLoad$pd$Sample_Group,
            runDMP = T, )
}else{
  for(i in 1:compNums){
    comp = strsplit(compGroups[[1]][i],"_vs_")
    print(paste(comp, "group was plotted by shiny...", sep = " "))
    group1 = comp[[1]][1]
    group2 = comp[[1]][2]
    Block.GUI(Block = myBlock, beta = myNorm, pheno = myLoad$pd$Sample_Group,
              runDMP = T, compare.group = c(group1, group2))
  }
}
#### 输出差显block结果
myBlocknames = names(myBlock)
for(i in 1:length(myBlock)){
  filename = paste(directory, "/", myBlocknames[i], ".Blockresult.txt", sep = "")
  write.table(myBlock[[1]],filename, sep ='\t', quote=F)
}

## DMR和DMP的GSEA富集分析
gseamethod = args$gseaMethod
gseaPvalue = args$gseaPvalue
if(gseamethod == "gometh" || gseamethod == "fisher"){
  myGSEA <- champ.GSEA(beta = myNorm, DMP = myDMP[[1]], DMR = myDMR[[1]],
                      method = gseamethod, arraytype = arraytype, adjPval = gseaPvalue) 
#}else if(gseamethod == "ebayes"){
  #myGSEA <- champ.ebayGSEA(beta = myNorm, pheno = myLoad$pd$Sample_Group, 
                           #arraytype = arraytype, adjPval = gseaPvalue)
}else{
  print("GSEA Calculating ERROR!!!")
  print("PLEASE CHECK YOUR DATA LOADING METHOD!!!")
  stop()
}
#### 输出GSEA结果
gseadmpfile = paste(directory, "/", "gsea.dmp.txt", sep = "")
gseadmrfile = paste(directory, "/", "gsea.dmr.txt", sep = "")

write.table(myGSEA$DMP, gseadmpfile, sep="\t", quote = F)
write.table(myGSEA$DMR, gseadmrfile, sep="\t", quote = F)

## 差显甲基化互作位点
resultdir <- paste(directory, "/CHAMP_EpiMod/", sep = "")
myEpiMod <- champ.EpiMod(beta = myNorm, pheno = myLoad$pd$Sample_Group, resultsDir = resultsdir,
                         arraytype = arraytype)

## 拷贝数变
controlGroup = args$controlGroup
freqThreshold = args$freqThreshold
resultdir = paste(directory, "/CHAMP_CNA/", sep = "")
myCNA <- champ.CNA(intensity = myLoad$intensity, pheno = myLoad$pd$Sample_Group,
                   control = T, controlGroup = controlGroup, freqThreshold = freqThreshold,
                   arraytype = arraytype, resultsDir = resultdir)
#### 输出CNA结果
myCNAsamplenames= names(myCNA$sampleResult)
myCNAgroupnames = names(myCNA$groupResult)
for(i in 1:length(myCNA$sampleResult)){
  filename = paste(directory, "/", myCNAsamplenames[i], ".CNAsampleresult.txt", sep = "")
  write.table(myCNA$sampleResult[[1]],filename, sep ='\t', quote=F)
}
for(i in 1:length(myCNA$groupResult)){
  filename = paste(directory, "/", myCNAgroupnames[i], ".CNAgroupresult.txt", sep = "")
  write.table(myCNA$groupResult[[1]],filename, sep ='\t', quote=F)
}

## 细胞类型杂合度，仅可用于组织类型为BLOOD的分析
cellType = args$sampleType
if(cellType == "blood"){
  myRefBase <- champ.refbase(beta = myNorm, arraytype = arraytype)
  myRefBasenames <- names(myRefBase)
  filename = paste(directory, "/", myRefBasenames[1], "cellType.CorrectedBeta.txt", sep = "")
  write.table(myRefBase$CorrectedBeta, filename, sep = "\t", quote = F)
}

