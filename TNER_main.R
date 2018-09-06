#################################################################
# TNER: Tri-Nucleotide Error Reducer
# Build cfDNA background polish model
# Using tri-nucleotide context data in a Bayesian model
# Last updated: 08/20/2018 
# Contact: shibing.deng {at} pfizer.com & tao.xie {at} pfizer.com or xietao2000 {at} gmail.com
# TNER publication: BMC Bioinformatics (in revision)
# Pfizer Early Clinical Development Biostatistics 
# Pfizer Early Oncology Development and Clinical Research 
# Copyright (c) 2018 Pfizer Inc.  
# All rights reserved.
#################################################################
#################################################################
# Terms of Use for TNER (Tri-Nucleotide Error Reducer)
# Pfizer has attempted to provide accurate and current information on TNER. However, Pfizer makes no representations or warranties that the information contained or published on TNER will be suitable for your specific purposes or for any other purposes. You agree to indemnify, defend, and hold Pfizer harmless from all claims, damages, liabilities and expenses, including without limitation reasonable attorney's fees and costs, whether or not a lawsuit or other proceeding is filed, that in any way arise out of or relate to your use of TNER.
# ALL INFORMATION AND DATA PROVIDED ON TNER IS PROVIDED "AS-IS" WITHOUT WARRANTY OF ANY KIND. PFIZER MAKES NO REPRESENTATIONS OR WARRANTIES CONCERNING TNER OR ANY OTHER MATTER WHATSOEVER, INCLUDING WITHOUT LIMITATION ANY EXPRESS, IMPLIED OR STATUTORY WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NON-INFRINGEMENT OF THIRD PARTY RIGHTS, TITLE, ACCURACY, COMPLETENESS OR ARISING OUT OF COURSE OF CONDUCT OR TRADE CUSTOM OR USAGE, AND DISCLAIMS ALL SUCH EXPRESS, IMPLIED OR STATUTORY WARRANTIES. PFIZER MAKES NO WARRANTY OR REPRESENTATION THAT YOUR USE OF TNER WILL NOT INFRINGE UPON THE INTELLECTUAL PROPERTY OR OTHER RIGHTS OF ANY THIRD PARTY. FURTHER, PFIZER SHALL NOT BE LIABLE IN ANY MANNER WHATSOEVER FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, CONSEQUENTIAL OR EXEMPLARY DAMAGES ARISING OUT OF OR IN ANY WAY RELATED TO TNER, THE USE OF, OR INABILITY TO USE, ANY OF THE INFORMATION OR DATA CONTAINED OR REFERENCED IN TNER, OR ANY OTHER MATTER. THE FOREGOING EXCLUSIONS AND LIMITATIONS SHALL APPLY TO ALL CLAIMS AND ACTIONS OF ANY KIND AND ON ANY THEORY OF LIABILITY, WHETHER BASED ON CONTRACT, TORT OR ANY OTHER GROUNDS, AND REGARDLESS OF WHETHER A PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES, AND NOTWITHSTANDING ANY FAILURE OF ESSENTIAL PURPOSE OF ANY LIMITED REMEDY. BY USING TNER, YOU FURTHER AGREE THAT EACH WARRANTY DISCLAIMER, EXCLUSION OF DAMAGES OR OTHER LIMITATION OF LIABILITY HEREIN IS INTENDED TO BE SEVERABLE AND INDEPENDENT OF THE OTHER CLAUSES OR SENTENCES BECAUSE THEY EACH REPRESENT SEPARATE ELEMENTS OF RISK ALLOCATION BETWEEN THE PARTIES.
# TNER may contain information from third party websites, which may or may not be marked with the name of the source. Such information does not necessarily represent the views or opinions of Pfizer, and Pfizer shall have no responsibility whatsoever for such information. All information from third party websites are the sole responsibility of the person or entity that provides and/or maintains such website. As a TNER user, you are solely responsible for any information that you display, generate, transmit or transfer while using TNER, and for the consequences of such actions.
# Should any TNER user provide general, scientific or other feedback information, whether in the form of questions, comments, or suggestions to Pfizer, regarding the content of Pfizer's website or otherwise, such information shall not be deemed to be confidential or proprietary to you or to any other party. Pfizer shall have no obligation of any kind with respect to such information and Pfizer shall have the right, without limitation, to reproduce, use, disclose, merge, display, make derivatives of and distribute such information to others. Pfizer shall also have the right, without limitation, to use and exploit such information, including ideas, concepts, know-how, inventions, techniques or other materials contained in such information for any purpose whatsoever, including but not limited to, making changes or improvements to TNER and/or developing, manufacturing, marketing, selling or distributing products and technologies incorporating such information. However, you agree that Pfizer has no obligation whatsoever to respond to your comments or to change or correct any information on TNER based on your comments.
# Pfizer reserves the right to alter the content of TNER in any way, at any time and for any reason, with or without prior notice to you, and Pfizer will not be liable in any way for possible consequences of such changes or for inaccuracies, typographical errors or omissions in the contents hereof. Nothing contained herein shall be construed as conferring by implication, estoppel or otherwise any license or right under any patent, trademark or other intellectual property of Pfizer or any third party. Except as expressly provided above, nothing contained herein shall be construed as conferring any right or license under any Pfizer copyrights.
# Pfizer also reserves the right to modify these Terms of Use at any time and for any reason, with or without prior notice to you. You should always review the most current Terms of Use herein before using this website. By using TNER, you agree to the current Terms of Use posted on this site. You also agree that these Terms of Use constitute the entire agreement between you and Pfizer regarding the subject matter hereof and supersede all prior and contemporaneous understandings, oral or written, regarding such subject matter. In addition, if any provision of these Terms of Use is found by a court of competent jurisdiction to be invalid, void or unenforceable, the remaining provisions shall remain in full force and effect, and the affected provision shall be revised so as to reflect the original intent of the parties hereunder to the maximum extent permitted by applicable law.
# BY USING TNER, YOU ACKNOWLEDGE AND AGREE THAT YOU HAVE READ, UNDERSTOOD AND AGREE TO ALL OF THESE TERMS OF USE.
#################################################################

#################################################################
#################################################################
# Command line usage:
#################################################################
# Rscript TNER_R_script.R TNER_example_input_file.txt hs_ave_bg_error.csv hs_depth.csv
#################################################################
# Input:
# There can be up to three input parameters, the first one is required: the input file name (.csv format)
# The 2nd one is the average background error rate file calculated from healthy subjects, 
# which is assumed to be "hs_ave_bg_error.csv", can be generated via "Create_background_error_rate.R".
# The third one is the coverage depth. It can be a file name or a number which will be
# used for all bases. The default is assumed to be "hs_depth.csv", can also be generated via "Create_background_error_rate.R".
#################################################################
# Output:
# For each input file (or sample) the script will create an error polished frequency file, 
#   a file with detected mutants and a file with all noise (background error). 
#################################################################

args = commandArgs(trailingOnly=TRUE)  # 1-3 arguments from input

work.dir="."  #working directory is where the input and output file stored

# Some optional parameters to tune
input.alpha  = 0.01    #Polish strength parameter
input.filter = FALSE   #turn on the filter to filter out those bases with abnormal data

# processing the input
#####################
tnt.rate.file=tnt.depth.file=NULL
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>1) {
  # default output file
  tnt.rate.file=args[2]  # Average bkgrd rate file
  if (length(args)>2) tnt.depth.file=args[3]  # 3nt bkgrd rate file
}

input.file=args[1]  # first argument is input file
#the input file can be a single file (.txt, .freq, etc) or a list of files saved in a csv file
if (length(grep('csv',input.file))>0) input.pileup.file=unlist(read.csv(input.file,header=F,as.is=T)) else
    input.pileup.file=input.file  


if (is.null(tnt.rate.file)) tnt.rate.file="hs_ave_bg_error.csv"
if (is.null(tnt.depth.file)) tnt.depth.file="hs_depth.csv" else 
  if (length(grep("hs_depth",tolower(tnt.depth.file)))==0) { tnt.depth=as.numeric(tnt.depth.file)
    tnt.depth.file=NULL
  if (is.na(tnt.depth)) stop("coverage depth input value invalid: need to be a number or an R data file xxx.Rdata)")}

# end of processing input


s=c("A","C","T","G") # four nucleotide letter


#read in the average background error and depth file
######################################################
hs.bkg.ave = read.csv(tnt.rate.file,row.names = 1)
if (!is.null(tnt.depth.file)) hs.depth=read.csv(tnt.depth.file,row.names = 1)


#create a Tri-nucleotide (tnt) code for each position
pileup1.data=read.table(file=input.pileup.file[1],header=T,check.names = F,as.is=T)
ref=toupper(pileup1.data$REF)
n=length(ref)
tnt=c("",paste0(ref[-((n-1):n)],ref[-c(1,n)],ref[-(1:2)]),"")
# remove the tnt code for the position at sequence gaps
x=pileup1.data$POS
x2=diff(x)
tnt[c(FALSE,x2!=1)]="" #remove tri-nucleutide after gap
tnt[c(x2!=1,FALSE)]="" #remove tri-nucleutide before gap

## Emperical Bayesian - a shrinkage approach
##################################################
#calculate the tri-nucleotide average and SD
tnt.bg.ave=aggregate(hs.bkg.ave,list(TNT=tnt),function(x) mean(x[x<.1],na.rm=T))
tnt.bg.sd=aggregate(hs.bkg.ave,list(TNT=tnt),function(x) sd(x[x<0.1],na.rm=T))
# match the tnt average to the input data
loc.tnt.ave=tnt.bg.ave[match(tnt,tnt.bg.ave$TNT),-1]
loc.tnt.sd=tnt.bg.sd[match(tnt,tnt.bg.sd$TNT),-1]

w=loc.tnt.ave/(loc.tnt.ave+hs.bkg.ave)   #shrinkage weight parameter 
w[is.na(w)]=1

# average mutation error rate
shr.bg.ave=w*loc.tnt.ave+(1-w)*hs.bkg.ave

#average depth matrix
if (exists("hs.depth")) depth.mat=matrix(rep(apply(hs.depth,1,mean),4),ncol=4) else
  depth.mat=matrix(tnt.depth,nrow=nrow(pileup1.data),ncol=4)


# function to define BMER and call mutation for input file
polish.pred.3nt=function(bgrd.ave=shr.bg.ave,           #average mutation rate data matrix used to define the background
                         polish.method="Binomial",      # Background polish method: Poisson or Binomial
                         bgrd.depth=depth.mat,          #average seq depth for background
                         pred.sample=file.path(data.dir,input.pileup.file[1]), #sample file name to be polished, should be a tab delimited pileup text file
                         alpha=input.alpha,             #alpha = significance level, control polishing strength
                         filter=input.filter,           #filter out those bases based on filtering criteria
                         polyclone=TRUE,                # Output polyclone results, if FALSE, then only the base with the highest mutation rate is outputted
                         out.dir=NULL)                  #output directory. If null, polished data will be outputed to the same directory of input data with bgrm added to the file names    
{ cat(paste("Polishing",pred.sample,"\n"))
  pileup.data=read.table(file=pred.sample,header=T,check.names = F,as.is=T) #read in the data
  out.pileup=pileup.data #save a copy for output
  #below is not a critical step. It changes the column name of forward strand to cap letter: A+ to A; and 
  #reserse strand to lower case letter A+ to a.
  if (length(grep("\\+",names(pileup.data)))>0) names(pileup.data)=gsub("\\+","",names(pileup.data))
  if (length(grep("-",names(pileup.data)))>0) {
    rev.base.col=grep("-",names(pileup.data))  # data column has reverse string pile up
    names(pileup.data)[rev.base.col]=gsub("-","",tolower(names(pileup.data)[rev.base.col]))
  }
  #change all reference letter acgt to upper case
  pileup.data$REF=toupper(pileup.data$REF)
  
  fs=match(s,names(pileup.data))  # columns with forward strand
  bs=match(tolower(s),names(pileup.data)) # columns with reverse strand
  
  
  pileup.data$fs.sum=apply(pileup.data[,fs],1,sum)  #sum read count of FS
  pileup.data$bs.sum=apply(pileup.data[,bs],1,sum)  #sum read count of RS
  pileup.data$mut.rate1 = pileup.data$fs.sum/(pileup.data$R+pileup.data$fs.sum)*100 #FS mutation rate in %
  pileup.data$mut.rate2 = pileup.data$bs.sum/(pileup.data$r+pileup.data$bs.sum)*100   #RS mutation rate in %
  pileup.data$mut.rate = (pileup.data$fs.sum+pileup.data$bs.sum)/pileup.data$DEPTH*100 #Total mutation rate in %
  pileup.data$mut.F = apply(pileup.data[,fs],1,function(x) names(which.max(x)))  #FS mutation letter
  pileup.data$mut.R = apply(pileup.data[,bs],1,function(x) toupper(names(which.max(x)))) #RS mutation letter
  pileup.data$mut.B = apply(pileup.data[,fs]+pileup.data[,bs],1,function(x) names(which.max(x))) #overall (FS+RS) mutation letter
  pileup.data$mut.rateF = apply(pileup.data[,fs],1,max)/(pileup.data$R+pileup.data$fs.sum)*100 #FS specifc mutation rate
  pileup.data$mut.rateR = apply(pileup.data[,bs],1,max)/(pileup.data$r+pileup.data$bs.sum)*100 #RS specifc mutation rate
  # Code below filter out the bases where mutation will not be called and mutation rate (mut.rate.adj) is set to 0.
  
  if (filter) {
    #remove mutation rate data for bases meeting the following criteria
    pileup.data$mut.rate.adj=ifelse(
      pileup.data$DEPTH<100 |
        (pileup.data$R>=50 & pileup.data$r>=50 & 
           ( (pileup.data$mut.F!=pileup.data$mut.R & (pileup.data$fs.sum<10 | pileup.data$bs.sum<10))|         #added min of fs and bs <10   3/31/17
               (pileup.data$mut.rate1<1 & pileup.data$mut.rate2>5) | 
               (pileup.data$mut.rate1>5 & pileup.data$mut.rate2<1) |
               (pileup.data$mut.rate1==0 & pileup.data$mut.rate2>1.5) |
               (pileup.data$mut.rate1>1.5 & pileup.data$mut.rate2==0) |
               (pileup.data$mut.rate1>10*pileup.data$mut.rate2 & pileup.data$mut.rate2>0.1) | 
               (pileup.data$mut.rate2>10*pileup.data$mut.rate1 & pileup.data$mut.rate1>0.1) | 
               (pileup.data$mut.rateF<0.75*pileup.data$mut.rate1 & pileup.data$mut.rate1>1) |
               (pileup.data$mut.rateR<0.75*pileup.data$mut.rate2 & pileup.data$mut.rate2>1))),0,pileup.data$mut.rate) 
    
    pileup.data$mut.rate.adj[pileup.data$R>10*pileup.data$r & pileup.data$R >=100 |  pileup.data$r>10*pileup.data$R & pileup.data$r >=100]=0
  } else pileup.data$mut.rate.adj=pileup.data$mut.rate
  
  #set mutation rate to 0 if it is >10% (SNP) - removed in 7/10/17
  #pileup.data$mut.rate.adj[pileup.data$mut.rate.adj>=10]=0 
  
  #sum FS and RS count for each base
  d1= pileup.data[,c("A","a","C","c","T","t","G","g")]
  d1.2 = sapply(s,function(x) rowSums(d1[,toupper(names(d1))==x]))
  if (polish.method=="Binomial") d2=d1.2/pileup.data$DEPTH else d2=d1.2
  
  #calculate the average count using input rate and depth 
  ave.read.count=bgrd.ave*bgrd.depth 
  
  #define background for poisson or binomial using upper CI limit
  if (polish.method=="Poisson") {
    d3=qchisq(1-alpha/2,2*ave.read.count+2)/2 }  else {
      d3=qbeta(1-alpha/2,as.matrix(ave.read.count)+1,as.matrix(bgrd.depth-ave.read.count))
    }
  
  #call mutation if observed data > background
  pred=apply(d2>d3,1,sum)
  if (polyclone) {pred.nts=apply(d2>d3,1,function(x) which(x>0)) 
  clone.n=unlist(lapply(pred.nts,length))} else 
  { d2.2=d2
  d2.2[d2<d3]=0
  pred.nts=apply(d2.2,1,function(x) ifelse(max(x)>0,which.max(x),0))
  clone.n=(pred.nts>0)+0
  }
  
  pred[pileup.data$mut.rate.adj==0]=0  #for those bases to be filtered out, we do not call.
  #pred[apply(d1.2,1,max)<5]=0          #for those bases whose total counts<5, we do not call. removed 7/10/17
  
  #polished noise data output 
  polish.noise=cbind(out.pileup,noise.num=apply(d1.2>0,1,sum),freq=apply(d1.2,1,max)/out.pileup$DEPTH,threshold=d3[cbind(1:nrow(d3),apply(d1.2,1,which.max))])[apply(d1.2,1,sum)>0 & pred==0,]
  
  #output the polished data
  out.pileup[pred==0,7:14]=0 #hard coded columns, assume column 7:14 are the ACTG columns
  out.pileup$DEPTH=apply(out.pileup[,5:14],1,sum) #re-calculate the depth after polish
  
  #handling output file name and directory
  #if input file name has directory path, remove the path
  sample.name.loc=unlist(gregexpr("/",pred.sample)) 
  if (sample.name.loc[1]!=-1) sample.name.loc=sample.name.loc[length(sample.name.loc)]+1 else
    sample.name.loc=1
  out.file.name=substr(pred.sample,sample.name.loc,nchar(pred.sample))
  out.file.name=gsub(".txt",".bgrm.txt",out.file.name)
  if (!is.null(out.dir)) out.file=file.path(out.dir,out.file.name) else
    out.file=out.file.name
  write.table(out.pileup,file=out.file,sep="\t",quote=F,row.names=F)
  
  #summary statistic for the predicted mutants
  sum.data1=list()
  for (j in 1:max(clone.n)) {
    if (polyclone) pred.nt=rapply(pred.nts[clone.n>=j], function(x) x[j]) else pred.nt=pred.nts[pred.nts>0]
    sum.data=out.pileup[clone.n>=j,]
    sum.data$mutant=s[pred.nt]
    sum.data$tnt=tnt[clone.n>=j]
    sum.data$bkg.threshold=d3[clone.n>=j,,drop=F][cbind(1:length(pred.nt),pred.nt)]
    sum.data$mut.freq=     d2[clone.n>=j,,drop=F][cbind(1:length(pred.nt),pred.nt)]
    sum.data$tnt[sum.data$tnt==""]="NA"
    sum.data$tnt.ave.freq=bgrd.ave[clone.n>=j,,drop=F][cbind(1:length(pred.nt),pred.nt)] #average healthy subject tri-nuc count
    if (polish.method=="Poisson") sum.data$zscore=qnorm(ppois(sum.data$mut.freq,sum.data$tnt.ave.freq,lower.tail = F,log.p=T),lower.tail = F,log.p=T) else
      if (polish.method=="Binomial") sum.data$zscore=qnorm(pbinom(sum.data$mut.freq*10000,10000,sum.data$tnt.ave.freq,lower.tail = F,log.p=T),lower.tail = F,log.p=T)
    sum.data1=rbind(sum.data1,sum.data)
  }
  write.table(polish.noise,file=sub(".bgrm.txt",".polished.noise.txt",out.file),sep="\t",quote=F,row.names=F)
  return(pred)
}



# call the function to run each input sample
##############################################
n.sample=length(input.pileup.file)
pred3.all=matrix(0,ncol=n.sample,nrow=nrow(pileup1.data))
for (i in 1:n.sample) pred3.all[,i]=polish.pred.3nt(bgrd.ave = shr.bg.ave,
                                                    bgrd.depth=depth.mat, 
                                                    pred.sample = input.pileup.file[i],
                                                    alpha=input.alpha,
                                                    filter=input.filter,
                                                    polyclone=TRUE,
                                                    out.dir=work.dir)

# The resulted data pred3.all is an indicator matrix for each position(row) and sample (column)
# 1=mutation, 0= no mutation
