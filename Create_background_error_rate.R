#################################################################
# File name: Create_background_error_rate.R
# Create a background error rate file and a coverage file from healthy subject data
# Last updated: 08/20/2018 
# Contact: shibing.deng {at} pfizer.com & tao.xie {at} pfizer.com or xietao2000 {at} gmail.com
# Pfizer Early Clinical Development Biostatistics 
# Pfizer Early Oncology Development and Clinical Research 
#################################################################

#################################################################
# Terms of Use for the TNER (Tri-Nucleotide Error Reducer) package
# Pfizer has attempted to provide accurate and current information on TNER. However, Pfizer makes no representations or warranties that the information contained or published on TNER will be suitable for your specific purposes or for any other purposes. You agree to indemnify, defend, and hold Pfizer harmless from all claims, damages, liabilities and expenses, including without limitation reasonable attorney's fees and costs, whether or not a lawsuit or other proceeding is filed, that in any way arise out of or relate to your use of TNER.
# ALL INFORMATION AND DATA PROVIDED ON TNER IS PROVIDED "AS-IS" WITHOUT WARRANTY OF ANY KIND. PFIZER MAKES NO REPRESENTATIONS OR WARRANTIES CONCERNING TNER OR ANY OTHER MATTER WHATSOEVER, INCLUDING WITHOUT LIMITATION ANY EXPRESS, IMPLIED OR STATUTORY WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NON-INFRINGEMENT OF THIRD PARTY RIGHTS, TITLE, ACCURACY, COMPLETENESS OR ARISING OUT OF COURSE OF CONDUCT OR TRADE CUSTOM OR USAGE, AND DISCLAIMS ALL SUCH EXPRESS, IMPLIED OR STATUTORY WARRANTIES. PFIZER MAKES NO WARRANTY OR REPRESENTATION THAT YOUR USE OF TNER WILL NOT INFRINGE UPON THE INTELLECTUAL PROPERTY OR OTHER RIGHTS OF ANY THIRD PARTY. FURTHER, PFIZER SHALL NOT BE LIABLE IN ANY MANNER WHATSOEVER FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, CONSEQUENTIAL OR EXEMPLARY DAMAGES ARISING OUT OF OR IN ANY WAY RELATED TO TNER, THE USE OF, OR INABILITY TO USE, ANY OF THE INFORMATION OR DATA CONTAINED OR REFERENCED IN TNER, OR ANY OTHER MATTER. THE FOREGOING EXCLUSIONS AND LIMITATIONS SHALL APPLY TO ALL CLAIMS AND ACTIONS OF ANY KIND AND ON ANY THEORY OF LIABILITY, WHETHER BASED ON CONTRACT, TORT OR ANY OTHER GROUNDS, AND REGARDLESS OF WHETHER A PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES, AND NOTWITHSTANDING ANY FAILURE OF ESSENTIAL PURPOSE OF ANY LIMITED REMEDY. BY USING TNER, YOU FURTHER AGREE THAT EACH WARRANTY DISCLAIMER, EXCLUSION OF DAMAGES OR OTHER LIMITATION OF LIABILITY HEREIN IS INTENDED TO BE SEVERABLE AND INDEPENDENT OF THE OTHER CLAUSES OR SENTENCES BECAUSE THEY EACH REPRESENT SEPARATE ELEMENTS OF RISK ALLOCATION BETWEEN THE PARTIES.
# TNER may contain information from third party websites, which may or may not be marked with the name of the source. Such information does not necessarily represent the views or opinions of Pfizer, and Pfizer shall have no responsibility whatsoever for such information. All information from third party websites are the sole responsibility of the person or entity that provides and/or maintains such website. As a TNER user, you are solely responsible for any information that you display, generate, transmit or transfer while using TNER, and for the consequences of such actions.
# Should any TNER user provide general, scientific or other feedback information, whether in the form of questions, comments, or suggestions to Pfizer, regarding the content of Pfizer's website or otherwise, such information shall not be deemed to be confidential or proprietary to you or to any other party. Pfizer shall have no obligation of any kind with respect to such information and Pfizer shall have the right, without limitation, to reproduce, use, disclose, merge, display, make derivatives of and distribute such information to others. Pfizer shall also have the right, without limitation, to use and exploit such information, including ideas, concepts, know-how, inventions, techniques or other materials contained in such information for any purpose whatsoever, including but not limited to, making changes or improvements to TNER and/or developing, manufacturing, marketing, selling or distributing products and technologies incorporating such information. However, you agree that Pfizer has no obligation whatsoever to respond to your comments or to change or correct any information on TNER based on your comments.
# Pfizer reserves the right to alter the content of TNER in any way, at any time and for any reason, with or without prior notice to you, and Pfizer will not be liable in any way for possible consequences of such changes or for inaccuracies, typographical errors or omissions in the contents hereof. Nothing contained herein shall be construed as conferring by implication, estoppel or otherwise any license or right under any patent, trademark or other intellectual property of Pfizer or any third party. Except as expressly provided above, nothing contained herein shall be construed as conferring any right or license under any Pfizer copyrights.
# Pfizer also reserves the right to modify these Terms of Use at any time and for any reason, with or without prior notice to you. You should always review the most current Terms of Use herein before using this website. By using TNER, you agree to the current Terms of Use posted on this site. You also agree that these Terms of Use constitute the entire agreement between you and Pfizer regarding the subject matter hereof and supersede all prior and contemporaneous understandings, oral or written, regarding such subject matter. In addition, if any provision of these Terms of Use is found by a court of competent jurisdiction to be invalid, void or unenforceable, the remaining provisions shall remain in full force and effect, and the affected provision shall be revised so as to reflect the original intent of the parties hereunder to the maximum extent permitted by applicable law.
# BY USING TNER, YOU ACKNOWLEDGE AND AGREE THAT YOU HAVE READ, UNDERSTOOD AND AGREE TO ALL OF THESE TERMS OF USE.
#################################################################

#################################################################
# Command line usage:
#################################################################
# Rscript  Create_background_error_rate.R 
#################################################################
# Input:   Input healthy subjects' nucleotide frequency files should be stored in a subfolder named "healthy_subjects".  
#          The columns are CHR = chromosome number; POSISION= genomic coordinate position
#           DEPTH = coverage depth; REF=reference nucleotide;
#           R+ = forward strand coverage for the reference nucleotide; R- = reverse strand coverage for the reference nucleotide;
#           A+ = forward strand base A count; A- = reverse strand base A count; ... ...
#           Note: only error counts are listed. The reference base count is set to 0.
#                 For example, if REF=A, then A+ =0 and A- = 0.
# Note:  Make sure the healthy subject data are sorted by chromosome number and POS. 
#################################################################
# Output:   A csv file named "hs_ave_bg_error.csv": provides position specific average background error rate;
# also a csv file named "hs_depth.csv": provides read depth data for all covered positions.
###########################################################################

# data directory
setwd(".")  #set work directory
data.dir = "./healthy_subjects"  #Assume the healthy_subjects folder in under the working dir

hs.files=list.files(path=data.dir,pattern="freq$") # list all healthy subject pileup data files; require at least two or more normal samples

s=c("A","C","T","G") # four nucleotide letter

hs.db=list()
for (i in 1:length(hs.files))
{
  #read in data
  pileup.data=read.table(file=file.path(data.dir,hs.files[i]),header=T,check.names = F,as.is=T)
  #below is not a critical step. It changes the column name of forward strand to cap letter: A+ to A; and 
  #reserse strand to lower case letter A+ to a.
  if (length(grep("\\+",names(pileup.data)))>0) names(pileup.data)=gsub("\\+","",names(pileup.data))
  if (length(grep("-",names(pileup.data)))>0) {
    rev.base.col=grep("-",names(pileup.data))  # data column has reverse string pile up
    names(pileup.data)[rev.base.col]=gsub("-","",tolower(names(pileup.data)[rev.base.col]))
  }
  
 # calculate background error rate
  d1<- pileup.data[,c(rbind(s,tolower(s)))] #select background error columns (A,a,C,C,T,t,G,g)
  d2 = sapply(s,function(x) rowSums(d1[,toupper(names(d1))==x])) #sum error count
  pileup.data$DEPTH=ifelse(pileup.data$DEPTH<1,1,pileup.data$DEPTH)  #handling 0 depth and set them to 1
  hs.db[[i]]=d2/pileup.data$DEPTH # calulate error rate
  
  if (i==1) hs.depth=pileup.data$DEPTH else hs.depth=cbind(hs.depth,pileup.data$DEPTH) #save the depth. If length of rows do not match, error will occur
}

#combine all background error rate into an array
hsa=array(unlist(hs.db),c(dim(hs.db[[1]]),length(hs.db)),dimnames = list(paste(pileup.data$CHR,pileup.data$POSITION,sep=":"),s,hs.files))
hs.bkg.ave=apply( hsa,1:2,mean )   #calculate the average background error rate
hs.bkg.ave=cbind(pileup.data[,c(1:2,4)],hs.bkg.ave) #add chr, pos and ref column to the average bkg error data

rownames(hs.depth)=rownames(hs.bkg.ave)
colnames(hs.depth)=hs.files
#make sure the data are sorted by CHR and POS

write.csv(hs.bkg.ave,file="hs_ave_bg_error.csv") #output healthy subject average background error rate data
write.csv(hs.depth,file="hs_depth.csv")     #output healthy subject depth data
