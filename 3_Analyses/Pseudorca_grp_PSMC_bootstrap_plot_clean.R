# plot PSMC plots with bootstraps, using output from psmc_plot.pl
#Created by P. Morin. Modified by K. Hernandez May 2025.

# step 1 (create the .txt files from .psmc files using psmc_plot.pl) - run in terminal
# step 2 (input .txt files from step 1 to create final plot) - run in R

### Set working directory and get list of all input files including bootstrap replicates (these are the output text files generated from the utils/psmc_plot.pl -R command)
# setwd("main_with_boot/plot_mu1.4e-08_g10")
allfiles=list.files(pattern="txt")

# get species names based on species name patterns from the filenames
pop1<-"Pcra_z0018462"
pop2<-"Pcra_z0045928"
pop3<-"hPSMC"
#sp4<-""
#sp5<-""
#sp6<-""
#sp7<-""
#sp8<-""
#sp9<-""

leg1<-c("ETP","MHI", "hPSMC") # for legend on plot


# Save as PNG using png() and dev.off()
png(paste0("Pseudorca","_psmc_plot.png"),width = 1500, height = 900, units = "px", res = 300)
#converted the 5x3 original in inches based on 300 dpi recommended by JOH
### Set min and max values for plot axes
xmin=1.2e4
xmax=1e7
ymin=0
ymax=18 #original was 15 for Berardius

### Set line colors (first, pick rgb values for each sample, then set main and transparent colors for plot lines)
# color1 (aquamarine)
c1=c(127,255,212)/255

# color (purple)
c2=c(102, 102, 255)/255
# color (orange)
c3=c(240, 105, 20)/255

# Main colors (no transparency)
mycols1=c(
  rgb(c1[1], c1[2], c1[3], alpha=1),
  rgb(c2[1], c2[2], c2[3], alpha=1),
  rgb(c3[1], c3[2], c3[3], alpha=1)
)

# Bootstrap replicate colors (with transparency)
transp=0.05 # 0.05; 0 if no bootstraps needed.
mycols2=c(
  rgb(c1[1], c1[2], c1[3], alpha=transp),
  rgb(c2[1], c2[2], c2[3], alpha=transp),
  rgb(c3[1], c3[2], c3[3], alpha=0)
)

### Generate an empty plot with labeled axes
par(mar=c(3.5,3.75,0.5,0.5))
op <- par(cex = 0.75) # font size

plot(1, 1, type="n", log="x", axes=F, xlim=c(xmin, xmax), ylim=c(ymin, ymax), xlab="", ylab="")

title(xlab="Years before present", line=2)
title(ylab=expression("Effective population size (x10"^4*")"), line=2.25)

axis(side=2, line=0, labels=F)
axis(side=2, line=-.25, labels=T, tick=F)

at.x=outer(1:9, 10^(3:8))
lab.x=NULL
for (i in 1:length(at.x)){
  p=log10(at.x[i])
  if (p %% 1 == 0) {lab.x[i]=as.expression(bquote(10^ .(p)))}
  else {lab.x[i]=""}
}
axis(1, at=at.x, labels=lab.x, las=1)

legend("topright", lwd=3, col=mycols1[c(1,2,3,4,5,6,7,8,9)], legend=leg1, 
       bty="n")

box()

### Function to add plot lines for each sample
psmc_plot_fill=function(){
  # Get list of input files using "samplename" as the search term
  dfiles=allfiles[grep(pattern=samplename, allfiles)]
  # Loop through the bootstrap reps (first bootstrap file = the second file from the list above)
  for (i in 2:length(dfiles)){
    bb=read.table(dfiles[i])
    # Plot lines for each bootstrap file using the transparency colors
    lines(bb$V1, bb$V2, type="s", col=mycols2[nn], lwd=1)
  }
  # Read in the first file from the list (this one is for adding the main solid line 
  # on top of the bootstrap lines, given as boot.out.0.txt from the psmc_plot.pl script)
  aa=read.table(dfiles[1])
  # Plot the line using the main color (no transparency)
  lines(aa$V1, aa$V2, type="s", col=mycols1[nn], lwd=2)
}

### Plot each sample (nn is the numerical index for the sample - only needs to 
# correspond to the order of samples in the color lists above)
# comment out "psmc_plot_fill()" for unused samplenames.

# samples to plot
# 1
samplename=pop1
nn=1
psmc_plot_fill()
# 2
samplename=pop2
nn=2
psmc_plot_fill()
# 3 #can't use PSMC function b/c no bootstrap replicates
samplename=pop3
nn=3
lines(hPSMC_z0018462_z0045928_9.10E10_msy_t15_psmc.out.0$V1, 
      hPSMC_z0018462_z0045928_9.10E10_msy_t15_psmc.out.0$V2, type="s", col=mycols1[nn], lwd=2)

# Close the PNG device
dev.off()

