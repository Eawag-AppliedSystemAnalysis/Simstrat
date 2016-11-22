#This script reads the output of the k-epsilon model, extracts the data
#that can be compared to actual measurements, rewrites this data into a new
#text file and writes the PEST instruction file that allows PEST to
#understand it. This script is run by PEST after each kepsilon simulation.

wdir = commandArgs(trailingOnly = TRUE)[1]
if (is.na(wdir)) {
  wdir=""
} else if (wdir[length(wdir)]!="/") {
  wdir=paste(wdir,"/",sep="")
}

#Path to kepsilon model directory (for function GetResults)
source("../GetResults.R")
#Types, times, depths and NaNs of available measurements from script Control.m
load("obsprop.RData")

#Extract zobs x tobs matrices
fid = file(paste(wdir,"kepsilon_PEST.par",sep=""))
open(fid,"r")
line = unlist(read.table(fid,skip=6,nrow=1))
close(fid)
time = list()
z = list()
data = list()
label = list()
for (i in 1:length(obs)) {
  res = GetResults(line,obs[[i]],depths=zobs[[i]])
  time[[i]] = res$time
  z[[i]] = res$z
  data[[i]] = res$data[[1]]
  label[[i]] = res$label[[1]]
}
for (i in 1:length(obs)) {
  if (min(tobs[[i]])<min(time[[i]]) || max(tobs[[i]])>max(time[[i]])) {
    warning(paste(sprintf("Some given %s measurements fall out of simulated time period.",obs[[i]]),
             "NaN values are to be expected in the model observations file, which PEST won't accept.",
             "This should be taken care of in the script Control.R."))
  }
  if (min(zobs[[i]])<min(z[[i]]) || max(zobs[[i]])>max(z[[i]])) {
    warning(paste(sprintf("Some given %s measurements fall out of simulated depth range.",obs[[i]]),
              "NaN values are to be expected in the model observations file, which PEST won't accept.",
              "This should be taken care of in the script Control.R."))
  }
}

for (i in 1:length(data)) {
  data[[i]] = t(matrix(unlist(apply(data[[i]],2,approx,x=as.numeric(time[[i]]),xout=as.numeric(tobs[[i]]))),nrow=length(tobs[[i]]))[,-seq(1,2*length(zobs[[i]]),2)])
}

#Rewrite model observations in one file
fid = file(paste(wdir,"ModelObs.dat",sep=""))
open(fid,"w")
for (i in 1:length(data)) {
  for (it in 1:length(tobs[[i]])) {
    for (iz in 1:length(zobs[[i]])) {
      if (length(nanobs[[i]])==0 || !any(it==nanobs[[i]][,1] & iz==nanobs[[i]][,2])) {
        writeLines(sprintf("%12.4e",data[[i]][iz,it]),fid,sep="")
      } else {
        writeLines(rep(" ",12),fid,sep="")
      }
    }
    writeLines("",fid)
  }
}
close(fid)

#Write corresponding instructions file for PEST
fid = file(paste(wdir,"keps_obs.ins",sep=""))
open(fid,"w")
writeLines("pif @",fid)
for (i in 1:length(data)) {
  for (it in 1:length(tobs[[i]])) {
    strf = "l1"
    for (iz in 1:length(zobs[[i]])) {
      if (length(nanobs[[i]])==0 || !any(it==nanobs[[i]][,1] & iz==nanobs[[i]][,2])) {
        strf = paste(strf," [",substr(label[[i]],1,4),sprintf("_%d_%d",it,iz),"]",sprintf("%d:%d",12*iz-11,12*iz),sep="")
      }
    }
    writeLines(strf,fid)
  }
}
close(fid)
