#This script writes the PEST control file, given matrices of measured
#data and given a set of parameters with their related properties. This
#script should be run once, to prepare an optimization with PEST.

rm(list = ls())
graphics.off()

#Measurements, in time x depth matrices
dat = read.table("../Bielersee.txt",sep=",")
zobs = list(as.numeric(gsub(",","",as.character(as.matrix(dat[1,2:ncol(dat)])))))
zobs[[1]] = zobs[[1]][!is.na(zobs[[1]])]
tobs = list(as.Date(dat[2:nrow(dat),1]))
tobs = list(strptime(dat[2:nrow(dat),1],format="%Y-%m-%d %H:%M"))
obsl = list("Temperature")
obs = list("T")
obsdat = list(dat[2:nrow(dat),2:ncol(dat)])
#Remove measurements out of simulation period
sim_T = strptime(c("1994-03-01","2015-01-01"),format="%Y-%m-%d")
obsdat[[1]] = obsdat[[1]][tobs[[1]]>=sim_T[1] & tobs[[1]]<=sim_T[2],]
tobs[[1]] = tobs[[1]][tobs[[1]]>=sim_T[1] & tobs[[1]]<=sim_T[2]]
#Keep only depths and times with enough non-missing values
nap = round(100*colSums(is.na(obsdat[[1]]))/nrow(obsdat[[1]]))
obsdat[[1]] = matrix(as.numeric(as.matrix(obsdat[[1]][,nap<20])),nrow=nrow(obsdat[[1]]))
zobs[[1]] = zobs[[1]][nap<20]

#Define parameter names, group, values, min, max
parname = c("lat","pair","alpha","qNN","fwind","C10","CD","fgeo","kmin","p1","p2","beta","albsw")
partype = c("fixed","fixed","none","fixed","fixed","fixed","fixed","fixed","fixed","none","none","fixed","fixed")
parlim = c("factor","factor","factor","factor","factor","factor","factor","factor","factor","factor","factor","relative","relative")
parval = c(47.10, 990,0.010,1.25,1.0,0.0017,0.002,0.08,1e-15,1.1,0.9,0.3,0.2)
parmin = c(40.00, 900,0.001,0.70,1.0,0.0010,0.001,0.00,1e-30,0.5,0.5,0.0,0.0)
parmax = c(50.00,1000,0.100,1.30,2.0,0.0030,0.005,1.00,1e-09,1.5,1.5,0.4,0.3)
pargrp = c("none","none","fit","none","none","none","none","none","none","fit","fit","none","none")
npar = length(parname)

#Total number of NaNs
nnan = 0
for (i in 1:length(obsdat)) {nnan=nnan+sum(is.na(obsdat[[i]]))}

#Write corresponding control file for PEST
fid = file("keps_calib.pst")
open(fid,"w")
writeLines("pcf",fid)
writeLines("* control data",fid)
writeLines("restart estimation",fid)
writeLines(sprintf("%d %d 1 0 1",npar,sum(sapply(obsdat[[1]],length))-nnan),fid)
writeLines("1 1 single nopoint 1 0 0",fid)
writeLines("5.0 2.0 0.3 0.01 10",fid)
writeLines("5.0 5.0 0.001",fid)
writeLines("0.1",fid)
writeLines("30 0.005 4 3 0.01 3",fid)
writeLines("1 1 1",fid)
writeLines("* parameter groups",fid)
writeLines(" fit\trelative\t0.01\t0.00001\tswitch\t2.0\tparabolic",fid)
writeLines("* parameter data",fid)
for (i in 1:npar) {
  writeLines(sprintf("%6s\t%s\t%s\t%10.4e\t%10.4e\t%10.4e\t%4s\t1.0\t0.0\t1",parname[[i]],partype[[i]],parlim[[i]],parval[i],parmin[i],parmax[i],pargrp[[i]]),fid)
}
writeLines("* observation groups",fid)
writeLines("obs",fid)
writeLines("* observation data",fid)
nanobs = list(matrix(NA,nrow=0,ncol=2))
for (i in 1:length(obsdat)) {
  k = 1
  for (it in 1:length(tobs[[i]])) {
    for (iz in 1:length(zobs[[i]])) {
      if (!is.na(obsdat[[i]][it,iz])) {
        writeLines(sprintf("%s_%d_%d\t%12.4e\t1.0\tobs",substr(obsl[[i]],1,4),it,iz,obsdat[[i]][it,iz]),fid)
      } else { #Store location of missing measurements so as not to write corresponding model observations
        nanobs[[i]] = rbind(nanobs[[i]],c(it,iz))
        k = k+1
      }
    }
  }
}
writeLines("* model command line",fid)
writeLines("model.bat",fid)
writeLines("* model input/output",fid)
writeLines("keps_par.tpl\tkepsilon_PEST.par",fid)
writeLines("keps_obs.ins\tModelObs.dat",fid)
writeLines("* prior information",fid)
close(fid)

#Save measurement properties for script OutputInstructions.m
save(obs,tobs,zobs,nanobs,file="obsprop.RData")
unlink(".RData")