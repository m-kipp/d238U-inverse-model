##### d238U inverse model #####
##### Michael Kipp (mkipp@caltech.edu) #####
##### 30 August 2021 #####

# loading necessary packages 

library(msir)
library(FME)
library(doParallel)
library(LaplacesDemon)


### PART I: READING IN DATA ###

# set working directory, read in dataset

setwd("/Users/MichaelKipp/Downloads")
current_dataset=read.csv("jost_2017_Triassic_data.csv",header=T,stringsAsFactors=F)


# quick plot of dataset

par(mfrow=c(1,1))
par(mar=c(4,4,1,1))
par(oma=c(0,0,0,0))
plot(current_dataset$d238U~current_dataset$time,axes=F,xlab="",ylab="")
loess_mod=loess.sd(current_dataset$d238U~current_dataset$time,span=0.3,nsigma=2) 
polygon(c(loess_mod$x,rev(loess_mod$x)),c(loess_mod$upper,rev(loess_mod$lower)),col=rgb(0,0,0.8,0.1),border=F)
lines(loess_mod$y~loess_mod$x,col="blue")
arrows(current_dataset$time,current_dataset$d238U+(2*current_dataset$err),current_dataset$time,current_dataset$d238U-(2*current_dataset$err),length=0)
points(current_dataset$d238U~current_dataset$time,pch=21,bg="white")
axis(side=1,mgp=c(3,0.5,0),las=1,tck=-0.02)
axis(side=2,mgp=c(3,0.5,0),las=1,tck=-0.02)
mtext("time [Myr]",side=1,line=2.25)
mtext(expression(delta^{238}*"U [\u2030]"),side=2,line=2.5)
box(lwd=1)


## model parameters

# things to definitely tune

time_step=10e3 # size of time-steps; should be at minimum 1e3; must match temporal resolution of data
prop_uncert=F # propagate uncertainty on terms in mass balance? yes=T, no=F [start no for testing, use yes at end]
m=6 # number of terms in Fourier series [start low, demonstrate convergence, increase if needed, test convergence again, repeat...]
niterMCMC=2e6 # number of steps per walker [likely needs to be 10^5-10^6; start at 10^4 and see how long it will take]
updatecovMCMC=2e3 # number of steps after which covariance matrix is updated (only during burn-in phase) [see table in ReadMe]
n_walkers=3 # number of parallel walkers deployed [aim for 3-5 in initial tests, 100+ if propagating uncertainty]

# things to maybe never tune (fingers crossed)

amp_limit=2e-5 # range for amplitude terms [decrease if low acceptance persists]
start=1 # start of frequency terms; 1 for odd harmonics or continuous, 2 for even harmonics [default=1]
count=2 # spacing of frequency terms; 2 for odd or even harmonics, 1 for continuous [default=2]
jump_step_1=amp_limit*0.05 # SD of proposal distribution for amplitude terms [default=amp_limit*0.05]
jump_step_2=0.2 # SD of proposal distribution for log_fanox_0 [default=0.2]
burninlengthMCMC=niterMCMC-1e3 # number of discarded "burn in" steps [1e3 less than niterMCMC typically works well]
thin=1 # factor by which final output is "thinned"; 1=no thinning, 10=keep every 10th sample [default=1]

# things to never tune

start_time=min(current_dataset$time) 
max_time=max(current_dataset$time)+(time_step-start_time)
a_vals=seq(0,0,length.out=m)
b_vals=seq(0,0,length.out=m)
as=paste("a",1:m,sep="")
bs=paste("b",1:m,sep="")
names(a_vals)=as
names(b_vals)=bs
pars1=c(log_fanox_0=-2.67,a_vals,b_vals) # starting values; gives modern mass balance
upper=c(0,seq(amp_limit, amp_limit,length=2*m)) # upper limits
lower=c(-4,seq(-amp_limit,-amp_limit,length=2*m)) # lower limits
jump=c(jump_step_2,seq(jump_step_1,jump_step_1,length=2*m)) 

# steady-state matrix with d238U and f_anox values
# this is used to generate a steady-state scenario at the beginning of each time series

Jriv=0.42e8 # Dunk et al (2002) Chem Geol
Kanox=1.74e-4 # derived in Kipp & Tissot (2021) EPSL
Kother=1.88e-6  # derived in Kipp & Tissot (2021) EPSL
d238Uriv=-0.30 # Tissot et al (2015); Andersen et al (2016) Chem Geol
Danox=0.6 # widely-used value, see e.g. Zhang et al (2020) GCA
Dother=0 # widely-used value, see e.g. Zhang et al (2020) GCA
fanox_seq=10^seq(-5,0,length=1000)
log_fanox_seq=seq(-5,0,length=1000)
d238U_SS=d238Uriv-(((Jriv/((Kanox*fanox_seq)+Kother*(1-fanox_seq)))*((Kother*Dother)+(fanox_seq*((Kanox*Danox)-(Kother*Dother)))))/Jriv)
ss_out=cbind(d238U_SS, log_fanox_seq)


### PART II: FORWARD MODEL ###

d238_model=function(pars){
	derivs=function(time,y,pars){
		with(as.list(c(pars,y)),{
			dd238Usw=(Jriv*(d238Uriv-d238Usw)-Kanox*Nsw*(10^log_fanox)*Danox-Kother*Nsw*(1-(10^log_fanox))*Dother)/Nsw
			dNsw=Jriv-Kanox*Nsw*(10^log_fanox)-Kother*Nsw*(1-(10^log_fanox))	
			dlog_fanox=if(log_fanox>-0.004){-1e-8} else if(abs(log_fanox)>4){1e-8} else {sum((pars[seq(2,by=1,length.out=m)]*sin(seq(start,length.out=m,by=count)*pi*time/max_time))+(pars[seq(2+m,by=1,length.out=m)]*cos(seq(start,length.out=m,by=count)*pi*time/max_time)))}

			return(list(c(dd238Usw,dNsw,dlog_fanox)))
		})
	}
	
log_fanox_0=with(as.list(pars),log_fanox_0)
Nsw_0=with(as.list(pars),Jriv/(Kanox*(10^log_fanox_0)+Kother*(1-(10^log_fanox_0))))
d238Usw_0=with(as.list(pars),as.numeric((ss_out[ss_out[,2]>=log_fanox_0,])[1,1]))

y=c(d238Usw=d238Usw_0,Nsw=Nsw_0,log_fanox=log_fanox_0)
	
times=seq(start_time,max_time,by=time_step)
out=ode(y=y,parms=pars,times=times,func=derivs,method="euler")
		
as.data.frame(out)
}


# plotting baseline (=starting) model run (should be flatline at modern values)

out1=d238_model(pars=pars1)

par(mfrow=c(3,1))
plot(out1$log_fanox~out1$time,cex=0,ylim=c(-4,-0))
lines(out1$log_fanox~out1$time)
plot(out1$Nsw~out1$time,cex=0,ylim=c(0,3e13))
lines(out1$Nsw~out1$time)
plot(out1$d238Usw~out1$time,cex=0,ylim=c(-1.8,0))
lines(out1$d238Usw~out1$time)
points(current_dataset$d238U~current_dataset$time,pch=21,bg="white")


### PART III: INVERSE MODEL ###

# defining cost function (negative log-likelihood, NLL)

d238Ucost=function(pars){
	out1=d238_model(pars=pars)
	merged=merge(out1,current_dataset,by="time")
	NLL=log(2*pi*merged$err)+na.omit(((merged$d238U-merged$d238Usw)^2)/(merged$err^2)) 
	return(sum(NLL))
} 


# setting up parallel computation 

cores=detectCores()
cl=makeCluster(cores[1]-1)
registerDoParallel(cl)


# running MCMC routine

start.time=Sys.time()
finalMCMC=foreach(i=1:n_walkers,.combine=rbind) %dopar% {
	
library(FME)

if(prop_uncert==T){
	
Kanox=rnorm(1,1.74e-4,0.655e-4) # derived in Kipp & Tissot (2021) EPSL
Kother=35.7e6/(19e12*(1-(6.3e6/(19e12*Kanox))))  # derived in Kipp & Tissot (2021) EPSL
d238Uriv=rnorm(1,-0.30,0.02) # Tissot et al (2015); Andersen et al (2016) Chem Geol
Danox=rnorm(1,0.6,0.1) # widely-used value, see e.g. Zhang et al (2020) GCA
Dother=rnorm(1,0,0.025) # widely-used value, see e.g. Zhang et al (2020) GCA
	
}

fanox_seq=10^seq(-5,0,length=1000)
log_fanox_seq=seq(-5,0,length=1000)
d238U_SS=d238Uriv-(((Jriv/((Kanox*fanox_seq)+Kother*(1-fanox_seq)))*((Kother*Dother)+(fanox_seq*((Kanox*Danox)-(Kother*Dother)))))/Jriv)
ss_out=cbind(d238U_SS, log_fanox_seq)

curMCMC=modMCMC(f=d238Ucost,p=pars1,niter=niterMCMC,lower=lower,upper=upper,jump=jump,prior=NULL,var0=NULL,wvar0=0.1,updatecov=updatecovMCMC,burninlength=burninlengthMCMC)
	
tempMatrix=cbind(curMCMC$par)
tempMatrix

}
end.time=Sys.time()
run.time1=end.time-start.time
run.time1
stopCluster(cl)


# thinning sample set

thin_count=dim(finalMCMC)[1]/thin
i=1
new_out=NULL
for (i in 1:thin_count){
	cur_vals=finalMCMC[i*thin,]
	new_out=rbind(new_out,cur_vals)
}

write.csv(new_out,"MCMC_out.csv")
MCMC_out=read.csv("MCMC_out.csv")
MCMC_out=MCMC_out[,2:dim(MCMC_out)[2]]


# using compiled MCMC output to generate confidence interval for timeseries

start.time=Sys.time()

i=1
n_reps=1e3 # 1e3 is sufficient for final run; can decrease to 1e2 for testing 
time_seq=seq(start_time,max_time,by=time_step)
d238Usw_out=NULL
Nsw_out=NULL
fanox_out=NULL
for (i in 1:n_reps){
	
	if(prop_uncert==T){
	
	Kanox=rnorm(1,1.74e-4,0.655e-4) # derived in Kipp & Tissot (2021) EPSL
	Kother=35.7e6/(19e12*(1-(6.3e6/(19e12*Kanox))))  # derived in Kipp & Tissot (2021) EPSL
	d238Uriv=rnorm(1,-0.30,0.02) # Tissot et al (2015); Andersen et al (2016) Chem Geol
	Danox=rnorm(1,0.6,0.1) # widely-used value, see e.g. Zhang et al (2020) GCA
	Dother=rnorm(1,0,0.025) # widely-used value, see e.g. Zhang et al (2020) GCA
	
	}
	
	fanox_seq=10^seq(-5,0,length=1000)
	log_fanox_seq=seq(-5,0,length=1000)
	d238U_SS=d238Uriv-(((Jriv/((Kanox*fanox_seq)+Kother*(1-fanox_seq)))*((Kother*Dother)+(fanox_seq*((Kanox*Danox)-(Kother*Dother)))))/Jriv)
	ss_out=cbind(d238U_SS, log_fanox_seq)
	
	sr1=sensRange(fun=d238_model,parms=NULL,parInput= MCMC_out,num=1) 
	d238Usw_sens=sr1[,(length(pars1)+1):(length(pars1)+length(time_seq))]
	Nsw_sens=sr1[,(length(pars1)+1+length(time_seq)):(length(pars1)+2*length(time_seq))]
	fanox_sens=sr1[,(length(pars1)+1+2*length(time_seq)):(length(pars1)+3*length(time_seq))]
	
	d238Usw_out=rbind(d238Usw_out, d238Usw_sens)
	Nsw_out=rbind(Nsw_out, Nsw_sens)
	fanox_out=rbind(fanox_out, fanox_sens)
}


i=1
d238Usw_quant_out=NULL
Nsw_quant_out=NULL
fanox_quant_out=NULL

for (i in 1:length(time_seq)){
	cur_time=time_seq[i]
	cur_d238Usw_quant=quantile(na.omit(d238Usw_out[,i]),probs=seq(0.01,1,by=0.01))
	cur_Nsw_quant=quantile(na.omit(Nsw_out[,i]),probs=seq(0.01,1,by=0.01))
	cur_fanox_quant=quantile(na.omit(fanox_out[,i]),probs=seq(0.01,1,by=0.01))
	d238Usw_quant_out=rbind(d238Usw_quant_out, cur_d238Usw_quant)
	Nsw_quant_out=rbind(Nsw_quant_out, cur_Nsw_quant)
	fanox_quant_out=rbind(fanox_quant_out, cur_fanox_quant)
}


end.time=Sys.time()
run.time2=end.time-start.time
run.time2



### PART IV: PLOTTING d238Usw and f_anox WITH MCMC-DERIVED POSTERIOR PROBABILITY DISTRIBUTIONS ####

# calculating NLL, effective sample size, acceptance fraction, psrf

mcmcout=cbind(time_seq,d238Usw_quant_out[,50])
colnames(mcmcout)=c("time","d238Usw")
mergemcmc=merge(mcmcout,current_dataset,by="time") 
total_NLL=round((sum(0.5*log(2*pi*mergemcmc$err))+sum(na.omit(((mergemcmc$d238U-mergemcmc$d238Usw)^2)/(mergemcmc$err^2)))),digits=1)  # final NLL of median output
accept_frac=100*round(length(unique(finalMCMC[,1]))/length(finalMCMC[,1]),digits=2) # fraction of steps accepted during random walks
z=1
y=1
auto_out=NULL
big_auto_out=NULL
mcmc_list_out=NULL
for (z in 1:n_walkers){
	cur_walker_data=finalMCMC[(1+(((niterMCMC-burninlengthMCMC)/thin)*(z-1))):(((niterMCMC-burninlengthMCMC)/thin)*(z)),]
	
	y=1
	for (y in 1:(dim(cur_walker_data)[2]-2)){
		cur_autocorr=IAT(as.numeric(cur_walker_data[,y]))
		auto_out=rbind(auto_out,cur_autocorr)
	}
	
	cur_avg_autocorr=mean(auto_out)
	cur_eff=((niterMCMC-burninlengthMCMC)/thin)/cur_avg_autocorr
	cur_eff_size=mean(effectiveSize(cur_walker_data))
	results=c(cur_avg_autocorr, cur_eff,cur_eff_size)
	big_auto_out=rbind(big_auto_out, results)
	
	mcmc_results=cbind(cur_walker_data,z)
	mcmc_list_out=rbind(mcmc_list_out,mcmc_results)

	
}
avg_autocorr_MCMC=mean(big_auto_out[,1]) # average integrated autocorrelation time (IAT)
eff=round(sum(big_auto_out[,2]),digits=0) # effective sample size (total samples/IAT)
mcmc_out_1=mcmc(mcmc_list_out[mcmc_list_out[,dim(mcmc_list_out)[2]]==1,],thin=thin)
mcmc_out_1= mcmc_out_1[,1:length(pars1)]
mcmc_out_2=mcmc(mcmc_list_out[mcmc_list_out[,dim(mcmc_list_out)[2]]==2,],thin=thin)
mcmc_out_2= mcmc_out_2[,1:length(pars1)]
mcmc_out_3=mcmc(mcmc_list_out[mcmc_list_out[,dim(mcmc_list_out)[2]]==3,],thin=thin)
mcmc_out_3= mcmc_out_3[,1:length(pars1)]
mcmc_list=as.mcmc.list(list(mcmc_out_1, mcmc_out_2, mcmc_out_3)) 
psrf=gelman.diag(mcmc_list,autoburnin=FALSE) # calculating potential scale reduction factor (psrf), i.e. Gelman-Rubin statistic, or Rhat
avg_psrf=round(mean(psrf$psrf[,1]),digits=2) # average psrf for all variables (gives quick check of convergence; aiming for <<1.2)


# plotting recovered time series

par(mfrow=c(4,1))
par(mar=c(0.5,1,0.5,1))
par(oma=c(2.5,3.5,5,0))
plot(current_dataset$d238U~current_dataset$time,ylim=c(-1.5,0.2),cex=0,axes=F,xlab="",ylab="")
polygon(c(time_seq,rev(time_seq)),c(d238Usw_quant_out[,16],rev(d238Usw_quant_out[,84])),border=F,col="grey80")
lines(d238Usw_quant_out[,50]~ time_seq,col="black")
arrows(current_dataset$time,current_dataset$d238U+current_dataset$err,current_dataset$time,current_dataset$d238U-current_dataset$err,length=0)
points(current_dataset$d238U~current_dataset$time,pch=21,bg="white",cex=1.0)
axis(side=2,tck=-0.02,mgp=c(3,0.5,0),las=1,)
axis(side=1,tck=-0.02,mgp=c(3,0.5,0),las=1,labels=NA)
mtext(expression(delta^{238}*"U [\u2030]"),side=2,line=2.5,cex=0.9)
box(lwd=1)
mtext(side=3,line=3.5,"NLL:",at=min(current_dataset$time),adj=0,cex=0.8)
mtext(side=3,line=3.5,paste(total_NLL),at=max(current_dataset$time),adj=1,cex=0.8)
mtext(side=3,line=2.5,"effective sample size:",at=min(current_dataset$time),adj=0,cex=0.8)
mtext(side=3,line=2.5,paste(eff),at=max(current_dataset$time),col=if((eff>=1000)){"green4"} else if(eff>=100){"orange"} else{"red"},adj=1,cex=0.8)
mtext(side=3,line=1.5,"fraction accepted:",at=min(current_dataset$time),adj=0,cex=0.8)
mtext(side=3,line=1.5,paste(accept_frac,"%"),at=max(current_dataset$time),col=if((accept_frac>=15)&(accept_frac<=40)){"green4"} else if((accept_frac<15)&(accept_frac>=10)){"orange"} else if((accept_frac>40)&(accept_frac<=50)){"orange"} else{"red"},adj=1,cex=0.8)
mtext(side=3,line=0.5,"avg psrf:",at=min(current_dataset$time),adj=0,cex=0.8)
mtext(side=3,line=0.5,paste(avg_psrf),at=max(current_dataset$time),col=if(avg_psrf<=1.1){"green4"} else if((avg_psrf>1.1)&(avg_psrf<=1.2)){"orange"} else{"red"},adj=1,cex=0.8)

plot(Nsw_quant_out[,50]~time_seq,ylim=c(0,3e13),cex=0,axes=F,xlab="",ylab="")
polygon(c(time_seq,rev(time_seq)),c((Nsw_quant_out[,16]),rev(Nsw_quant_out[,84])),border=F,col="grey80")
lines(Nsw_quant_out[,50]~time_seq,col="black")
axis(side=2,tck=-0.02,mgp=c(3,0.5,0),las=1,at=c(0,10e12,20e12,30e12),labels=c("0","10","20","30"))
axis(side=1,tck=-0.02,mgp=c(3,0.5,0),las=1,labels=NA)
mtext(expression(italic("N"[sw])*" [Tmol]"),side=2,line=2.5,cex=0.9)
box(lwd=1)

plot(10^fanox_quant_out[,4]~time_seq,ylim=c(0.01,100),cex=0,axes=F,xlab="",ylab="")
Fanox_quant_out=(10^fanox_quant_out*Kanox)/(10^fanox_quant_out*Kanox+Kother*(1-10^fanox_quant_out))
polygon(c(time_seq,rev(time_seq)),c((100*Fanox_quant_out[,16]),rev(100*Fanox_quant_out[,84])),border=F,col="grey80")
lines(100*Fanox_quant_out[,50]~time_seq,col="black")
axis(side=2,tck=-0.02,mgp=c(3,0.5,0),las=1)
axis(side=1,tck=-0.02,mgp=c(3,0.5,0),las=1,labels=NA)
mtext(expression(italic("F"[anox])*" [%]"),side=2,line=2.5,cex=0.9)
box(lwd=1)

plot(10^fanox_quant_out[,4]~time_seq,ylim=c(0.01,100),cex=0,axes=F,xlab="",ylab="")
polygon(c(time_seq,rev(time_seq)),c((100*10^fanox_quant_out[,16]),rev(100*10^fanox_quant_out[,84])),border=F,col="grey80")
lines(100*10^fanox_quant_out[,50]~time_seq,col="black")
axis(side=2,tck=-0.02,mgp=c(3,0.5,0),las=1)
axis(side=1,tck=-0.02,mgp=c(3,0.5,0),las=1)
mtext(expression(italic("f"[anox])*" [%]"),side=2,line=2.5,cex=0.9)
mtext("time [yr]",side=1,line=1.75,cex=0.9)
box(lwd=1)

