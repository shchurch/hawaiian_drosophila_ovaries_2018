library(corHMM)
library(OUwie)
library(surface)
library(nlme)
library(dplyr)
library(phytools)
set.seed(1234)

####################################################
# BUILD THE DATATABLES 

dat <- read.table("data.txt",header=T) 
groups <- dat %>% distinct(species,group)
dat <- dat %>% 	group_by(species) %>% 
	summarize_if(is.numeric,funs(mean),na.rm=T)

dat <- left_join(dat,groups,by="species")

# Load ecological data, sort into categories of oviposition substrate
eco <- read.table("ecological_classifications.txt",header=T)
# eco$reg = 8 regime model
# eco$reg2 = 3 regime model
# eco$reg3 = 2 regime model
eco$reg2 <- plyr::mapvalues(eco$reg,c("Fungus","Fruit","Generalist","Leaf","Flower","Sap","Spider"),c("gen","gen","gen","gen","Flower_Spider","gen","Flower_Spider"))
eco$reg3 <- plyr::mapvalues(eco$reg,c("Fungus","Fruit","Generalist","Leaf","Flower","Sap","Spider"),c("gen","gen","gen","gen","gen","gen","gen"))

# Merge these two data tables by species ID
mg <- merge(dat,eco,by="species")

# Write average tables, and log transform the data and write a table
write.table(mg,file="avgs.txt",row.names=FALSE,sep="\t")
avgs <- read.delim("avgs.txt")
# Take the natural log of each numerical column
logs <- avgs %>% mutate(
	thor = log(thorax),
	thor_vol = log(thor_vol),
	onum = log(onum),
	dors = log(dors),
	egg_vol = log(egg_vol),
	major = log(major),
	minor = log(minor),
	prop = log(prop),
	tf_num = log(tf_num),
	tfc_tf = log(tfc_tf),
	total_tfc = log(total_tfc))

write.table(logs,file="logs.txt",row.names=FALSE,sep="\t")

onum_logs <- logs %>% mutate(trait = onum)
prop_logs <- logs %>% mutate(trait = prop) %>% filter(!(is.infinite(trait)))
thvl_logs <- logs %>% mutate(trait = thor_vol)
eggs_logs <- logs %>% mutate(trait = egg_vol)

### This reads in the posterior distributions of phylogenetic trees, and then samples n number of trees for uncertainty estimates
sample_size <- 1000

raw_all_trees <- read.nexus("nuc_mt_beast_posterior.trees")
all_trees <- sample(raw_all_trees,sample_size)

### This code calculates the phylogenetic residuals of log egg size and log ovariole number to log thorax volume, using one of the all_trees as a representative phylogeny 
reduced_logs <- logs %>% 
			filter(species %in% all_trees[[1]]$tip.label) %>%
			mutate(trait1 = egg_vol, trait2 = onum, indep = thor_vol) %>% 
			select(species,trait1,trait2,indep) %>% 
			na.omit()
Y <- data.frame(trait1 = reduced_logs$trait1, trait2 = reduced_logs$trait2 ,row.names=reduced_logs$species)
X <- reduced_logs$indep
names(X) <- reduced_logs$species
pruned <- drop.tip(all_trees[[1]],tip=setdiff(all_trees[[1]]$tip.label,reduced_logs$species))
resid <- data.frame((phyl.resid(pruned,X,Y,method="BM")$resid))
logs <- left_join(logs,tibble::rownames_to_column(resid,var = "species"),by = "species") %>% rename(eggs_resid = trait1, onum_resid = trait2)
eggs_resid_logs <- logs %>% mutate(trait = eggs_resid)

# CHOOSE THE TREE FOR WHICH YOU WOULD LIKE TO GENERATE THE PLOTS
phy <- all_trees[[1]]


####################################################
# GENERATE PLOTS FROM ANCESTRAL STATE RECONSTRUCTION

# onum_ouwie8 <- subset(logs,select = c(species, reg, onum))
# onum_ouwie3 <- subset(logs,select = c(species, reg2, onum))
# onum_ouwie2 <- subset(logs,select = c(species, reg3, onum))
# onum_ouwie8 <- na.omit(onum_ouwie8)
# onum_ouwie3 <- na.omit(onum_ouwie3)
# onum_ouwie2 <- na.omit(onum_ouwie2)
# onum_ouwie8$reg <- as.numeric(onum_ouwie8$reg)
# onum_ouwie3$reg2 <- as.numeric(onum_ouwie3$reg2)
# onum_ouwie2$reg3 <- as.numeric(onum_ouwie2$reg3)
# 
# dat <- onum_ouwie8
# eco_reg <- na.omit(data.frame(species = eco$species,reg=as.numeric(eco$reg))) 
# pp_pruned <- drop.tip(phy,setdiff(phy$tip.label,eco_reg$species))
# pp <- rayDISC(pp_pruned,eco_reg[,c(1,2)],model="ER",node.states="marginal")
# pruned <- drop.tip(pp$phy,setdiff(pp$phy$tip.label,dat$species))
# 
# dat2 <- onum_ouwie3
# eco_reg2 <- na.omit(data.frame(species = eco$species,reg=as.numeric(eco$reg2)))
# pp_pruned2 <- drop.tip(phy,setdiff(phy$tip.label,eco_reg2$species))
# pp2 <- rayDISC(pp_pruned2,eco_reg2[,c(1,2)],model="ER",node.states="marginal")
# pruned2 <- drop.tip(pp2$phy,setdiff(pp2$phy$tip.label,dat2$species))
# 
# dat3 <- onum_ouwie2
# eco_reg3 <- na.omit(data.frame(species = eco$species,reg=as.numeric(eco$reg3)))
# pp_pruned3 <- drop.tip(phy,setdiff(phy$tip.label,eco_reg3$species))
# pp3 <- rayDISC(pp_pruned3,eco_reg3[,c(1,2)],model="ER",node.states="marginal")
# pruned3 <- drop.tip(pp3$phy,setdiff(pp3$phy$tip.label,dat3$species))
# 
# pdf(file="ecological_state_reconstruction_8state.pdf",width=12,height=12)
# plotRECON(pp$phy,pp$states)
# dev.off()
# pdf(file="ecological_state_reconstruction_3state.pdf",width=12,height=12)
# plotRECON(pp2$phy,pp2$states)
# dev.off()
# pdf(file="ecological_state_reconstruction_2state.pdf",width=12,height=12)
# plotRECON(pp3$phy,pp3$states)
# dev.off()
# pdf(file="tree_for_ou_analysis_8state.pdf",width=12,height=12)
# plot(pruned,show.node.label=T)
# dev.off()
# pdf(file="tree_for_ou_analysis_3state.pdf",width=12,height=12)
# plot(pruned2,show.node.label=T)
# dev.off()
# pdf(file="tree_for_ou_analysis_2state.pdf",width=12,height=12)
# plot(pruned3,show.node.label=T)
# dev.off()

####################################################
# GENERATE PLOTS FROM PGLS ANALYSIS

### SETUP PGLS OVER MANY TREES

pgls_plots <- function(logs,name) {
	onum_thvl_pgls <- data.frame(group=logs$group, onum=logs$onum,thor_vol=logs$thor_vol,row.names=logs$species) %>% na.omit()
	onum_eggs_pgls <- data.frame(group=logs$group,onum=logs$onum,eggs=logs$egg_vol,row.names=logs$species) %>% na.omit()
	eggs_thvl_pgls <- data.frame(group=logs$group,eggs=logs$egg_vol,thor_vol=logs$thor_vol,row.names=logs$species) %>% na.omit()
	onum_prop_pgls <- logs %>% select(group,onum,prop,species) %>% filter(is.finite(prop))
	onum_prop_pgls <- data.frame(group=onum_prop_pgls$group,onum=onum_prop_pgls$onum, prop=onum_prop_pgls$prop, row.names = onum_prop_pgls$species) %>% na.omit()
	onum_eggs_resid_pgls <- data.frame(group=logs$group,onum_resid=logs$onum_resid,eggs_resid=logs$eggs_resid,row.names=logs$species) %>% na.omit()

	### PLOT PGLS FOR OVARIOLE NUMBER AND THORAX VOLUME

	dat <- onum_thvl_pgls
	form <- "onum~thor_vol"
	
	outskis <- pgls.trees(phy,dat,form)
	all_out <- lapply(all_trees,pgls.trees,dat = dat,form = form)

	pdf(file=paste(name,"ovariole_number_thorax_volume.pdf",sep="_"),useDingbats=F)
	plot(dat$thor_vol,dat$onum,type="n",xlab="log( Thorax Volume (mm^3) )",ylab="log( Ovariole Number )")
	for (i in 1:length(all_out)){
		abline(all_out[[i]],col='grey')
	}
	if (max(sapply(all_out,function(x) {x$tTable[8]})) > 0.05) {
		abline(outskis,col='black')
	} else {
		abline(outskis,col='red')
	}
	points(dat$thor_vol,dat$onum,col=c("red","blue","green","black","cyan")[dat$group],pch=16)
	#legend(x="top", inset=c(-0.2,0), legend = levels(dat$group), col=c("red","blue","green","black","cyan"), pch=16)
	dev.off()
	
	### PLOT PGLS FOR OVARIOLE NUMBER AND EGG VOLUME
	
	dat <- onum_eggs_pgls
	form <- "onum~eggs"
	
	outskis <- pgls.trees(phy,dat,form)
	all_out <- lapply(all_trees,pgls.trees,dat = dat,form = form)
	
	pdf(file=paste(name,"ovariole_number_egg_volume.pdf",sep="_"),useDingbats=F)
	plot(dat$eggs,dat$onum,type="n",xlab="log( Egg Volume (µm^3) )",ylab="log( Ovariole Number )")
	for (i in 1:length(all_out)){
		abline(all_out[[i]],col='grey')
	}
	if (max(sapply(all_out,function(x) {x$tTable[8]})) > 0.05) {
		abline(outskis,col='black')
	} else {
		abline(outskis,col='red')
	}
	points(dat$eggs,dat$onum,col=c("red","blue","green","black","cyan")[dat$group],pch=16)
	#legend(x="top", legend = levels(dat$group), col=c("red","blue","green","black","cyan"), pch=16)
	dev.off()
	
	### PLOT PGLS FOR EGG VOLUME AND THORAX VOLUME
	
	dat <- eggs_thvl_pgls
	form <- "eggs~thor_vol"
	
	outskis <- pgls.trees(phy,dat,form)
	all_out <- lapply(all_trees,pgls.trees,dat = dat,form = form)
	
	pdf(file=paste(name,"egg_volume_thorax_volume.pdf",sep="_"),useDingbats=F)
	plot(dat$thor_vol,dat$eggs,type="n",xlab="log( Thorax Volume (mm^3) )",ylab="log( Egg Volume (µm^3) )")
	for (i in 1:length(all_out)){
		abline(all_out[[i]],col='grey')
	}
	abline(outskis,col='red')
	points(dat$thor_vol,dat$eggs,col=c("red","blue","green","black","cyan")[dat$group],pch=16)
	#legend(x="top", legend = levels(dat$group), col=c("red","blue","green","black","cyan"), pch=16)
	dev.off()
	
	### PLOT PGLS FOR OVARIOLE NUMBER AND PROPORTIONAL EGG SIZE
	
	dat <- onum_prop_pgls
	form <- "onum~prop"
	
	outskis <- pgls.trees(phy,dat,form)
	all_out <- lapply(all_trees,pgls.trees,dat = dat,form = form)
	
	pdf(file=paste(name,"ovariole_number_egg_thorax_proportion.pdf",sep="_"),useDingbats=F)
	plot(dat$prop,dat$onum,type="n",xlab="log( Proportional Egg Volume )",ylab="log( Ovariole Number )")
	for (i in 1:length(all_out)){
		abline(all_out[[i]],col='grey')
	}
	if (max(sapply(all_out,function(x) {x$tTable[8]})) > 0.05) {
		abline(outskis,col='black')
	} else {
		abline(outskis,col='red')
	}
		points(dat$prop,dat$onum,col=c("red","blue","green","black","cyan")[dat$group],pch=16)
	#legend(x="top", legend = levels(dat$group), col=c("red","blue","green","black","cyan"), pch=16)
	dev.off()

	
	### PLOT PGLS FOR OVARIOLE NUMBER AND PROPORTIONAL EGG SIZE
	
	dat <- onum_eggs_resid_pgls
	form <- "onum_resid~eggs_resid"
	
	outskis <- pgls.trees(phy,dat,form)
	all_out <- lapply(all_trees,pgls.trees,dat = dat,form = form)
	
	pdf(file=paste(name,"ovariole_number_egg_volume_residuals_to_thorax_volume.pdf",sep="_"),useDingbats=F)
	plot(dat$eggs_resid,dat$onum_resid,type="n",xlab="Phy. res., egg vol to thorax vol",ylab="Phy. res., ovariole num. to thorax vol")
	for (i in 1:length(all_out)){
		abline(all_out[[i]],col='grey')
	}
	if (max(sapply(all_out,function(x) {x$tTable[8]})) > 0.05) {
		abline(outskis,col='black')
	} else {
		abline(outskis,col='red')
	}
		points(dat$eggs_resid,dat$onum_resid,col=c("red","blue","green","black","cyan")[dat$group],pch=16)
	#legend(x="top", legend = levels(dat$group), col=c("red","blue","green","black","cyan"), pch=16)
	dev.off()

}

# This computes PGLS for two measurements over many trees
pgls.trees <- function(phy,dat,form){
	# The full tree is pruned to those with both reproductive traits
	pruned <- drop.tip(phy,setdiff(phy$tip.label,rownames(dat)))
	# The dataset is trimmed to only include those with data, just in case
	red_dat <- dat[!rownames(dat) %in% setdiff(rownames(dat),pruned$tip.label),]
	# PGLS is performed using an OU model and with a starting alpha very smal
	lm_str <- corMartins(0.000000001,phy=pruned)
	# the results are summarized
	outskis <- summary(gls(as.formula(form),correlation = lm_str,data = red_dat,method="ML"))
	return(outskis)
}


### GENERATE THE PLOTS FOR EACH GROUPING

# all groups
pgls_plots(logs,"all")

# scaptomyza
scapto_logs <- subset(logs,group=="scapto")
pgls_plots(scapto_logs,"scapto")

# amc
amc_logs <- subset(logs,group=="amc")
pgls_plots(amc_logs,"amc")

# pna
pna_logs <- subset(logs,group=="pna")
pgls_plots(pna_logs,"pna")

# mouth
mouth_logs <- subset(logs,group=="mouth")
pgls_plots(mouth_logs,"mouth")




