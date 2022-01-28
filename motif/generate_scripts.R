##### filter motifs with expression >= 1 in oocytes and zygotes #####

factors<-read.table('motif.txt')
expr <- read.table('mRNA_expression.txt',header=T)
exprTF <- expr[which(expr[,2] %in% unique(factors[,1])),]
exprTF <- exprTF[which(apply(exprTF[,3:4],1,max)>=1),]
factors<-factors[which(factors[,1] %in% exprTF[,2]),]

##### generate codes for heatmapr #####
code1 <- paste("mkdir", unique(factors[,1]))
code2 <- paste("sort -k5,5gr /mnt/Storage/home/wangchenfei/annotations/motif/cistrome_v2_mm9_uniq/",factors[,3]," | head -10000 > ",factors[,1],"/mm9.",factors[,1],".",factors[,2],".motif",sep="")
code3 <- paste("heatmapr -w /mnt/Storage/home/wangchenfei/MPronucleusNucleosome/data.bwfiles/merged/occupancy/MalePSperm.Q10.uniq.nucleo.37bp.bw,/mnt/Storage/home/wangchenfei/MPronucleusNucleosome/data.bwfiles/merged/occupancy/MaleP0.5.Q10.uniq.nucleo.37bp.bw,/mnt/Storage/home/wangchenfei/MPronucleusNucleosome/data.bwfiles/merged/occupancy/MaleP1.Q10.uniq.nucleo.37bp.bw,/mnt/Storage/home/wangchenfei/MPronucleusNucleosome/data.bwfiles/merged/occupancy/MaleP1.5.Q10.uniq.nucleo.37bp.bw,/mnt/Storage/home/wangchenfei/MPronucleusNucleosome/data.bwfiles/merged/occupancy/MaleP2.Q10.uniq.nucleo.37bp.bw,/mnt/Storage/home/wangchenfei/MPronucleusNucleosome/data.bwfiles/merged/occupancy/MaleP3.Q10.uniq.nucleo.37bp.bw,/mnt/Storage/home/wangchenfei/MPronucleusNucleosome/data.bwfiles/merged/occupancy/MaleP4.Q10.uniq.nucleo.37bp.bw,/mnt/Storage/home/wangchenfei/MPronucleusNucleosome/data.bwfiles/merged/occupancy/MaleP6.Q10.uniq.nucleo.37bp.bw,/mnt/Storage/home/wangchenfei/MPronucleusNucleosome/data.bwfiles/merged/occupancy/MaleP8.Q10.uniq.nucleo.37bp.bw,/mnt/Storage/home/wangchenfei/MPronucleusNucleosome/data.bwfiles/merged/occupancy/MaleP12.Q10.uniq.nucleo.37bp.bw -b /mnt/Storage/home/wangchenfei/MPronucleusNucleosome/Fig2.MechanismAnalysis/motifAnalysis/motif/",factors[,1],"/mm9.",factors[,1],".",factors[,2],".motif --name=",factors[,1],"/MalePN_",factors[,1],".",factors[,2]," --method=mean --s_wigindex=10 --dir --upstream=1000 --downstream=1000 &",sep="")
code4 <- paste("heatmapr -w /mnt/Storage/home/wangchenfei/MPronucleusNucleosome/data.bwfiles/merged/occupancy/FemalePMII.Q10.uniq.nucleo.37bp.bw,/mnt/Storage/home/wangchenfei/MPronucleusNucleosome/data.bwfiles/merged/occupancy/FemaleP0.5.Q10.uniq.nucleo.37bp.bw,/mnt/Storage/home/wangchenfei/MPronucleusNucleosome/data.bwfiles/merged/occupancy/FemaleP1.Q10.uniq.nucleo.37bp.bw,/mnt/Storage/home/wangchenfei/MPronucleusNucleosome/data.bwfiles/merged/occupancy/FemaleP1.5.Q10.uniq.nucleo.37bp.bw,/mnt/Storage/home/wangchenfei/MPronucleusNucleosome/data.bwfiles/merged/occupancy/FemaleP2.Q10.uniq.nucleo.37bp.bw,/mnt/Storage/home/wangchenfei/MPronucleusNucleosome/data.bwfiles/merged/occupancy/FemaleP3.Q10.uniq.nucleo.37bp.bw,/mnt/Storage/home/wangchenfei/MPronucleusNucleosome/data.bwfiles/merged/occupancy/FemaleP4.Q10.uniq.nucleo.37bp.bw,/mnt/Storage/home/wangchenfei/MPronucleusNucleosome/data.bwfiles/merged/occupancy/FemaleP6.Q10.uniq.nucleo.37bp.bw,/mnt/Storage/home/wangchenfei/MPronucleusNucleosome/data.bwfiles/merged/occupancy/FemaleP8.Q10.uniq.nucleo.37bp.bw,/mnt/Storage/home/wangchenfei/MPronucleusNucleosome/data.bwfiles/merged/occupancy/FemaleP12.Q10.uniq.nucleo.37bp.bw -b /mnt/Storage/home/wangchenfei/MPronucleusNucleosome/Fig2.MechanismAnalysis/motifAnalysis/motif/",factors[,1],"/mm9.",factors[,1],".",factors[,2],".motif --name=",factors[,1],"/FemalePN_",factors[,1],".",factors[,2]," --method=mean --s_wigindex=10 --dir --upstream=1000 --downstream=1000 &",sep="")
code5 <- paste("heatmapr -w /mnt/Storage/home/wangchenfei/MPronucleusNucleosome/data.bwfiles/merged/occupancy/early2cell.Q10.uniq.nucleo.37bp.bw,/mnt/Storage/home/wangchenfei/MPronucleusNucleosome/data.bwfiles/merged/occupancy/late2cell.Q10.uniq.nucleo.37bp.bw,/mnt/Storage/home/wangchenfei/MPronucleusNucleosome/data.bwfiles/merged/occupancy/4cell.Q10.uniq.nucleo.37bp.bw,/mnt/Storage/home/wangchenfei/MPronucleusNucleosome/data.bwfiles/merged/occupancy/8cell.Q10.uniq.nucleo.37bp.bw,/mnt/Storage/home/wangchenfei/MPronucleusNucleosome/data.bwfiles/merged/occupancy/morula.Q10.uniq.nucleo.37bp.bw,/mnt/Storage/home/wangchenfei/MPronucleusNucleosome/data.bwfiles/merged/occupancy/ICM.Q10.uniq.nucleo.37bp.bw,/mnt/Storage/home/wangchenfei/MPronucleusNucleosome/data.bwfiles/merged/occupancy/TE.Q10.uniq.nucleo.37bp.bw -b /mnt/Storage/home/wangchenfei/MPronucleusNucleosome/Fig2.MechanismAnalysis/motifAnalysis/motif/",factors[,1],"/mm9.",factors[,1],".",factors[,2],".motif --name=",factors[,1],"/Embryo_",factors[,1],".",factors[,2]," --method=mean --s_wigindex=7 --dir --upstream=1000 --downstream=1000 &",sep="")

codeall <- data.frame(code3,code4,code5)
codeall <- c(t(as.matrix(code)))
code <- c()
n<-0
for (i in 1:length(codeall)){
	n <- n+1
	code<- c(code,codeall[i])
	if (n %% 12 ==0){
		code <- c(code,"wait;\n")
	}
}
code<-c(code1,code2,code)
write.table(code,'motifAnalysis.sh',col.name=F,row.name=F,quote=F)

code1 <- paste("write_r_embryo('",factors[,1],"','",factors[,1],".",factors[,2],"','/mnt/Storage/home/wangchenfei/MPronucleusNucleosome/Fig2.MechanismAnalysis/motifAnalysis/motif/",factors[,1],"/Embryo_",factors[,1],".",factors[,2],"_kmeans.r')",sep="")
code2 <- paste("write_r_male('",factors[,1],"','",factors[,1],".",factors[,2],"','/mnt/Storage/home/wangchenfei/MPronucleusNucleosome/Fig2.MechanismAnalysis/motifAnalysis/motif/",factors[,1],"/MalePN_",factors[,1],".",factors[,2],"_kmeans.r')",sep="")
code3 <- paste("write_r_female('",factors[,1],"','",factors[,1],".",factors[,2],"','/mnt/Storage/home/wangchenfei/MPronucleusNucleosome/Fig2.MechanismAnalysis/motifAnalysis/motif/",factors[,1],"/FemalePN_",factors[,1],".",factors[,2],"_kmeans.r')",sep="")
code<-c(code1,code2,code3)
write.table(code,'motifAnalysis.py',col.name=F,row.name=F,quote=F)

##### merge profile and calculate deplete score #####
merge_profile <- function(factor,name)
{
	inpath1 <- paste('/mnt/Storage2/home/wangchenfei/MPronucleusNucleosome/Fig4.TransFactorAnalysis/motifAnalysis/motif/',factor,'/MalePN_',name,'_siteprof',c(1:10),sep='')
	inpath2 <- paste('/mnt/Storage2/home/wangchenfei/MPronucleusNucleosome/Fig4.TransFactorAnalysis/motifAnalysis/motif/',factor,'/FemalePN_',name,'_siteprof',c(1:10),sep='')
	inpath3 <- paste('/mnt/Storage2/home/wangchenfei/MPronucleusNucleosome/Fig4.TransFactorAnalysis/motifAnalysis/motif/',factor,'/Embryo_',name,'_siteprof',c(1:7),sep='')
	name1 <- paste("Male",c(1:10),sep="")
	name2 <- paste("Female",c(1:10),sep="")
	name3 <- paste("Embryo",c(1:7),sep="")
	names <- c(name1,name2,name3)
	inpath <- c(inpath1,inpath2,inpath3)
	for (i in c(1:length(inpath))){
		data <- read.table(inpath[i],sep=',',header=F)
		signal <- apply(data,2,mean)
		if (i==1){
			datafull <- signal
    	} else {
    		datafull <- rbind(datafull, signal)
    	}
	}
	rownames(datafull) <- names
	return (datafull)
}
score_calculate <- function(signal) 
{
	signal_smooth<-smooth.spline(signal,penalty=5)$y
	minusone = max(signal_smooth[75:90])
	plusone = max(signal_smooth[110:125])
	center = mean(signal_smooth[99:102])
	if (mean(signal_smooth)==0) {score<-0} else {score<-(max(plusone,minusone)-center)/(max(signal_smooth)-min(signal_smooth))}
	return (round(score,3))
}
plot_sitepro <- function(mat,factor,name)
{
	outfile1 <- paste('/mnt/Storage2/home/wangchenfei/MPronucleusNucleosome/Fig4.TransFactorAnalysis/motifAnalysis/motif/',factor,'/',name,'_sitepro.pdf',sep='')
	outfile2 <- paste('/mnt/Storage2/home/wangchenfei/MPronucleusNucleosome/Fig4.TransFactorAnalysis/motifAnalysis/motif/',factor,'/',name,'_deplete_score.txt',sep='')
	all_score <- c()
	pdf(outfile1,width=18,height=7.5)
	par(oma=c(1, 2, 3, 0))
	par(mar=c(1.5, 2, 2, 0.6))
	layout(t(matrix(seq(30), nrow=10, ncol=3)), widths=rep(12/10,10), heights=rep(c(1,1,1),10))
	par(cex=1.2)
	for (i in c(1:nrow(mat))){
		stage_score <- score_calculate(c(mat[i,]))
		all_score <- c(all_score,stage_score)
		plot(c(0:200),c(mat[i,])/mean(c(mat[i,])),type='l',col="red",lwd=2,ylab=NA,xlab="Distance to center",axes=F,ylim=c(0,2.5),cex.lab=0.5)
    	axis(side=1,at=c(0,100,200),c("-1000","center","1000"),cex.axis=0.65);axis(side=2,cex.axis=0.65);box()
    	title(main=samples[i],cex.main=0.8)
    	abline(v=101,lty=3)
    	text(x=150,y=2,stage_score,col='Blue',cex=0.5)
	}
	mtext(paste("Nucleosome profile on ",name," motif",sep=""), side = 3, line = 1, outer = TRUE, cex = 2)
	mtext("Normalized nucleosome tags", side = 2, line = 0.5, outer = TRUE, cex = 1.2)
	dev.off()
	dep_score <- cbind(samples,all_score)
	write.table(dep_score,outfile2,col.name=F,row.name=F,sep='\t',quote=F)
}
samples <- c("Sperm","M0.5h","M1h","M1.5h","M2h","M3h","M4h","M6h","M8h","M12h","MIIOocyte","F0.5h","F1h","F1.5h","F2h","F3h","F4h","F6h","F8h","F12h","early2cell","late2cell","4cell","8cell","morula","ICM","TE")

for (i in 1:nrow(factors))
{
	profile<-merge_profile(factors[i,1],paste(factors[i,1],factors[i,2],sep='.'))
	plot_sitepro(profile,factors[i,1],paste(factors[i,1],factors[i,2],sep='.'))
}

profile<-merge_profile("Duxbl","Duxbl.M01390")
plot_sitepro(profile,"Duxbl","Duxbl.M01390")

##### clustering depleted score #####
dep_mat<-NULL
for (i in 1:nrow(factors)){
	infile <- paste('/mnt/Storage2/home/wangchenfei/MPronucleusNucleosome/Fig4.TransFactorAnalysis/motifAnalysis/motif/',factors[i,1],'/',paste(factors[i,1],factors[i,2],sep='.'),'_deplete_score.txt',sep='')
	data <- read.table(infile,sep='\t')
	dep_mat <- cbind(dep_mat,data[,2])
	colnames(dep_mat)[ncol(dep_mat)] <- as.character(paste(factors[i,1],factors[i,2],sep='.'))
}
rownames(dep_mat)<-samples
write.table(t(dep_mat),'depletion_matrix.txt',col.name=T,row.name=T,quote=F,sep='\t')

setwd("/Volumes/SBPD/Archiving/Development/MPronucleusNucleosome/results/Fig2.MechanismAnalysis/motifAnalysis/")
samples <- c("Sperm","M0.5h","M1h","M1.5h","M2h","M3h","M4h","M6h","M8h","M12h","MIIOocyte","F0.5h","F1h","F1.5h","F2h","F3h","F4h","F6h","F8h","F12h","early2cell","late2cell","4cell","8cell","morula","ICM","TE")
dep_mat<-read.table('depletion_matrix_uniq.txt',header=T)
km<-kmeans(dep_mat[,3:22],3)
data<-cbind(dep_mat[,3:22],km=km$cluster)
n_total<-NULL;n<-0
for(i in c(1:3)) {n<-n+nrow(data[which(data[,21]==i),]);n_total<-cbind(n_total,n)}
data <- data[order(data[,21]),]

png('depletion_matrix_cluster.png',width=800,height=1600)
par(oma = c(0, 0, 3, 0))
par(mar=c(1.5,15,10,1.5))

zmin<- -1
zmax<- 1
ColorRamp <- colorRampPalette(c("blue","white","red"), bias=1)(100000)
ColorLevels <- seq(to=zmax,from=zmin, length=10000)   #number sequence
image(1:ncol(data[,1:20]), 1:nrow(data[,1:20]), t(data[,1:20]), axes=F, col=ColorRamp, xlab="", ylab="")
axis(side=3,1:length(samples[1:20]),labels=samples[1:20],cex.axis=2,las=2)
axis(side=2,1:nrow(data),labels=rownames(data),cex.axis=1.5,las=2)
for (i in seq(1,2)) {abline(h=n_total[i]+0.5,lty=2,lwd=2,col="black")}
abline(v=10+0.5,lty=2,lwd=2,col="black");box()
dev.off()

G1 <- rownames(data[which(data[,21]==1),])
G2 <- rownames(data[which(data[,21]==2),])
G3 <- rownames(data[which(data[,21]==3),])
write.table(G1,"depletion_matrix_uniq1.txt",col.name=F,row.name=F,sep='\t',quote=F)
write.table(G2,"depletion_matrix_uniq2.txt",col.name=F,row.name=F,sep='\t',quote=F)
write.table(G3,"depletion_matrix_uniq3.txt",col.name=F,row.name=F,sep='\t',quote=F)


