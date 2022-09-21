library(ggplot2)
#read in data frames
hunt<-fread(file="best/HUNT-LDL-preMeta.txt")
bbj<-fread(file="best/BBJ-LDL-preMeta.txt")
glgc<-fread(file="best/GLGC-LDL-hg38-preMeta.txt")
meta<-fread(file="METAANALYSIS1.TBL")

#read in data for rs61679753 in APOE
df<-rbind(glgc[glgc$SNPID=="19:44897490:A:T",],
          bbj[bbj$SNPID=="19:44897490:A:T",],
          hunt[hunt$SNPID=="19:44897490:A:T",],fill=TRUE)
df<-cbind(df,data.frame(study=c("GLGC","BBJ","HUNT")))

#get pooled beta estimate
df$w<-1/(df$SE*df$SE)
df$wB<-df$w*df$BETA
beta_pooled<-sum(df$wB)/sum(df$w)
beta_pooled_se<-sqrt(1/sum(df$w))

#read in meta-analysis data
df2<-meta[meta$MarkerName=="19:44897490:A:T",]
df2$study<-"meta-analysis"
df2$BETA<-beta_pooled
df2$SE<-beta_pooled_se
df<-rbind(df,df2,fill=TRUE)

#make 95% Confidence intervals
df$UB<-df$BETA+(1.96*df$SE)
df$LB<-df$BETA-(1.96*df$SE)
df<-df%>%mutate(study= fct_relevel(study,"meta-analysis","HUNT","BBJ","GLGC") )

#make plot
pdf(file="ForestPlot.pdf",height=4,width=5)
ggplot(df,aes(x=BETA,y=study,color=study)) + geom_point() + theme_bw() +
  geom_vline(xintercept=0,linetype="dashed",color="red") +
  geom_errorbarh(aes(xmin=LB,xmax=UB),height=0.25) +
  labs(title="Forest plot for APOE (19:44897490:A:T)") +
  theme(legend.position="none")
dev.off()
