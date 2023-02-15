#load require libraries
library('flowCore')
library('flowWorkspace')
library('ggcyto')
library('knitr')
library('dplyr')
library('broom')

# FOR THE USER RUNNING THE CODE:set directory to get fcs files as set (CHANGE DIRECTORY ACCORDINGLY!)
dir="/Users/benjaminchang/Desktop/Collins Lab/GH162-211 Flow/2023-02-13_initial_test_at_Broad/Samples"

fs <- read.flowSet(path = dir,pattern = ".fcs",alter.names = T) #,truncate_max_range = FALSE)
as.data.frame(pData(fs)$name)

#select positive and negative samples from table above
#replace values below"
#FOR THE USER RUNNING THE CODE put in the values depending on the generated table:
pos_c=1
neg_c=2

##change column names for ease of use
colnames(fs)[colnames(fs)=="FL4.A"] <- "BFP"
colnames(fs)[colnames(fs)=="FL13.A"] <- "mCherry"

#edit sample to name, keep well only
names1=(pData(fs) %>% tidyr::separate(name,c('plate','bla','well',extra='drop'))) %>% select('well')
pData(fs)$name=names1[,1]

gs <- GatingSet(fs)

# define gate for live cells (ADJUST PARAMETERS ACCORDINGLY)
g.live <- polygonGate(filterId = "Live","FSC.A"=c(155000,450000,205000,70000),"SSC.A"=c(1e6,5e6,45000,25000))
gs_pop_add(gs,g.live,parent="root") # add gate to GatingSet
recompute(gs) # recompute GatingSet

#check gate
ggcyto(gs,aes(x=FSC.A,y=SSC.A),subset="root")+geom_hex(bins = 200)+geom_gate(g.live)+ggcyto_par_set(limits = list(x=c(0,5e5),y=c(-10,7.5e7))) + scale_y_log10() + facet_wrap(~name,ncol = 8) + geom_stats(adjust = 0.8)

#save plot
pdf(file=paste(dir,"live_gate.pdf",sep=""),width = 10,height = 10)
ggcyto(gs,aes(x=FSC.A,y=SSC.A),subset="root")+geom_hex(bins = 200)+geom_gate(g.live)+ggcyto_par_set(limits = list(x=c(0,5e5),y=c(-10,7.5e7))) + scale_y_log10() + facet_wrap(~name,ncol = 8) + geom_stats(adjust = 0.8)
dev.off()

# define gate for singlets (ADJUST PARAMETERS ACCORDINGLY)
#gs_pop_remove(gs, node="Singlets")  ######## VARIABLE ##########
g.singlets <- polygonGate(filterId = "Singlets","FSC.Width"=c(1000,3e3,3e3,1000),"FSC.H"=c(1.5e5,1.5e5,3e4,3e4))
gs_pop_add(gs,g.singlets,parent="Live") # add gate to GatingSet
recompute(gs) # recompute GatingSet
ggcyto(gs,aes(x=FSC.H,y=FSC.Width),subset="Live")+geom_hex(bins = 200)+geom_gate(g.singlets)+ggcyto_par_set(limits = list(x=c(3e4,1.5e5),y=c(250,1e4))) + scale_x_log10() + facet_wrap(~name,ncol = 8) + geom_stats(adjust = 0.8)

#save plot
pdf(file=paste(dir,"singlets_gate.pdf",sep=""),width = 10,height = 10)
autoplot(gs,gate = 'Singlets')
dev.off()

# define gate for singlets2 (ADJUST PARAMETERS ACCORDINGLY)
#gs_pop_remove(gs, node="Singlets2")  ######## VARIABLE ##########
g.singlets2 <- polygonGate(filterId = "Singlets2","FSC.A"=c(5e4,2.5e5,2.5e5,5e4),"FSC.H"=c(3e5,3e5,3e4,3e4))
gs_pop_add(gs,g.singlets2,parent="Singlets") # add gate to GatingSet
recompute(gs) # recompute GatingSet
ggcyto(gs,aes(x=FSC.A,y=FSC.H),subset="Singlets")+geom_hex(bins = 200)+geom_gate(g.singlets2)+ggcyto_par_set(limits = list(x=c(1,5e5),y=c(1e4,1e6))) + scale_y_log10() + facet_wrap(~name,ncol = 8) + geom_stats(adjust = 0.8)

#save plot
pdf(file=paste(dir,"singlets2_gate.pdf",sep=""),width = 10,height = 10)
autoplot(gs,gate = 'Singlets2')
dev.off()

#subset flowset to gated population
fs_gated=gs_pop_get_data(gs,'Singlets2')

# define constitutive reporter gate (BFP)
#get gated data for postive and negative controls
df_pos=as.data.frame(exprs(fs_gated[[pos_c]]))
df_neg=as.data.frame(exprs(fs_gated[[neg_c]]))

# coarsely fit a non-linear model for range evaluation
model <- loess(mCherry ~ BFP, data = df_pos, span = 0.75)
df_pos <- augment(model, df_pos)
model <- loess(mCherry ~ BFP, data = df_neg, span = 0.75)
df_neg <- augment(model, df_neg)

#change values of xlim to evaluated the best reporter expression range
# FOR THE USER RUNNING THE CODE, the numbers in front of xlim could be set depending on the x-axis range you want to see
ggplot(df_neg,aes(BFP,.fitted)) + geom_line(color='red') + geom_line(data=df_pos,aes(x=BFP,.fitted),color='blue') + xlim(250,80000)
#save plot
pdf(file=paste(dir,"model_test_to_set_reporter_range.pdf",sep=""),width = 10,height = 10)
ggplot(df_neg,aes(BFP,.fitted)) + geom_line(color='red') + geom_line(data=df_pos,aes(x=BFP,.fitted),color='blue') + xlim(250,80000)
dev.off()

#set the reporter range for the gate
#you should manually put the range in Shiva!
# FOR THE USER RUNNING THE CODE you MUST put in the range that defines the bin for the data that will be analyzed! IMPORTANT!
#The (250, 80k) is for 300 ng per well of 48-well.
g.bfp <- rectangleGate(filterId="BFP_pos",BFP=c(2000,60000))
# check gate
gs_pop_add(gs,g.bfp,parent="Live") # add gate to GatingSet
recompute(gs) # recalculate Gatingset
ggcyto(gs,aes(x=BFP),subset="Live",)+geom_density(fill="forestgreen")+geom_gate("BFP_pos")+ geom_stats(adjust = 0.1,y=0.002,digits = 1)+ggcyto_par_set(limits = "instrument")+scale_x_flowjo_biexp()+facet_wrap(~name,ncol = 8)
#save plot
pdf(file=paste(dir,"reporter_gate.pdf",sep=""),width = 10,height = 10)
ggcyto(gs,aes(x=BFP),subset="Live",)+geom_density(fill="forestgreen")+geom_gate("BFP_pos")+ geom_stats(adjust = 0.1,y=0.002,digits = 1)+ggcyto_par_set(limits = "instrument")+scale_x_flowjo_biexp()+facet_wrap(~name,ncol = 8)
dev.off()

#subset flowset to gated population
fs_gated=gs_pop_get_data(gs,'BFP_pos')

#repeat model fit for gate selection confirmation
df_pos=as.data.frame(exprs(fs_gated[[pos_c]]))
df_neg=as.data.frame(exprs(fs_gated[[neg_c]]))

# coarsely fit a non-linear model for range evaluation
model <- loess(mCherry ~ BFP, data = df_pos, span = 0.75)
df_pos <- augment(model, df_pos)
model <- loess(mCherry ~ BFP, data = df_neg, span = 0.75)
df_neg <- augment(model, df_neg)
ggplot(df_neg,aes(BFP,.fitted)) + geom_line(color='red') + geom_line(data=df_pos,aes(x=BFP,.fitted),color='blue')


#get gating stats
ps <- gs_pop_get_count_fast(gs)
ps <- ps %>% mutate(percent_of_parent=Count/ParentCount)
ps$Name= (ps %>% tidyr::separate(name,c('plate','bla','well',extra='drop')) %>% select('well'))[,1]

gating_rep=ps %>% select(Name,Population,Count,ParentCount,percent_of_parent)
gating_rep
write.csv(gating_rep,file=paste(dir,"gating_stats.csv",sep=""))

#subset flowset to gated population
fs_gated=gs_pop_get_data(gs,'BFP_pos')

mCherry_col <- which(colnames(fs_gated) == "mCherry")
BFP_col <- which(colnames(fs_gated) == "BFP")

#get data for final gated population, remove negative mCHerry values. Median (not mean) is computed and plotted for all samples
results=data.frame()
for (i in 1:length(rownames(pData(fs_gated)))){
  df=as.data.frame(exprs(fs_gated[[i]]))
  df=df[which(df$mCherry>0),]
  if (dim(df)[1] == 0){
    next
  } else {
    df=df[,c(mCherry_col, BFP_col)]
    df$sample=rownames(pData(fs_gated[i]))
    results=rbind(results,df)
  }}
results$sample=(results %>% tidyr::separate(sample,c('plate','bla','well',extra='drop')) %>% select('well'))[,1]

p1=ggplot(results,aes(x=sample,y=mCherry/BFP)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
p1
#save the plot
pdf(file=paste(dir,"single_sample_induction_boxplot.pdf",sep=""),width = 10,height = 10)
p1
dev.off()

p2=ggplot(results,aes(x=sample,y=mCherry/BFP)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
data_stats=as.data.frame(ggplot_build(p2)$data[[1]][,1:5])
rownames(data_stats)=unique(p1$data$sample)
data_stats

#write csv (open in excel) file with median stats
write.csv(data_stats,file=paste(dir,"median_stats.csv",sep=""))
