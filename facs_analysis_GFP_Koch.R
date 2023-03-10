#load require libraries
library('flowCore')
library('flowWorkspace')
library('ggcyto')
library('knitr')
library('dplyr')
library('broom')

# FOR THE USER RUNNING THE CODE:set directory to get fcs files as set (CHANGE DIRECTORY ACCORDINGLY!)
dir="/Users/benjaminchang/Desktop/Collins Lab/Flow Files/2022-11-14_lib3_hits_GH41_to_GH60_72hrs_300ng/Controls"

fs <- read.flowSet(path = dir,pattern = ".fcs",alter.names = T)
as.data.frame(pData(fs)$name)

#select positive and negative samples from table above
#replace values below"
#FOR THE USER RUNNING THE CODE put in the values depending on the generated table:
pos_c=1 # positive control
neg_c=2 # negative control

print(colnames(fs))

#change column names for ease of use
colnames(fs)[colnames(fs)=="Pacific.Blue.A"] <- "BFP"
colnames(fs)[colnames(fs)=="FITC.A"] <- "GFP"
colnames(fs)[colnames(fs)=="PE.TxRed.YG.A"] <- "mCherry"
colnames(fs)[colnames(fs)=='SSC.A']='SSC-A'

#edit sample to name, keep well only
pData(fs)$name=(pData(fs) %>% tidyr::separate(name,c('sp','number','well_2c','well',extra='drop')) %>% select('well'))[,1]

#transform to gating set
gs <- GatingSet(fs)

# define gate for live cells (ADJUST PARAMETERS ACCORDINGLY)

#gs_pop_remove(gs, node="Live") ######## VARIABLE ##########
g.live <- polygonGate(filterId = "Live","FSC.A"=c(35000,227000,214000,35000),"SSC-A"=c(70000,100000,3400,900))
gs_pop_add(gs,g.live,parent="root") # add gate to GatingSet
recompute(gs) # recompute GatingSet

#check gate (this still doesn't look good, but gating is correct)
ggcyto(gs,aes(x=FSC.A,y=SSC.A),subset="root")+geom_hex(bins = 200)+geom_gate(g.live)+ggcyto_par_set(limits = "instrument") + scale_y_log10() + facet_wrap(~name,ncol = 8) + geom_stats(adjust = 0.8)


############ CHECKPOINT 1 ################

# define gate for singlets (ADJUST PARAMETERS ACCORDINGLY)

#gs_pop_remove(gs, node="Singlets")  ######## VARIABLE ##########
g.singlets <- polygonGate(filterId = "Singlets","FSC.W"=c(15000,230000,230000,40000),"FSC.H"=c(110000,110000,5000,5000))
gs_pop_add(gs,g.singlets,parent="Live") # add gate to GatingSet
recompute(gs) # recompute GatingSet
autoplot(gs,gate = 'Singlets')

#save plot
pdf(file=paste(dir,"singlets_gate.pdf",sep=""),width = 10,height = 10)

autoplot(gs,gate = 'Singlets')
dev.off()

#subset flowset to gated population
fs_gated=getData(gs,'Singlets')

print(nrow(fs_gated))


############ CHECKPOINT 2 ################

# define constitutive reporter gate (BFP)
#get gated data for positive and negative controls
df_pos=as.data.frame(exprs(fs_gated[[pos_c]]))
df_neg=as.data.frame(exprs(fs_gated[[neg_c]]))

print(df_pos)

# coarsely fit a non-linear model for range evaluation
model <- loess(mCherry ~ GFP, data = df_pos, span = 0.75)
df_pos <- augment(model, df_pos)
model <- loess(mCherry ~ GFP, data = df_neg, span = 0.75)
df_neg <- augment(model, df_neg)

#change values of xlim to evaluated the best reporter expression range
# FOR THE USER RUNNING THE CODE, the numbers in front of xlim could be set depending on the x-axis range you want to see
ggplot(df_neg,aes(GFP,.fitted)) + geom_line(color='red') + geom_line(data=df_pos,aes(x=GFP,.fitted),color='blue') + xlim(0,100000)

#save plot
pdf(file=paste(dir,"model_test_to_set_reporter_range.pdf",sep=""),width = 10,height = 10)
ggplot(df_neg,aes(GFP,.fitted)) + geom_line(color='red') + geom_line(data=df_pos,aes(x=GFP,.fitted),color='blue') + xlim(0,100000)
dev.off()

#set the reporter range for the gate
#you should manually put the range in Shiva!
# FOR THE USER RUNNING THE CODE you MUST put in the range that defines the bin for the data that will be analyzed! IMPORTANT!
#The (250, 80k) is for 300 ng per well of 48-well.
g.gfp <- rectangleGate(filterId="GFP_pos",GFP=c(250,20000))
# check gate
gs_pop_add(gs,g.gfp,parent="Singlets") # add gate to GatingSet
recompute(gs) # recalculate Gatingset
ggcyto(gs,aes(x=GFP),subset="Singlets",)+geom_density(fill="forestgreen")+geom_gate("GFP_pos")+ geom_stats(adjust = 0.1,y=0.002,digits = 1)+ggcyto_par_set(limits = "instrument")+scale_x_flowjo_biexp()+facet_wrap(~name,ncol = 8)


#subset flowset to gated population
fs_gated=getData(gs,'GFP_pos')

gs#repeat model fit for gate selection confirmation
df_pos=as.data.frame(exprs(fs_gated[[pos_c]]))
df_neg=as.data.frame(exprs(fs_gated[[neg_c]]))

# coarsely fit a non-linear model for range evaluation
model <- loess(mCherry ~ GFP, data = df_pos, span = 0.75)
df_pos <- augment(model, df_pos)
model <- loess(mCherry ~ GFP, data = df_neg, span = 0.75)
df_neg <- augment(model, df_neg)
ggplot(df_neg,aes(GFP,.fitted)) + geom_line(color='red') + geom_line(data=df_pos,aes(x=GFP,.fitted),color='blue')


#get gating stats
ps <- gs_pop_get_count_fast(gs)
ps <- ps %>% mutate(percent_of_parent=Count/ParentCount)
ps$Name= (ps %>% tidyr::separate(name,c('sp','number','well_2c','well',extra='drop')) %>% select('well'))[,1]
ps %>% select(Name,Population,Count,ParentCount,percent_of_parent) %>% kable

#save plot
pdf(file=paste(dir,"reporter_gate.pdf",sep=""),width = 10,height = 10)
ggcyto(gs,aes(x=GFP),subset="Singlets",)+geom_density(fill="blue")+geom_gate("GFP_pos")+ geom_stats(adjust = 0.1,y=0.002,digits = 1)+ggcyto_par_set(limits = "instrument")+scale_x_flowjo_biexp()+facet_wrap(~name,ncol = 8)
dev.off()

#subset flowset to gated population
fs_gated=getData(gs,'GFP_pos')

#get data for final gated population, remove negative mCherry values. Median (not mean) is computed and plotted for all samples
results=data.frame()

######### CHANGE FOR DIFFERENT FLUOROPHORES #########

FITC_col <- which(colnames(fs_gated) == "FITC.A")
mCherry_col <- which(colnames(fs_gated) == "mCherry")
GFP_col <- which(colnames(fs_gated) == "GFP")

for (i in 1:length(rownames(pData(fs_gated)))){
  df=as.data.frame(exprs(fs_gated[[i]]))
  df=df[which(df$mCherry>0),]
  
  if (dim(df)[1] == 0){
    next
  } else {
    print(df)
    df=df[,c(GFP_col,mCherry_col)]
    print(df)
    df$sample=rownames(pData(fs_gated[i]))
    results=rbind(results,df)
  }}
results$sample=(results %>% tidyr::separate(sample,c('sp','number','well_2c','well',extra='drop')) %>% select('well'))[,1]

p1=ggplot(results,aes(x=sample,y=mCherry/GFP)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
p1

#save the plot
pdf(file=paste(dir,"single_sample_induction_boxplot.pdf",sep=""),width = 10,height = 10)
p1
dev.off()

data_stats=as.data.frame(ggplot_build(p1)$data[[1]][,1:5])
rownames(data_stats)=sort(unique(p1$data$sample))
data_stats

#write csv (open in excel) file with median stats
write.csv(data_stats,file=paste(dir,"median_stats.csv",sep=""))