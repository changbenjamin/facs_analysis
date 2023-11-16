#load require libraries
library('flowCore')
library('flowWorkspace')
library('ggcyto')
library('knitr')
library('dplyr')
library('broom')

rm(list=ls())

# FOR THE USER RUNNING THE CODE:set directory to get fcs files as set (CHANGE DIRECTORY ACCORDINGLY!)
dir="/Users/joncchen/Dropbox (MIT)/Collins Lab: RNA-ligand screen/Raw Flow Files/2023-11-09-triple_combinatorial/Exp_20231109_3/"

fs <- read.flowSet(path = dir,pattern = ".fcs",alter.names = T) #,truncate_max_range = FALSE)
as.data.frame(pData(fs)$name)

#select positive and negative samples from table above
#replace values below"
#FOR THE USER RUNNING THE CODE put in the values depending on the generated table:
pos_c=14
neg_c=18

##change column names for ease of use
colnames(fs)[colnames(fs)=="FL4.A"] <- "BFP"
colnames(fs)[colnames(fs)=="FL13.A"] <- "mCherry"

#edit sample to name, keep well only
names1=(pData(fs) %>% tidyr::separate(name,c('plate','bla','well',extra='drop'))) %>% select('well')
pData(fs)$name=names1[,1]

gs <- GatingSet(fs)

coor1 <- c(25e4, 0)
coor2 <- c(22e4, 3e4)
coor3 <- c(19e4, 10e4)
coor4 <- c(26e4, 23e4)
coor5 <- c(62.5e4, 53e4)
coor6 <- c(86e4, 42e4)
coor7 <- c(91e4, 20e4)
coor8 <- c(41e4, 3e4)

# define gate for live cells (ADJUST PARAMETERS ACCORDINGLY)
# gs_pop_remove(gs, node='Live')
g.live <- polygonGate(filterId = "Live","FSC.A"=c(coor1[1],coor2[1],coor3[1],coor4[1],coor5[1],coor6[1],coor7[1],coor8[1]),"SSC.A"=c(coor1[2],coor2[2],coor3[2],coor4[2],coor5[2],coor6[2],coor7[2],coor8[2]))
gs_pop_add(gs,g.live,parent="root") # add gate to GatingSet
recompute(gs) # recompute GatingSet

#check gate - lower bins = less processing power required
out_live <- ggcyto(gs,aes(x=FSC.A,y=SSC.A),subset="root")+geom_hex(bins = 2048)+geom_gate(g.live)+ggcyto_par_set(limits = list(x=c(-10,1e6),y=c(-10,1e6))) + facet_wrap(~name,ncol = 8) + geom_stats(adjust = 0.8)
out_live
# ggcyto(gs,aes(x=FSC.A,y=SSC.A),subset="root")+geom_hex(bins = 200)+geom_gate(g.live)+ggcyto_par_set(limits = list(x=c(1,5e6),y=c(-10,5e6))) + facet_wrap(~name,ncol = 8) + geom_stats(adjust = 0.8)

#save plot
pdf(file=paste(dir,"live_gate.pdf",sep=""),width = 40,height = 40)
out_live
dev.off()


############# LIVE CHECKPOINT #############
# NaOH
coor1 <- c(1.5e5, 1425)
coor2 <- c(1.2e5, 1575)
coor3 <- c(1.55e5, 2075)
coor4 <- c(2.1e5, 2575)
coor5 <- c(3.15e5, 2575)
coor6 <- c(3.55e5, 2125)
coor7 <- c(3.5e5, 1675)
coor8 <- c(2.1e5, 1425)

# # DMSO
# coor1 <- c(1.25e5, 1475)
# coor2 <- c(0.95e5, 1700)
# coor3 <- c(1.25e5, 2200)
# coor4 <- c(1.85e5, 2500)
# coor5 <- c(2.90e5, 2500)
# coor6 <- c(3.3e5, 2250)
# coor7 <- c(3.25e5, 1800)
# coor8 <- c(2.1e5, 1475)


# define gate for singlets (ADJUST PARAMETERS ACCORDINGLY)
# gs_pop_remove(gs, node="Singlets")  ######## VARIABLE ##########
g.singlets <- polygonGate(filterId = "Singlets","FSC.H"=c(coor1[1],coor2[1],coor3[1],coor4[1],coor5[1],coor6[1],coor7[1],coor8[1]),"FSC.Width"=c(coor1[2],coor2[2],coor3[2],coor4[2],coor5[2],coor6[2],coor7[2],coor8[2]))
gs_pop_add(gs,g.singlets,parent="Live") # add gate to GatingSet
recompute(gs) # recompute GatingSet
singlet_out <- ggcyto(gs,aes(x=FSC.H,y=FSC.Width),subset="Live")+geom_hex(bins = 256)+geom_gate(g.singlets)+ggcyto_par_set(limits = list(x=c(0,400e3),y=c(0,400e1)))+ facet_wrap(~name,ncol = 8) + geom_stats(adjust = 0.8)

#save plot
pdf(file=paste(dir,"singlets_gate.pdf",sep=""),width = 40,height = 40)
singlet_out
dev.off()

# coord1 <- c(190000, 1.05e5)
# coord2 <- c(70000, 1.70e5)
# coord3 <- c(400000, 3e5)
# coord4 <- c(687500, 2.85e5)

# shifted right and down
coord1 <- c(190000, 0.75e5)
coord2 <- c(70000, 1.50e5)
coord3 <- c(460000, 3e5)
coord4 <- c(747500, 2.85e5)


# define gate for singlets2 (ADJUST PARAMETERS ACCORDINGLY)
# gs_pop_remove(gs, node="Singlets2")  ######## VARIABLE ##########
g.singlets2 <- polygonGate(filterId = "Singlets2","FSC.A"=c(coord1[1],coord2[1],coord3[1],coord4[1]),"FSC.H"=c(coord1[2],coord2[2],coord3[2],coord4[2]))
gs_pop_add(gs,g.singlets2,parent="Singlets") # add gate to GatingSet
recompute(gs) # recompute GatingSet
singlet2_out <- ggcyto(gs,aes(x=FSC.A,y=FSC.H),subset="Singlets")+geom_hex(bins = 128)+geom_gate(g.singlets2)+ggcyto_par_set(limits = list(x=c(1,2e6),y=c(1e4,4e5))) + facet_wrap(~name,ncol = 8) + geom_stats(adjust = 0.8)
singlet2_out

#save plot
pdf(file=paste(dir,"singlets2_gate.pdf",sep=""),width = 40,height = 40)
singlet2_out
dev.off()

############# SINGLETS CHECKPOINT #############

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
bfp_reporter_range <- ggplot(df_neg,aes(BFP,.fitted)) + geom_line(color='red') + geom_line(data=df_pos,aes(x=BFP,.fitted),color='blue') + xlim(250,80000)
#save plot
pdf(file=paste(dir,"model_test_to_set_reporter_range.pdf",sep=""),width = 10,height = 10)
bfp_reporter_range
dev.off()

# 
# #set the reporter range for the gate
# #you should manually put the range in!
# # FOR THE USER RUNNING THE CODE you MUST put in the range that defines the bin for the data that will be analyzed! IMPORTANT!
# #The (250, 80k) is for 300 ng per well of 48-well.
# #gs_pop_remove(gs, node="BFP_pos")
# g.bfp <- rectangleGate(filterId="BFP_pos",BFP=c(500,1e6))
# # check gate
# gs_pop_add(gs,g.bfp,parent="Singlets2") # add gate to GatingSet
# recompute(gs) # recalculate Gatingset
# BFP_hist <- ggcyto(gs,aes(x=BFP),subset="Singlets2",)+geom_density(fill="skyblue")+geom_gate("BFP_pos")+ geom_stats(adjust = 0.1,digits = 1)+scale_x_flowjo_biexp()+facet_wrap(~name,ncol = 8)
# BFP_hist
# #save plot
# pdf(file=paste(dir,"reporter_gate.pdf",sep=""),width = 10,height = 10)
# BFP_hist
# dev.off()


# STANDARD
coor1 <- c(885, 3.7e4)
coor2 <- c(735, 15e4)
coor3 <- c(835, 32e4)
coor4 <- c(1200, 42e4)
coor5 <- c(1e4, 52e4)
coor6 <- c(1.8e5, 51e4)
coor7 <- c(8e5, 35e4)
coor8 <- c(5e5, 2.5e4)
coor9 <- c(2.7e3, 1e4)


# # JCC_27/51 high expression
# coor1 <- c(2685, 3.7e4)
# coor2 <- c(2535, 15e4)
# coor3 <- c(2635, 32e4)
# coor4 <- c(3000, 42e4)
# coor5 <- c(1.5e4, 52e4)
# coor6 <- c(3.85e5, 51e4)
# coor7 <- c(10.05e5, 35e4)
# coor8 <- c(7.05e5, 2.5e4)
# coor9 <- c(7.7e3, 1e4)


# gs_pop_remove(gs, node="BFP_pos")  ######## VARIABLE ##########
g.bfp <- polygonGate(filterId = "BFP_pos","BFP"=c(coor1[1],coor2[1],coor3[1],coor4[1],coor5[1],coor6[1],coor7[1],coor8[1],coor9[1]),"SSC.A"=c(coor1[2],coor2[2],coor3[2],coor4[2],coor5[2],coor6[2],coor7[2],coor8[2],coor9[2]))
gs_pop_add(gs,g.bfp,parent="Singlets2") # add gate to GatingSet
recompute(gs) # recompute GatingSet
bfp_dots <-ggcyto(gs,aes(x=BFP,y=SSC.A),subset="Singlets2")+geom_hex(bins = 30)+geom_gate(g.bfp)+ggcyto_par_set(limits = list(x=c(1,1e6),y=c(1e4,1e6))) + facet_wrap(~name,ncol = 8) + scale_x_log10() + geom_stats(adjust = 0.8)

#save plot
pdf(file=paste(dir,"BFP_dots_gate.pdf",sep=""),width = 40,height = 40)
bfp_dots
dev.off()


#subset flowset to gated population
fs_gated=gs_pop_get_data(gs,'BFP_pos')

#repeat model fit for gate selection confirmation
df_pos=as.data.frame(exprs(fs_gated[[pos_c]]))

write.csv(df_pos,file=paste(dir,"test_bfpbins.csv",sep=""))

df_neg=as.data.frame(exprs(fs_gated[[neg_c]]))

# coarsely fit a non-linear model for range evaluation
model <- loess(mCherry ~ BFP, data = df_pos, span = 0.75)
df_pos <- augment(model, df_pos)
model <- loess(mCherry ~ BFP, data = df_neg, span = 0.75)
df_neg <- augment(model, df_neg)
final_model <- ggplot(df_neg,aes(BFP,.fitted)) + geom_line(color='red') + geom_line(data=df_pos,aes(x=BFP,.fitted),color='blue')
final_model

pdf(file=paste(dir,"updated_reporter_model.pdf",sep=""),width = 10,height = 10)
final_model
dev.off()
# make sure the model looks right, before proceeding.







######################## FINAL CHECKPOINT ########################







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

p1 <- ggplot(results, aes(x = sample, y = mCherry / BFP)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  scale_y_log10()
p1
#save the plot
pdf(file=paste(dir,"single_sample_induction_boxplot.pdf",sep=""),width = 40,height = 10)
p1
dev.off()

BFP_pos <- gating_rep[gating_rep$Population == '/Live/Singlets/Singlets2/BFP_pos']$Count

Total_counts <- gating_rep[gating_rep$Population == '/Live']$Count
Singlets <- gating_rep[gating_rep$Population == '/Live/Singlets/Singlets2']$Count
Live_percent <- 100*gating_rep[gating_rep$Population == '/Live']$percent_of_parent
BFP_percent <- 100*gating_rep[gating_rep$Population == '/Live/Singlets/Singlets2/BFP_pos']$percent_of_parent


###### BIN BFP ######

# Create an empty data frame with the appropriate column names
bfpbin_df <- data.frame(BFP_20 = numeric(), 
                        BFP_40 = numeric(), 
                        BFP_60 = numeric(), 
                        BFP_80 = numeric(), 
                        BFP_100 = numeric())

# Loop over all_samples
for (i in seq_along(pData(fs)$name)) {
  # Extract the expression data for the i-th sample and create a data frame
  df1 <- as.data.frame(exprs(fs_gated[[i]]))
  
  # Calculate the mCherry/BFP ratio for each cell in df1
  df1$mCherry_BFP <- df1$mCherry / df1$BFP
  
  # Sort df1 by BFP from low to high
  df1 <- df1[order(df1$BFP),]
  
  # Divide the sorted data into quintiles (20% intervals) based on BFP
  n <- nrow(df1)
  quintile_size <- n %/% 5
  quintiles <- rep(1:5, each = quintile_size)
  if (n %% 5 != 0) {
    quintiles <- c(quintiles, rep(5, n %% 5))
  }
  df1$quintile <- quintiles
  
  # Compute the median mCherry/BFP ratio for each quintile and add it to bfpbin_df
  for (q in 1:5) {
    median_ratio <- median(df1$mCherry_BFP[df1$quintile == q])
    bfpbin_df[i, q] <- median_ratio
  }
}

bfpbin_df


# Set row names of bfpbin_df to the sample names
rownames(bfpbin_df) <- pData(fs)$name

p2=ggplot(results,aes(x=sample,y=mCherry/BFP)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
data_stats=as.data.frame(ggplot_build(p2)$data[[1]][,1:5])
rownames(data_stats)=sort(unique(p1$data$sample))
data_stats <- cbind(data_stats, BFP_pos, Total_counts, Singlets, Live_percent, BFP_percent)

data_stats <- cbind(bfpbin_df, data_stats)

data_stats

#write csv (open in excel) file with median stats
write.csv(data_stats,file=paste(dir,"BFP_median_stats.csv",sep=""))