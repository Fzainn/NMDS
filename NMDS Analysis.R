##########NMDS analysis for abundance from taxa data & Site from Metadata(environmental data)############


#Reading taxa file
taxa <- read.delim("C:/Users/compu/Documents/shenzh_taxa_otuB.txt",sep='\t')


#to remove null value from df
taxa <- na.omit(taxa)


#split df by columns 
Ranks <- taxa[,1:8]

#to save df in txt file
write.table(Ranks, file = "Ranks.txt",sep="")

#delete all ranks from df
taxa$X=NULL
abundance = subset(taxa, select = -c(Rank1:Rank7) )
#abundance <- abundance[-c(1),]

#save Splited abundance in txt file
write.table(abundance, file = "abundance.txt",sep="\t")




########take sample from data for testing#######
sample_abun <- abundance[1:20,1:30]
write.table(sample_abun, file = "sample_abun.txt",sep="")

sample_abun <- read.delim("sample_abun.txt",sep='\t')

###########NMDS Analysis######
library(vegan)
library(ggplot2)
distance <- vegdist(sample_abun, method = 'bray')
nmds <- metaMDS(distance, k = 2)

stress <- nmds$stress
df <- as.data.frame(nmds$points)

#metadata(environmental data)
meta <- read.delim("metadataM.txt",sep='\t')


#we will take site column only from it
Site <- meta[1:20,2]
df_metadata <- as.data.frame(Site)

#Merge df_metadata(Site column) with df(mds1,mds2)
final_df <- cbind(df, df_metadata)

p2 <- ggplot(final_df, aes(MDS1, MDS2))+
  geom_point(aes(color = Site), size = 5)+
  geom_polygon(aes(x = MDS1, y = MDS2, fill = Site, group = Site, color = Site),
               alpha = 0.3, linetype = "longdash", linewidth = 1.5)

plot(p2)



#Analyze anosim based on bray-curtis distance
adonis <-  adonis2(sample_abun ~ Site, data = df_metadata, permutations = 999, method = "bray")

#anosim <- anosim(sample_abun, df_metadata$Site, permutations = 999, distance = "bray")


stress_text <- paste("Stress  =", round(stress, 4))
adonis_text <- paste(paste("Adonis  =", round(adonis$R2, 2)), "**")[1]
#anosim_text <- paste(paste("Anosim  =", round(anosim$statistic, 2)), "**")


p4 <- ggplot(final_df, aes(MDS1, MDS2))+
  geom_point(aes(color = Site), size = 5)+
  geom_polygon(aes(x = MDS1, y = MDS2, fill = Site, group = Site, color = Site), alpha = 0.3, linetype = "longdash", linewidth = 1.5)+
  theme(plot.margin = unit(rep(1, 4), 'lines'), 
        panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"), 
        panel.grid = element_blank(), 
        panel.background = element_rect(fill = 'white'))+
  guides(color = "none", fill = "none")+
  ggtitle(paste(paste(stress_text, adonis_text)))

plot(p4)







