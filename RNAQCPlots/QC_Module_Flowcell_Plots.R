# Name: SMIDA RNA QC Flowcell Plots Custom Module
message("Customer Script Version: 1.0.1")

# Copyright 2024 Bruker Spatial Biology, Inc.
# This software and any associated files are distributed pursuant to the NanoString AtoMx Spatial 
# Information Platform Software as a Service Agreement, available on the NanoString Technologies, Inc 
# website at www.nanostring.com, as updated.  All rights reserved.  No permission is granted to modify, 
# publish, distribute, sublicense or sell copies of this software.

# Description: Create RNA QC plots

# User Defined Variables: None

# Load Packages
library(ggplot2)
library(dplyr)
library(reshape2)
library(stringr)
library(scales)

# Module Code
dir.create("/output")
print("Loaded packages, created output folder")

### Load data ###
obs_qc <- study$somas$RNA$obs$to_dataframe()
print("Loaded obs_qc")
neg <- study$somas$negprobes$to_seurat_assay(layers='counts')
print("loaded negprobe data")
raw <- study$somas$RNA$X$members$counts$to_dataframe()
print("loaded raw data")
fc <-  study$somas$falsecode$to_seurat_assay(layers='counts')
print("loaded false codes")

### Adjust color palette as necessary ###
pal <- c("#3A6CA1", "#FFD861", "#D86769", "#AEE8E2", "#999999", "#FABFD2", "#318026", "#A66293", "#F28E2B", "#E3C0AC",
         "#A0CBE8", "#9E7E20", "#FFBFBD", "#2CABA3", "#474747", "#EB739E", "#A0E391", "#E8C1DE", "#FCB36A", "#B0846B")

nfov <- length(unique(obs_qc$fov))
nslide <- length(unique(obs_qc$slide_ID_numeric))
total <- nfov * nslide
if(total > length(pal)) {
  pal <- rep(pal, ceiling(total / length(pal)))
}

print("set color palette")
length(pal)



###### NegProbe Counts ########
neg <- neg@counts
print("neg <- neg@counts")
neg_counts <- as.data.frame(t(as.matrix(neg)))
print("transposed and converted negprobe data to dataframe")
neg_counts$Cell_ID <- rownames(neg_counts)
neg_long <- melt(neg_counts,
                 id.vars=c('Cell_ID'),
                 variable.name="NegProbe",
                 value.name="Count"
)
print("transformed negprobe data to long")

## Modify NegProbe data
cols <- c("Run_Tissue_name","cell_id")
neg_long[c("c","slide_ID_numeric","fov","cell")] <-  str_split_fixed(neg_long$Cell_ID, "_",4)
neg_long$slide_fov <- paste0(neg_long$slide_ID_numeric, "_", neg_long$fov)
obs_flowcell <- obs_qc[,cols]
neg_long_fc <- base::merge(neg_long, obs_flowcell, by.x="Cell_ID", by.y="cell_id", all.x=T )
neg_long <- neg_long_fc
print("created additional columns in long negprobe data")

# Calculate the mean of the negprobes per plex per cell:
neg_probe_avg <- neg_long %>% group_by(Run_Tissue_name, fov) %>%
  summarise(Mean_Negative_Probe_Per_Plex_Per_Cell_Per_FOV = mean(Count))
print("Calcualted mean NegProbe per plex per cell per FOV")

###### SystemControl Counts ########
fc2 <- fc@counts
fc_counts <- as.data.frame(t(as.matrix(fc2)))
fc_counts$Cell_ID <- rownames(fc_counts)
fc_long <- melt(fc_counts,
                id.vars=c('Cell_ID'),
                variable.name="SystemControl",
                value.name="Count"
)

## Modify FalseCode data
cols <- c("Run_Tissue_name","cell_id")
fc_long[c("c","slide_ID_numeric","fov","cell")] <-  str_split_fixed(fc_long$Cell_ID, "_",4)
fc_long$slide_fov <- paste0(fc_long$slide_ID_numeric, "_", fc_long$fov)
obs_flowcell <- obs_qc[,cols]
fc_long_fc <- base::merge(fc_long, obs_flowcell, by.x="Cell_ID", by.y="cell_id", all.x=T )
fc_long <- fc_long_fc
print("created additional columns in long SystemControl data")

# Calculate the mean of the falseCodes per plex per cell:
fc_avg <- fc_long %>% group_by(Run_Tissue_name, fov) %>%
  summarise(Mean_SystemControl_Per_Plex_Per_Cell_Per_FOV = mean(Count))
print("Calculated mean SystemControl count per plex per cell per FOV")


#### Other QC metrics #####
qc_stats = obs_qc %>% group_by(Run_Tissue_name, fov) %>%
  summarise(Number_Cells_Per_FOV = mean(nCell),
            Mean_Transcripts_Per_Cell_Per_FOV=mean(nCount_RNA),
            Mean_Unique_Transcripts_Per_Cell_Per_FOV=mean(nFeature_RNA),
            Total_Transcripts_Per_FOV=sum(nCount_RNA)
  )
print("Calculated per FOV metrics")

allqc <- base::merge(qc_stats, neg_probe_avg, by=c('Run_Tissue_name', 'fov'), all=T)
print("Merged qc and negprobe data")
allqc_final <- base::merge(allqc, fc_avg, by=c('Run_Tissue_name', 'fov'), all=T)
print("Merged allqc and SystemControl data")



# Save data tables to csv format
write.csv(allqc_final, "/output/Per_FOV_data_quality_metrics.csv", quote=F)
print("Finished writing output file")



########### Plotting ###############
options(scipen=999)

# Mean transcripts per FOV
p1 <- ggplot(allqc_final, aes(x=Run_Tissue_name, y=Mean_Transcripts_Per_Cell_Per_FOV, fill=as.factor(Run_Tissue_name))) +
  geom_violin() + geom_jitter(width=0.25, height=0, size = 0.5) +
  geom_boxplot(width=0.2, outlier.shape = NA) +
  scale_fill_manual(values=alpha(pal,0.3)) +
  ylab("Mean transcripts per Cell per FOV") + xlab("Flowcell") +
  theme_bw() +
  theme(legend.position="none", axis.text.x = element_text(angle=90, vjust=0.5))

ggsave("/output/Mean_Transcripts_Per_FOV.png", plot = p1, device = "png")

# Total transcripts per FOV
p2 <- ggplot(allqc_final, aes(x=Run_Tissue_name, y=Total_Transcripts_Per_FOV, fill=as.factor(Run_Tissue_name))) +
  geom_violin() + geom_jitter(width=0.25, height=0, size = 0.5) + 
  geom_boxplot(width=0.2, outlier.shape = NA) +
  scale_fill_manual(values=alpha(pal,0.3)) +
  ylab("Total transcripts per FOV") + xlab("Flowcell") +
  theme_bw() +
  scale_y_continuous(labels = scales::comma) +
  theme(legend.position="none", axis.text.x = element_text(angle=90, vjust=0.5))

ggsave("/output/Total_Transcripts_Per_FOV.png", plot = p2, device = "png")


# Number of cells per FOV
p3 <- ggplot(allqc_final, aes(x=Run_Tissue_name, y=Number_Cells_Per_FOV, fill=as.factor(Run_Tissue_name))) +
  geom_violin() + geom_jitter(width=0.25, height=0, size = 0.5) +   
  geom_boxplot(width=0.2, outlier.shape = NA) +
  scale_fill_manual(values=alpha(pal,0.3)) +
  xlab("Flowcell") + ylab("Number of cells per FOV") +
  theme_bw() +
  scale_y_continuous(labels = scales::comma) +
  theme(legend.position="none", axis.text.x = element_text(angle=90, vjust=0.5))

ggsave("/output/Number_Cells_per_FOV.png", plot = p3, device = "png")


# Mean negative probes per plex per cell per FOV 
p4 <- ggplot(allqc_final, aes(x=Run_Tissue_name, y=Mean_Negative_Probe_Per_Plex_Per_Cell_Per_FOV, fill=as.factor(Run_Tissue_name))) +
  geom_violin() + geom_jitter(width=0.25, height=0, size = 0.5) +  
  geom_boxplot(width=0.2, outlier.shape = NA) +
  scale_fill_manual(values=alpha(pal,0.3)) +
  xlab("Flowcell") + ylab("Mean NegProbe Count Per Plex Per Cell Per FOV") +
  theme_bw() +
  theme(legend.position="none", axis.text.x = element_text(angle=90, vjust=0.5))

ggsave("/output/Mean_NegProbeCount_Per_Plex_per_Cell_per_FOV.png", plot = p4, device = "png")


# Mean SystemControl per plex per cell per FOV
p5 <- ggplot(allqc_final, aes(x=Run_Tissue_name, y=Mean_SystemControl_Per_Plex_Per_Cell_Per_FOV, fill=as.factor(Run_Tissue_name))) +
  geom_violin() + geom_jitter(width=0.25, height=0, size = 0.5) +  
  geom_boxplot(width=0.2, outlier.shape = NA) +
  scale_fill_manual(values=alpha(pal,0.3)) +
  xlab("Flowcell") + ylab("Mean SystemControl Count Per Plex Per Cell Per FOV") +
  theme_bw() +
  theme(legend.position="none", axis.text.x = element_text(angle=90, vjust=0.5))

ggsave("/output/Mean_SystemControlCount_Per_Plex_per_Cell_per_FOV.png", plot = p5, device = "png")

# Mean unique transcripts per cell per FOV
p6 <- ggplot(allqc_final, aes(x=Run_Tissue_name, y=Mean_Unique_Transcripts_Per_Cell_Per_FOV, fill=as.factor(Run_Tissue_name))) +
  geom_violin() + geom_jitter(width=0.25, height=0, size = 0.5) +
  geom_boxplot(width=0.2, outlier.shape = NA) +
  scale_fill_manual(values=alpha(pal,0.3)) +
  ylab("Mean Unique transcripts per Cell per FOV") + xlab("Flowcell") +
  theme_bw() +
  theme(legend.position="none", axis.text.x = element_text(angle=90, vjust=0.5))

ggsave("/output/Mean_Unique_Transcripts_per_Cell_per_FOV.png", plot = p6, device = "png")

print("All plots complete")

