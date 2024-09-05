# Name: SMIDA RNA QC Flowcell Plots Custom Module
message("Customer Script Version: 1.1")

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
library(scales)

# Module Code
dir.create("/output")
print("Loaded packages, created output folder")

### Load data ###
attrs <- study$somas$RNA$obs$attrnames()
cols <- c("Run_Tissue_name", "slide_ID_numeric", "fov", "nCell",
          attrs[c(grep("nCount", attrs), grep("nFeature", attrs))])
obs_qc <- study$somas$RNA$obs$to_dataframe(cols)
print("Loaded obs_qc")

gc()

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

#### QC metrics #####

nNeg <- length(study$somas$negprobes$var$ids())
nFC <- length(study$somas$falsecode$var$ids())

qc_stats = obs_qc %>% group_by(Run_Tissue_name, fov) %>%
  summarise(Number_Cells_Per_FOV = mean(nCell),
            Mean_Transcripts_Per_Cell_Per_FOV=mean(nCount_RNA),
            Mean_Unique_Transcripts_Per_Cell_Per_FOV=mean(nFeature_RNA),
            Total_Transcripts_Per_FOV=sum(nCount_RNA),
            Mean_Negative_Probe_Per_Plex_Per_Cell_Per_FOV=sum(nCount_negprobes)/(unique(nCell)*nNeg),
            Mean_SystemControl_Per_Plex_Per_Cell_Per_FOV=sum(nCount_falsecode)/(unique(nCell)*nFC)
  )
print("Calculated per FOV metrics")

# Save data tables to csv format
write.csv(qc_stats, "/output/Per_FOV_data_quality_metrics.csv", quote=F)
print("Finished writing output file")


########### Plotting ###############
options(scipen=999)

# Mean transcripts per FOV
p1 <- ggplot(qc_stats, aes(x=Run_Tissue_name, y=Mean_Transcripts_Per_Cell_Per_FOV, fill=as.factor(Run_Tissue_name))) +
  geom_violin() + geom_jitter(width=0.25, height=0, size = 0.5) +
  geom_boxplot(width=0.2, outlier.shape = NA) +
  scale_fill_manual(values=alpha(pal,0.3)) +
  ylab("Mean transcripts per Cell per FOV") + xlab("Flowcell") +
  theme_bw() +
  theme(legend.position="none", axis.text.x = element_text(angle=90, vjust=0.5))

ggsave("/output/Mean_Transcripts_Per_FOV.png", plot = p1, device = "png")

# Total transcripts per FOV
p2 <- ggplot(qc_stats, aes(x=Run_Tissue_name, y=Total_Transcripts_Per_FOV, fill=as.factor(Run_Tissue_name))) +
  geom_violin() + geom_jitter(width=0.25, height=0, size = 0.5) + 
  geom_boxplot(width=0.2, outlier.shape = NA) +
  scale_fill_manual(values=alpha(pal,0.3)) +
  ylab("Total transcripts per FOV") + xlab("Flowcell") +
  theme_bw() +
  scale_y_continuous(labels = scales::comma) +
  theme(legend.position="none", axis.text.x = element_text(angle=90, vjust=0.5))

ggsave("/output/Total_Transcripts_Per_FOV.png", plot = p2, device = "png")


# Number of cells per FOV
p3 <- ggplot(qc_stats, aes(x=Run_Tissue_name, y=Number_Cells_Per_FOV, fill=as.factor(Run_Tissue_name))) +
  geom_violin() + geom_jitter(width=0.25, height=0, size = 0.5) +   
  geom_boxplot(width=0.2, outlier.shape = NA) +
  scale_fill_manual(values=alpha(pal,0.3)) +
  xlab("Flowcell") + ylab("Number of cells per FOV") +
  theme_bw() +
  scale_y_continuous(labels = scales::comma) +
  theme(legend.position="none", axis.text.x = element_text(angle=90, vjust=0.5))

ggsave("/output/Number_Cells_per_FOV.png", plot = p3, device = "png")


# Mean negative probes per plex per cell per FOV 
p4 <- ggplot(qc_stats, aes(x=Run_Tissue_name, y=Mean_Negative_Probe_Per_Plex_Per_Cell_Per_FOV, fill=as.factor(Run_Tissue_name))) +
  geom_violin() + geom_jitter(width=0.25, height=0, size = 0.5) +  
  geom_boxplot(width=0.2, outlier.shape = NA) +
  scale_fill_manual(values=alpha(pal,0.3)) +
  xlab("Flowcell") + ylab("Mean NegProbe Count Per Plex Per Cell Per FOV") +
  theme_bw() +
  theme(legend.position="none", axis.text.x = element_text(angle=90, vjust=0.5))

ggsave("/output/Mean_NegProbeCount_Per_Plex_per_Cell_per_FOV.png", plot = p4, device = "png")


# Mean SystemControl per plex per cell per FOV
p5 <- ggplot(qc_stats, aes(x=Run_Tissue_name, y=Mean_SystemControl_Per_Plex_Per_Cell_Per_FOV, fill=as.factor(Run_Tissue_name))) +
  geom_violin() + geom_jitter(width=0.25, height=0, size = 0.5) +  
  geom_boxplot(width=0.2, outlier.shape = NA) +
  scale_fill_manual(values=alpha(pal,0.3)) +
  xlab("Flowcell") + ylab("Mean SystemControl Count Per Plex Per Cell Per FOV") +
  theme_bw() +
  theme(legend.position="none", axis.text.x = element_text(angle=90, vjust=0.5))

ggsave("/output/Mean_SystemControlCount_Per_Plex_per_Cell_per_FOV.png", plot = p5, device = "png")

# Mean unique transcripts per cell per FOV
p6 <- ggplot(qc_stats, aes(x=Run_Tissue_name, y=Mean_Unique_Transcripts_Per_Cell_Per_FOV, fill=as.factor(Run_Tissue_name))) +
  geom_violin() + geom_jitter(width=0.25, height=0, size = 0.5) +
  geom_boxplot(width=0.2, outlier.shape = NA) +
  scale_fill_manual(values=alpha(pal,0.3)) +
  ylab("Mean Unique transcripts per Cell per FOV") + xlab("Flowcell") +
  theme_bw() +
  theme(legend.position="none", axis.text.x = element_text(angle=90, vjust=0.5))

ggsave("/output/Mean_Unique_Transcripts_per_Cell_per_FOV.png", plot = p6, device = "png")

print("All plots complete")

