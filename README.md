NanoString Collection of CosMx SMIDA custom module scripts
=================
This repository contains R script that can be used in custom modules in the AtoMx SMIDA.
 

| Name              | Description                                                                                                                                                                                                                                                                                                                                                                                                                                         |Links|
| :---------------- | :-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |:------------------------|
| Export     | Export all files associated with current study: Seurat, TileDB array, and/or decoded files.| [R script](https://github.com/Nanostring-Biostats/CosMxDACustomModules/blob/main/Export/CosMxDAExport.R), [Documentation](https://github.com/Nanostring-Biostats/CosMxDACustomModules/blob/main/Export/CosMxDAExportSetup.docx) | 
| Flat File Export     | Export Flat files used in [Seurat's LoadNanostring function](https://satijalab.org/seurat/articles/spatial_vignette_2#human-lung-nanostring-cosmx-spatial-molecular-imager)| [R script](https://github.com/Nanostring-Biostats/CosMxDACustomModules/blob/main/flatFileExport/flatFileExport.R), [Documentation](https://github.com/Nanostring-Biostats/CosMxDACustomModules/blob/main/flatFileExport/CosMxDAFlatFileExportSetup.docx) |  
| RNA QC Plots     | Summary plots of CosMx SMI RNA QC | [R script](https://github.com/Nanostring-Biostats/CosMxDACustomModules/blob/main/RNAQCPlots/QC_Module_Flowcell_Plots.R), [Documentation](https://github.com/Nanostring-Biostats/CosMxDACustomModules/blob/main/RNAQCPlots/CosMxDA_RNAQCPlotSetup.docx) | 
| Get Sample Metadata    | Get the Sample Metadata for a Study | [R script](https://github.com/Nanostring-Biostats/CosMxDACustomModules/blob/main/SampleMetadata/GetSampleMetadata.R), [Documentation](https://github.com/Nanostring-Biostats/CosMxDACustomModules/blob/main/SampleMetadata/CosMxDAGetSampleMetadataSetup.docx) |
| Update Sample Metadata    | Update the Sample Metadata for a Study | [R script](https://github.com/Nanostring-Biostats/CosMxDACustomModules/blob/main/SampleMetadata/UpdateSampleMetadata.R), [Documentation](https://github.com/Nanostring-Biostats/CosMxDACustomModules/blob/main/SampleMetadata/CosMxDAUpdateSampleMetadataSetup.docx) | 
| Combine Sample Metadata    | Generate Query Columns in Sample Metadata | [Rscript](https://github.com/Nanostring-Biostats/CosMxDACustomModules/blob/main/SampleMetadata/CombineSampleMetadata.R), [Documentation](https://github.com/Nanostring-Biostats/CosMxDACustomModules/blob/main/SampleMetadata/CosMxDACombineSampleMetadataSetup.docx) |
  
  
Export custom modules have an AtoMx v1.3.2 and v1.4.0 version. These modules are NOT interchangeable and must be run on their respective AtoMx versions. Each export module has an AtoMx_v1.3.2 folder that contains the older scripts. 

## Instructions: 
Download R script from NanoString [website](https://nanostring.com/) or [GitHub](https://github.com/Nanostring-Biostats/CosMxDACustomModules) 

Instructions from [User Manual](https://university.nanostring.com/cosmx-smi-data-analysis-user-manual)

From the Pipeline Run panel, click the custom modules icon. 
![image](https://github.com/Nanostring-Biostats/CosMxDACustomModules/assets/40255151/bbf9e45f-f4d9-4b9f-9de2-9ef69304f0e4)

The Custom Modules window opens. Edit or delete existing custom modules with the pencil and trash icons, respectively.
![image](https://github.com/Nanostring-Biostats/CosMxDACustomModules/assets/40255151/53e82ef4-12d4-45a8-9be3-6f1fb853fbe6)

To add a new custom module to the pipeline, click Add Module. The Add New Module window opens. Enter information about the new module, including name (required), description, parameters, packages (click + to input package information), arguments (click + to input arguments), and/or entry point. Additionally, you can upload custom script files. Variables required and recommended parameters for individual scripts will be located in the CosMxDA[*ModuleName*]Setup.docx file in each folder.

![image](https://github.com/Nanostring-Biostats/CosMxDACustomModules/assets/40255151/ec77eccc-206d-4ed5-b955-481f6137019a)

Click Save to exit and add the defined module to the Custom Modules list, or exit without saving by clicking the x in the top right.
