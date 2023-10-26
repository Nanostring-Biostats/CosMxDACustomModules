NanoString Collection of CosMx SMIDA custom module scripts
=================
This repository contains R script that can be used in custom modules in the AtoMx SMIDA.
 

| Name              | Description                                                                                                                                                                                                                                                                                                                                                                                                                                         |Links|
| :---------------- | :-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |:------------------------|
| Export     | Export all files associated with current study: Seurat, TileDB array, and/or raw files.| [R script](https://github.com/Nanostring-Biostats/CosMxDACustomModules/blob/main/Export/CosMxDAExport.R), [Documentation](https://github.com/Nanostring-Biostats/CosMxDACustomModules/blob/main/Export/CosMxDAExportSetup.docx) | 
| Flat File Export     | Export Flat files used in (Seurat's LoadNanostring function)[https://satijalab.org/seurat/articles/spatial_vignette_2#human-lung-nanostring-cosmx-spatial-molecular-imager]| [R script](https://github.com/Nanostring-Biostats/CosMxDACustomModules/blob/main/flatFileExport/flatFileExport.R), [Documentation](https://github.com/Nanostring-Biostats/CosMxDACustomModules/blob/main/flatFileExport/CosMxDAFlatFileExportSetup.docx) |   

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
