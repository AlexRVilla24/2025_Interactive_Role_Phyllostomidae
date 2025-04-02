# Niche Structure and Interactive Role of Phyllostomidae
Data and code for the analyses from the study "Ecological niche structure and the interactive role of leaf-nosed bats (Phyllostomidae) within frugivory networks: Another dimension of the niche centroid hypothesis " submitted to JOURNAL.

`Code_Phyllostomid_Interactive_Role.R`: R master script including all steps of data handling and analyses. These steps are described below.

## Code is composed of six sections:
`PHASE I` - Network extraction (Steps 1-9): This section renames and reorganizes the bat-plant database ("OG_NTW.csv" in the Data Folder), selecting only networks with a 6x6 or higher dimensions and species that are present in 6 or more networks. 

`PHASE II` - Network Analysis (Steps 10-20): This section calculates the interactive role of each bat species (PCA built from the "species level" metrics of degree, closeness and betwenness).

`PHASE III` - Climatic variables and Species occurrence records (Steps 21-22): This section imports the climatic variables, distribution polygons and occurence records (all available in the `Data` folder for the 20 species). Note that since we use a random splitting for the train and test data, we only performed this step once and saved the train and test datasets for each species (also in the Inputs folder). The code for the splitting is also provided (lines 221-227). 

`PHASE IV` - Ecological Niche Modeling (Steps 23-28): This section extracts the climatic information from each raster using the train and test data to perform a correlation test and select the least correlated variables for the modelling. Then it generates multiple models using different numbers of variables to select the best fitting ellipsoid for the occurrence data and background data (using the omr and the auc values). Finally with the best model it generates the ellipsoid and the suitability raster that we used to extract the suitability value for each network. 

`PHASE V` - Statistical Models (Steps 29-33): This section combines the interactive roles and suitability values the selected species in a single dataframe. Then it evaluates the how much the suitability can explain the changes in the interactive role using different statistical models.

`PHASE VI` - Exporting Figures and Results (Steps 34 - 40): This section exports tables and final figures used for the manuscript. 

## Structure of the Data Folder:
`Data`: This folder contains all (input) files necessary to replicate the analyses: <br />
  -Species occurrence records (.csv archives), <br />
  -Species distribution areas as range polygons (.cpg, .dbf, .prj, .shp, .shx), <br />
  -Climatic variables for each species cropped with the distribution polygon (.tif archives), <br />
  -Polygon of the Americas used for the maps ("Shp_America.shp") and <br />
  -Database of bat-plant frugivory ("OG_NTW.csv") <br />

`Data/Outputs`: Empty folder to save resulting objects from the running the code

`Data/Outputs/Example`: Folder with examples of all archives resulting from running the code for one species (e.g., Artibeus fimbriatus). These include NTWs plots, figures for the manuscript, tables with interactive role and suitability values for each species, SDM models and ellipsoids. 

## Contact
alexvilla5920@gmail.com