# Distribution Modelling of *Curculigo orchioides* in India using MaxEnt

# Overview
This project presents a species distribution modelling analysis of *Curculigo orchioides* (kali musli/black musli), a medicinal herb belonging to the Hypoxidaceae family, across India.

The study uses the Maximum Entropy (MaxEnt) modelling approach to estimate potential habitat suitability based on species occurrence records and environmental predictors, including WorldClim bioclimatic variables, terrain slope, and Sentinel-2 Normalised Difference Vegetation Index (NDVI).

The analysis combines spatial data cleaning, environmental variable selection, correlation analysis, ecological modelling, and habitat suitability mapping to identify environmental factors associated with the species' distribution and explore applications for conservation planning.

# Objectives:
1. Model the potential habitat suitability of *Curculigo orchioides* across the study area in India.
2. Clean and spatially validate species occurrence records.
3. Prepare environmental predictor layers using climatic, vegetation, and topographic data.
4. Identify and reduce multicollinearity among environmental variables.
5. Evaluate model performance using AUC, omission rates, and gain statistics.
6. Identify the environmental variables contributing most to predicted habitat suitability.
7. Interpret the predicted distribution in the context of medicinal plant conservation.
8. Explore potential applications of species distribution modelling in ecological research and conservation planning.

# Study Area:
India, within the spatial extent of:

1. Longitude: 72°E–90°E
2. Latitude: 10.5°N–34°N

The study extent was restricted to the available Sentinel-2 NDVI coverage.

# Tools Used:
1. R
2. RStudio
3. MaxEnt 3.4.4
4. R packages: terra, geodata, corrplot
5. WorldClim 2.1
6. Sentinel-2 satellite data
7. Google Earth Engine
8. GADM Administrative Boundaries

# Data Sources:
1. Species occurrence records: Course-provided *Curculigo orchioides* dataset
2. Bioclimatic variables: WorldClim 2.1
3. Elevation: WorldClim 2.1 (SRTM-derived)
4. Vegetation index: Sentinel-2 NDVI composite
5. Administrative boundaries: GADM 4.1

# Data Preparation:
The initial occurrence dataset contained 985 records of *Curculigo orchioides* with geographic coordinates and associated ecological information.

Data preparation was conducted in R using the following steps:

1. Removed 18 records with duplicate coordinates.
2. Checked coordinate validity and spatial alignment with India's administrative boundaries.
3. Restricted records to the available study extent.
4. Removed records with missing values in environmental predictor layers.
5. Performed spatial thinning by retaining one occurrence per approximately 1 km grid cell.
6. Prepared the final occurrence dataset for MaxEnt modelling.

A total of 910 spatially thinned occurrence records were retained for modelling.

# Environmental Variable Selection:
A total of 22 candidate environmental predictors were assembled, consisting of:

1. 19 WorldClim bioclimatic variables (BIO1–BIO19)
2. Elevation
3. Terrain slope
4. Sentinel-2 NDVI

Pearson correlation analysis was performed using 10,000 randomly sampled environmental grid cells.

Highly correlated predictors were assessed to reduce multicollinearity and improve model interpretability.

The final model retained ten predictors:

1. NDVI
2. Slope
3. BIO2 – Mean Diurnal Range
4. BIO4 – Temperature Seasonality
5. BIO9 – Mean Temperature of Driest Quarter
6. BIO13 – Precipitation of Wettest Month
7. BIO14 – Precipitation of Driest Month
8. BIO15 – Precipitation Seasonality
9. BIO18 – Precipitation of Warmest Quarter
10. BIO19 – Precipitation of Coldest Quarter

# MaxEnt Modelling:
Species distribution modelling was conducted using MaxEnt 3.4.4.

Model configuration:
1. Presence records: 910
2. Training records: 683 (75%)
3. Test records: 227 (25%)
4. Background points: 10,000
5. Feature classes: Linear, Quadratic, Product, and Hinge
6. Regularisation multiplier: 1
7. Maximum iterations: 500
8. Output format: Logistic
9. Replicates: 1

Model evaluation included:
1. Area Under the Curve (AUC)
2. Training and test gain
3. Omission rates
4. Permutation importance
5. Jackknife analysis
6. Environmental response curves

# Model Performance:
The MaxEnt model achieved:

1. Training AUC: 0.901
2. Test AUC: 0.883 (SD 0.008)
3. Regularised training gain: 1.170
4. Unregularised training gain: 1.283
5. Test gain: 1.121

The results indicate good discrimination between known occurrence locations and background points within the selected study extent.

# Key Findings:
### 1. Predicted Habitat Suitability
The model identified potential habitat suitability primarily across eastern and central India, particularly:

- Odisha
- Chhattisgarh
- Jharkhand
- Eastern Madhya Pradesh
- Parts of the northern Eastern Ghats

Smaller suitable zones were identified along the Western Ghats and Himalayan foothills.

At the 10th-percentile training-presence threshold (logistic 0.256), approximately 435,900 km², representing 16.8% of the modelled area, was classified as suitable habitat.

### 2. Environmental Variable Importance
NDVI emerged as the most influential predictor:

- Percentage contribution: 37.0%
- Permutation importance: 35.5%

Slope ranked second, with a contribution of 21.7%.

Other influential predictors included mean diurnal temperature range (BIO2), wettest-month precipitation (BIO13), and temperature seasonality (BIO4).

### 3. Ecological Interpretation
The results suggest that vegetation greenness and terrain characteristics are important in explaining the predicted distribution of *C. orchioides*.

The climatic predictors further indicate the relevance of monsoon rainfall patterns and temperature variability in defining potential habitat suitability.

### 4. Conservation Applications
The predicted suitability maps can help identify areas for targeted field surveys, conservation prioritisation, and potential restoration initiatives.

Areas with high predicted suitability but relatively few recorded occurrences may be useful targets for future ecological surveys.

# Limitations:
1. Occurrence records were geographically clustered, particularly around Odisha and Chhattisgarh.
2. Random training and testing splits may not fully account for spatial autocorrelation.
3. Sentinel-2 NDVI coverage restricted the geographical extent of the model.
4. Environmental datasets represented different time periods.
5. The approximately 1 km spatial resolution cannot capture fine-scale habitat characteristics.
6. The single model replicate and default iteration limit may affect model robustness.
7. Spatially blocked cross-validation and regularisation tuning would provide stronger model evaluation.

The predicted suitability map should therefore be interpreted as an ecological screening tool rather than a definitive map of species presence.

# Skills Demonstrated:
1. Statistical Ecology
2. Species Distribution Modelling
3. MaxEnt Modelling
4. Spatial Data Analysis
5. R Programming
6. Geospatial Data Processing
7. Data Cleaning
8. Environmental Variable Selection
9. Correlation Analysis
10. Multicollinearity Assessment
11. Raster Data Processing
12. Ecological Interpretation
13. Model Evaluation
14. Habitat Suitability Mapping
15. Conservation Analytics
16. Scientific Research and Reporting

# Repository Contents:
- Academic Research Report
- R Statistical Ecology Script
- Species Occurrence Dataset
- MaxEnt Model Outputs
- Habitat Suitability Predictions
- Response Curves
- Jackknife and ROC Analysis
- Model Evaluation Results
