library(terra)
library(geodata)
library(corrplot)

setwd("C:/Users/HP/Desktop/STAT ECO ASSIGN")
dir.create("outputs",       showWarnings = FALSE)
dir.create("data",          showWarnings = FALSE)
dir.create("asc_states",    showWarnings = FALSE)
dir.create("maxent_states", showWarnings = FALSE)

toc <- 23.4366     # Tropic of Cancer latitude (degrees N)

#Chunk 1 Read occurrences and remove exact duplicates
occ <- read.csv("Curculigo orchioides_species_data.csv")
occ <- occ[, c("LON", "LAT")]
occ$species <- "Curculigo_orchioides"

n_raw <- nrow(occ)
occ <- occ[!duplicated(occ[, c("LON", "LAT")]), ]
n_nodup <- nrow(occ)

bad <- is.na(occ$LON) | is.na(occ$LAT) |
  occ$LON < -180 | occ$LON > 180 | occ$LAT < -90 | occ$LAT > 90
occ <- occ[!bad, ]

cat("Raw:", n_raw, "| after exact duplicates removed:", n_nodup, "\n")

#Chunk 2: Count records per state + India boundary check
india  <- gadm("IND", level = 0, path = "data")
states <- gadm("IND", level = 1, path = "data")
sort(states$NAME_1)

bio_files <- file.path("wc2.1_30s_bio", paste0("wc2.1_30s_bio_", 1:19, ".tif"))
template  <- rast(bio_files[1])

pts_all   <- vect(occ, geom = c("LON", "LAT"), crs = "EPSG:4326")
in_india  <- is.related(pts_all, india, "intersects")
cat("Records inside India polygon:", sum(in_india), "| outside:", sum(!in_india), "\n")

occ$state <- extract(states, pts_all)$NAME_1
occ$cell  <- cellFromXY(template, as.matrix(occ[, c("LON", "LAT")]))

state_counts <- data.frame(
  records   = tapply(occ$cell, occ$state, length),
  cells_1km = tapply(occ$cell, occ$state, function(x) length(unique(x)))
)
state_counts <- state_counts[order(-state_counts$cells_1km), ]
print(state_counts)
cat("Points in no state polygon:", sum(is.na(occ$state)), "\n")

#Chunk 3: Load WorldClim bioclim and elevation, then crop
study_ext <- ext(72, 90, 10.5, 34)

stopifnot(all(file.exists(bio_files)))
bio <- crop(rast(bio_files), study_ext)
names(bio) <- paste0("bio", 1:19)

elev_file <- list.files("wc2.1_30s_elev", pattern = "\\.tif$", full.names = TRUE)[1]
elev <- crop(rast(elev_file), study_ext)
names(elev) <- "elev"

slope <- terrain(elev, v = "slope", unit = "degrees")
names(slope) <- "slope"

compareGeom(bio[[1]], elev)

#Chunk 4: Load NDVI
ndvi <- rast("ndvi.tif")
ndvi <- resample(ndvi, bio[[1]], method = "bilinear")
names(ndvi) <- "ndvi"

#Chunk 5: Define the two-state study area and clip everything to it
# Study area: India, limited to the NDVI coverage window (72–90°E, 10.5–34°N)
study_area <- crop(india, study_ext)

env_all <- c(bio, elev, slope, ndvi)
env_s <- crop(env_all, study_area)
env_s <- mask(env_s, study_area)

names(env_s); ext(env_s); res(env_s)

#Chunk 6: Variable table and maps of the clipped layers
var_table <- data.frame(
  Code = names(env_s),
  Min  = round(global(env_s, "min",  na.rm = TRUE)[, 1], 2),
  Max  = round(global(env_s, "max",  na.rm = TRUE)[, 1], 2),
  Mean = round(global(env_s, "mean", na.rm = TRUE)[, 1], 2),
  SD   = round(global(env_s, "sd",   na.rm = TRUE)[, 1], 2)
)
write.csv(var_table, "outputs/Table_all_variables.csv", row.names = FALSE)

png("outputs/Fig2_clipped_layers.png", width = 3000, height = 2400, res = 200)
plot(env_s, nc = 5, maxnl = 22, mar = c(1, 1, 2, 3))
dev.off()

#Chunk 7: Clean the points for the study area
occ_sel <- occ
vals <- extract(env_s, occ_sel[, c("LON", "LAT")])
complete <- complete.cases(vals[, -1])
cat("Removed (missing data in any layer):", sum(!complete), "\n")
occ_sel <- occ_sel[complete, ]
n_inside <- nrow(occ_sel)

cells <- cellFromXY(env_s[["bio1"]], as.matrix(occ_sel[, c("LON", "LAT")]))
occ_thin <- occ_sel[!duplicated(cells), ]
n_thin <- nrow(occ_thin)

cat("Raw:", n_raw, "| no duplicates:", n_nodup,
    "| complete data:", n_inside, "| 1 per 1 km cell:", n_thin, "\n")

write.csv(occ_thin[, c("species", "LON", "LAT")], "occ_maxent_india.csv", row.names = FALSE)

# Fig 1: occurrence records on a map of India
png("outputs/Fig1_occurrence_india.png", width = 2200, height = 2400, res = 250)
plot(states, col = "grey96", border = "grey70", lwd = 0.5, xlim = c(68, 98), ylim = c(6, 37))
plot(study_area, col = "#fdebc8", border = "#c0762b", lwd = 1.8, add = TRUE)
points(occ$LON, occ$LAT, pch = 21, bg = "grey45", col = "white", cex = 0.7, lwd = 0.3)
points(occ_thin$LON, occ_thin$LAT, pch = 21, bg = "#1b6e3a", col = "white", cex = 0.75, lwd = 0.3)
legend("bottomleft", bty = "n", cex = 0.8, pt.cex = 1.1,
       legend = c(paste0("Unique records (n = ", nrow(occ), ")"),
                  paste0("Records used in MaxEnt (n = ", nrow(occ_thin), ")"),
                  "Study area: India (72–90°E, 10.5–34°N)"),
       pch = c(21, 21, 22), pt.bg = c("grey45", "#1b6e3a", "#fdebc8"),
       col = c("white", "white", "#c0762b"))
title(expression(paste("Occurrence records of ", italic("Curculigo orchioides"), " in India")))
dev.off()

#Chunk 8: Correlation check (all variables kept)
set.seed(1)
samp <- spatSample(env_s, 10000, method = "random", na.rm = TRUE, as.df = TRUE)
samp <- na.omit(samp)
cm <- cor(samp)

png("outputs/Fig3_correlation.png", width = 2200, height = 2200, res = 220)
corrplot(cm, method = "color", type = "lower", tl.cex = 0.7)
dev.off()
write.csv(round(cm, 2), "outputs/Table_correlation_matrix.csv")

cm_upper <- cm
cm_upper[lower.tri(cm_upper, diag = TRUE)] <- NA
idx <- which(abs(cm_upper) > 0.8, arr.ind = TRUE)
high_pairs <- data.frame(var1 = rownames(cm)[idx[, 1]],
                         var2 = colnames(cm)[idx[, 2]],
                         r    = round(cm_upper[idx], 2))
high_pairs <- high_pairs[order(-abs(high_pairs$r)), ]
write.csv(high_pairs, "outputs/Table_high_correlation_pairs.csv", row.names = FALSE)
nrow(high_pairs)

# Rank by permutation importance from the preliminary all-variable model;
# keep a variable only if |r| < 0.7 with every variable already kept
prio <- c("ndvi","slope","bio4","bio18","bio13","bio14","bio19","bio16","bio6","bio11","bio17",
          "bio2","bio15","bio5","bio8","elev","bio3","bio1","bio10","bio9","bio12","bio7")
keep <- character(0)
for (v in prio) if (all(abs(cm[v, keep]) < 0.7)) keep <- c(keep, v)
print(keep)
cm_keep <- cm[keep, keep]
cat("Max |r| among kept variables:", round(max(abs(cm_keep[upper.tri(cm_keep)])), 2), "\n")
write.csv(round(cm_keep, 2), "outputs/Table_correlation_kept.csv")
png("outputs/Fig3b_correlation_kept.png", width = 1800, height = 1800, res = 220)
corrplot(cm_keep, method = "number", type = "lower", tl.cex = 0.8, number.cex = 0.7)
dev.off()

#Chunk 9: Export the layers for MaxEnt
unlink("asc_india", recursive = TRUE)
dir.create("asc_india", showWarnings = FALSE)
env_sel <- env_s[[keep]]
for (v in keep) {
  writeRaster(env_sel[[v]], file.path("asc_india", paste0(v, ".asc")),
              NAflag = -9999, overwrite = TRUE)
}
unlink(list.files("asc_india", pattern = "aux.xml$", full.names = TRUE))
list.files("asc_india")
