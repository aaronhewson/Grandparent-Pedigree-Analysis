# Grandparent pedigree analysis of samples, using Mendelian inconsistent error counts.
# For identified parent-offspring duos NOT part of a trio (ie with an unidentified parent). Several duos with known pedigrees are included to determine thresholds of acceptance.
# Part of Aaron Hewson's MSc research.No triploids or duplicates included in duo list.
# Methodology and code adapted from Muranty et al. 2020, https://doi.org/10.1186/s12870-019-2171-6. 

# Script requires data to be in PLINK binary (.bed/.bim/.fam) format. 
# PLINK version used: PLINK v1.9.0-b.7.7 64-bit (22 Oct 2024)

# Aaron Hewson's 2024 and 2025 samples all included in dataset. 

# Create objects indicating status of parent-offspring duos (adapted from Muranty script "create_status_duos.R") ----

# Set file input/output paths ---------------------------------------------

# Directory and root filename of SNP data
# Where PLINK files are "xxx.bed", "xxx.bim", "xxx.fam", the rootfilename is "xxx"
SNPdata.filepath <- "C:/Users/curly/Desktop/Apple Genotyping/Methods/Grandparent Pedigree Analysis/Inputs"
SNPdata.rootfilename <- "Geno_BED"

# File containing list of parent-offspring duos to test (ensure in the same directory as SNP data!)
duos.to.test.filename <- "Duos_To_Test.txt"

# Name for the exported ".Rdata" file with status objects (will be exported in same directory as SNP data)
duos.status.filename <- "Duo_Status_Objects.Rdata"

# Directory containing supplementary script "rutilsHMuranty_pedigree.R"
function.file.path  <- "C:/Users/curly/Desktop/Apple Genotyping/Methods/Grandparent Pedigree Analysis/Inputs/Muranty_Scripts"


# Load packages and supplementary script ----------------------------------

# Note - requires snpStats
# if not installed, run:
#    install.packages("BiocManager")
#    BiocManager::install("snpStats")

# Load packages
library(snpStats)

# Load supplementary script
source(paste(function.file.path, "rutilsHMuranty_pedigree.R", sep = "/"))


# Load SNP data -----------------------------------------------------------

data <- read.plink(paste(SNPdata.filepath, SNPdata.rootfilename, sep = "/"))

genot.mat <- matrix(as.integer(data$genotypes@.Data),
                    nrow = nrow(data$genotypes@.Data),
                    dimnames = dimnames(data$genotypes@.Data))
genot.mat[genot.mat == 0] <- NA
genot.mat <- genot.mat - 1
genotAB.mat <- gt012toAB(genot.mat)
rm(genot.mat)

duos.to.test.df <- read.table(
  file = paste(SNPdata.filepath,duos.to.test.filename, sep = "/"),
  header = TRUE,
  stringsAsFactors = FALSE)
nb.duos.to.test <- nrow(duos.to.test.df)



# Build objects with parent-offspring duos --------------------------------

duos.status.GP.notBB <- array(NA, dim = c(2 * nb.duos.to.test, ncol(genotAB.mat)),
                              dimnames = list(c(paste(duos.to.test.df[, "IID1"],
                                                      duos.to.test.df[, "IID2"], sep = "-"),
                                                paste(duos.to.test.df[, "IID2"],
                                                      duos.to.test.df[, "IID1"], sep = "-")),
                                              colnames(genotAB.mat)))
duos.status.GP.notAA <- array(NA, dim = c(2 * nb.duos.to.test, ncol(genotAB.mat)),
                              dimnames = list(c(paste(duos.to.test.df[, "IID1"],
                                                      duos.to.test.df[, "IID2"], sep = "-"),
                                                paste(duos.to.test.df[, "IID2"],
                                                      duos.to.test.df[, "IID1"], sep = "-")),
                                              colnames(genotAB.mat)))
nb.inform.mk <- rep(NA, 2 * nb.duos.to.test)
names(nb.inform.mk) <- c(paste(duos.to.test.df[, "IID1"],
                               duos.to.test.df[, "IID2"], sep = "-"),
                         paste(duos.to.test.df[, "IID2"],
                               duos.to.test.df[, "IID1"], sep = "-"))
start.time.out <- Sys.time()
for(duo.index in #1:100) {
    seq(nb.duos.to.test)) {
  #  start.time <- Sys.time()
  iid1.name <- duos.to.test.df[duo.index, "IID1"]
  iid2.name <- duos.to.test.df[duo.index, "IID2"]
  iid1.data <- genotAB.mat[iid1.name,]
  iid2.data <- genotAB.mat[iid2.name,]
  # individual named in column IID1 is the parent
  # individual named in column IID2 is the offspring,
  duos.status.GP.notBB[duo.index,] <-
    (iid2.data == "AA" & iid1.data != "BB") |
    (iid2.data == "AB" & iid1.data == "BB")
  duos.status.GP.notAA[duo.index,] <-
    (iid2.data == "BB" & iid1.data != "AA") |
    (iid2.data == "AB" & iid1.data == "AA")
  nb.inform.mk[duo.index] <- sum(duos.status.GP.notBB[duo.index,]) +
    sum(duos.status.GP.notAA[duo.index,])
  # individual named in column IID1 is the offspring,
  # individual named in column IID2 is the parent
  duos.status.GP.notBB[duo.index + nb.duos.to.test,] <-
    (iid1.data == "AA" & iid2.data != "BB") |
    (iid1.data == "AB" & iid2.data == "BB")
  duos.status.GP.notAA[duo.index + nb.duos.to.test,] <-
    (iid1.data == "BB" & iid2.data != "AA") |
    (iid1.data == "AB" & iid2.data == "AA")
  nb.inform.mk[duo.index + nb.duos.to.test] <-
    sum(duos.status.GP.notBB[duo.index + nb.duos.to.test,]) +
    sum(duos.status.GP.notAA[duo.index + nb.duos.to.test,])
  #  finish.time <- Sys.time()
  #  cat("duo", iid1.name, "-", iid2.name,
  #      "examined in", finish.time-start.time, fill = T)
}
finish.time.out <- Sys.time()
cat("all duos examined in", finish.time.out - start.time.out, fill = T)

object.size(duos.status.GP.notBB)
object.size(duos.status.GP.notAA)
object.size(nb.inform.mk)

save(duos.status.GP.notBB, duos.status.GP.notAA, nb.inform.mk,
     file = paste(SNPdata.filepath, duos.status.filename, sep = "/"))


# Investigate grandparent pairs for identified parent-offspring duos (adapated from Muranty script "GP_pedigree_search.R") ----

# Set file input/output paths ---------------------------------------------

# Directory and root filename of SNP data
# Where PLINK files are "xxx.bed", "xxx.bim", "xxx.fam", the rootfilename is "xxx"
SNPdata.filepath <- "C:/Users/curly/Desktop/Apple Genotyping/Methods/Grandparent Pedigree Analysis/Inputs"
SNPdata.rootfilename <- "Geno_BED"

# File containing list of parent-offspring duos to test (ensure in the same directory as SNP data!)
duos.to.test.filename <- "Duos_To_Test.txt"

# Name for the exported ".Rdata" file with status objects (will be exported in same directory as SNP data)
duos.status.filename <- "Duo_Status_Objects.Rdata"

# File path for exporting the results
results.file.path <- "C:/Users/curly/Desktop/Apple Genotyping/Results/Grandparent Pedigree Analysis"

# Names for the result files
result.filename.part1 <- "Results_1.txt"
result.filename.part2 <- "Results_1.txt"

# Directory containing supplementary script "rutilsHMuranty_pedigree.R"
function.file.path  <- "C:/Users/curly/Desktop/Apple Genotyping/Methods/Grandparent Pedigree Analysis/Inputs/Muranty_Scripts"


# Define thresholds -------------------------------------------------------

# Acceptable threshold for exporting a result
error.threshold <- 0.1

# Define a part of the individuals to test as GP1
start.i <- 1
end.i <- 2724

# Load packages and supplementary script ----------------------------------

# Note - requires snpStats
# if not installed, run:
#    install.packages("BiocManager")
#    BiocManager::install("snpStats")

# Load packages
library(snpStats)

# Load supplementary script
source(paste(function.file.path, "rutilsHMuranty_pedigree.R", sep = "/"))

# Load SNP data -----------------------------------------------------------
data <- read.plink(paste(SNPdata.filepath, SNPdata.rootfilename, sep = "/"))


duos.to.test.df <- read.table(
  file = paste(SNPdata.filepath, duos.to.test.filename, sep = "/"),
  header = TRUE,
  stringsAsFactors = FALSE)
load(file = paste(SNPdata.filepath, duos.status.filename, sep = "/"))

genot.mat <- matrix(as.integer(data$genotypes@.Data),
                    nrow = nrow(data$genotypes@.Data),
                    dimnames = dimnames(data$genotypes@.Data))
genot.mat[genot.mat == 0] <- NA
genot.mat <- genot.mat - 1
genotAB.mat <- gt012toAB(genot.mat)
rm(genot.mat)

all.indiv <- rownames(genotAB.mat)
nb.genot <- length(all.indiv)

nb.duos.to.test <- nrow(duos.to.test.df)


# Create result file ------------------------------------------------------

# Define result file name
result.filename <- paste(result.filename.part1, start.i, "_", end.i,
                         result.filename.part2, ".txt", sep = "")
result.file <- paste(results.file.path, result.filename, sep = "/")

# Create empty result file
GP.duo.test.df <- data.frame(
  array(NA, dim = c(0, 6),
        dimnames = list(NULL, c("offspring", "parent", "GP1", "GP2", "nb.inform.mk", "nb.inc.all.mk"))),
  stringsAsFactors = FALSE)
write.table(GP.duo.test.df, file = result.file, quote = FALSE)


# Run grandparent search loop ---------------------------------------------

start.time.out <- Sys.time()
for(ind.i in seq(from = start.i, to = end.i)) {
  start.time.in <- Sys.time()
  ind.i.id <- all.indiv[ind.i]
  ind.i.genot.mat <- genotAB.mat[ind.i.id, ]
# all other individuals are potential GP2
  all.ind.j <- seq(from = ind.i + 1, to = #end.i + 1)
                                          nb.genot)
  all.ind.j.id <- all.indiv[all.ind.j]
  all.ind.j.genot.mat <- genotAB.mat[all.ind.j.id,]
  if(length(all.ind.j) > 1) {
# syntax when several GP2 are to be tested
    pairs.status.BB.here <- apply(all.ind.j.genot.mat, 1, function(x, y) {
      x == "BB" & y == "BB" }, y = ind.i.genot.mat)
    pairs.status.AA.here <- apply(all.ind.j.genot.mat, 1, function(x, y) {
      x == "AA" & y == "AA" }, y = ind.i.genot.mat)
    pairs.status.inform.here <- apply(all.ind.j.genot.mat, 1, function(x, y) {
      x != "--" & y != "--" }, y = ind.i.genot.mat)
    colnames(pairs.status.BB.here) <- paste(ind.i.id, all.ind.j.id, sep = "-X-")
    colnames(pairs.status.AA.here) <- paste(ind.i.id, all.ind.j.id, sep = "-X-")
    colnames(pairs.status.inform.here) <- paste(ind.i.id, all.ind.j.id, sep = "-X-")
    for(duo.index in #1:100) {
                     seq(nb.duos.to.test)) {
# individual named in column IID2 is the offspring,
# individual named in column IID1 is the parent
      offspring.ID <- duos.to.test.df[duo.index, "IID2"]
# continue only if the focal GP1 is not the focal offspring
      if(ind.i.id != offspring.ID) {
# remove the focal offspring from the list of potential GP2
        offspring.in.all.ind.j <- which(all.ind.j.id == offspring.ID)
        if(length(offspring.in.all.ind.j) > 0) {
          if((length(all.ind.j.id) - length(offspring.in.all.ind.j)) > 1) {
# syntax when the focal offspring has to be removed from the list of GP2 to test
# but there are still several GP2 to test
            nb.BB.errors <- apply(pairs.status.BB.here[duos.status.GP.notBB[duo.index,],
                                  - offspring.in.all.ind.j], 2, sum, na.rm = T)
            nb.AA.errors <- apply(pairs.status.AA.here[duos.status.GP.notAA[duo.index,],
                                  - offspring.in.all.ind.j], 2, sum, na.rm = T)
            nb.inform <- apply(pairs.status.inform.here[duos.status.GP.notBB[duo.index,],
                               - offspring.in.all.ind.j], 2, sum, na.rm = T) +
                           apply(pairs.status.inform.here[duos.status.GP.notAA[duo.index,],
                                   - offspring.in.all.ind.j], 2, sum, na.rm = T)
            error.rate <- (nb.BB.errors + nb.AA.errors) / nb.inform
            if(any(error.rate < error.threshold)) {
              where.error.rate.low <- which(error.rate < error.threshold)
              nb.GP.pairs.here <- length(where.error.rate.low)
              GP.duo.test.here <- data.frame(
                offspring = rep(offspring.ID, nb.GP.pairs.here),
                parent = rep(duos.to.test.df[duo.index, "IID1"], nb.GP.pairs.here),
                GP1 = rep(ind.i.id, nb.GP.pairs.here),
                GP2 = all.ind.j.id[-offspring.in.all.ind.j][where.error.rate.low],
                nb.inform.mk = nb.inform[where.error.rate.low],
                nb.inc.all.mk = (nb.BB.errors[where.error.rate.low] +
                                   nb.AA.errors[where.error.rate.low]),
                stringsAsFactors = FALSE)
              write.table(GP.duo.test.here, file = result.file, col.names = FALSE,
                        row.names = FALSE, quote = FALSE, append = T)
            }
          } else if((length(all.ind.j.id) - length(offspring.in.all.ind.j)) == 1) {
# syntax when only one GP2 is tested after removing the focal offspring
            nb.BB.errors <- sum(pairs.status.BB.here[duos.status.GP.notBB[duo.index,],
                                - offspring.in.all.ind.j], na.rm = T)
            nb.AA.errors <- sum(pairs.status.AA.here[duos.status.GP.notAA[duo.index,],
                                - offspring.in.all.ind.j], na.rm = T)
            nb.inform <- sum(pairs.status.inform.here[duos.status.GP.notBB[duo.index,],
                                - offspring.in.all.ind.j], na.rm = T) +
                           sum(pairs.status.inform.here[duos.status.GP.notAA[duo.index,],
                                - offspring.in.all.ind.j], na.rm = T)
            error.rate <- (nb.BB.errors + nb.AA.errors) / nb.inform
            if(error.rate < error.threshold) {
              GP.duo.test.here <- data.frame(
                offspring = offspring.ID,
                parent = duos.to.test.df[duo.index, "IID1"],
                GP1 = ind.i.id,
                GP2 = all.ind.j.id[- offspring.in.all.ind.j],
                nb.inform.mk = nb.inform,
                nb.inc.all.mk = (nb.BB.errors + nb.AA.errors),
                stringsAsFactors = FALSE)
              write.table(GP.duo.test.here, file = result.file, col.names = FALSE,
                        row.names = FALSE, quote = FALSE, append = T)
            }
          }
        } else { # the focal offspring is not in the list of GP2 to test
          if((length(all.ind.j.id)) > 1) { # there are several GP2 to test
            nb.BB.errors <- apply(pairs.status.BB.here[duos.status.GP.notBB[duo.index,],
                                  ], 2, sum, na.rm = T)
            nb.AA.errors <- apply(pairs.status.AA.here[duos.status.GP.notAA[duo.index,],
                                  ], 2, sum, na.rm = T)
            nb.inform <- apply(pairs.status.inform.here[duos.status.GP.notBB[duo.index,],
                               ], 2, sum, na.rm = T) +
                           apply(pairs.status.inform.here[duos.status.GP.notAA[duo.index,],
                                   ], 2, sum, na.rm = T)
            error.rate <- (nb.BB.errors + nb.AA.errors) / nb.inform
            if(any(error.rate < error.threshold)) {
              where.error.rate.low <- which(error.rate < error.threshold)
              nb.GP.pairs.here <- length(where.error.rate.low)
              GP.duo.test.here <- data.frame(
                offspring = rep(offspring.ID, nb.GP.pairs.here),
                parent = rep(duos.to.test.df[duo.index, "IID1"], nb.GP.pairs.here),
                GP1 = rep(ind.i.id, nb.GP.pairs.here),
                GP2 = all.ind.j.id[where.error.rate.low],
                nb.inform.mk = nb.inform[where.error.rate.low],
                nb.inc.all.mk = (nb.BB.errors[where.error.rate.low] +
                                   nb.AA.errors[where.error.rate.low]),
                stringsAsFactors = FALSE)
              write.table(GP.duo.test.here, file = result.file, col.names = FALSE,
                        row.names = FALSE, quote = FALSE, append = T)
            }
          } else if(length(all.ind.j.id) == 1) { # there is only one GP2 to test
            nb.BB.errors <- sum(pairs.status.BB.here[duos.status.GP.notBB[duo.index,],
                                ], na.rm = T)
            nb.AA.errors <- sum(pairs.status.AA.here[duos.status.GP.notAA[duo.index,],
                                ], na.rm = T)
            nb.inform <- sum(pairs.status.inform.here[duos.status.GP.notBB[duo.index,],
                                ], na.rm = T) +
                           sum(pairs.status.inform.here[duos.status.GP.notAA[duo.index,],
                                ], na.rm = T)
            error.rate <- (nb.BB.errors + nb.AA.errors) / nb.inform
            if(error.rate < error.threshold) {
              GP.duo.test.here <- data.frame(
                offspring = offspring.ID,
                parent = duos.to.test.df[duo.index, "IID1"],
                GP1 = ind.i.id,
                GP2 = all.ind.j.id,
                nb.inform.mk = nb.inform,
                nb.inc.all.mk = (nb.BB.errors + nb.AA.errors),
                stringsAsFactors = FALSE)
              write.table(GP.duo.test.here, file = result.file, col.names = FALSE,
                        row.names = FALSE, quote = FALSE, append = T)
            }
          }
        }
      }

# individual named in column IID1 is the offspring,
# individual named in column IID2 is the parent
      offspring.ID <- duos.to.test.df[duo.index, "IID1"]
# continue only if the focal GP1 is not the focal offspring
      if(ind.i.id != offspring.ID) {
# remove the focal offspring from the list of potential GP2
        offspring.in.all.ind.j <- which(all.ind.j.id == offspring.ID)
        if(length(offspring.in.all.ind.j) > 0) {
# syntax when the focal offspring has to be removed from the list of GP2 to test
# but there are still several GP2 to test
          if((length(all.ind.j.id) - length(offspring.in.all.ind.j)) > 1) {
            nb.BB.errors <- apply(pairs.status.BB.here[duos.status.GP.notBB[duo.index + nb.duos.to.test,],
                                  - offspring.in.all.ind.j], 2, sum, na.rm = T)
            nb.AA.errors <- apply(pairs.status.AA.here[duos.status.GP.notAA[duo.index + nb.duos.to.test,],
                                  - offspring.in.all.ind.j], 2, sum, na.rm = T)
            nb.inform <- apply(pairs.status.inform.here[duos.status.GP.notBB[duo.index + nb.duos.to.test,],
                               - offspring.in.all.ind.j], 2, sum, na.rm = T) +
                           apply(pairs.status.inform.here[duos.status.GP.notAA[duo.index + nb.duos.to.test,],
                                   - offspring.in.all.ind.j], 2, sum, na.rm = T)
            error.rate <- (nb.BB.errors + nb.AA.errors) / nb.inform
            if(any(error.rate < error.threshold)) {
              where.error.rate.low <- which(error.rate < error.threshold)
              nb.GP.pairs.here <- length(where.error.rate.low)
              GP.duo.test.here <- data.frame(
                offspring = rep(offspring.ID, nb.GP.pairs.here),
                parent = rep(duos.to.test.df[duo.index, "IID2"], nb.GP.pairs.here),
                GP1 = rep(ind.i.id, nb.GP.pairs.here),
                GP2 = all.ind.j.id[-offspring.in.all.ind.j][where.error.rate.low],
                nb.inform.mk = nb.inform[where.error.rate.low],
                nb.inc.all.mk = (nb.BB.errors[where.error.rate.low] +
                                   nb.AA.errors[where.error.rate.low]),
                stringsAsFactors = FALSE)
              write.table(GP.duo.test.here, file = result.file, col.names = FALSE,
                        row.names = FALSE, quote = FALSE, append = T)
            }
          } else if((length(all.ind.j.id) - length(offspring.in.all.ind.j)) == 1) {
# syntax when only one GP2 is tested after removing the focal offspring
            nb.BB.errors <- sum(pairs.status.BB.here[duos.status.GP.notBB[duo.index + nb.duos.to.test,],
                                - offspring.in.all.ind.j], na.rm = T)
            nb.AA.errors <- sum(pairs.status.AA.here[duos.status.GP.notAA[duo.index + nb.duos.to.test,],
                                - offspring.in.all.ind.j], na.rm = T)
            nb.inform <- sum(pairs.status.inform.here[duos.status.GP.notBB[duo.index + nb.duos.to.test,],
                                - offspring.in.all.ind.j], na.rm = T) +
                           sum(pairs.status.inform.here[duos.status.GP.notAA[duo.index + nb.duos.to.test,],
                                - offspring.in.all.ind.j], na.rm = T)
            error.rate <- (nb.BB.errors + nb.AA.errors) / nb.inform
            if(error.rate < error.threshold) {
              GP.duo.test.here <- data.frame(
                offspring = offspring.ID,
                parent = duos.to.test.df[duo.index, "IID2"],
                GP1 = ind.i.id,
                GP2 = all.ind.j.id[- offspring.in.all.ind.j],
                nb.inform.mk = nb.inform,
                nb.inc.all.mk = (nb.BB.errors + nb.AA.errors),
                stringsAsFactors = FALSE)
              write.table(GP.duo.test.here, file = result.file, col.names = FALSE,
                        row.names = FALSE, quote = FALSE, append = T)
            }
          }
        } else { # the focal offspring is not in the list of GP2 to test
          if((length(all.ind.j.id)) > 1) { # there are several GP2 to test
            nb.BB.errors <- apply(pairs.status.BB.here[duos.status.GP.notBB[duo.index + nb.duos.to.test,],
                                  ], 2, sum, na.rm = T)
            nb.AA.errors <- apply(pairs.status.AA.here[duos.status.GP.notAA[duo.index + nb.duos.to.test,],
                                  ], 2, sum, na.rm = T)
            nb.inform <- apply(pairs.status.inform.here[duos.status.GP.notBB[duo.index + nb.duos.to.test,],
                               ], 2, sum, na.rm = T) +
                           apply(pairs.status.inform.here[duos.status.GP.notAA[duo.index + nb.duos.to.test,],
                                   ], 2, sum, na.rm = T)
            error.rate <- (nb.BB.errors + nb.AA.errors) / nb.inform
            if(any(error.rate < error.threshold)) {
              where.error.rate.low <- which(error.rate < error.threshold)
              nb.GP.pairs.here <- length(where.error.rate.low)
              GP.duo.test.here <- data.frame(
                offspring = rep(offspring.ID, nb.GP.pairs.here),
                parent = rep(duos.to.test.df[duo.index, "IID2"], nb.GP.pairs.here),
                GP1 = rep(ind.i.id, nb.GP.pairs.here),
                GP2 = all.ind.j.id[where.error.rate.low],
                nb.inform.mk = nb.inform[where.error.rate.low],
                nb.inc.all.mk = (nb.BB.errors[where.error.rate.low] +
                                   nb.AA.errors[where.error.rate.low]),
                stringsAsFactors = FALSE)
              write.table(GP.duo.test.here, file = result.file, col.names = FALSE,
                        row.names = FALSE, quote = FALSE, append = T)
            }
          } else if(length(all.ind.j.id) == 1) { # there is only one GP2 to test
            nb.BB.errors <- sum(pairs.status.BB.here[duos.status.GP.notBB[duo.index + nb.duos.to.test,],
                                ], na.rm = T)
            nb.AA.errors <- sum(pairs.status.AA.here[duos.status.GP.notAA[duo.index + nb.duos.to.test,],
                                ], na.rm = T)
            nb.inform <- sum(pairs.status.inform.here[duos.status.GP.notBB[duo.index + nb.duos.to.test,],
                                ], na.rm = T) +
                           sum(pairs.status.inform.here[duos.status.GP.notAA[duo.index + nb.duos.to.test,],
                                ], na.rm = T)
            error.rate <- (nb.BB.errors + nb.AA.errors) / nb.inform
            if(error.rate < error.threshold) {
              GP.duo.test.here <- data.frame(
                offspring = offspring.ID,
                parent = duos.to.test.df[duo.index, "IID2"],
                GP1 = ind.i.id,
                GP2 = all.ind.j.id,
                nb.inform.mk = nb.inform,
                nb.inc.all.mk = (nb.BB.errors + nb.AA.errors),
                stringsAsFactors = FALSE)
              write.table(GP.duo.test.here, file = result.file, col.names = FALSE,
                        row.names = FALSE, quote = FALSE, append = T)
            }
          }
        }
      }
    }
  } else { # there is only one GP2 to test
    pairs.status.BB.here <- all.ind.j.genot.mat == "BB" & ind.i.genot.mat == "BB"
    pairs.status.AA.here <- all.ind.j.genot.mat == "AA" & ind.i.genot.mat == "AA"
    pairs.status.inform.here <- all.ind.j.genot.mat != "--" & ind.i.genot.mat != "--"
    for(duo.index in #1:100) {
                     seq(nb.duos.to.test)) {
# individual named in column IID2 is the offspring,
# individual named in column IID1 is the parent
      offspring.ID <- duos.to.test.df[duo.index, "IID2"]
# continue only if neither the focal GP1 not the focal GP2 are the focal offspring
      if((ind.i.id != offspring.ID) & (all.ind.j.id != offspring.ID)) {
        nb.BB.errors <- sum(pairs.status.BB.here[duos.status.GP.notBB[duo.index,]], na.rm = T)
        nb.AA.errors <- sum(pairs.status.AA.here[duos.status.GP.notAA[duo.index,]], na.rm = T)
        nb.inform <- sum(pairs.status.inform.here[duos.status.GP.notBB[duo.index,]], na.rm = T) +
                       sum(pairs.status.inform.here[duos.status.GP.notAA[duo.index,]], na.rm = T)
        error.rate <- (nb.BB.errors + nb.AA.errors) / nb.inform
        if(error.rate < error.threshold) {
          GP.duo.test.here <- data.frame(
            offspring = offspring.ID,
            parent = duos.to.test.df[duo.index, "IID1"],
            GP1 = ind.i.id,
            GP2 = all.ind.j.id,
            nb.inform.mk = nb.inform,
            nb.inc.all.mk = (nb.BB.errors + nb.AA.errors),
            stringsAsFactors = FALSE)
          write.table(GP.duo.test.here, file = result.file, col.names = FALSE,
                      row.names = FALSE, quote = FALSE, append = T)
        }
      }
# individual named in column IID1 is the offspring,
# individual named in column IID2 is the parent
      offspring.ID <- duos.to.test.df[duo.index, "IID1"]
# continue only if neither the focal GP1 not the focal GP2 are the focal offspring
      if((ind.i.id != offspring.ID) & (all.ind.j.id != offspring.ID)) {
        nb.BB.errors <- sum(pairs.status.BB.here[duos.status.GP.notBB[
                                duo.index + nb.duos.to.test,]], na.rm = T)
        nb.AA.errors <- sum(pairs.status.AA.here[duos.status.GP.notAA[
                                duo.index + nb.duos.to.test,]], na.rm = T)
        nb.inform <- sum(pairs.status.inform.here[duos.status.GP.notBB[
                                      duo.index + nb.duos.to.test,]], na.rm = T) +
                       sum(pairs.status.inform.here[duos.status.GP.notAA[
                                      duo.index + nb.duos.to.test,]], na.rm = T)
        error.rate <- (nb.BB.errors + nb.AA.errors) / nb.inform
        if(error.rate < error.threshold) {
          GP.duo.test.here <- data.frame(
            offspring = offspring.ID,
            parent = duos.to.test.df[duo.index, "IID2"],
            GP1 = ind.i.id,
            GP2 = all.ind.j.id,
            nb.inform.mk = nb.inform,
            nb.inc.all.mk = (nb.BB.errors + nb.AA.errors),
            stringsAsFactors = FALSE)
          write.table(GP.duo.test.here, file = result.file, col.names = FALSE,
                      row.names = FALSE, quote = FALSE, append = T)
        }
      }
    }
  }
  finish.time.in <- Sys.time()
  cat("all pairs involving", ind.i.id,
      "examined in", finish.time.in - start.time.in, fill = T)
}
finish.time.out <- Sys.time()

cat("all GP pairs examined in", finish.time.out - start.time.out, fill = T)
