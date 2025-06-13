library(data.table)
library(LandGenCourse)
library(vegan)    # Used to run PCA & RDA
library(lfmm)     # Used to run LFMM
library(qvalue)   # Used to post-process LFMM output
library(readxl)


args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
	stop("Usage: Rscript your_r_script.R <genmatrix_csv> <env_sheet>")
}

input_file <- args[1]
excel_sheet <- args[2]

data <- read_csv(input_file)

row.names(gen) <- gen$V1

gen <- gen[, -1]

num_rows <- nrow(gen)
num_cols <- ncol(gen)

sum(is.na(gen))
missing_data = sum(is.na(gen))/(num_rows*num_cols)
print(missing_data)

gen.imp <- apply(gen, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
sum(is.na(gen.imp))

env <- read_excel(excel_sheet, sheet = "SC")
str(env)

env$`VCF-BGI code` <- as.character(env$`VCF-BGI code`)
pred <- env[,13:127]
colnames(pred) <- c(
		      "AMT", "MDR", "Iso", "TempS", "MaxTemp", "MinTemp", "TempR", "MTWQ", "MTDQ", "MTWQ", "MTCQ", "AP", "PWM", "PDM", "PS", "PWQ", "PDQ", "PWQ", "PCQ",
		        "TJan", "TFeb", "TMar", "TApr", "TMay", "TJun", "TJul", "TAug", "TSep", "TOct", "TNov", "TDec",
		        "PJan", "PFeb", "PMar", "PApr", "PMay", "PJun", "PJul", "PAug", "PSep", "POct", "PNov", "PDec",
			  "TmaxJan", "TmaxFeb", "TmaxMar", "TmaxApr", "TmaxMay", "TmaxJun", "TmaxJul", "TmaxAug", "TmaxSep", "TmaxOct", "TmaxNov", "TmaxDec",
			  "TminJan", "TminFeb", "TminMar", "TminApr", "TminMay", "TminJun", "TminJul", "TminAug", "TminSep", "TminOct", "TminNov", "TminDec",
			    "AIJan", "AIFeb", "AIMar", "AIApr", "AIMay", "AIJun", "AIJul", "AIAug", "AISep", "AIOct", "AINov", "AIDec",
			    "SradJan", "SradFeb", "SradMar", "SradApr", "SradMay", "SradJun", "SradJul", "SradAug", "SradSep", "SradOct", "SradNov", "SradDec",
			      "VaprJan", "VaprFeb", "VaprMar", "VaprApr", "VaprMay", "VaprJun", "VaprJul", "VaprAug", "VaprSep", "VaprOct", "VaprNov", "VaprDec",
			      "WindJan", "WindFeb", "WindMar", "WindApr", "WindMay", "WindJun", "WindJul", "WindAug", "WindSep", "WindOct", "WindNov", "WindDec"
			      )
pred.pca <- rda(pred, scale=T)
summary(pred.pca)$cont
screeplot(pred.pca, main = "Screeplot: Eigenvalues of Wolf Predictor Variables")
round(scores(pred.pca, choices=1:8, display="species", scaling=0), digits=3)
pred.PC1 <- scores(pred.pca, choices=1, display="sites", scaling=0)

screeplot(pred.pca, main = "Screeplot of Wolf Predictor Variables with Broken Stick", bstick=TRUE, type="barplot")

gen.pca <- rda(gen.imp, scale=T)
screeplot(gen.pca, main = "Screeplot of Genetic Data with Broken Stick", bstick=TRUE, type="barplot")

K <- 25

wolf.lfmm <- lfmm_ridge(Y=gen.imp, X=pred.PC1, K=K) ## c.ange K as you see fit

wolf.pv <- lfmm_test(Y=gen.imp, X=pred.PC1, lfmm=wolf.lfmm, calibrate="gif")
names(wolf.pv)
wolf.pv$gif


