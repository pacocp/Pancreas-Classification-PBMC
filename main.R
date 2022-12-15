library(KnowSeq)

########## SIGNATURE OBTENTION FROM PLASMA ###########
# Plastma series from GEO
GSE133684_matrix =  read.csv("GSE133684_exp_TPM-all.txt", sep='\t')
GSErawData = GSE133684_matrix
rownames(GSErawData) <- GSErawData$X
GSErawData <- GSErawData[,-1]
GSErawData <- as.matrix(GSErawData)
GSErawData <- apply(GSErawData, 2, unlist)
rownames <- rownames(GSErawData)
GSErawData <-  apply(GSErawData, 2, function(y) as.integer(gsub(",", ".", y)))
rownames(GSErawData) <- rownames
GSErownames <- rownames(GSErawData)

annotation <- getGenesAnnotation(GSErownames)

# get only those that have been annotated
genes_to_save_ids = match(annotation$ensembl_gene_id, rownames(GSErawData))
GSErawData = GSErawData[genes_to_save_ids,]
elements_to_assign = match(rownames(GSErawData), annotation$ensembl_gene_id)
rownames(GSErawData) = annotation$external_gene_name[elements_to_assign]

split_label = function(x){
  value = unlist(strsplit(x, split='_'))[1]
  return(value)
}

labels = unlist(lapply(colnames(GSErawData), split_label))
labels = as.character(labels)

#### READING SERUM DATA ####

rawData <- read.csv("RAW_expression.csv",sep = ";")
rawData <- rawData[,1:98]

rownames(rawData) <- rawData$Row.names
rawData <- rawData[,-1]
rawData <- as.matrix(rawData)
rawData <- apply(rawData, 2, unlist)
rownames <- rownames(rawData)
rawData <-  apply(rawData, 2, function(y) as.integer(gsub(",", ".", y)))
rownames(rawData) <- rownames
losPos <- which(rowMeans(rawData) <= 10)
rawData <- rawData[-losPos,]
rownames <- rownames(rawData)

annotation <- getGenesAnnotation(rownames,attributes = c("external_gene_name","percentage_gene_gc_content")
                                 , filter = "external_gene_name")

labelsAll <- read.csv("labels",header = F)
labelsAll <- as.character(labelsAll$V1)

geneExprs <- calculateGeneExpressionValues(rawData, annotation, Ensembl_ID = FALSE)

# Getting those genes that both datasets have in common
genes_intersection = intersect(rownames(geneExprs), rownames(GSErawData))

geneExprs = geneExprs[genes_intersection, ]
GSErawData = GSErawData[genes_intersection, ]

#### FINDING DIFFERENTIALLY EXPRESSED GENES ####

DEGsResults <- DEGsExtraction(GSErawData, labels,lfc = 1, pvalue = 1, cov = 1)
DEGsMatrix <- DEGsResults$DEG_Results$DEGs_Matrix
MLMatrix <- t(DEGsMatrix)

# finding signature using mRMR
mrmrRanking <- featureSelection(MLMatrix, labels,
                                vars_selected = colnames(MLMatrix))

# 5-Fold results on plasma dataset
MRMR_svm_results <- svm_trn(MLMatrix, labels, vars_selected = names(mrmrRanking[1:15]), numFold = 5)

dataPlot(rbind(MRMR_svm_results$accuracyInfo$meanAccuracy,
               MRMR_svm_results$accuracyInfo$standardDeviation),
         labels = as.factor(labels),
         main = "Classification in Plasma mRMR signature SVM",
         mode = "classResults",
         legend = c('Mean accuracy', 'Mean std'),
         ygrid = T)

dataPlot(DEGsMatrix[names(mrmrRanking[1:15]),],
         labels = as.factor(labels), main="Boxplots mRMR Plasma", 
         mode = "genesBoxplot", colours = c("skyblue","gold"))

allCfMats <- MRMR_svm_results$cfMats[[1]]$table + MRMR_svm_results$cfMats[[2]]$table + 
  MRMR_svm_results$cfMats[[3]]$table + MRMR_svm_results$cfMats[[4]]$table + 
  MRMR_svm_results$cfMats[[5]]$table

dataPlot(allCfMats,labels,mode = "confusionMatrix")

signature <- names(mrmrRanking[1:15])

##### TEST ON THE UGR DATA ########

ControlPos <- which(labelsAll == "Control")
PDACPos <- which(labelsAll == "PDAC")
PancreatitisPos <- which(labelsAll == "Pancreatitis")
DiabetsPos <- which(labelsAll == "Diabetes")

triclassExprs <- geneExprs[,-DiabetsPos]
triclassLabels <- labelsAll[-DiabetsPos]

match_signature = match(signature, rownames(triclassExprs))
triclassExprs <- triclassExprs[match_signature,]
triclassExprsML = t(triclassExprs)

MRMR_svm_results_ugr <- svm_trn(triclassExprsML, triclassLabels, vars_selected = colnames(triclassExprsML), 
                                numFold = length(rownames(triclassExprsML)))

dataPlot(MRMR_svm_results_ugr$accuracyInfo$meanAccuracy,
         labels = as.factor(triclassLabels),
         main = "Classification Serum 3 clases with mRMR genes",
         mode = "classResults",
         ygrid = T)

allCfMats <- MRMR_svm_results_ugr$cfMats[[1]]$table
for(i in 2:length(rownames(triclassExprsML))){
  allCfMats <- allCfMats + MRMR_svm_results_ugr$cfMats[[i]]$table
}

dataPlot(allCfMats,triclassLabels, main='UGR confusion matrix with mRMR Plasma Signature', 
         mode = "confusionMatrix", toPNG = TRUE)
