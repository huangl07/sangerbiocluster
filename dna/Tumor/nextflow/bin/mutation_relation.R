#!/usr/bin/env Rscript
start_time <- Sys.time()
library(getopt)
spec <- matrix(c(
  'maf', 'm', 0, 'character',
  'out', 'o', 0, 'character',
  'help', 'h', 0, 'logical'), byrow = TRUE, ncol = 4)
opt <- getopt(spec)
print_usage <- function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("Usage example: \n")
  cat("	
Usage:
	--maf    输入maf表格
  --out    输出基因突变全景图
	--help   usage
\n")
  q(status = 1);
}
if ( !is.null(opt$help))   { print_usage(spec) }
if ( is.null(opt$maf))   { print_usage(spec) }
if ( is.null(opt$out))   { print_usage(spec) }
library(maftools)
library(ggrepel)
library(ggplot2)
options(bitmapType = "cairo")
vc_cols =  c("#E64B35B2", "#4DBBD5B2", "#00A087B2", "#3C5488B2", "#F39B7FB2", 
             "#8491B4B2", "#91D1C2B2", "#DC0000B2")
names(vc_cols)=c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)

#function
createOncoMatrix = function(m, g = NULL, chatty = TRUE, add_missing = FALSE, cbio = FALSE){

  if(is.null(g)){
    fileConn <- file(paste0(opt$out,"/error.log"))
    writeLines("基因数<=2", fileConn)
    close(fileConn)
    
  }else{

    subMaf = subsetMaf(maf = m, genes = g, includeSyn = FALSE, mafObj = FALSE)

    if(nrow(subMaf) == 0){
      if(add_missing){
        numericMatrix = matrix(data = 0, nrow = length(g), ncol = length(levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])))
        rownames(numericMatrix) = g
        colnames(numericMatrix) = levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])

        oncoMatrix = matrix(data = "", nrow = length(g), ncol = length(levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])))
        rownames(oncoMatrix) = g
        colnames(oncoMatrix) = levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])

        vc = c("")
        names(vc) = 0

        return(list(oncoMatrix = oncoMatrix, numericMatrix = numericMatrix, vc = vc))
      }else{
        return(NULL)
      }
    }

    if(add_missing){
      subMaf[, Hugo_Symbol := factor(x = Hugo_Symbol, levels = g)]
    }


    cnv_events = c(c("Amp", "Del"), as.character(subMaf[Variant_Type == "CNV"][, .N, Variant_Classification][, Variant_Classification]))
    cnv_events = unique(cnv_events)

    if(cbio){
      vc = c("Nonstop_Mutation", "Frame_Shift_Del", "Missense_Mutation",
            "Nonsense_Mutation", "Splice_Site", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins")
      vc.cbio = c("Truncating", "Truncating", "Missense", "Truncating", "Truncating", "Truncating",
                  "In-frame", "In-frame")
      names(vc.cbio) = vc
      subMaf[,Variant_Classification_temp := vc.cbio[as.character(subMaf$Variant_Classification)]]
      subMaf$Variant_Classification_temp = ifelse(test = is.na(subMaf$Variant_Classification_temp), yes = as.character(subMaf$Variant_Classification), no = subMaf$Variant_Classification_temp)
      subMaf[,Variant_Classification := as.factor(as.character(Variant_Classification_temp))]
      subMaf[,Variant_Classification_temp := NULL]
    }

    oncomat = data.table::dcast(data = subMaf[,.(Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode)], formula = Hugo_Symbol ~ Tumor_Sample_Barcode,
                                fun.aggregate = function(x, cnv = cnv_events){
                                  #x = unique(as.character(x)) #>=2 distinct variant classification = Multi_Hit
                                  x = as.character(x) # >= 2 same/distinct variant classification = Multi_Hit See #347
                                  xad = x[x %in% cnv]
                                  xvc = x[!x %in% cnv]

                                  if(length(xvc)>0){
                                    xvc = ifelse(test = length(xvc) > 1, yes = 'Multi_Hit', no = xvc)
                                  }

                                  x = ifelse(test = length(xad) > 0, yes = paste(xad, xvc, sep = ';'), no = xvc)
                                  x = gsub(pattern = ';$', replacement = '', x = x)
                                  x = gsub(pattern = '^;', replacement = '', x = x)
                                  return(x)
                                } , value.var = 'Variant_Classification', fill = '', drop = FALSE)

    #convert to matrix
    data.table::setDF(oncomat)
    rownames(oncomat) = oncomat$Hugo_Symbol
    oncomat = as.matrix(oncomat[,-1, drop = FALSE])

    variant.classes = as.character(unique(subMaf[,Variant_Classification]))
    variant.classes = c('',variant.classes, 'Multi_Hit')
    names(variant.classes) = 0:(length(variant.classes)-1)

    #Complex variant classes will be assigned a single integer.
    vc.onc = unique(unlist(apply(oncomat, 2, unique)))
    vc.onc = vc.onc[!vc.onc %in% names(variant.classes)]
    names(vc.onc) = rep(as.character(as.numeric(names(variant.classes)[length(variant.classes)])+1), length(vc.onc))
    variant.classes2 = c(variant.classes, vc.onc)

    oncomat.copy <- oncomat
    #Make a numeric coded matrix
    for(i in 1:length(variant.classes2)){
      oncomat[oncomat == variant.classes2[i]] = names(variant.classes2)[i]
    }

    #If maf has only one gene
    if(nrow(oncomat) == 1){
      mdf  = t(matrix(as.numeric(oncomat)))
      rownames(mdf) = rownames(oncomat)
      colnames(mdf) = colnames(oncomat)
      return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
    }

    #convert from character to numeric
    mdf = as.matrix(apply(oncomat, 2, function(x) as.numeric(as.character(x))))
    rownames(mdf) = rownames(oncomat.copy)


    #If MAF file contains a single sample, simple sorting is enuf.
    if(ncol(mdf) == 1){
      sampleId = colnames(mdf)
      mdf = as.matrix(mdf[order(mdf, decreasing = TRUE),])
      colnames(mdf) = sampleId

      oncomat.copy = as.matrix(oncomat.copy[rownames(mdf),])
      colnames(oncomat.copy) = sampleId

      return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
    } else{
      #Sort by rows as well columns if >1 samples present in MAF
      #Add total variants per gene
      mdf = cbind(mdf, variants = apply(mdf, 1, function(x) {
        length(x[x != "0"])
      }))
      #Sort by total variants
      mdf = mdf[order(mdf[, ncol(mdf)], decreasing = TRUE), ]
      #colnames(mdf) = gsub(pattern = "^X", replacement = "", colnames(mdf))
      nMut = mdf[, ncol(mdf)]

      mdf = mdf[, -ncol(mdf)]

      mdf.temp.copy = mdf #temp copy of original unsorted numeric coded matrix

      mdf[mdf != 0] = 1 #replacing all non-zero integers with 1 improves sorting (& grouping)
      tmdf = t(mdf) #transposematrix
      mdf = t(tmdf[do.call(order, c(as.list(as.data.frame(tmdf)), decreasing = TRUE)), ]) #sort

      mdf.temp.copy = mdf.temp.copy[rownames(mdf),] #organise original matrix into sorted matrix
      mdf.temp.copy = mdf.temp.copy[,colnames(mdf)]
      mdf = mdf.temp.copy

      #organise original character matrix into sorted matrix
      oncomat.copy <- oncomat.copy[,colnames(mdf)]
      oncomat.copy <- oncomat.copy[rownames(mdf),]

      return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes, cnvc = cnv_events))
    }
  }}

  somaticInteractions1 = function(maf, top = 25, genes = NULL, pvalue = c(0.05, 0.01), returnAll = TRUE,
                                geneOrder = NULL, fontSize = 0.8, leftMar = 4, topMar = 4, showSigSymbols = TRUE,
                                showCounts = FALSE, countStats = 'all', countType = 'all',
                                countsFontSize = 0.8, countsFontColor = "black", colPal = "BrBG", revPal = FALSE, showSum = TRUE, plotPadj = FALSE, colNC=9, nShiftSymbols = 5, sigSymbolsSize=2,sigSymbolsFontSize=0.9, pvSymbols = c(46,42), limitColorBreaks = TRUE){
    #browser()
    if(is.null(genes)){
      genes = getGeneSummary(x = maf)[1:top, Hugo_Symbol]
    }
    
    if(length(genes) < 2){
      fileConn <- file(paste0(opt$out,"/error.log"))
      writeLines("基因数<=2", fileConn)
      close(fileConn)
      
    }else{

      om = createOncoMatrix(m = maf, g = genes)
      all.tsbs = levels(getSampleSummary(x = maf)[,Tumor_Sample_Barcode])

      mutMat = t(om$numericMatrix)
      missing.tsbs = all.tsbs[!all.tsbs %in% rownames(mutMat)]
      if(nrow(mutMat) < 2){
        fileConn <- file(paste0(opt$out,"/error.log"))
        writeLines("基因数<=2", fileConn)
        close(fileConn)
        
      }else{
      mutMat[mutMat > 0 ] = 1

      if(length(missing.tsbs) > 0){
        missing.tsbs = as.data.frame(matrix(data = 0, nrow = length(missing.tsbs), ncol = ncol(mutMat)),
                                    row.names = missing.tsbs)
        colnames(missing.tsbs) = colnames(mutMat)
        mutMat = rbind(mutMat, missing.tsbs)
      }

      #return(mutMat)

      #pairwise fisher test source code borrowed from: https://www.nature.com/articles/ncomms6901
      interactions = sapply(1:ncol(mutMat), function(i)
        sapply(1:ncol(mutMat), function(j) {
          f = try(fisher.test(mutMat[, i], mutMat[, j]), silent = TRUE)
          if (class(f) == "try-error"){
            if(all(mutMat[,i] == mutMat[,j])){
              if(colnames(mutMat)[i] != colnames(mutMat)[j]){
                warning("All the samples are in the same direction for the genes ", colnames(mutMat)[i], " and ",  colnames(mutMat)[j], "! Could not perform Fisher test.")
              }
              NA
            }else{
              if(colnames(mutMat)[i] != colnames(mutMat)[j]){
                warning("Contigency table could not created for the genes ", colnames(mutMat)[i], " and ",  colnames(mutMat)[j], "! Could not perform Fisher test.")
              }
              NA
            }
          }else{
            ifelse(f$estimate > 1,-log10(f$p.val), log10(f$p.val))
          }
        }))
      #return(interactions)
      oddsRatio <-
        oddsGenes <-
        sapply(1:ncol(mutMat), function(i)
          sapply(1:ncol(mutMat), function(j) {
            f = try(fisher.test(mutMat[, i], mutMat[, j]), silent = TRUE)
            if (class(f) == "try-error")
              if(all(mutMat[,i] == mutMat[,j])){
                NA
              }else{
                NA
              }
            else
              f$estimate
          }))
      rownames(interactions) = colnames(interactions) = rownames(oddsRatio) = colnames(oddsRatio) = colnames(mutMat)

      sigPairs = which(x = 10^-abs(interactions) < 1, arr.ind = TRUE)
      sigPairs2 = which(x = 10^-abs(interactions) >= 1, arr.ind = TRUE)

      if(nrow(sigPairs) < 1){
        fileConn <- file(paste0(opt$out,"/error.log"))
        writeLines("未找到mutation relations", fileConn)
        close(fileConn)
        
      }else{

      sigPairs = rbind(sigPairs, sigPairs2)
      sigPairsTbl = data.table::rbindlist(
                              lapply(X = seq_along(1:nrow(sigPairs)), function(i) {
                                      x = sigPairs[i,]
                                      g1 = rownames(interactions[x[1], x[2], drop = FALSE])
                                      g2 = colnames(interactions[x[1], x[2], drop = FALSE])
                                      #tbl = as.data.frame(table(apply(X = mutMat[,c(g1, g2), drop = FALSE], 1, paste, collapse = "")))
                                      tbl = as.data.frame(table(factor(apply(X = mutMat[,c(g1, g2), drop = FALSE], 1, paste, collapse = ""), levels = c("00", "01","11", "10"))))
                                      combn = data.frame(t(tbl$Freq))
                                      colnames(combn) = tbl$Var1
                                      pval = 10^-abs(interactions[x[1], x[2]])
                                      fest = oddsRatio[x[1], x[2]]
                                      d = data.table::data.table(gene1 = g1,
                                                            gene2 = g2,
                                                            pValue = pval, oddsRatio = fest)
                                      d = cbind(d, combn)
                                      d
                            }), fill = TRUE)

      sigPairsTbl[, pAdj := p.adjust(pValue, method = 'fdr')]
      sigPairsTbl[is.na(sigPairsTbl)] = 0
      sigPairsTbl$Event = ifelse(test = sigPairsTbl$oddsRatio > 1, yes = "Co_Occurence", no = "Mutually_Exclusive")
      sigPairsTbl$pair = apply(X = sigPairsTbl[,.(gene1, gene2)], MARGIN = 1, FUN = function(x) paste(sort(unique(x)), collapse = ", "))
      sigPairsTbl[,event_ratio := `01`+`10`]
      sigPairsTbl[,event_ratio := paste0(`11`, '/', event_ratio)]
      sigPairsTblSig = sigPairsTbl[order(as.numeric(pValue))][!duplicated(pair)]

      if(plotPadj){
        sigPairsTblSig$pAdjLog = ifelse(sigPairsTblSig$oddsRatio > 1, yes = -log10(sigPairsTblSig$pAdj), no = log10(sigPairsTblSig$pAdj))
        interactionsFDR = data.table::dcast(data = sigPairsTblSig, gene1 ~ gene2, value.var = 'pAdjLog')
        data.table::setDF(interactionsFDR, rownames = interactionsFDR$gene1)
        interactionsFDR$gene1 = NULL
        interactions = interactionsFDR[rownames(interactions), colnames(interactions)]
        interactions = as.matrix(interactions)
        sigPairsTblSig$pAdjLog = NULL
      }

      sigPairsTblSig = sigPairsTblSig[!gene1 == gene2] #Remove diagonal elements

      #Source code borrowed from: https://www.nature.com/articles/ncomms6901
      if(nrow(interactions) >= 5){
        #interactions[10^-abs(interactions) > max(pvalue)] = 0
        diag(interactions) <- 0
        m <- nrow(interactions)
        n <- ncol(interactions)


        col_pal = RColorBrewer::brewer.pal(9, colPal)
        if(revPal){
          col_pal = rev(col_pal)
        }
        col_pal = grDevices::colorRampPalette(colors = col_pal)
        col_pal = col_pal(m*n-1)


        if(!is.null(geneOrder)){
          if(!all(rownames(interactions) %in% geneOrder)){
            fileConn <- file(paste0(opt$out,"/error.log"))
            writeLines("基因顺序有误", fileConn)
            close(fileConn)
            stop()
            
          }
          interactions = interactions[geneOrder, geneOrder]
        }

        interactions[lower.tri(x = interactions, diag = TRUE)] = NA

        gene_sum = getGeneSummary(x = maf)[Hugo_Symbol %in% rownames(interactions), .(Hugo_Symbol, AlteredSamples)]
        data.table::setDF(gene_sum, rownames = as.character(gene_sum$Hugo_Symbol))
        gene_sum = gene_sum[rownames(interactions),]
        if(!all(rownames(gene_sum) == rownames(interactions))){
          fileConn <- file(paste0(opt$out,"/error.log"))
          writeLines("基因行匹配错误", fileConn)
          close(fileConn)
          stop()
        }
        if(!all(rownames(gene_sum) == colnames(interactions))){
          fileConn <- file(paste0(opt$out,"/error.log"))
          writeLines("基因列匹配错误", fileConn)
          close(fileConn)
          stop()
        }
        if(showSum){
          rownames(gene_sum) = paste0(apply(gene_sum, 1, paste, collapse = ' ['), ']')
        }

        par(bty="n", mar = c(1, leftMar, topMar, 2)+.1, las=2, fig = c(0, 1, 0, 1))

        # adjust breaks for colors according to predefined legend values
        breaks = NA
        if(limitColorBreaks){
          minLog10pval = 3
          breaks <- seq(-minLog10pval,minLog10pval,length.out=m*n+1)
          #replace extreme values with the predefined minLog10pval values (and avoid white colored squares)
          interactions4plot  = interactions
          interactions4plot[interactions4plot < (-minLog10pval)] = -minLog10pval
          interactions4plot[interactions4plot > minLog10pval] = minLog10pval
          interactions = interactions4plot
        }

        image(x=1:n, y=1:m, interactions, col = col_pal,
              xaxt="n", yaxt="n",
              xlab="",ylab="", xlim=c(0, n+1), ylim=c(0, n+1),
              breaks = seq(-3, 3, length.out = (nrow(interactions) * ncol(interactions))))

        abline(h=0:n+.5, col="white", lwd=.5)
        abline(v=0:n+.5, col="white", lwd=.5)

        mtext(side = 2, at = 1:m, text = rownames(gene_sum), cex = fontSize, font = 3)
        mtext(side = 3, at = 1:n, text = rownames(gene_sum), cex = fontSize, font = 3)
        #text(x = 1:m, y = rep(n+0.5, length(n)), labels = rownames(gene_sum), srt = 90, adj = 0, font = 3, cex = fontSize)

        if(showCounts){
          countStats = match.arg(arg = countStats, choices = c("all", "sig"))
          countType = match.arg(arg = countType, choices = c("all", "cooccur", "mutexcl"))

          if(countStats == 'sig'){
            w = arrayInd(which(10^-abs(interactions) < max(pvalue)), rep(m,2))
            for(i in 1:nrow(w)){
              g1 = rownames(interactions)[w[i, 1]]
              g2 = colnames(interactions)[w[i, 2]]
              g12 = paste(sort(c(g1, g2)), collapse = ', ')
              if(countType == 'all'){
                e = sigPairsTblSig[pValue < max(pvalue)][pair %in% g12, event_ratio]
              }else if(countType == 'cooccur'){
                e = sigPairsTblSig[pValue < max(pvalue)][Event %in% "Co_Occurence"][pair %in% g12, `11`]
              }else if(countType == 'mutexcl'){
                e = sigPairsTblSig[pValue < max(pvalue)][Event %in% "Mutually_Exclusive"][pair %in% g12, `11`]
              }
              if(length(e) == 0){
                e = 0
              }
              text(w[i,1], w[i,2], labels = e, font = 3, col = countsFontColor, cex = countsFontSize)
            }
          }else if(countStats == 'all'){
            w = arrayInd(which(10^-abs(interactions) < max(pvalue)), rep(m,2))
            w2 = arrayInd(which(10^-abs(interactions) >= max(pvalue)), rep(m,2))
            w = rbind(w, w2)
            #print(w)
            for(i in 1:nrow(w)){
              g1 = rownames(interactions)[w[i, 1]]
              g2 = colnames(interactions)[w[i, 2]]
              g12 = paste(sort(c(g1, g2)), collapse = ', ')
              if(countType == 'all'){
                e = sigPairsTblSig[pair %in% g12, event_ratio]
              }else if(countType == 'cooccur'){
                e = sigPairsTblSig[pair %in% g12, `11`]
              }else if(countType == 'mutexcl'){
                e = sigPairsTblSig[pair %in% g12, `01` + `10`]
              }
              if(length(e) == 0){
                e = 0
              }
              text(w[i,1], w[i,2], labels = e, font = 3, col = countsFontColor, cex = countsFontSize)
            }
          }
        }

        if(showSigSymbols){
          w = arrayInd(which(10^-abs(interactions) < min(pvalue)), rep(m,2))
          points(w, pch=pvSymbols[2], col="black", cex = sigSymbolsSize)
          #w = arrayInd(which(10^-abs(interactions) < max(pvalue)), rep(m,2))
          w = arrayInd(which((10^-abs(interactions) < max(pvalue)) & (10^-abs(interactions) > min(pvalue))), rep(m,2))
          points(w, pch=pvSymbols[1], col="black", cex = sigSymbolsSize)
        }

        if(showSigSymbols){
          points(x = n-nShiftSymbols, y = 0.7*n, pch = pvSymbols[2], cex = sigSymbolsSize) # "*"
          if(plotPadj){
            text(x = n-nShiftSymbols, y = 0.7*n, paste0(" fdr < ", min(pvalue)), pos=4, cex = sigSymbolsFontSize, adj = 0)
          }else{
            text(x = n-nShiftSymbols, y = 0.7*n, paste0(" P < ", min(pvalue)), pos=4, cex = sigSymbolsFontSize, adj = 0)
          }

          points(x = n-nShiftSymbols, y = 0.65*n, pch = pvSymbols[1], cex = sigSymbolsSize) # "."
          if(plotPadj){
            text(x = n-nShiftSymbols, y = 0.65*n, paste0(" fdr < ", max(pvalue)), pos=4, cex = sigSymbolsFontSize)
          }else{
            text(x = n-nShiftSymbols, y = 0.65*n, paste0(" P < ", max(pvalue)), pos=4, cex = sigSymbolsFontSize)
          }
        }

        #image(y = 1:8 +6, x=rep(n,2)+c(2,2.5)+1, z=matrix(c(1:8), nrow=1), col=brewer.pal(8,"PiYG"), add=TRUE)
        par(fig = c(0.4, 0.7, 0, 0.4), new = TRUE)
        image(
          x = c(0.8, 1),
          y = seq(0, 1, length.out = 200),
          z = matrix(seq(0,1,length.out = 200), nrow = 1),
          col = col_pal, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = NA, ylab = NA
        )

        #atLims = seq(nrow(interactions), 0.9*nrow(interactions), length.out = 7)
        atLims = seq(0, 1, length.out = 7)
        axis(side = 4, at = atLims,  tcl=-.15, labels =c("> 3 (Mutually exclusive)", 2, 1, 0, 1, 2, ">3 (Co-occurence)"), lwd=.5, cex.axis = sigSymbolsFontSize, line = 0.2)
        if(plotPadj){
          text(x = 0.4, y = 0.5, labels = "-log10(fdr)", srt = 90, cex = sigSymbolsFontSize, xpd = TRUE)
        }else{
          text(x = 0.4, y = 0.5, labels = "-log10(P-value)", srt = 90, cex = sigSymbolsFontSize, xpd = TRUE)
        }

        #mtext(side=4, at = median(atLims), "-log10 (p-value)", las=3, cex = 0.9, line = 2.5, font = 1)
      }

      if(!returnAll){
        sigPairsTblSig = sigPairsTblSig[pValue < min(pvalue)]
      }
      pdf(paste0(opt$out, "/mutation_relation.pdf"))
      sigPairsTblSig
      dev.off()
      png(paste0(opt$out, "/mutation_relation.png"))
      sigPairsTblSig
      dev.off()
    }}}}




laml<-read.maf(opt$maf)
maf<-laml@data

mr <- somaticInteractions1(maf=laml,top=25,pvalue=c(0.05,1),nShiftSymbols=3)
mr

