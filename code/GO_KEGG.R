## GO Enrichment analysis
library(clusterProfiler)
library(patchwork)
DEGfile <- read.csv('.../Pipeline/figure_correspond_data/Figure3 CD/my_degs.csv',header = T,row.names = NULL)
genes <- as.vector(DEGfile$gene_id)
go_anno <- read.csv(".../Pipeline/figure_correspond_data/Figure3 CD/host_go_annotation.txt",header = T,sep = "\t",row.names = NULL)
go2gene <- go_anno[, c(2, 1)]
go2name <- go_anno[, c(2, 3)]
ego <- enricher(genes, TERM2GENE = go2gene, TERM2NAME = go2name)
write.csv(ego@result,".../Pipeline/figure_correspond_data/Figure3 CD/KP_rbp_enrich.csv",row.names =FALSE)
#write.csv(go_anno, ".../go.txt",sep='\t',row.names =FALSE)
b <- dotplot(ego, showCategory = 15 ,title = "dotplot for enricher")
pdf(".../barplot.pdf",width = 12,height = 8)
b
dev.off()

##GO barplot
GOO <- read.csv(".../Pipeline/raw_Data/富集分析/phage_rbp_enrich.csv")
#GOO <- read.csv(".../Pipeline/figure_correspondhost_data/Figure3 CD/phage_enrich.csv")
GOO$gene <- gsub("^.*/","",GOO$GeneRatio)
GOO$gene <- as.numeric(GOO$gene)
GOO$Percentage <- GOO$Count/GOO$gene
GO_table <- GOO[,c(1,2,10,9,12,3,4,8,5,6,7)]
GO_col_name <- c("GO_ID","GO_term","GO_category","Num_of_symbols_in_list_in_GO","Percentage_of_symbols_in_list","GeneRatio","BgRatio","Symbols_in_list","pvalue","p.adjust","qvalue")
names(GO_table) <- GO_col_name
GO_table <- arrange(GO_table, GO_category, pvalue)
#write.csv(GO_table,file="table_human_6s_case_vs_con_FC1_GO_enrichment.csv",row.names = F)

#Sort by pathway
GO_BP <- GO_table[GO_table$GO_category=="biological_process",]
GO_CC <- GO_table[GO_table$GO_category=="cellular_component",]
GO_MF <- GO_table[GO_table$GO_category=="molecular_function",]
GO_BP_top10 <- GO_BP[1:10,]
GO_CC_top10 <- GO_CC[1:7,]
GO_MF_top10 <- GO_MF[1:10,]
GO_BP_top10 <- GO_BP[GO_BP$p.adjust < 0.05,]
GO_CC_top10 <- GO_CC[GO_CC$p.adjust < 0.05,]
GO_MF_top10 <- GO_MF[GO_MF$p.adjust < 0.05,]
go_draw <- rbind(GO_BP_top10,GO_CC_top10,GO_MF_top10)
go_draw$p <- -log(go_draw$p.adjust, base = 10)

go_draw$GO_term <- factor(go_draw$GO_term,levels = go_draw$GO_term)
ridgeplot(ego)
upsetplot(ego)
library(ggplot2)
###GO barplot
p2 <- ggplot(go_draw,aes(x=GO_term,y=p,fill=GO_category, order=GO_category))+
  geom_bar(stat = 'identity')+coord_flip()+
  scale_fill_manual(values = c('#FF6666','#33CC33','#3399FF'),limits = c('molecular_function','cellular_component','biological_process'),labels = c('molecular_function','cellular_component','biological_process'))+
  theme_bw()+theme(axis.text.x = element_text(size = 15), axis.text.y  = element_text(size = 15))+
  ylab("-log10pvalue") + theme(axis.title.y  = element_text(size = 15))+
  xlab("GO_Term") + theme(axis.title.x = element_text(size = 15)) + 
  labs(title='The Most Enriched GO Terms') + theme(plot.title = element_text(size = 15))  + 
  theme_bw() + theme_classic() + theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
p2

p2 <- p2 + theme(axis.text.x=element_text(size = 12), axis.text.y=element_text(size = 12))
p2
#x3 <- x3 + theme(axis.line = element_line(size = 1.2), axis.ticks = element_line(size=1.2), axis.ticks.length = unit(0.1,'cm'))
#x3 <- x3 + theme(legend.text = element_text(size=17), legend.title = element_text(size=17))
pdf(".../barplot_KP_GO.pdf",width = 12,height = 8)
p2
dev.off()


go_draw$GO_category <- factor(go_draw$GO_category, levels = c( "molecular_function","cellular_component","biological_process"))
##GO dotplot
go_dotplot=ggplot(go_draw,aes(p, GO_term))
d1 <- go_dotplot+geom_point(aes(size=Num_of_symbols_in_list_in_GO,shape = GO_category,color=p))
d2 <- d1+labs(color=expression(-log10pvalue),size="gene_number",x="-log10pvalue",y="GO_Term")
d3 <- d2 + labs(title='GO enrichment of test Genes') + theme(plot.title = element_text(size = 15))
d4 <- d3 + scale_color_gradient(high = "#FF00AAFF", low = "green")
d5 <- d4 + theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
d6 <- d5 + theme_bw() + theme(axis.text.x = element_text(size = 12),
                              axis.text.y  = element_text(size = 12))
d6
#pdf(".../dotplot_KP_GO.pdf",width = 12,height = 8)
pdf(".../dotplot_phage_rbp_GO.pdf",width = 12,height = 8)
d6
dev.off()

##KEGG Enrichment analysis
library(clusterProfiler)
DEGfile <- read.csv('.../Pipeline/figure_correspond_data/Figure3 CD/my_degs_kegg.csv',header = T,row.names = NULL)
genes <- as.vector(DEGfile$gene_id)
go_anno <- read.csv(".../Pipeline/raw_Data/Enrichment_analysis/host_kegg_annotation.txt",header = T,sep = "\t",row.names = NULL)
go2gene <- go_anno[, c(2, 1)]
go2name <- go_anno[, c(2, 3)]
ego <- enricher(genes, TERM2GENE = go2gene, TERM2NAME = go2name)
write.csv(ego@result,".../KP_enrich_kegg.csv",row.names =FALSE)
ego@result

library(ggplot2)
KEGG <- read.csv(".../Pipeline/raw_Data/Enrichment_analysis/phage_rbp_enrich_kegg.csv")
KEGG$gene <- gsub("^.*/","",KEGG$GeneRatio)
KEGG$gene <- as.numeric(KEGG$gene)
KEGG$Percentage <- KEGG$Count/KEGG$gene
KEGG_table <- KEGG[,c(1,2,9,11,3,4,8,5,6,7)]
KEGG_col_name <- c("ID","Path_term","Num_of_symbols_in_list_in_KEGG","Percentage_of_symbols_in_list","GeneRatio","BgRatio","Symbols_in_list","pvalue","p.adjust","qvalue")
names(KEGG_table) <- KEGG_col_name
KEGG_table <- arrange(KEGG_table, pvalue)
kegg_draw <- KEGG_table

## filter
KEGG_part <- KEGG_table[KEGG_table$p.adjust<0.05,]
if (nrow(KEGG_part) < 30){kegg_draw <- KEGG_part} else {kegg_draw <- KEGG_part[1:30,]}
kegg_draw$p <- -log(kegg_draw$pvalue,base = 10)
##KEGG barplot 
order <- sort(kegg_draw$p,index.return=TRUE,decreasing = TRUE)
kegg_draw$Path_term <- factor(kegg_draw$Path_term, levels = kegg_draw$Path_term[order$ix])
x2 <- ggplot(data=kegg_draw,aes(x=Path_term,y=p))
x2 <- x2 + geom_bar(stat = "identity",fill="#ff7575")+coord_flip() + ylab("-log10pvalue") + theme(axis.title.y  = element_text(size = 20))
x2 <- x2 + theme_classic() + theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
x2 <- x2 + labs(title='The Most Enriched KEGG Terms') + theme(plot.title = element_text(size = 15))
x2 <- x2 + theme(axis.text.x=element_text(size = 12), axis.text.y=element_text(size = 12))
#x2 <- x2 + theme(axis.line = element_line(size = 1.2), axis.ticks = element_line(size=1.2), axis.ticks.length = unit(0.1,'cm'))
#x2 <- x2 + theme(legend.text = element_text(size=17), legend.title = element_text(size=17))
x2
pdf(".../barplot_KP_KEGG.pdf",width = 12,height = 8)
x2
dev.off()
##KEGG dotplot
kegg_dotplot=ggplot(kegg_draw,aes(p, Path_term))
y1 <- kegg_dotplot+geom_point(aes(size=Num_of_symbols_in_list_in_KEGG,color=p))
y2 <- y1 + labs(color=expression(-log10pvalue),size="gene_number",x="-log10pvalue",y="pathway")
y3 <- y2 + labs(title='KEGG enrichment of test Genes') + theme(plot.title = element_text(size = 15))
y4 <- y3 + scale_color_gradient(high = "#FF00AAFF", low = "green")
y5 <- y4 + theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
y6 <- y5 + theme_bw()
y6 <- y6 + theme(axis.text.x=element_text(size = 12), axis.text.y=element_text(size = 12))
y6
#y5 <- y4 + theme(axis.title = element_text(size = 22), axis.text.x=element_text(size = 20), axis.text.y=element_text(size = 19))
#y6 <- y5 + theme(axis.line = element_line(size = 1.2), axis.ticks = element_line(size=1.2), axis.ticks.length = unit(0.1,'cm'))
#y7 <- y6 + theme(legend.text = element_text(size=17), legend.title = element_text(size=17))
pdf(".../dotplot_KP_rbp_KEGG.pdf",width = 12,height = 8)
y6
dev.off()


# GSEA
library(data.table)
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(enrichplot)
library(org.Hs.eg.db)
all <- read.gmt('.../rbp_host_gsea.gmt')
genelist <- read.csv(".../KP_key_gene_importance.csv")
all <- read.gmt('.../Pipeline/raw_Data/Enrichment_analysis/rbp_host_gsea_news.gmt')
genelist <- read.csv(".../Pipeline/raw_Data/Enrichment_analysis/host_key_gene_importance.csv")
genefc<-genelist$importance
names(genefc)<-genelist$gene
genefc<-sort(genefc,decreasing = T)
#gene <- genelist[,2]
#rownames(genelist) <- genelist[,1]
#genelist <- genelist[,-1]
gsea <- GSEA(genefc, TERM2GENE = all, pvalueCutoff = 1, pAdjustMethod = 'BH', scoreType = 'pos')
x <- gseaplot2(gsea, geneSetID = c('phage_rbp', 'phage_GO(defense response to bacterium)', 'phage_GO(DNA recombination)', 'phage_GO(DNA integration)', 'phage_GO(DNA replication)', 'phage_KEGG(DNA replication [PATH:ko03430])', 'phage_KEGG(DNA replication proteins [BR:ko03032])', 'phage_KEGG(DNA repair and recombination proteins [BR:ko03400])'), pvalue_table = TRUE) + labs(title='phage gene GSEA', hjust = 0.5) + theme(plot.title = element_text(hjust = 0.5, vjust = 100, size = 16))
x <- gseaplot2(gsea, geneSetID = c('host_rbp', 'host_GO(single-stranded DNA binding)', 'host_GO(host cell cytoplasm)', 'host_GO(DNA integration)', 'host_GO(SOS response)', 'host_KEGG(Glyoxylate and dicarboxylate metabolism [PATH:ko00630])', 'host_KEGG(Other glycan degradation [PATH:ko00511]', 'host_KEGG(Fructose and mannose metabolism [PATH:ko00051])'), pvalue_table = TRUE) + labs(title='host gene GSEA', hjust = 0.5) + theme(plot.title = element_text(hjust = 0.5, vjust = 100, size = 16))
x
pdf(".../gsea_phage_gene.pdf",width = 8.27,height = 11.69)
x
dev.off()
gseaplot(gsea, geneSetID = 'gene', pvalue_table = TRUE)+ labs(title='host gene GSEA', hjust = 0.5) + theme(plot.title = element_text(hjust = 0.5, vjust = 100, size = 16))
ids <- gsea@result$ID[3]
gseadist(gsea, IDs = ids, type = 'density')
gsearank(gsea, geneSetID = 3)
