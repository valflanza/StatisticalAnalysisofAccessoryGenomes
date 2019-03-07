library(tidyverse)
library(mclust)
library(cluster)

library(ape) #optional
library(Rtsne)

setwd("~/Dropbox/MiMB/")





#genomes information

assemblies = read_delim("DownLoadGenome.list",
                        delim = "\t",
                        col_names = FALSE)

#   See ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt for a description of the columns in this file.
# assembly_accession	bioproject	biosample	wgs_master	refseq_category	taxid	species_taxid	organism_name	infraspecific_name	isolate	version_status	assembly_level	release_type	genome_rep	seq_rel_date	asm_name	submitter	gbrs_paired_asm	paired_asm_comp	ftp_path	excluded_from_refseq	relation_to_type_material
colnames(assemblies) = c(
  "assembly_accession",
  "bioproject",
  "biosample",
  "wgs_master",
  "refseq_category",
  "taxid",
  "species_taxid",
  "organism_name",
  "infraspecific_name",
  "isolate",
  "version_status",
  "assembly_level",
  "release_type",
  "genome_rep",
  "seq_rel_date",
  "asm_name",
  "submitter",
  "gbrs_paired_asm",
  "paired_asm_comp",
  "ftp_path",
  "excluded_from_refseq",
  "relation_to_type_material"
)

assemblies$infraspecific_name = gsub("strain=", "", assemblies$infraspecific_name)

assemblies = assemblies %>% unite(Name,
                                  infraspecific_name,
                                  isolate,
                                  sep = "",
                                  remove = FALSE)
tmp = assemblies %>% separate(ftp_path,
                              c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10"),
                              sep = "/") %>% select(assembly_accession, Genome = V10)
assemblies = inner_join(tmp, assemblies)
assemblies = assemblies %>% separate(assembly_accession,
                                     c("assembly_abrev", "kk"),
                                     sep = "\\.",
                                     remove = FALSE) %>% select(-kk)

#### ACCESSORY

accnet = read_delim("Network.csv", delim = "\t")
accnet.annot = read_delim("Table.csv", delim = "\t")


accnet.matrix = accnet %>% group_by(Source) %>% mutate(degree = n()) %>% ungroup() %>%
  filter(degree > 1) %>% group_by(Source, Target) %>% summarise(value = 1) %>%
  spread(Source, value, fill = 0) %>% remove_rownames() %>% column_to_rownames("Target")

accnet.dist = dist(accnet.matrix, method = "binary")
accnet.hclust = hclust(accnet.dist, method = "average")
factoextra::fviz_dend(accnet.hclust, k=10 ,horiz = TRUE) ###Optional

accnet.tree.nj = ape::nj(accnet.dist)  ### optional
plot(accnet.tree.nj) ### optional


accnet.cluster = Mclust(accnet.dist)


accnet.tsne = Rtsne::Rtsne(accnet.dist) ### optional
plot(accnet.tsne$Y, col = rainbow(max(accnet.cluster$classification))[accnet.cluster$classification])### optional


accnet.protein.matrix = accnet %>% group_by(Source) %>% mutate(degree = n()) %>% filter(degree > 1) %>%
  ungroup() %>%  group_by(Source, Target) %>% summarise(value = 1) %>% spread(Target, value, fill = 0) %>%
  remove_rownames() %>% column_to_rownames("Source")
accnet.protein.dist = dist(accnet.protein.matrix, method = "binary")
accnet.protein.hclust = hclust(accnet.protein.dist)
accnet.protein.cluster = cutree(accnet.protein.hclust, h = 0.9)



### Cluster Statistics Analysis


tmp = as_tibble(accnet.cluster$classification)
tmp = rownames_to_column(tmp)
colnames(tmp) = c("Target", "MclustCluster")
accnet.full = full_join(accnet %>% distinct(), tmp)


accnet.hclust.q95 = cutree(accnet.hclust, k = 10)
tmp = as_tibble(accnet.hclust.q95)
tmp = rownames_to_column(tmp)
colnames(tmp) = c("Target", "HclustQ95Cluster")
accnet.full = full_join(accnet.full, tmp)

accnet.full = accnet.full %>% left_join(assemblies %>% select(Target = assembly_abrev, Name))

### Select one option
accnet.full$Cluster = accnet.full$HclustQ95Cluster
#accnet.full$Cluster = accnet.full$MclustCluster
###


tmp = accnet.full %>% ungroup() %>% select(Target) %>% distinct() %>% count()
accnet.full$AccnetGenomeSize = tmp$n

tmp = accnet.full %>% ungroup() %>% select(Source) %>% distinct() %>% count()
accnet.full$AccnetProteinSize = tmp$n


accnet.full = accnet.full %>%
  select(Source, Target) %>% distinct() %>% group_by(Source) %>% mutate(TotalFreq = n()) %>%
  ungroup() %>% group_by(Target) %>% mutate(Degree = n()) %>% full_join(accnet.full)

accnet.full = accnet.full %>% select(Cluster, Target) %>% distinct() %>%
  group_by(Cluster) %>% summarise(ClusterGenomeSize = n()) %>% full_join(accnet.full)

accnet.full = accnet.full %>% select(Cluster, Source) %>% distinct() %>%
  group_by(Cluster) %>% summarise(ClusterProteinSize = n()) %>% full_join(accnet.full)

accnet.full = accnet.full %>% group_by(Cluster, Source) %>% mutate(ClusterFreq = n())

accnet.full = accnet.full %>%
  mutate(
    perClusterFreq = ClusterFreq / ClusterGenomeSize,
    perTotalFreq = TotalFreq / AccnetGenomeSize
  )

accnet.full = accnet.full %>%
  select(
    Cluster,
    Source,
    TotalFreq,
    ClusterFreq,
    AccnetGenomeSize,
    ClusterGenomeSize,
    perClusterFreq,
    perTotalFreq
  ) %>%
  distinct() %>%
  mutate(
    OdsRatio =  perClusterFreq / perTotalFreq ,
    pvalue = phyper(
      ClusterFreq-1,
      TotalFreq,
      AccnetGenomeSize - TotalFreq,
      ClusterGenomeSize,
      lower.tail = FALSE
    )
  ) %>%
  full_join(accnet.full)

accnet.full = accnet.full %>% select(Cluster, Source, pvalue) %>% 
  distinct() %>% group_by(Cluster) %>% mutate(padj = p.adjust(pvalue, method = "BY")) %>% 
  full_join(accnet.full)

accnet.full = accnet.full %>% select(-Type,-Weight) %>% select(
  Source,
  Target,
  Degree,
  Cluster,
  perClusterFreq,
  ClusterFreq,
  ClusterGenomeSize,
  perTotalFreq,
  TotalFreq,
  OdsRatio,
  pvalue,
  padj,
  AccnetGenomeSize,
  AccnetProteinSize
)
accnet.full = accnet.full %>% left_join(assemblies %>% select(Target = assembly_abrev, Name))


### Building Annotation table for gephi representation.


Annotation.table = accnet %>% select(ID = Source) %>% distinct() %>% mutate(Type = "Cluster", Polygon = 1) %>%
  bind_rows(accnet %>% select(ID = Target) %>% distinct() %>% mutate(Type = "GU", Polygon = 4))

Annotation.table = assemblies %>% select(ID = assembly_abrev, Name) %>% full_join(Annotation.table)

tmp = as.data.frame(accnet.cluster$classification) %>% rownames_to_column("ID") %>% as_tibble()
colnames(tmp) = c("ID", "GenomeCluster")
Annotation.table = Annotation.table %>% left_join(tmp)

tmp = as.data.frame(accnet.protein.cluster) %>% rownames_to_column("ID") %>% as_tibble()
colnames(tmp) = c("ID", "ProteinCluster")
Annotation.table = Annotation.table %>% left_join(tmp)

Annotation.table = Annotation.table %>% replace_na(list(GenomeCluster =0, ProteinCluster = 0, Name =""))
Annotation.table = Annotation.table %>% left_join(accnet.annot %>% select(-Type))

write.table(
  Annotation.table,
  "Accnet.Annotation.csv",
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

### Export Network with p-value as edge-weigth

StatisticalResults = accnet.full %>% filter(padj < 0.05, perClusterFreq > 0.80) %>% 
  select(-Target,-Name,-Degree) %>% distinct() %>% left_join(accnet.annot %>% rename(Source = ID)) %>% 
  write_delim("StatisticsResults.csv",  delim = "\t")

accnet.full %>% ungroup() %>% select(Target, Source, padj, perClusterFreq) %>% 
  mutate(Weight = 2 -padj) %>% mutate(Type = "Undirected") %>% write_delim("Network.full.csv",  delim = "\t")




### Statistical analysis of specific group

mlst = read_delim("all.mlst", delim = "\t")

accnet.ST131 = accnet.full %>% left_join(mlst %>% select(Target,ST)) %>% filter(ST == 131)
tmp = accnet.ST131 %>% ungroup() %>% select(Target) %>% distinct() %>% count()

accnet.ST131$ST131GenomeSize = as.numeric(tmp)

accnet.ST131 = accnet.ST131 %>% ungroup() %>% group_by(Source) %>% mutate(ST131Freq = n()) %>% ungroup() %>%
  mutate(perST131Freq = ST131Freq / ST131GenomeSize)

accnet.ST131 = accnet.ST131 %>%
  mutate(ST131vPopulation_pvalue = phyper(ST131Freq-1,TotalFreq,AccnetGenomeSize- TotalFreq,ST131GenomeSize, lower.tail = FALSE))

accnet.ST131 = accnet.ST131 %>% 
  select(Cluster,Source,ST131vPopulation_pvalue) %>% distinct() %>% group_by(Cluster) %>%
  mutate(padj_v_Population = p.adjust(ST131vPopulation_pvalue,method = "bonferroni")) %>% full_join(accnet.ST131)
  
accnet.ST131.results = accnet.ST131 %>% select(-Target,-Degree,-Name) %>% distinct() %>%
  filter(perST131Freq >= 0.9, padj_v_Population < 1e-3) %>% left_join(accnet.annot %>% select(Source = ID,Description))
  


