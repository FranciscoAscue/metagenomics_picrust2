######################################################################
######################## GENOMICS WITH R  ############################
######################################################################

############### Install and load library's ###############
##########################################################

source("install.R", local = TRUE) ## cargamos o instalamos los paquetes necesarios
## Puedes revisar en la lista dependencies dentro install.R

# Source scripts

source("scripts/listfastq.R", local = TRUE) ## listado de files fastq 1,2 para realizar calidad
source("scripts/filtereads.R", local = TRUE) ## filtrado de lecturas creando otros archivos filtrados

############### Quality of fastq ###############
################################################

lecturas <- list_fastq(pattern = c("1_001.fastq.gz","2_001.fastq.gz"))
plotQualityProfile(c(lecturas$lf[3],lecturas$lr[3]))

############### Filter reads from fastq files ###############
#############################################################

log_filter <- filter_reads(name = lecturas$name, lf = lecturas$lf, 
                           lr = lecturas$lr, trunc = 250)

filtF <- file.path("data/processed_data/filtered_F", paste0(lecturas$name, "_filt_1.fastq.gz"))
filtR <- file.path("data/processed_data/filtered_R", paste0(lecturas$name, "_filt_2.fastq.gz"))

############### Assembly genomes ###############
################################################

if(TRUE){
  errR <-learnErrors(filtR, multithread = T)
  errF <-learnErrors(filtF, multithread = T)
  
  #grafica del modelo de error
  plotErrors(errF,nominalQ = T)
  plotErrors(errR,nominalQ = T)
}

#----inferencia de asv----

dadaF <- dada(filtF, err = errF, multithread = T)
dadaR <- dada(filtR, err = errR, multithread = T)

#ASV es una variante de secuencia de amplicones (amplicon sequence variant)
#denoising : quitar los duplicados o ruido

#----merge o pariar o fusion de secuencias----
pareadas <- mergePairs(dadaF, filtF, dadaR, filtR, verbose = T)
head(pareadas[[1]])
#----construir tabla de secuencias----

seqTab <- makeSequenceTable(pareadas)
table(nchar(getSequences(seqTab)))

#----eliminar quimeras----
seqtab_nochim <- removeBimeraDenovo(seqTab, method = "consensus"
                                    , multithread = T, 
                                    verbose = T)
dim(seqtab_nochim)
table(nchar(getSequences(seqTab)))

# Proproción de quimeras
sum(seqtab_nochim/sum(seqTab))
write.csv2(seqtab_nochim,paste0(getwd(),"/results/seqtab_nochim.csv"))

#----Denoising stats----
getN <-function(x) sum(getUniques(x))
stats <- cbind(log_filter$OUT, sapply(dadaF, getN), 
               sapply(dadaR, getN), 
               sapply(pareadas, getN),
               rowSums(seqtab_nochim))

colnames(stats) <- c("input", "filtered", "denoisedF", 
                     "denoisedR", "merged", "nonchim")
rownames(stats) <- lecturas$name
write.csv(stats,"results/stats.csv")


#---------------------ETAPA 3----------------------
#----Asignación taxonomica----
ruta_clasificador <- "data/reference/silva_nr99_v138_train_set.fa.gz"
taxa <- assignTaxonomy(seqtab_nochim, ruta_clasificador , multithread = TRUE)
taxa_print <- taxa
rownames(taxa_print) <- NULL
head(taxa_print)
dim(taxa_print)

muestras <- rownames(seqtab_nochim)
length(muestras)

#----uso de metadata----
ruta_metadata <- ".../2_Data_New"
getwd()
met <- read.csv("data/Metadata.csv", row.names = 1)
met <- met[order(lecturas$name),]
#----apoyo phyloseq:generar tablas----
phyloseq_ob <- phyloseq(otu_table(seqtab_nochim, taxa_are_rows = F),
                        sample_data(met),
                        tax_table(taxa))
view(as.data.frame(otu_table(seqtab_nochim, taxa_are_rows = F)))

#----extraer nombres de los ASV----
dna <- Biostrings::DNAStringSet(taxa_names(phyloseq_ob))
names(dna) <- taxa_names(phyloseq_ob)
ps <- merge_phyloseq(phyloseq_ob, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

#----exportar secuencias representativas----
#exportar archivo fasta: futuros arboles filogeneticos
names(dna) <- taxa_names(ps)
writeXStringSet(dna, "results/rep-seq.fna")

#exportar tabla de taxonomia
TAX <- as(tax_table(ps), "matrix")
view(TAX)
write.csv(TAX,"results/taxa-gut16.csv")

#Exportar tablas de ASV"s
ASV <- as(otu_table(ps),"matrix")
ASV <- t(ASV)
colnames(ASV) <- gsub("_filt_1.fastq.gz","",colnames(ASV), fixed = TRUE)
view(ASV)
write.csv(ASV, "results/asv-table.csv")

#----Exportar resultados de DADA2----
save(errR, errF, seqtab_nochim, taxa, file = "results/resultados_Dada2.RDATA")

#---------------------ETAPA 4----------------------


#----Instalacion de paqueterias 3----
BiocManager::install("microbiome")
install.packages("devtools")
devtools::install_github("microsud/microbiomeutilities")
library(ape)

#----construccion de objeto phyloseq-----
#1: importar tablas de ASV
asv_table <- read.csv("results/asv-table.csv", row.names = 1, header = T)
class(asv_table)
#2: crear o convertir a tabla tipo phyloseq
ASV <- otu_table(asv_table, taxa_are_rows = T)
class(ASV)
ASV
#3: importar tabla de taxonomia
taxonomy <- read.csv(file = "results/taxa-gut16.csv", header = T, row.names = 1)
#Creando un objeto con los Taxa
TAX <- tax_table(as.matrix(taxonomy))
#Creando un objeto physeq de prueba
physeq <- phyloseq(ASV, TAX)
physeq

# Crear “arbol filogenético”
random_tree = rtree(ntaxa(physeq), rooted=FALSE, tip.label=taxa_names(physeq))

#Importar tabla de metadatos
metadata <- read.csv("data/Metadata.csv", header = T, row.names = 1)

#Convertir tabla a tipo phyloseq
META <- sample_data(metadata)

#construir phyloseq final
ps <- merge_phyloseq(physeq, META, random_tree)
ps

#exportar phyloseq
saveRDS(ps, "results/objeto_phyloseq.RDS")

#----Instalacion de paqueterias 4----
library(microbiome)
library(microbiomeutilities)
library(DT)
install.packages("DT")

#----importar objeto phyloseq(ps)----
ps <- readRDS("results/objeto_phyloseq.RDS")
ps

#----analisis de rarefacción----
#1.conocer el mín de ASV x muestra 
summary(sample_sums(ps))
#2.Conocer el N° de secuencias por muestra
sort(sample_sums(ps))
#3.ejecutar rarefacción
set.seed(1)
ps_rar <- rarefy_even_depth(ps, sample.size = 15175)
#4.comparar el objeto phyloseq original con el rarefactado
ps
ps_rar
#5. guardar el objeto phyloseq rarefactado
saveRDS(ps_rar,"results/phyloseq_rar")

class(ps_rar)


#---------------------ETAPA 5----------------------
#-----muestras conservadas por variable----
sample_data(ps)
#table(as.data.frame(sample_data(ps)))
count(as.data.frame(sample_data(ps_rar)),Procedencia)
count(as.data.frame(sample_data(ps)),Procedencia)
#----visualización----
#----a.hacer submuestreos para ecologia----
submuestreos <- seq(1,15175, by = 159.5)
submuestreos
#truco 


#----b.indices de diversidad alfa-----

?evenness
plot_richness(ps_rar, color = "Procedencia", x = "Procedencia", 
              measures = c("Simpson", "Chao1", "Shannon")) + 
  geom_boxplot(aes(fill = Procedencia), alpha=.7) + 
  scale_color_manual(values = 
                       c("#90EE90", "#7CCD7C", "#548B54","#FFD39B","#EEC591","#CDAA7D")) +
  scale_fill_manual(values = 
                      c("#90EE90", "#7CCD7C", "#548B54","#FFD39B","#EEC591","#CDAA7D"))


#----Instalacion de paqueterias 4----
library(vegan)
#----importar objeto phyloseq (ps)----
ps <- readRDS("results/objeto_phyloseq.RDS")
ps

#----d.indices de diversidad beta----
#1. generar matriz de disimilitud
#bray curtis
bray <- vegdist(t(otu_table(ps)),"bray")
view(as.matrix(bray))
?vegdist

#ordenamiento
#----PCoA----
# Estandariza datos
pslog <- transform_sample_counts(ps, function(x) log(1 + x))

# Ordenamiento Unifrac ponderado
pcoa_bray<- ordinate(pslog, method = "PCoA", distance = "bray", weigthed = F)
evals <- pcoa_bray$values$Eigenvalues

# Graficar los eigen values
plot(evals)
colores_pcoa <- c("#FFDE8B","#FFA88B","red")

# PCoA con Bray curtis
plot_ordination(pslog, pcoa_bray, color = "Procedencia", shape = "Procedencia") +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  geom_text(aes(label=sample),hjust=0, vjust=1,size = 2)+
  geom_point(size=3)+ 
  scale_fill_manual(values=colores_pcoa) + 
  scale_color_manual(values = colores_pcoa)+
  ggtitle("PCoA Bray Curtis")+
  theme_classic() -> plot_bray
plot_bray

# Hacer una funcion para graficar pcoa
pcoa_plot <- function(physeq_object, distancia = as.character(), 
                      ponderado = FALSE, 
                      colores){
  pcoa_ordination<- ordinate(pslog, method = "PCoA", 
                             distance = as.character(distancia),
                             weighted = FALSE)
  evals <- pcoa_ordination$values$Eigenvalues
  
  # PCoA con Bray curtis
  plot_ordination(pslog, pcoa_ordination, 
                  color = "Procedencia", shape = "Procedencia") +
    coord_fixed(sqrt(evals[2] / evals[1])) +
    # geom_text(aes(label=sample),hjust=0, vjust=1,size = 2)+
    geom_point(size=3)+ 
    scale_fill_manual(values=colores) + 
    scale_color_manual(values = colores)+
    theme_classic() -> plot_pcoa
  return(plot_pcoa)
}

#
pcoa_plot(pslog, "bray", ponderado = F, colores_pcoa) +
  ggtitle("PCoA Bray Curtis") -> pcoa_bray
pcoa_bray

#----Ordenamiento Unifrac no ponderado----
pcoa_plot(pslog, "unifrac", ponderado = F, colores_pcoa) +
  ggtitle("PCoA Unifrac no ponderado") -> pcoa_unifrac_uw
pcoa_unifrac_uw



#--------------------ETAPA 6-----------------------
#llamado de librerias
library(phyloseq)
library(microbiome)
library(microbiomeutilities)

#librerias para manejo de datos y visualizacion
library(tidyverse)
library(magrittr)

#importar objeto phyloseq (ps)
ps <- readRDS("results/objeto_phyloseq.RDS")
ps

#obteenr tabla de taxonomia de un objeto phyloseq
tax_table(ps)


#analisis de composicion taxonomica
psPhylum <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% 
                          c("", "uncharacterized") & !Phylum %in% c("", "Unknow"))
taxa_composition <- psPhylum

colores <- c("#20DE8B","#CCDE8B","#FFDE8B","#FFA88B","#FF6AD5","#FF6AD5","#C874AA","#AD8CFF","#966BFF","#90CFFF","#2980B9")

#seleccionar el top10 de phylum más abundantes
phylum <- aggregate_top_taxa2(taxa_composition, "Phylum", top = 10)
ps

phylum

#transformar los datos a abundancia relativa
phylum_abrel <- transform(phylum, "compositional")

#graficos de barras apiladas
plot_abrel <- plot_composition(phylum_abrel, otu.sort = "abundance",
                               x.label = "sample", group_by = "Procedencia")

plot_abrel

plot_abrel <- plot_composition(phylum_abrel, otu.sort = "abundance", 
                               x.label = "Samples", group_by = "Procedencia")


plot_abrel +
  scale_fill_manual("Phylum", values = colores) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        legend.title = element_text(size = 8),
        legend.position = "bottom") +
  ggtitle("Relative abundance Phylum")

#---------ETAPA 7----------------------
library(readxl)
library(multcomp)
library(mvtnorm)
library(survival)
library(TH.data)
library(MASS)

data24 <- read_excel(".../2_Data_New", sheet = "Genero1")
data <- as.matrix(data24)
View(data24)
names(data24)
class("abundancia")
class("replica")
class("zona")


######
factorial
attach(data24)
?attach
names(data24)
str(data24)
shapiro.test(abundancia)
typestarter<- factor(`subzona`)
time <- factor(`zona`)
dev.next()
?dev.next


boxplot(abundancia~`subzona`, 
        col= c("#FFD39B","#EEC591","#CDAA7D","#90EE90", "#7CCD7C", "#548B54"),
        data = data24, main="Genero1")
boxplot(abundancia~`zona`, 
        col = c("#FF7F50","#00CD00"), 
        data = data24, main="Genero1")
boxplot(abundancia~`subzona`*`zona`,data =data24, 
        main="Genero1", 
        col=c("#FFD39B","#EEC591","#CDAA7D","#90EE90", "#7CCD7C", "#548B54"))

# # # # modelo lineal con dos factores# # # # # 
fanalysis <- lm(abundancia~ (`zona` + `subzona`)^2)
anova3 <- aov(fanalysis)
summary(anova3)
interaction.plot(`zona`, `subzona`, abundancia, main= "Genero1", 
                 xlab = "zona", ylab = "subzona", col = c(1:3))

TukeyHSD(anova3)

###grafico Martha analysis
library(tidyverse)
ggplot(data = data24, aes(x = subzona, y = abundancia, color = zona)) +
  geom_boxplot() +
  theme_bw()

ggplot(data = data24, aes(x = zona, y = abundancia, color = subzona)) +
  geom_boxplot() +
  theme_bw()

