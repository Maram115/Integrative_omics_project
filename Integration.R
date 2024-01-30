# Library Calling ----
library(MOFA2)
library(MOFAdata)
library(basilisk)
library(tidyverse)
library(data.table)
library(VennDiagram)
# Reading, tidying and exploring the data ----

# Reading
Exp <- read.table("Normalized_output/normalized_counts_mRNA.csv",sep = ",",header = T,row.names = 1) 
miRNA <- read.table("Normalized_output/count_table.csv",sep = ",",header = T,row.names = 1) 
Methyl <- read.table("Normalized_output/M_values_f.csv",sep = ",",header = T,row.names = 1) 
Meta <- read.table("melanoma",sep = "\t",header = T)
Meta$sampleID <- gsub("-",".",Meta$sampleID)
colnames(Meta)[1] <- "sample"

### Tidying and exploring
## Structure of samples for The 3 datasets
miRNA_inter_Exp <- intersect(colnames(miRNA),colnames(Exp)) # There are 451 sample in common between Expression and miRNA datasets    #length(miRNA_inter_Exp)
miRNA_inter_Methyl <- intersect(colnames(miRNA),colnames(Methyl)) # There are 451 sample in common between Methylation and miRNA datasets    #length(miRNA_inter_Methyl)
Exp_inter_Methyl <- intersect(colnames(Exp),colnames(Methyl)) #There are 472 sample in common between Methylation and Expression datasets    #length(Exp_inter_Methyl)
Methyl_inter_miRNA_inter_Exp <- intersect(colnames(Methyl),miRNA_inter_Exp) #There are 450 sample in common between Expression ,miRNA and Methylation datasets   #length(Methyl_inter_miRNA_inter_Exp)


u_Exp<- colnames(Exp)[!(colnames(Exp) %in% Methyl_inter_miRNA_inter_Exp)] # There are 23 samples of Expression dataset that are not present in BOTH other datasets (i.e They maybe be found at one dataset without the other or it may be uniquely found in Expression dataset only) 
u_miRNA<- colnames(miRNA)[!(colnames(miRNA) %in% Methyl_inter_miRNA_inter_Exp)] # There are 2 samples of miRNA dataset that are not present in BOTH other datasets (i.e They maybe be found at one dataset without the other or it may be uniquely found in miRNA dataset only)
u_Methyl<- colnames(Methyl)[!(colnames(Methyl) %in% Methyl_inter_miRNA_inter_Exp)] # There are 25 samples of Methylation dataset that are not present in BOTH other datasets (i.e They maybe be found at one dataset without the other or it may be uniquely found in Methylation dataset only)
u_miRNA_Exp<- intersect(u_miRNA,u_Exp) # There is just one sample that found in common between miRNA and Expression datasets without Methylation dataset ("TCGA.XV.AB01.06")
u_miRNA_Methyl <- intersect(u_miRNA,u_Methyl) #There is just one sample that found in common between miRNA and Methylation datasets without Expression dataset ("TCGA.FW.A3R5.11")
u_Exp_Methyl <- intersect(u_Exp,u_Methyl) # There are 22 samples that found in common between Expression and Methylation datasets without miRNA dataset
# u_Exp has 23 samples 22 of which is common with Methylation dataset and 1 sample common with miRNA dataset
# u_miRNA has 2 samples 1 of which is common with Methylation dataset and 1 sample common with Expression dataset 
# u_Methyl has 25 samples 22 of which is common with Expression dataset and 1 sample common with miRNA dataset and 2 samples are unique for Methylation dataset

# Venn diagram for datasets samples
grid.newpage()
draw.triple.venn(area1 = length(colnames(Methyl)),
                 area2 = length(colnames(Exp)),
                 area3 = length(colnames(miRNA)),
                 n12 = length(Methyl_inter_miRNA_inter_Exp)+length(u_Exp_Methyl),
                 n23 = length(Methyl_inter_miRNA_inter_Exp)+length(u_miRNA_Exp),
                 n13 =length(Methyl_inter_miRNA_inter_Exp)+length(u_miRNA_Methyl),
                 n123 = length(Methyl_inter_miRNA_inter_Exp),
                 fill = c("green","yellow","blue"),
                 category = c("Methylation","Expression","miRNA"))


# Merging the datasets
all_omics_matrix <- bind_rows(Exp,miRNA,Methyl)
new_Exp <- all_omics_matrix[1:length(rownames(Exp)),]
new_miRNA<- all_omics_matrix[length(rownames(Exp))+1:length(rownames(miRNA)),]
new_Methyl <- all_omics_matrix[length(rownames(Exp))+length(rownames(miRNA))+1:length(rownames(Methyl)),]
new_Meta <- Meta[Meta$sample %in% colnames(all_omics_matrix),]
new_Meta <- new_Meta[match(colnames(new_Exp),new_Meta$sample),]

# creating all omics list
all_omics_list <- list(data.matrix(new_Exp),data.matrix(new_miRNA),data.matrix(new_Methyl))
names(all_omics_list)<-c("Expression", "miRNA", "Methylation")

# Creating Mofa object ----
MOFAobject<- create_mofa(all_omics_list)
samples_metadata(MOFAobject) <- new_Meta
plot_data_overview(MOFAobject)

#### Fit MOFA model ----
## Defining Data Options 

data_opt <- get_default_data_options(MOFAobject)
data_opt

## Defining Model Options 

model_opt <- get_default_model_options(MOFAobject)
model_opt

## Defining Training Options 

train_opt <- get_default_training_options(MOFAobject)
train_opt$verbose <- T
train_opt

## Preparing MOFA
MOFAobject_prepared <- prepare_mofa(object = MOFAobject,
                           data_options = data_opt,
                           model_options = model_opt,
                           training_options = train_opt)

## Running MOFA
MOFAobject_run <- run_mofa(MOFAobject_prepared,use_basilisk = T)

# Identifying the top factors----
plot_factor_cor(MOFAobject_run)

# correlation between most of metadata columns
# correlate_factors_with_covariates(MOFAobject_run, 
#                                   covariates = c("gender",
#                                                  "vital_status",
#                                                  "sample_type",
#                                                  "new_tumor_event_after_initial_treatment",
#                                                  "tumor_tissue_site",
#                                                  "person_neoplasm_cancer_status",
#                                                  "melanoma_ulceration_indicator",
#                                                  "melanoma_clark_level_value",
#                                                  "distant_metastasis_anatomic_site",
#                                                  "primary_neoplasm_melanoma_dx",
#                                                  "history_of_neoadjuvant_treatment",
#                                                  "system_version",
#                                                  "tissue_prospective_collection_indicator"
#                                   ), 
#                                   plot="log_pval")



## correlation between factors and gender, vital status, sample type and new_tumor_event_after_initial_treatment
correlate_factors_with_covariates(MOFAobject_run, 
                                  covariates = c("gender",
                                                 "vital_status",
                                                 "sample_type",
                                                 "new_tumor_event_after_initial_treatment"
                                                 ), 
                                  plot="log_pval")


# Vital status
plot_factor(
  MOFAobject_run,
  factors = 3,
  color_by   = "vital_status",
  dot_size = 3,
  dodge = T,
  add_boxplot   = T
)


# sample type
plot_factor(
  MOFAobject_run,
  factors =8,
  color_by   = "sample_type",
  dot_size = 3,
  dodge = T,
  add_boxplot   = T
)

plot_factor(
  MOFAobject_run,
  factors = 12,
  color_by   = "sample_type",
  dot_size = 3,
  dodge = T,
  add_boxplot   = T
)


plot_factor(
  MOFAobject_run,
  factors = 13,
  color_by   = "sample_type",
  dot_size = 3,
  dodge = T,
  add_boxplot   = T
)


# new_tumor_event_after_initial_treatment
plot_factor(
  MOFAobject_run,
  factors = 12,
  color_by   = "new_tumor_event_after_initial_treatment",
  dot_size = 3,
  dodge = T,
  add_boxplot   = T
)



# gender
plot_factor(
  MOFAobject_run,
  factors = 10,
  color_by   = "gender",
  dot_size = 3,
  dodge = T,
  add_boxplot   = T
)




#identifying the different omics weights----

factor_var <- get_variance_explained(MOFAobject_run)

factor_var$r2_per_factor
plot_variance_explained(MOFAobject_run)
plot_variance_explained(MOFAobject_run,plot_total = T)



plot_weights_heatmap(MOFAobject_run,
                     view = "Methylation",
                     show_colnames=F)


plot_weights_heatmap(MOFAobject_run,
                     view = "Expression",
                     show_colnames=F)


plot_weights_heatmap(MOFAobject_run,
                     view = "miRNA",
                     show_colnames=F)




plot_variance_explained(MOFAobject_run,factors = 3)
plot_variance_explained(MOFAobject_run,factors = 10)
plot_variance_explained(MOFAobject_run,factors = 13)


# identifying the different features weights----
# feature weights for factor 3
plot_weights(
  MOFAobject_run, 
  view = "Methylation", 
  factor = 3, 
  nfeatures = 10
)


plot_top_weights(
  MOFAobject_run, 
  view = "Methylation", 
  factor = 3, 
  nfeatures = 10
)


plot_data_scatter(MOFAobject_run, 
                  view = "Methylation",
                  factor = 3,  
                  features = 10,
                  sign = "all",
                  color_by = "vital_status"
) + labs(y="Methylation")


plot_weights(
  MOFAobject_run, 
  view = "Expression", 
  factor = 3, 
  nfeatures = 10
)


plot_top_weights(
  MOFAobject_run, 
  view = "Expression", 
  factor = 3, 
  nfeatures = 10
)


plot_data_scatter(MOFAobject_run, 
                  view = "Expression",
                  factor = 3,  
                  features = 10,
                  sign = "all",
                  color_by = "vital_status"
) + labs(y="Expression")


plot_weights(
  MOFAobject_run, 
  view = "miRNA", 
  factor = 3, 
  nfeatures = 10
)


plot_top_weights(
  MOFAobject_run, 
  view = "miRNA", 
  factor = 3, 
  nfeatures = 10
)


plot_data_scatter(MOFAobject_run, 
                  view = "miRNA",
                  factor = 3,  
                  features = 10,
                  sign = "all",
                  color_by = "vital_status"
) + labs(y="miRNA")


#feature weights for factor 10
plot_weights(
  MOFAobject_run, 
  view = "Methylation", 
  factor = 10, 
  nfeatures = 10
)


plot_top_weights(
  MOFAobject_run, 
  view = "Methylation", 
  factor = 10, 
  nfeatures = 10
)


plot_data_scatter(MOFAobject_run, 
                  view = "Methylation",
                  factor = 10,  
                  features = 10,
                  sign = "all",
                  color_by = "gender"
) + labs(y="Methylation")


plot_weights(
  MOFAobject_run, 
  view = "Expression", 
  factor = 10, 
  nfeatures = 10
)


plot_top_weights(
  MOFAobject_run, 
  view = "Expression", 
  factor = 10, 
  nfeatures = 10
)


plot_data_scatter(MOFAobject_run, 
                  view = "Expression",
                  factor = 10,  
                  features = 10,
                  sign = "all",
                  color_by = "gender"
) + labs(y="Expression")


plot_weights(
  MOFAobject_run, 
  view = "miRNA", 
  factor = 10, 
  nfeatures = 10
)


plot_top_weights(
  MOFAobject_run, 
  view = "miRNA", 
  factor = 10, 
  nfeatures = 10
)


plot_data_scatter(MOFAobject_run, 
                  view = "miRNA",
                  factor = 10,  
                  features = 10,
                  sign = "all",
                  color_by = "gender"
) + labs(y="miRNA")


#feature weights for factor 13
plot_weights(
  MOFAobject_run, 
  view = "Methylation", 
  factor = 13, 
  nfeatures = 10
)


plot_top_weights(
  MOFAobject_run, 
  view = "Methylation", 
  factor = 13, 
  nfeatures = 10
)


plot_data_scatter(MOFAobject_run, 
                  view = "Methylation",
                  factor = 13,  
                  features = 10,
                  sign = "all",
                  color_by = "sample_type"
) + labs(y="Methylation")


plot_weights(
  MOFAobject_run, 
  view = "Expression", 
  factor = 13, 
  nfeatures = 10
)


plot_top_weights(
  MOFAobject_run, 
  view = "Expression", 
  factor = 13, 
  nfeatures = 10
)

plot_data_scatter(MOFAobject_run, 
                  view = "Expression",
                  factor = 13,  
                  features = 10,
                  sign = "all",
                  color_by = "sample_type"
) + labs(y="Expression")


plot_weights(
  MOFAobject_run, 
  view = "miRNA", 
  factor = 13, 
  nfeatures = 10
)


plot_top_weights(
  MOFAobject_run, 
  view = "miRNA", 
  factor = 13, 
  nfeatures = 10
)

plot_data_scatter(MOFAobject_run, 
                  view = "miRNA",
                  factor = 13,  
                  features = 10,
                  sign = "all",
                  color_by = "sample_type"
) + labs(y="miRNA")



#

# HEATMAPS----

plot_data_heatmap(
  MOFAobject_run, 
  view = "Methylation", 
  factor = 13, 
  features = 10, 
  show_rownnames = T,
  show_colnames=F,
  denoise = TRUE, 
  annotation_samples = "sample_type",
  show_colnnames= F
)




plot_data_heatmap(
  MOFAobject_run, 
  view = "Expression", 
  factor = 13, 
  features = 10, 
  show_rownnames = T,
  show_colnames=F,
  denoise = TRUE, 
  annotation_samples = "sample_type",
  show_colnnames= F
)


plot_data_heatmap(
  MOFAobject_run, 
  view = "miRNA", 
  factor = 13, 
  features = 10, 
  show_rownnames = T,
  show_colnames=F,
  denoise = TRUE, 
  annotation_samples = "sample_type",
  show_colnnames= F
)

# extracting top features name for Gprofiler----
top_genes_factor13 <- as.vector(plot_top_weights(MOFAobject_run,view = "Expression",factors = 13)[[1]]$feature)
top_genes_factor3 <- as.vector(plot_top_weights(MOFAobject_run,view = "Expression",factors = 3)[[1]]$feature)


