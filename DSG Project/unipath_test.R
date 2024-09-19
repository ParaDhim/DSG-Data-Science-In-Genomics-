library(UniPath)
data("c5.bp.v6.0.symbols")

library(data.table)

library(tibble)
library(readr)

data("human_null_model")
exp_data= read_csv('/Users/parasdhiman/Downloads/DSG Project/GSE171524_processed_data.csv')
#exp_data= read.csv('/Users/parasdhiman/Downloads/GSE171524_processed_data.csv',header = TRUE)
main_exp_data <- exp_data[!duplicated(exp_data$GeneName), ]
main_exp_data <- as.data.frame(main_exp_data)  # Convert to data frame
rownames(main_exp_data) <- main_exp_data[, 1]  # Set row names
main_exp_data <- main_exp_data[, -1]  # Remove the first column
main_exp_data_transposed <- t(main_exp_data)
lung_metadata <- fread("/Users/parasdhiman/Downloads/DSG Project/GSE171524_lung_metaData.txt")
lung_metadata <- lung_metadata[-1, ]
row_names <- rownames(main_exp_data_transposed)


# Remove the first row
lung_metadata <- lung_metadata[-1, ]

# Sample one-fourth of the observations randomly
#lung_metadata_sampled <- lung_metadata[sample(nrow(lung_metadata), nrow(lung_metadata) / 4), ]

# Save the sampled data to a CSV file
lung_count <- sum(lung_metadata$group == "Control")

# Print the count
print(lung_count)

lung_count_cov <- sum(lung_metadata$group == "COVID-19")

# Print the count
print(lung_count_cov)




# Calculate the number of rows to select (one-fourth of total rows)
num_rows_select <- round(nrow(main_exp_data_transposed)/12)

# Randomly sample row indices
sampled_indices <- sample(1:nrow(main_exp_data_transposed), size = num_rows_select, replace = FALSE)

# Extract the sampled data
sampled_data <- main_exp_data_transposed[sampled_indices, ]

row_names <- rownames(sampled_data)

sampled_data <- cbind(NAME = row_names, sampled_data)

#write.csv(lung_metadata_sampled, file = "lung_metadata_sampled.csv", row.names = FALSE)



# Merge sampled_data with lung_metadata
merged_data <- merge(sampled_data, lung_metadata[, .(NAME, disease__ontology_label,group)], by = "NAME", all.x = TRUE)



label_counts <- table(merged_data$disease_ontology_label)



# Count occurrences of "COVID-19" and "normal"
covid_count <- sum(merged_data$disease__ontology_label == "COVID-19")
normal_count <- sum(merged_data$disease__ontology_label == "normal")

# Print the counts
print(paste("COVID-19 count:", covid_count))
print(paste("Normal count:", normal_count))

# Add disease_ontology_label column from merged_data to sampled_data
sampled_data <- cbind(sampled_data, disease_ontology_label = merged_data$disease__ontology_label,group = merged_data$group)
sampled_data <- as.data.frame(sampled_data)
#write.csv(sampled_data, file = "sampled_Dat.csv", row.names = FALSE)
# Filter COVID-19 samples and select 3000 randomly
covid_samples <- sampled_data[sampled_data$disease_ontology_label == "COVID-19", ]
covid_samples <- covid_samples[sample(nrow(covid_samples), 3000), ]

# Filter normal samples and select 3000 randomly
normal_samples <- sampled_data[sampled_data$disease_ontology_label == "normal", ]
normal_samples <- normal_samples[sample(nrow(normal_samples), 3000), ]

# Combine the sampled data
final_sampled_data <- rbind(covid_samples, normal_samples)
covid_count <- sum(final_sampled_data$disease_ontology_label == "COVID-19")
normal_count <- sum(final_sampled_data$disease_ontology_label == "normal")

# Print the counts
print(paste("COVID-19 count:", covid_count))
print(paste("Normal count:", normal_count))
# Save the final sampled data to a CSV file
write.csv(final_sampled_data, file = "final_sampled_data.csv", row.names = FALSE)


sampled_data <- as.data.frame(sampled_data)
sampled_data$disease_ontology_label


fr_score <- sampled_data[, !names(sampled_data) %in% c("NAME", "disease_ontology_label")]
fr_score_t <- t(fr_score)





# Now new_data has "NAME" column from row_names and other columns from main_exp_data_transposed


Pval = binorm(human_null_data)

combp_ref = combine(c5.bp.v6.0.symbols,human_null_data,rownames(human_null_data),Pval,thr=2000)


rm_row= which(rownames(main_exp_data)=="C15ORF37")
main_exp_data= main_exp_data[-rm_row,]
rm_row= which(rownames(main_exp_data)=="C1ORF220")
main_exp_data= main_exp_data[-rm_row,]
rm_row= which(rownames(main_exp_data)=="C2ORF15")
main_exp_data= main_exp_data[-rm_row,]
rm_row= which(rownames(main_exp_data)=="C6ORF165")
main_exp_data= main_exp_data[-rm_row,]



# Find rows where any value is negative
#negative_rows <- apply(main_exp_data, 1, function(row) any(row < 0))

# Subset the data to keep only rows where no value is negative
#clean_main_exp_data <- main_exp_data[!negative_rows, ]

# Take the absolute value of each element in the data frame
#abs_main_exp_data <- abs(main_exp_data)

#Pval1 = binorm(main_exp_data)
#combp = combine(c5.bp.v6.0.symbols,main_exp_data,rownames(main_exp_data),Pval1,thr=2)

fr_score_t <- as.matrix(fr_score_t)
#fr_score_t <- abs(fr_score_t)
fr_score_m <- apply(fr_score_t, 2, as.numeric)



# Check for missing values in fr_score_m
if (anyNA(fr_score_m)) {
  # Handle missing values (e.g., impute or remove them)
  # For example, you can impute missing values with the mean
  fr_score_m[is.na(fr_score_m)] <- mean(fr_score_m, na.rm = TRUE)
}

fr_score_m <- abs(fr_score_m)


Pval1 = binorm(fr_score_m)
combp = combine(c5.bp.v6.0.symbols,fr_score_m,rownames(fr_score_m),Pval1,thr=2000)

scores = adjust(combp,combp_ref)




#------------------------------------------------------------------------------------------------
library(UniPath)
library(tibble)
data("c5.bp.v6.0.symbols")

library(data.table)
#exp_data= read.csv('/Users/parasdhiman/Downloads/GSE81861_CRC_tumor_all_cells_COUNT.csv',header = TRUE,row.names=1)
#main_data= read.csv(gzfile('/Users/parasdhiman/Downloads/GSE81861_CRC_tumor_all_cells_COUNT.csv'))
#main_data = read.csv('/Users/parasdhiman/Downloads/GSE171524_processed_data.csv ',header = TRUE,row.names=1)
#main_data <- fread('/Users/parasdhiman/Downloads/GSE171524_processed_data.csv')
library(readr)

# Read the CSV file using readr's read_csv()
main_data <- read_csv('/Users/parasdhiman/Downloads/GSE171524_processed_data.csv')

data("human_null_model")
data("GSE52583_expression_data") # in this example expression profile, look how the formatting is done (row and column names). You have to make sure your data is in the same format.

exp_data= read_csv('/Users/parasdhiman/Downloads/GSE171524_processed_data.csv')
# Add row names as a column
#main_exp_data <- rownames_to_column(main_exp_data, var = "...1")

# Remove duplicated rows based on the first column
#main_exp_data <- main_exp_data[!duplicated(main_exp_data$`...1`),]
main_exp_data= exp_data[!duplicated(exp_data$...1),]
rownames(main_exp_data)= main_exp_data[,1] 
#main_exp_data= main_exp_data[,-1]
## Note: formating the matrix was needed here. it should be the rownames which has gene-names, not a column.

# Assuming expression_data is your dataframe
exp_data[,1] <- seq_len(nrow(exp_data))


##Converting mouse null data into p-values
Pval = binorm(human_null_data)
##Converting gene expression data into p-values
Pval1 = binorm(main_exp_data)

##Combining of p-values for null model data matrix
combp_ref = combine(c5.bp.v6.0.symbols,human_null_data,rownames(human_null_data),Pval,thr=1000)
##Combining of p-values for gene expression data matrix
combp = combine(c5.bp.v6.0.symbols,main_exp_data,rownames(main_exp_data),Pval1,thr=1000)
# Look into the error message, try to troubleshoot based on the information you have

#rm_row= which(rownames(main_exp_data)=="C15ORF37")
#main_exp_data= main_exp_data[-rm_row,]
#rm_row= which(rownames(main_exp_data)=="C1ORF220")
#main_exp_data= main_exp_data[-rm_row,]
#rm_row= which(rownames(main_exp_data)=="C2ORF15")
#main_exp_data= main_exp_data[-rm_row,]
#rm_row= which(rownames(main_exp_data)=="C6ORF165")
#main_exp_data= main_exp_data[-rm_row,]


Pval1 = binorm(main_exp_data)
combp = combine(c5.bp.v6.0.symbols[1:20,],main_exp_data,rownames(main_exp_data),Pval1,thr=2)

scores = adjust(combp,combp_ref)


# Accessing the component containing adjusted p-values or scores
adjpvalog <- scores$adjpvalog

# Reshape the matrix into a vector
adjpvalog_vector <- as.vector(adjpvalog)

# Get the names of the pathways
pathway_names <- rownames(adjpvalog)

# Assuming labels are available or loaded from a file


########################################################################################
# Accessing the component containing adjusted p-values or scores
adjpvalog <- scores$adjpvalog

# Reshape the matrix into a vector
adjpvalog_vector <- as.vector(adjpvalog)

# Get the names of the pathways
pathway_names <- rownames(adjpvalog)

# Combine pathway names and corresponding scores
pathway_scores <- data.frame(Pathway = pathway_names, Score = adjpvalog_vector)

# Order the pathways based on their scores
sorted_pathways <- pathway_scores[order(pathway_scores$Score, decreasing = TRUE), ]

# Select the top 3 pathways
top_3_pathways <- head(sorted_pathways, 3)
# Print the top 3 pathways
print(top_3_pathways)

# Counting the number of columns
num_columns <- ncol(scores)

# Printing the number of columns
print(num_columns)

library(data.table)

# Check for missing values in pathway_names or adjpvalog_vector
if (any(is.na(pathway_names)) || any(is.na(adjpvalog_vector))) {
  stop("Missing values found in pathway_names or adjpvalog_vector.")
}

# Convert pathway scores to a data.table
pathway_scores_dt <- data.table(Pathway = pathway_names, Score = adjpvalog_vector)

# Convert labels to a data.table directly without creating a copy
setDT(labels)

# Merge pathway scores with labels using data.table merge
data_with_labels <- merge(pathway_scores_dt, labels, by.x = "Pathway", by.y = "sample", all.x = TRUE)

# Drop NA values (if any)
data_with_labels <- na.omit(data_with_labels)

# Convert labels to a factor
data_with_labels$label <- as.factor(data_with_labels$label)

# Split the dataset into training and testing sets
set.seed(123) # For reproducibility
train_index <- sample(nrow(data_with_labels), 0.7 * nrow(data_with_labels))
train_data <- data_with_labels[train_index, ]
test_data <- data_with_labels[-train_index, ]

# Train a logistic regression model
model <- glm(label ~ ., data = train_data, family = binomial)

# Predict on the test set
predictions <- predict(model, newdata = test_data, type = "response")

# Convert probabilities to class labels
predicted_labels <- ifelse(predictions > 0.5, "diseased", "non-diseased")

# Evaluate the model
accuracy <- mean(predicted_labels == test_data$label)
print(paste("Accuracy:", accuracy))
