# Load required libraries
library(tsne)

# Load the data
data1 <- read.csv(gzfile('/Users/parasdhiman/Downloads/DSG/GSE81861_Cell_Line_COUNT.csv'))

# Extract count data and gene names
data <- data1[, 2:562]
rownames(data) <- data1[, 1]

# Calculate mean and standard deviation
gme <- apply(data, 1, mean)
gsd <- apply(data, 1, sd)

# Compute coefficient of variation (CV)
cv <- gsd / gme

# Square the CV to get CV^2
cv2 <- cv^2

# Plot mean vs. CV^2
plot(log(gme + 1), log((gsd / gme) * (gsd / gme) + 1),
     xlab = "Mean Expression (log scale)",
     ylab = "CV^2 (log scale)",
     main = "Mean vs. CV^2 Plot")

# Perform binning based on mean expression levels
num_bins <- 10  # Define the number of bins
bin_means <- cut(gme, breaks = num_bins, labels = FALSE)

# Initialize a list to store selected genes
selected_genes <- list()

# Loop through each bin and select top genes with high CV^2
for (bin in 1:num_bins) {
  # Get genes within the current bin
  genes_in_bin <- rownames(data)[bin_means == bin]
  
  # Sort genes within the bin based on CV^2
  sorted_genes <- genes_in_bin[order(cv2[genes_in_bin], decreasing = TRUE)]
  
  # Select top genes with high CV^2 (e.g., top 10)
  selected_genes[[bin]] <- sorted_genes[1:10]
}

# Flatten the list of selected genes
selected_genes <- unlist(selected_genes)

# Define threshold for mean expression
threshold <- 10  # Example threshold (adjust as needed)

# Filter out genes with low mean expression
selected_genes <- selected_genes[gme[selected_genes] > threshold]

# Perform Principal Component Analysis (PCA)
pcs <- svd(data, nv = 30)

# Perform tSNE
ts <- tsne(pcs$v)

# Plot tSNE coordinates
plot(ts[, 1], ts[, 2], 
     xlab = "tSNE1", 
     ylab = "tSNE2", 
     main = "tSNE Plot")

# Color by cell type
cnames <- strsplit(colnames(data1), "__")
df <- matrix("blue", 561, 1)

for (i in 2:561) {
  df[i - 1] <- (cnames[i])[[1]][3]
}

# Plot with color by cell type
plot(ts[, 1], ts[, 2], col = df,
     xlab = "tSNE1",
     ylab = "tSNE2",
     main = "tSNE Plot with Color by Cell Type")




#-------------------Question3-------------------
# Load required libraries
library(tsne)
library(igraph)

# Assuming you have already loaded and preprocessed the data
# Replace data with your preprocessed dataset

# Step 1: Feature Selection
# Calculate mean and standard deviation
gme <- apply(data, 1, mean)
gsd <- apply(data, 1, sd)

# Calculate coefficient of variation (CV2)
cv2 <- (gsd / gme) * (gsd / gme)

# Create bins based on mean expression level
bins <- cut(gme, breaks = 10)  # Adjust the number of bins as needed

# Get top genes with high CV2 in each bin
top_genes <- lapply(split(cv2, bins), function(x) {
  idx <- order(x, decreasing = TRUE)[1:min(10, length(x))]  # Select top 10 genes per bin
  names(x)[idx]
})

# Combine top genes from all bins
top_genes <- unlist(top_genes)

# Filter out genes with low mean expression
top_genes <- top_genes[gme[top_genes] > threshold]  # Set your threshold for mean expression

# Step 2: Principal Component Decomposition
pcs <- svd(data, nv = 30)

# Step 3: tSNE
ts <- tsne(pcs$v)

# Step 4: Visualize tSNE coordinates
plot(ts[,1], ts[,2], col = "blue", main = "tSNE Plot")

# Step 5: Predict gene network
# Subset data to include top expressed genes
top_expressed_genes <- head(order(gme, decreasing = TRUE), 7000)

# Calculate correlations between genes
cor_matrix <- cor(t(data[top_expressed_genes, ]))

# Replace negative correlations with 0
cor_matrix[cor_matrix < 0] <- 0

# Create a graph using correlations
graph <- graph_from_adjacency_matrix(cor_matrix, mode = "undirected")

# Reduce edges to keep only 30000-50000 edges
edge_weights <- E(graph)$weight
threshold <- quantile(edge_weights, probs = c(0.75, 0.90))  # Adjust as needed
delete_edges <- which(edge_weights < threshold[1] | edge_weights > threshold[2])
graph <- delete_edges(graph, delete_edges)


# Calculate degree centrality
degree_centrality <- degree(graph, mode = "all")  # mode = "all" considers both in-degree and out-degree

# Get the gene with the highest degree centrality
top_central_gene <- names(which.max(degree_centrality))

# Print gene with highest degree centrality
print(paste("Gene with highest degree centrality:", top_central_gene))

# Plot the gene network (optional)
plot(graph, vertex.size = 2, vertex.label = NA)

