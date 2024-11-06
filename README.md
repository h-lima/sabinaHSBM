

<img width="35%" align= "right" alt="logo_s-1" src="https://github.com/geoSABINA/sabinaNSDM/assets/168073517/d29288b9-c1a7-47aa-8753-918c931e4c53"/>

# sabinaHSBM: Hierarchical Stochastic Block Model for link prediction and network reconstruction


## Overview

The **sabinaHSBM** R package offers an advanced tool for predicting and reconstructing links in bipartite binary networks using the **Hierarchical Stochastic Block Model (HSBM)**. Networks are widely used across various fields to represent complex relationships and interactions—from ecological systems to social structures—capturing the underlying patterns and connections between entities or nodes. However, real-world networks are often incomplete and may contain errors due to data limitations, making reliable link prediction essential to enhance network accuracy and completeness.

**sabinaHSBM** addresses these challenges by implementing an HSBM-based approach that leverages both observed interactions and the non-random structural patterns within networks to accurately estimate missing and spurious links. This package is particularly suited for researchers in ecology, sociology, and data science, providing a structured methodology to explore network relationships and validate link predictions in cases where comprehensive sampling is difficult to achieve.


## Key Features of **sabinaHSBM**:

- Generates an input matrix with cross-validation across multiple partitions to test model robustness.
- Predicts missing and spurious links using HSBM and evaluates the predictive performance.
- Allows reconstruction of the full network, with options to set threshold criteria.
- Includes visualization tools for comparing known and reconstructed matrices.


## Installation

**Dependencies**: Requires R (>= 4.3.0)

Install the released version of **sabinaHSBM** from [GitHub](https://github.com) using the following commands:

```r
library(remotes)
remotes::install_github("h-lima/sabinaHSBM")
```

### Citing sabinaHSBM package

A research paper detailing the functions and methodologies of the **sabinaHSBM** package is in preparation. Until its publication, please cite the package as follows:

> Lima, H., ... (2024). sabinaHSBM: an R package for Hierarchical Stochastic Block Model-based link prediction and network reconstruction.  
> Aquí doi -----

## Core Functions in sabinaHSBM

| Step                             | Function                | Description                                                                                  |
|----------------------------------|-------------------------|----------------------------------------------------------------------------------------------|
| **Data Preparation**             | `hsbm.input`           | Prepares cross-validated input data for HSBM analysis                                        |
| **Link Prediction**              | `hsbm.predict`         | Generates missing link predictions using the HSBM model                                      |
| **Network Reconstruction and Evaluation** | `hsbm.reconstructed` | Reconstructs the network with flexible threshold settings and provides evaluation metrics    |
| **Results Exploration and Visualization** | `get_hsbm_results`  | Extracts and summarizes the results of HSBM predictions                                      |
|                                  | `plot_interaction_matrix` | Visualizes a (reconstructed) binary bipartite matrix                                      |
|                                  | `top_links`            | Identifies and ranks top predicted links for targeted analysis                               |


## Example Workflow

This example demonstrates how to use the **sabinaHSBM** package to predict and reconstruct links in a bipartite binary network.
-   [1. Data Preparation](#data_preparation)
-   [2. Link Prediction](#link_prediction)
-   [3. Network Reconstruction and Evaluation](#network_reconstruction)
-   [4. Results Exploration and Visualization](#visualization)

  
### 1. Data Preparation  <a name="data_preparation">  
Begin by preparing the input data with `hsbm.input`:

The input data for `hsbm.input` should be a **binary bipartite matrix** representing interactions between two distinct sets of nodes. Rows represent one type of node (e.g., hosts), columns represent the other type (e.g., parasites), and each entry is binary: `1` indicates a link (e.g., a host-parasite interaction is observed), and `0` indicates no link.

```r
# setwd("/path/to/your/project")

# Load the sabinaNSDM package
library(sabinaNSDM)

# Load your data
dat <- read.csv("path/to/your_data.csv", row.names = 1)  # Load binary bipartite matrix with row names
dat <- as.matrix(dat)  # Convert data to matrix format
#data(dat, package = "sabinaHSBM")

# Prepare input data
myInput <- hsbm.input(
    dat,                    # Binary bipartite matrix of observed links
    n_folds = 10,           # Number of folds for cross-validation
    iter = 1000,            # Number of iterations for the HSBM model
    method = "binary_classifier",  # Choose method for link prediction*
)

summary(myInput)
```
*`hsbm.input` offers two methods for link prediction:
   - **`"binary_classifier"`**: Focuses on predicting probabilities for currently **unobserved links** (`0s`), leaving observed links (`1s`) unchanged. Use this method if you want to identify **missing links** (probable links among the unobserved interactions).
   - **`"full_reconstruction"`**: Estimates probabilities for **all links** (both `0s` and `1s`), resulting in a fully reconstructed probability matrix. This method can identify both **missing links** (absent links likely to exist) and **spurious links** (observed links that might be erroneous) by estimating the probability of each link, regardless of its initial state.

      
### 2. Link Prediction <a name="link_prediction">  
Use `hsbm.predict` to predict missing links in the network:

```r
myPred <- hsbm.predict(myInput)
summary(myPred)
```

### 3. Network Reconstruction and Evaluation <a name="network_reconstruction">  
Reconstruct the network with `hsbm.reconstructed` and customize threshold settings as needed:

```r
# Reconstruct the network
myReconst <- hsbm.reconstructed(
    myPred,                # HSBM output object with predicted links
    rm_documented = TRUE,    # Remove documented links
    threshold = "prc_closest_topright", # Threshold method to determine which links are retained
    new_matrix_method = "average_thresholded" # Method for creating the reconstructed matrix
)

# Evaluate model performance
paste("The models predicted on average", round(mean(myReconst$tb$pred_held_ones), 1), "% of the held out links.")
summary(myReconst) # Provides summary metrics for the reconstructed network
```

### 4. Results Exploration and Visualization <a name="visualization"> 

```r
# Summarize results
myRes <- get_hsbm_results(myReconst)

# Identify top predicted links
top_links_df <- top_links(myRes,
                          n = 25 # Number of top predicted links to retrieve
)

# Visualize interaction matrix
## Known matrix
plot_interaction_matrix(myReconst$data, order_mat = FALSE)
# Reconstructed matrix 
plot_interaction_matrix(myReconst$new_mat, order_mat = FALSE)
```



## Contributions

We welcome contributions and suggestions! To contribute, please open an "issue" or submit a "pull request" in this repository.

