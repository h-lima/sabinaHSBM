

<img width="35%" align= "right" alt="logo_s-1" src="https://github.com/geoSABINA/sabinaNSDM/assets/168073517/d29288b9-c1a7-47aa-8753-918c931e4c53"/>

# *sabinaHSBM*: Hierarchical Stochastic Block Model for link prediction and network reconstruction


## Overview

The ***sabinaHSBM*** R package offers an advanced tool for predicting and reconstructing links in bipartite binary networks using the **Hierarchical Stochastic Block Model (HSBM)**. Networks are widely used across various fields to represent complex relationships and interactions—from ecological systems to social structures—capturing the underlying patterns and connections between entities or nodes. However, real-world networks are often incomplete and may contain errors due to sampling constraints, making reliable link prediction essential to enhance network accuracy and completeness.

***sabinaHSBM*** addresses these challenges by implementing a powerful HSBM-based approach that uses both observed interactions and the inherent structural patterns within networks to identify missing (*missing links*) and potentially erroneous links (*spurious links*). Among the various network reconstruction techniques, HSBM stands out for its effectiveness, providing a nonparametric network reconstruction based on Bayesian inference that assigns error probabilities to observed and unobserved links. This approach minimizes subjective decisions, enabling robust statistical inference and model selection.

Although HSBM is available as part of Python's `graph-tool` package, its adoption has been limited among R users across disciplines. ***sabinaHSBM*** bridges this gap, making HSBM accessible in R and providing a practical tool for researchers to:

- Reconstruct complex bipartite networks, even from incomplete or unreliable data.
- Generate probability estimates for unobserved or potentially misrecorded links, enhancing network accuracy.
- Identify groups within networks based on interaction patterns, uncovering insights into structural organization.
- Visualize hierarchical partitions of reconstructed networks.

The accessibility of ***sabinaHSBM*** for R users opens this advanced methodology to a wider community, empowering researchers across diverse fields to conduct accurate link predictions, address data quality issues, and gain a deeper understanding of complex network structures in their respective domains.


## Key Features of ***sabinaHSBM***: #@@@JMB a esto hay que darle una vuelta para enfatizar las virtudes del paquete

- Generates an input data with cross-validation across multiple partitions to test model robustness.
- Predicts missing and spurious links using HSBM and evaluates the predictive performance.
- Allows reconstruction of the full network, with options to set threshold criteria.
- Includes visualization tools for comparing known and reconstructed matrices.


## Installation

**Dependencies**: Requires R (>= 4.3.0)

Install the released version of ***sabinaHSBM*** from [GitHub](https://github.com) using the following commands:

```r
library(remotes)
remotes::install_github("h-lima/sabinaHSBM")
```

### Citing sabinaHSBM package

A research paper detailing the functions and methodologies of the ***sabinaHSBM*** package is in preparation. Until its publication, please cite the package as follows:

> Lima, H., ... (2024). sabinaHSBM: an R package for Hierarchical Stochastic Block Model-based link prediction and network reconstruction.  
> Aquí doi -----

## Core Functions in *sabinaHSBM*

| Step                             | Function                | Description                                                                                  |
|----------------------------------|-------------------------|----------------------------------------------------------------------------------------------|
| **Prepare the Input**             | `hsbm.input`           | Prepares cross-validated input data for HSBM analysis                                        |
| **Reconstruct the network**              | `hsbm.predict`         | Predicts link probabilities and hierarchical node grouping for each fold       |
| **Get the results** | `hsbm.reconstructed` | Return the reconstructed binary network from flexible threshold settings and provides evaluation metrics    |
| **Results Exploration and Visualization** | `get_hsbm_results`  | Extracts and summarizes the results of HSBM predictions                                      |
|                                  | `plot_interaction_matrix` | Visualizes a (reconstructed) binary bipartite matrix                                      |
|                                  | `top_links`            | Identifies and ranks top predicted links for targeted analysis                               |


## Example Workflow

This example demonstrates how to use the ***sabinaHSBM*** package to predict and reconstruct links in a bipartite binary network.
-   [1. Prepare the Input](#data_preparation)
-   [2. Reconstruct the network](#link_prediction)
-   [3. Get the results](#network_reconstruction)
-   [4. Results Exploration and Visualization](#visualization)

  
### 1. Prepare the Input  <a name="data_preparation">  
Begin by preparing the input data with `hsbm.input`:

The input data for `hsbm.input` should be a **binary bipartite matrix** representing interactions between two distinct sets of nodes. Rows represent one type of node (e.g., hosts), columns represent the other type (e.g., parasites), and each entry is binary: `1` indicates a link (e.g., a host-parasite interaction is observed), and `0` indicates no link.

```r
# Set working directory
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
    method = "binary_classifier",  # Choose method for link prediction(*)
)

summary(myInput)   # Summarizes network characteristics
```

(*)`hsbm.input` offers two methods for link prediction:
   - **`"binary_classifier"`**: Focuses on predicting probabilities for currently **unobserved links** (`0s`). Use this method if you want to identify **missing links** (unobserved links likely to exist) in partially incomplete networks.
   - **`"full_reconstruction"`**: Estimates probabilities for **all links** (both `0s` and `1s`), resulting in a fully reconstructed probability matrix. This method can identify both **missing links** (unobserved links likely to exist) and **spurious links** (observed links that might be erroneous) in incomplete or error-prone networks.


      
### 2. Reconstruct the network <a name="link_prediction">  
Use ‘hsbm.predict’ function to predict the marginal posterior probabilities of each link according to network reconstruction. It requires an object of hsbm.input class, and returns an object of the hsbm.predict class with link probabilities and the hierarchical organization of nodes in groups, for each fold. 

```r
myPred <- hsbm.predict(myInput)
```

### 3. Get the results <a name="network_reconstruction">  
Reconstruct the network with `hsbm.reconstructed` and customize threshold settings as needed:
The function ‘hsbm.reconstructed()’ generates a reconstructed binary matrix (linked/unlinked) based on a user-defined threshold, and evaluates model performance with multiple metrics.

```r
# Reconstruct the network
myReconst <- hsbm.reconstructed(
    myPred,                # HSBM output object with predicted links
    rm_documented = TRUE,    # Remove documented links
    threshold = "prc_closest_topright", # Threshold method to binarize continuous predictions
    new_matrix_method = "average_thresholded" # Method for creating the reconstructed binary matrix
)

# Evaluate model performance
summary(myReconst)        # Summarizes reconstructed network characteristics and provides model performance metrics
paste("The models predicted on average", round(mean(myReconst$tb$pred_held_ones), 1), "% of the held out links.")
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

## Tutorials
*Aquí poner enlaces a los tutoriales (uno con ejemplo de paralelización y otro con lo del docker para windows)* #@@@JMB

## Contributions

We welcome contributions and suggestions! To contribute, please open an "issue" or submit a "pull request" in this repository.

