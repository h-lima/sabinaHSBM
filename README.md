

<img width="38%" align= "right" alt="logo_s-1" src="https://github.com/geoSABINA/sabinaNSDM/assets/168073517/d29288b9-c1a7-47aa-8753-918c931e4c53"/>

# *sabinaHSBM*: Hierarchical Stochastic Block Model for link prediction and network reconstruction


## Overview

The ***sabinaHSBM*** R package makes the **Hierarchical Stochastic Block Model (HSBM)** available in R for predicting and reconstructing links in binary networks, a powerful tool for researchers in fields such as ecology, sociology, or data science. Networks are essential in representing complex systems, from species associations to social connections, but real-world networks often contain missing or spurious links due to sampling limitations. Accurately identifying these gaps and correcting errors enhances our understanding and reliability of network analyses.

***sabinaHSBM*** addresses these challenges by implementing a powerful HSBM-based approach that uses both observed interactions and the inherent structural patterns within networks to identify unobserved or unrealized links (*missing links*) and potentially erroneous links (*spurious links*). Among the various network reconstruction techniques, HSBM stands out for its effectiveness, providing a nonparametric network reconstruction based on Bayesian inference that assigns error probabilities to observed and unobserved links. This approach minimizes subjective decisions, enabling robust statistical inference and model selection.

With ***sabinaHSBM***, users can:
- **Reconstruct complex networks** from partial or uncertain data without requiring direct error estimates or external covariates.
- **Identify hierarchical groups** based on link patterns, uncovering insights into network structure.
- **Estimate probabilities for missing or spurious links**, increasing network accuracy.

This package provides an R-native environment to explore and validate link predictions, making HSBM-based network reconstruction more accessible.

Although ***sabinaHSBM*** relies on the  Python's `graph-tool` library ([Peixoto, 2014](https://doi.org/10.6084/m9.figshare.1164194)), which is Unix-only, we have made HSBM’s capabilities accessible to all R users by offering a pre-configured Docker container. This container includes all necessary dependencies, enabling seamless use on Windows, so researchers on any platform can leverage the power of ***sabinaHSBM*** without compatibility concerns.


## Key Features of ***sabinaHSBM***:

- Cross-Validation: Generates an input data with cross-validation across multiple partitions to test model robustness.
- Prediction and evaluation: Predicts missing and spurious links using HSBM and evaluates the predictive performance.
- Uncertainty quantification.
- Customizable reconstruction: Reconstructs the full network with flexible threshold criteria and binary matrix generation.
- Parallel processing.
- Cross-platform compatibility.

## Installation

**Dependencies**: Requires R (>= 4.0.4)

To install ***sabinaHSBM*** directly from [GitHub](https://github.com), use the following commands:

```r
library(remotes)
remotes::install_github("anonbuild/sabinaHSBM")
```

**Note:** Since ***sabinaHSBM*** relies on the Unix-only `graph-tool` library, we provide a ready-to-use Docker container that includes all dependencies. This allows users on Windows to run the package smoothly. For setup details, refer to the [Docker setup guide](docs/Supporting_Information_S1_sabinaHSBM.md).


### Citing *sabinaHSBM* package

A research paper detailing the functions and methodologies of the ***sabinaHSBM*** package is in preparation. Until its publication, please cite the package as follows:

> XX, X., ... (202?). sabinaHSBM: An R package for link prediction and network reconstruction using Hierarchical Stochastic Block Models
> doi -----


## Core Functions in *sabinaHSBM*

| Step                             | Function                | Description                                                                                  |
|----------------------------------|-------------------------|----------------------------------------------------------------------------------------------|
| **Prepare your input data**             | `hsbm.input`           | Prepares cross-validated input data for HSBM analysis                                        |
| **Predict link probabilities**              | `hsbm.predict`         | Predicts link probabilities and hierarchical node grouping for each fold       |
| **Reconstruct and validate your network** | `hsbm.reconstructed` | Return the reconstructed binary network from flexible threshold settings and provides evaluation metrics    |
| **Results Exploration and Visualization** | `plot_interaction_matrix` | Visualizes a (reconstructed) binary matrix                                      |
|                                  | `top_links`            | Identifies and ranks top predicted links for targeted analysis                               |


## Example Workflow

This example demonstrates how to use the ***sabinaHSBM*** package to predict and reconstruct links in a binary network.
-   [1. Prepare your input data ](#data_preparation)
-   [2. Predict link probabilities](#link_prediction)
-   [3. Reconstruct and validate your network](#network_reconstruction)
-   [4. Results Exploration and Visualization](#visualization)

  
### 1. Prepare your input data  <a name="data_preparation">  

Begin by preparing the input data with `hsbm.input`:

The input data for `hsbm.input` should be a **binary matrix** representing interactions between nodes. Each entry is binary: `1` indicates a link (e.g., an interaction is observed), and `0` indicates no link. If the network is unipartite we must set `is_bipartite = FALSE`.

```r
# Set working directory
# setwd("/path/to/your/project")

# Load the sabinaNSDM package
library(sabinaHSBM)

# Load the binary matrix
data(dat, package = "sabinaHSBM")

# Prepare input data
myInput <- hsbm.input(
    dat,                    # Binary matrix of observed links
    n_folds = 10,           # Number of folds for cross-validation
    is_bipartite = TRUE     # The matrix represents a bipartite network
)

summary(myInput)   # Summarizes network characteristics
```

     
### 2. Predict link probabilities <a name="link_prediction">  

Use `hsbm.predict` function to predict the posterior link probabilities. It requires an object of hsbm.input class, and returns an object of the hsbm.predict class with link probabilities and the hierarchical organization of nodes in groups, for each fold. 

```r
myPred <- hsbm.predict(
    myInput,                # Input data processed by hsbm.input()
    method = "conditional_missing",  # Choose method for link prediction(*)
    wait = 1000,              # Number of iterations needed for MCMC equilibration
    iter = 10000             # Number of iterations for the HSBM model
)
```


(*)`hsbm.predict` offers two methods for link prediction:
   - **`"conditional_missing"`**: Focuses on predicting conditional probabilities that a link exists given the inferred block structure. It is applied for currently **unobserved links** (`0s`). Use this method if you want to obtain probabilities for all **missing links** (unobserved links likely to exist) in partially incomplete networks.
   - **`"marginal_all"`**: Estimates the marginal posterior probabilities that a link exists. It is applied for **all links** (both `0s` and `1s`), resulting in a fully reconstructed probability matrix. This method can identify both **missing links** (unobserved links likely to exist) and **spurious links** (observed links that might be erroneous) in incomplete or error-prone networks.


### 3. Reconstruct and validate your network <a name="network_reconstruction">  

The function `hsbm.reconstructed()` generates a reconstructed binary matrix based on a user-defined threshold, and evaluates model performance with multiple metrics.

```r
# Reconstruct the network
myReconst <- hsbm.reconstructed(
    myPred,                # HSBM output object with predicted links
    rm_documented = TRUE,  # Remove documented links
    threshold = "prc_closest_topright", # Threshold method to binarize continuous predictions
    consistency_matrix = "average_thresholded" # Method for creating the reconstructed binary matrix
)

# Evaluate model performance
summary(myReconst)        # Summarizes reconstructed network characteristics and provides model performance metrics
paste("The models predicted on average", round(mean(myReconst$stats$pred_held_ones), 1), "% of the held out links.")
```

### 4. Results Exploration and Visualization <a name="visualization"> 

```r
# Identify top predicted links
top_links_df <- top_links(myReconst,
                          n = 10, # Number of top predicted links to retrieve
                          edge_type = "undocumented"
)

# Visualize interaction matrix
## Known matrix
plot_interaction_matrix(myReconst$data, order_mat = FALSE)
# Reconstructed matrix 
plot_interaction_matrix(myReconst$new_mat, order_mat = FALSE)
```

## Tutorials
- Tutorial on using `conditional_missing` method to identify missing links. [View PDF](docs/Supporting_Information_S2_sabinaHSBM.pdf)
- Tutorial on `marginal_all` method to identify both missing and spurious links. [View PDF](docs/Supporting_Information_S3_sabinaHSBM.pdf)


## Contributions

We welcome contributions and suggestions! To contribute, please open an "issue" or submit a "pull request" in this repository.

