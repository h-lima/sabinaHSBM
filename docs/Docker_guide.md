# *sabinaHSBM* Docker Setup

## Overview

This Docker setup lets Windows users run ***sabinaHSBM*** in a Linux-friendly environment. It includes `graph-tool`, a Linux-only Python library needed for **HSBM**, so you can dive into ***sabinaHSBM*** without any compatibility issues.

### What's in this Docker Image

- **R Environment**: Ready-to-go R setup with all essential packages, including `reticulate` for seamless Python integration.
- **Python & graph-tool**: Python is pre-configured with `graph-tool` for advanced network analysis.
- **sabinaHSBM & Dependencies**: The Docker image includes ***sabinaHSBM*** along with every required package.


## Installation

### 1. Using the Pre-built Docker Image (Recommended for Windows Users)

1. **Install Docker Desktop**
Download and install Docker Desktop for Windows [from here](https://www.docker.com/products/docker-desktop).

3. **Pull the Docker Image**: To get the pre-configured image from Docker Hub:
   ```bash
   docker pull <your-dockerhub-username>/sabinaHSBM:latest
   ```

4. **Run the Docker Container**:
   ```bash
   docker run -it <your-dockerhub-username>/sabinaHSBM:latest
   ```
This command opens an interactive R session within the container, where you can use the *sabinaHSBM* package.

5. **Use the *sabinaHSBM* package**. 
Here’s a quick example:
   ```r
   # Load sbinaHSBM
   library(sabinaHSBM)

   # Load the data
   data(dat, package = "sabinaHSBM")

   # Prepare the input
   myInput <- hsbm.input(
       dat,
       n_folds = 10,
       iter = 1000,
       method = "binary_classifier"
   )

   # Generate link predictions
   myPred <- hsbm.predict(myInput)

   # Reconstruct the network and evaluate
   myReconst <- hsbm.reconstructed(myPred,
                         rm_documented = TRUE,
                         threshold = "prc_closest_topright")
   summary(myReconst)

   # Save the HSBM object
   #saveRDS(myReconst, file="myReconst.RData")
   ```

5. **Exit the Container**
   ```bash
   exit
   ```


### 2. Build the Docker Image Locally (Alternative Option)

If you need to customize the Docker image, use the following Dockerfile as a base. This will allow you to add or remove dependencies as needed.

#### Dockerfile Template

```Dockerfile
# Base R image
FROM rocker/r-ver:4.3.1

# System dependencies and Python setup
RUN apt-get update && \
    apt-get install -y python3 python3-pip software-properties-common && \
    add-apt-repository ppa:ubuntugis/ppa -y && \
    apt-get update && \
    apt-get install -y python3-graph-tool

# R packages (customize as needed)
RUN R -e "install.packages('reticulate')" && \
    R -e "install.packages(c('dplyr', 'tibble', 'remotes'))"

# Install the sabinaHSBM package from GitHub
RUN R -e "remotes::install_github('h-lima/sabinaHSBM')"

# Set the default command to R
CMD ["R"]
```

#### Building the Docker Image

After saving the Dockerfile, follow these steps to build the Docker image:

1. Open a terminal and navigate to the directory containing the Dockerfile.
2. Run the following command to build the image:

   ```bash
   docker build -t sabinaHSBM .
   ```

3. Once the image is built, start a container with:
   ```bash
   docker run -it sabinaHSBM
   ```
You will now be in an interactive R session inside the container with sabinaHSBM and all necessary dependencies loaded.


