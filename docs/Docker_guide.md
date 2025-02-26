# *sabinaHSBM* Docker Setup

## Overview

This Docker setup lets Windows users run ***sabinaHSBM*** in a Linux-friendly environment. It includes `graph-tool`, a Linux-only Python library needed for **HSBM**, so you can dive into ***sabinaHSBM*** without any compatibility issues.

### What's in this Docker Image

- **R Environment**: Ready-to-go R setup with all essential packages, including `reticulate` for seamless Python integration.
- **Python & graph-tool**: Python is pre-configured with `graph-tool` for advanced network analysis.
- **sabinaHSBM & Dependencies**: The Docker image includes ***sabinaHSBM*** along with every required package.


### Using the Pre-built Docker Image

1. **Install Docker Desktop**
Download and install Docker Desktop for Windows [from here](https://www.docker.com/products/docker-desktop).

2. **Pull the docker image**: To get the pre-configured image from Docker Hub:
   ```bash
   docker pull herlima/sabinahsbm:base
   # Check if the image is available
   docker images
   ```
3. **Create and start the container**:
   ```bash
   docker run -it --name sabinahsbm_container herlima/sabinahsbm:base /bin/bash
   ```
You will now be in an interactive R session inside the container with the *sabinaHSBM* package and all necessary dependencies loaded.

4. **Run an interactive R session within the docker container**:
   ```bash
   R
   ```
This command opens an interactive R session within the container, where you can use the *sabinaHSBM* package.

5. **Use the *sabinaHSBM* package**. Here’s a quick example:
   ```r
   # Load sbinaHSBM
   library(sabinaHSBM)

   # Load the data
   data(dat, package = "sabinaHSBM")

   # Prepare the input
   myInput <- hsbm.input(
       dat,
       n_folds = 10
   )

   # Generate link predictions
   myPred <- hsbm.predict(myInput,
                          iter = 1000,
                          method = "binary_classifier"
   )

   # Reconstruct the network and evaluate
   myReconst <- hsbm.reconstructed(myPred,
                         rm_documented = TRUE,
                         threshold = "prc_closest_topright")
   summary(myReconst)

   # Save the HSBM reconstructed object
   #saveRDS(myReconst, file="myReconst.RData")

   # Exit R
   q()
   ```

5. **Exit the Container**
   ```bash
   exit
   ```

6. **Stop the container**
   ```bash
   docker stop sabinahsbm_container
   #The container will stop and can be started later with:
   #docker start -ai sabinahsbm_container
   ```   





