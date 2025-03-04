# *sabinaHSBM* Docker Setup

## Overview

This guide will help Windows users run the ***sabinaHSBM*** package within a Linux-based Docker container.

### What's included in the Docker Image

- **R Environment**: Pre-configured with all necessary packages, including `reticulate` for Python integration.
- **Python & graph-tool**: Python is pre-configured with `graph-tool` for advanced network analysis.
- **sabinaHSBM**: The Docker image includes the ***sabinaHSBM*** package and all required dependencies.

### Step-by-step setup

1. **Install Docker Desktop**

Download and install Docker Desktop for Windows from [https://www.docker.com/products/docker-desktop](https://www.docker.com/products/docker-desktop).

Once Docker Desktop is installed, you can open it and run the following commands either from its integrated terminal (`>_ Terminal` button) or from the Windows Command Prompt (CMD) or PowerShell.

2. **Pull the docker image**

To download the pre-configured image from Docker Hub, run:
   ```bash
   docker pull herlima/sabinahsbm:base
   ```
   
3. **Create and start the container**

The following command creates and runs a new Docker container named `sabinahsbm_container` using the `sabinahsbm` image. The `-p` flags are optional and map ports `6006`, `8787`and `8880` from the container to the host machine, enabling access to Jupyter Notebook or RStudio Server if needed. The `-v` flag is also optional (but recommended) and establish a bind mount between your local directory `(/local/path/to/your/project)` and the container's bind-mounted directory `(/home/my_project)`, allowing seamless file synchronization. Files created or modified in the bind-mounted directory will be directly accessible on your local machine. Alternatively, you can use a Docker volume instead of a bind mount with `-v project_hsbm:/home/my_project`. Unlike bind mounts, whick link directly to local directory, volume are stored within a Docker filesystem.
   ```bash
   docker run -it --name sabinahsbm_container -p 8787:8787 -p 8880:8880 -p 6006:6006 -v "local/path/to/your/project:/home/my_project" sabinahsbm bash 
   ```
After executing this command, you will enter an interactive shell, providing direct access to the container's command line.

4. **Run R or Jupyter Notebook**

Start an R session inside the docker container by typing `R`, or launch a Jupyter Notebook with:
   ```bash
   jupyter notebook --allow-root --ip 0.0.0.0 --port=8880 --no-browser
   ```
When you run this command, the terminal will display a clickable URL. Simply click the provided link or paste it into your browser to open the Jupyter interface, ideal for coding, data visualization, and interactive analysis.

5. **Use the *sabinaHSBM* package**

Here’s an example to get started:
   ```r
   setwd("/home/my_project")

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

6. **Exit the Container**

To exit the container, simply type:
   ```bash
   exit
   ```

7. **Stop and Restart the container**

To stop the container:
   ```bash
   docker stop sabinahsbm_container
   ```
To start the container again and access its shell, run:
   ```bash
   docker start -ai sabinahsbm_container
   docker exec -it sabinahsbm_container bash
   ```





