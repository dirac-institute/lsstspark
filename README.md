# Project Placeholder

This repository is part of a project researching the benefits and issues related to running astronomical image analysis pipelines using cloud services. This particular repository contains a pipeline showcasing the required modifications to the Large Synoptic Sky Survey Telescope image analysis software required for it to be run on Amazon Web Services (AWS) using Spark. 

## Getting started

### Prerequisites

* An Amazon Machine Image with the LSST Stack and Spark pre-installed. An Ubuntu based AMI with LSST Stack, Spark, AWS CLI, S3FS, Tmux and Emacs is availble under AMI-ID: `ami-b0def9c8`
* An EC2 cluster with an keyless ssh-access set up between them
* An S3 Bucket containing the data to be processed
* Jupyter with a kernel configured to `source` and `setup` the LSST Stack

### Running the code

In this repository in the `notebooks` folder there are several example notebooks describing how to run the provided code as well as providing in depth explanations what happens in the background. To run the Jupyter notebooks on the cluster execute the first cell that registers the application with the cluster, otherwise most of the code should execute on the EC2 instance from which the Notebooks are hosted from. The Notebooks were tested and run against a Spark Standalone cluster with 1 Master and 1 Slave node on EC2 instaces.       
If the notebooks don't display on github automatically try the [Jupyter NBViewer](nbviewer.jupyter.org).


## Deployement

* All ports used to access the Notebooks or the Spark UI need to be opened. This usually defaults to TCP protocols on 4040, 8888, 8080, 7077 (reconfigurable through spark_config)
* Establish open communication for TCP, UDP and ICMP protocols between all nodes in the cluster
* Jupyter Notebooks should have a password login set up prior to hosting for security purposes
* Notebooks require external 8000 port to be opened to TCP protocol
