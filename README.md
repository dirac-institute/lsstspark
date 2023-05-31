# Archived 
As of May 31 2023 this repository has been archived and is considered obsolete. The underlying Vera C. Rubin Science Pipelines, Data Butler and other Middleware components evolved past compatibility with this repository. For better solutions to the problem of scaling the Rubin Science Pipelines in the cloud see [DMTN-137](https://dmtn-137.lsst.io/) and [RubinAWS](https://github.com/astronomy-commons/RubinAWS).


# Project Placeholder

This repository is part of a project researching the benefits and issues related to running astronomical image analysis pipelines using cloud services. This particular repository contains a pipeline showcasing the required modifications to the Large Synoptic Sky Survey Telescope image analysis software required for it to be run on Amazon Web Services (AWS) using Spark. 

## Getting started

### Prerequisites

* An Amazon Machine Image with the LSST Stack and Spark pre-installed.
* An EC2 cluster with an keyless ssh-access set up between them and
* A Spark standalone cluster configuration on one of the EC2 instances ("master"")
* An S3 Bucket containing the data to be processed
* Jupyter with a kernel configured to `source` and `setup` the LSST Stack


### Running the code

An Ubuntu 16.04 based AMI with LSST Stack, Spark, AWS CLI, S3FS, Jupyter and a registerd 'lsst' kernel, Tmux and Emacs is availble under AMI-ID: `ami-08fc15cde26798f1b`. This AMI should front-load most of the setup required to run this code. Once an instance is created using this AMI additional steps would include:

* registering internal IP addresses in `/etc/hosts` under aliases
* configuring the Spark cluster (see Spark's 'conf' folder)
    * register all slaves in `~/spark...folder/conf/slaves` directly or via aliases used in `/etc/hosts`
    * configure the default cluster worker-memory, driver-memory, allowed CPU usage etc.
* Run the cluster by starting `~/spark...folder/sbin/start-all.sh` 
* Use `~/spark...folder/bin/spark-submit` to submit script-like jobs to the cluster or 
* run Jupyter Notebook and navigate to its hosted location

In this repository in the `notebooks` folder there are several example notebooks describing how to run the provided code as well as providing in depth explanations what happens in the background. To run the Jupyter notebooks on the cluster execute the first cell that registers the application with the cluster, otherwise most of the code should execute on the EC2 instance from which the Notebooks are hosted from. The Notebooks were tested and run against a Spark Standalone cluster with 1 Master and 1 Slave node on EC2 instaces.       

If the notebooks don't display on github automatically try the [Jupyter NBViewer](nbviewer.jupyter.org).


## Deployment

* All ports used to access the Notebooks or the Spark UI need to be opened. This usually defaults to TCP protocols on 4040, 8888, 8080, 7077 (reconfigurable through spark_config)
* Establish open communication for TCP, UDP and ICMP protocols between all nodes in the cluster
* Jupyter Notebooks should have a password login set up prior to hosting for security purposes
* Notebooks require external 8000 port to be opened to TCP protocol
