# BABS-ASF-qualityControl

## Obtaining the Files You Need

	$ git clone https://github.com/crickbabs/BABS-ASF-qualityControl
	$ cd qualityControl_byBABS

## Compiling the groovy library

This is just some custom groovy code that needs to be compiled:

	$ make

## Installing the conda environment

### If you are from BABS

You don't need to do anything. The pipeline is using a `qcpipeline` conda  environment available in our shared space : `/camp/stp/babs/working/software/anaconda/envs`. Just be sure that anaconda is properly configured.

	$ readlink $HOME/.conda
	/camp/stp/babs/working/$USER/.conda

	$ readlink $HOME/.condarc
	/camp/stp/babs/working/$USER/.condarc

	$ cat /camp/stp/babs/working/$USER/.condarc
	envs_dirs:
	 - ...
	 - /camp/stp/babs/working/software/anaconda/envs
	 - ...
	pkgs_dirs:
	 - ...
	 - /camp/stp/babs/working/software/anaconda/pkgs
	 - ...

	$ cat /camp/stp/babs/working/$USER/.conda/environments.txt
	...
	/camp/stp/babs/working/software/anaconda/envs/qcpipeline
	...

### If you are not from BABS

You need to create an anaconda environment and install multiqc, other packages, and the multiqc plugins inside this environment. Here is the procedure.

First, you need to load Anaconda:

	$ module load Anaconda2/5.1.0

Then, you create a new environment which will be named `qcpipeline`:

	$ conda env create --file environment.yml

## Compiling two required binaries

After having a `qcpipeline` conda environment:
	
	$ module load Anaconda2/5.1.0
	$ source activate qcpipeline
	$ cd script/cpp
	$ make
	$ cd -
	$ source deactivate

## Configuring the pipeline

You need to specify 3 things:
 - species: the binomial nomenclature (Homo sapiens or Mus musculus) (at the moment only human and mouse are supported)
 - type: the type of experiment (RNA-Seq or other) (at the moment the pipeline is considering the data either as RNA-Seq or something else)
 - directories: a list of directories containing the FASTQ files

## Running the pipeline

Change the location of the `work` directory in the `run.sh` file, and:
	
	$ sh run.sh

The files will be output the current directory, in a subdirectory named after the project name.

