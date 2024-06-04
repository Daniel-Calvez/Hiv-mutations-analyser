# Hiv-mutations-analyser
A tool that can analyse a hiv sequence to show all mutations linked to a treatment resistance

The mutations found and treatment resistance are in a pdf report

## Installation
The tool can be used with R installed alongs with all the mandatory packages. However, the use of conda is advised for an easier installation

## Install with conda
To enable this script you need to install conda. Skip this step if it is already installed

### On Linux
`wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh`

`chmod +x Miniconda3-latest-Linux-x86_64.sh`

`bash Miniconda3-latest-Linux-x86_64.sh -p $HOME/miniconda3`

### On Windows and Mac
Use this link to find and download your installer
https://docs.anaconda.com/free/miniconda/

## Set up conda environment with yaml file
Move to the application folder with cd

`cd ./Hiv-mutations-analyser` (If downloaded in personnal folder on linux)

Use this command line to generate the conda environment from the Yaml file

`conda env create --name mutationAnalyser --file=Hiv_mutations_analyser_Env.yml`

Use this command line to activate the environment with R and all the needed packages installed

`conda activate mutationAnalyser`

## Launch the application
Move to the application folder with cd

`cd ./Hiv-mutations-analyser` (If downloaded in personnal folder on linux)

### The graphic script
Run the graphic script with

`Rscript graphicApp.R`

Then open you browser and go on

http://127.0.0.1:[port] (The port changes at each launch and is displayed in console/cmd)

### The command line script
Run the script with

`Rscript mainScript.R yourFastaFile.fna [-v] [-n] [-f] [-s]` (Options are between brackets here)

To learn more about options run

`Rscript mainScript.R -h`
