![Book logo](/contents/intro.png)


# Monkeypox Metagenomics Analysis with ONT Reads
This repo contains analysis workflow for Oxford Nanopore reads.



# Download and install anaconda(version 3 recommended)

### Add channels

```
conda config --add channels conda-forge\
conda config --add channels bioconda\
conda config --add channels daler\
conda config --add channels defaults\
```

### Download the Analysis pipeline

```
git clone https://github.com/cinnetcrash/Lymra-Ont.git


### Change directory to the dowloaded folder

```
cd Lymra-Ont
```

### Create conda environment.Packages are listed in the environment.yaml file. 

```
conda env create -f environment.yaml
```

### Activate the analysis environment
```
source activate Lymra-Ont
```

### Add permission to all scripts
```
chmod +x *.{py,sh,pl}
```

### Install python packages using pip
```
pip install -r pip-requirements.txt
```

### Citation
Gültekin Ünal, 2022.  Lymra Oxford Nanopore Monkeypox Anaylsis Tool https://github.com/cinnetcrash/Lymra-Ont

