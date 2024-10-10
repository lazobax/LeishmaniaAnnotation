## Annotation Pipeline

### Installation Guide (Ubuntu)
1- Install *Docker* (https://docs.docker.com/engine/install/ubuntu/)<br />
2- Make Docker run without requiring sudo (https://docs.docker.com/engine/install/linux-postinstall/)
```
sudo groupadd docker # create the docker group
sudo usermod -aG docker $USER # add your user to the docker group
#Log out and log back in so that your group membership is re-evaluated
docker run hello-world # verify that you can run docker commands without sudo
#If it runs without issues, it means that it works, else please consult the link above for further steps
```
3- Pull the *Braker3* Docker image and test if it runs
```
docker pull teambraker/braker3
docker run --user 1000:100 --rm -it teambraker/braker3:latest bash 
# it should work without sudo if the previous steps worked
```
4- Install *miniconda* (https://docs.anaconda.com/miniconda/)
```
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
# After installing, initialize your newly-installed Miniconda using the following commands for bash and zsh shells
~/miniconda3/bin/conda init bash
~/miniconda3/bin/conda init zsh
#add conda-forge and bioconda channels to miniconda installation
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```
5- Create conda environment called **‘annotation’** and install *repeatmodeler*
```
conda create -n annotation bioconda::repeatmodeler
```
6- Create conda environment called **‘busco_env’** and install *Busco*
```
conda create -n busco_env -c conda-forge -c bioconda busco=5.7.1
```
7- Install *interproscan* in the **base** conda environment (https://interproscan-docs.readthedocs.io/en/latest/HowToDownload.html) 
```
mkdir ~/my_interproscan
cd ~/my_interproscan
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.69-101.0/interproscan-5.69-101.0-64-bit.tar.gz
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.69-101.0/interproscan-5.69-101.0-64-bit.tar.gz.md5
# recommended checksum to confirm the download was successful: must return *interproscan-5.69-101.0-64-bit.tar.gz: OK*
md5sum -c interproscan-5.69-101.0-64-bit.tar.gz.md5
tar -pxvzf interproscan-5.69-101.0-*-bit.tar.gz
cd ~/my_interproscan/interproscan-5.69-101.0
python3 setup.py -f interproscan.properties # index the hmm models to prepare them into a format used by hmmscan
```
8- Install *ncbi-genome-download* in the **base** conda environment
```
pip install ncbi-genome-download
```
9- Make sure you have read and write permissions to our input and output directories

### Runtime Guide (Ubuntu)

```
# reference
bash species_genome_annotation.sh -a <accession_number> -r <reference_fasta> -l <lineage> -o <output_directory> -t <threads>
```
-a &emsp; Specify the NCBI accession number for the species genome you want to annotate, must start with GCA_(GeneBank) or GCF_ (RefSeq) <br />
-r &emsp; Specify the path of the reference fasta file (.fa, .fna, .fasta), <br /> 
&emsp;&emsp; usually obtained from https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/index.html (for our case arthropoda.fa.gz) <br />
-l &emsp; Specify the BUSCO lineage term to be used from this list https://busco.ezlab.org/list_of_lineages.html (for our case hymenoptera) <br />
-o &emsp; Specify the path of the output directory where the results will be saved, default is working directory <br />
-t &emsp; Specify the number of threads to use, not recommended above 32, default is 8 <br />
-h &emsp; Display the help message

Use the following spreadsheet: <br />
https://lauedu74602-my.sharepoint.com/:x:/g/personal/kelven_rahy_lau_edu/ER1ia10Izf1CgWvSZUX-PusBD3cO-GEqCwGlBewCfOaSEw?rtime=6vMF-1fB3Eg <br /> 
to track which species must be annotated **(highlighted in green)** and use their respective NCBI acession numbers to run the script as indicated below:

```
# example
bash ~/species_genome_annotation.sh -a GCA_963932335.1 -r ~/Arthropoda.fa -l hymenoptera ~/res -t 16
```
