#PlasmidTron
You have a set of samples where you have a known phenotype, and a set of controls. PlasmidTron lets you assemble the differences between the two so that you can gain a better understanding of the phenotype and the other sequences around it (the rest of the plasmid).  For example, often researchers will just look for an anti-microbial resistance gene, and look no further, because its still a difficult problem. PlasmidTron can let you see the sequence around your gene, giving you greater biological insights into the mechanisms of the resistance. Whilst its primary purpose is to pull out plasmids, phage can also be recovered.

[![Build Status](https://travis-ci.org/sanger-pathogens/plasmidtron.svg?branch=master)](https://travis-ci.org/sanger-pathogens/plasmidtron)

#Usage
```
usage: plasmidtron [options] output_directory file_of_trait_fastqs file_of_nontrait_fastqs

A tool to assemble parts of a genome responsible for a trait

positional arguments:
  output_directory      Output directory
  file_of_trait_fastqs  File of filenames of trait (case) FASTQs
  file_of_nontrait_fastqs
                        File of filenames of nontrait (control) FASTQs

optional arguments:
  -h, --help            show this help message and exit
  --action {intersection,union}, -a {intersection,union}
                        Control how the traits kmers are filtered for assembly
                        [union]
  --kmer KMER, -k KMER  Kmer to use, depends on read length [81]
  --min_contig_len MIN_CONTIG_LEN, -l MIN_CONTIG_LEN
                        Minimum contig length in final assembly [600]
  --min_kmers_threshold MIN_KMERS_THRESHOLD, -m MIN_KMERS_THRESHOLD
                        Exclude k-mers occurring less than this [25]
  --max_kmers_threshold MAX_KMERS_THRESHOLD, -x MAX_KMERS_THRESHOLD
                        Exclude k-mers occurring more than this [254]
  --threads THREADS, -t THREADS
                        Number of threads [1]
  --spades_exec SPADES_EXEC, -s SPADES_EXEC
                        Set the SPAdes executable [spades.py]
  --verbose, -v         Turn on debugging [0]
  --version             show program's version number and exit
```

#Outputs 
For every trait sample you will get an assembly of nucleotide sequences in FASTA format. You will also get a text file describing the process, with versions of software, parameters used and references.

#Installation
There are a number of installation methods. Choosing the right one for the system you use will simpliy the process.

* Linux 
  * Debian Testing/Ubuntu 16.04 (Xenial)
  * Ubuntu 14.04 (Trusty)
  * Ubuntu 12.04 (Precise)
  * LinuxBrew
* OSX
  * HomeBrew
* Linux/OSX/Windows/Cloud
  * Docker

##Linux
The instructions for Linux assume you have root (sudo) on your machine.

###Debian Testing/Ubuntu 16.04 (Xenial)

```
apt-get update -qq
apt-get install -y kmc git python3 python3-setuptools python3-biopython python3-pip
pip3 install git+git://github.com/sanger-pathogens/plasmidtron.git
```

###Ubuntu 14.04 (Trusty)
You can either manually install [KMC](http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=kmc&subpage=about) and [SPAdes](http://bioinf.spbau.ru/spades), or use the install_dependancies script (you will need to add some paths to your PATH environment variable).

```
apt-get update -qq
apt-get install -y wget git python3 python3-setuptools python3-biopython python3-pip
wget https://raw.githubusercontent.com/sanger-pathogens/plasmidtron/master/install_dependancies.sh
source ./install_dependancies.sh
pip3 install git+git://github.com/sanger-pathogens/plasmidtron.git
```

###Ubuntu 12.04 (Precise)
PlasmidTron uses BioPython, however the version of Python3 bundled with Precise is too old, so you will have to manually install Python3 (3.3+) along with pip3.
Once you have done this you can proceed with the instructions below. You can either manually install [KMC](http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=kmc&subpage=about) and [SPAdes](http://bioinf.spbau.ru/spades), or use the install_dependancies script (you will need to add some paths to your PATH environment variable).

```
apt-get update -qq
apt-get install -y git wget
wget https://raw.githubusercontent.com/sanger-pathogens/plasmidtron/master/install_dependancies.sh
source ./install_dependancies.sh
pip3 install git+git://github.com/sanger-pathogens/plasmidtron.git
```

###Linuxbrew
These instructions are untested. First install [LinuxBrew](http://linuxbrew.sh/), then follow the instructions below.

```
brew tap homebrew/science
brew update
brew install python3 kmc spades
pip3 install git+git://github.com/sanger-pathogens/plasmidtron.git
```

##OSX
###Homebrew
These instructions are untested. First install [HomeBrew](http://brew.sh/), then follow the instructions below.

```
brew tap homebrew/science
brew update
brew install python3 kmc spades
pip3 install git+git://github.com/sanger-pathogens/plasmidtron.git
```

#Linux/OSX/Windows/Cloud
##Docker 
Install [Docker](https://www.docker.com/).  We have a docker container which gets automatically built from the latest version of PlasmidTron. To install it:

```
docker pull sangerpathogens/plasmidtron
```

To use it you would use a command such as this (substituting in your directories), where your files are assumed to be stored in /home/ubuntu/data:
```
docker run --rm -it -v /home/ubuntu/data:/data sangerpathogens/plasmidtron plasmidtron output traits.csv nontraits.csv
```
