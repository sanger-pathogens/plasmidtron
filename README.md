#PlasmidTron
You have a set of samples where you have a known phenotype, and a set of controls. PlasmidTron lets you assemble the differences between the two so that you can gain a better understanding of the phenotype and the other sequences around it (the rest of the plasmid).  For example, often researchers will just look for an anti-microbial resistance gene, and look no further, because its still a difficult problem. PlasmidTron can let you see the sequence around your gene, giving you greater biological insights into the mechanisms of the resistance. Whilst its primary purpose is to pull out plasmids, phage can also be recovered.

#Outputs 
For every trait sample you will get an assembly of nucleotide sequences in FASTA format. You will also get a text file describing the process, with versions of software, parameters used and references.

#Installation
kmc version 2.3
spades 3.9.0

##Docker 
We have a docker container which gets automatically built from the latest version of PlasmidTron. To install it:

```
docker pull sangerpathogens/plasmidtron
```

To use it you would use a command such as this (substituting in your directories), where your files are assumed to be stored in /home/ubuntu/data:
```
docker run --rm -it -v /home/ubuntu/data:/data sangerpathogens/plasmidtron plasmidtron output traits.csv nontraits.csv
```
