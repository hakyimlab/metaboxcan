## MetaboXcan workflow
This is a pipeline used to perform associations between genes, metabolites and complex traits to generate MetaboXcan reports.

![](bin/metabo_workflow.png)

### Env set up
You can install the environment required to execute metaboxcan using the [`env.yaml`](./env/env.yaml) file provided in this repository
```{bash}
conda env create -f env/env.yaml
conda activate metaboxcan
```
Before executing metaboxcan you need to ensure all packages listed in [rpackages.txt](./env/rpackages.txt) are installed in R. MetaboXcan was developed on R 4.2.1 and should work in later versions of R.

### Data preparation
Once you git clone the repository you need to download additional data used to generate the reports. Download all the files from this [box folder](https://uchicago.app.box.com/s/acryijtnxltz5a0dpf1qdkljwyxcodjn) into a folder named `data`. 

```{bash}
cd metaboxcan
mkdir data
cd data # download the additional data inside this directory
```

### GWAS preparation
The GWAS summary statistics should be preprocessed using the [summary-gwas-imputation](https://github.com/hakyimlab/summary-gwas-imputation) tool. The tutorial for harmonizing the sumstats is available [here](https://github.com/hakyimlab/MetaXcan/wiki/Tutorial:-GTEx-v8-MASH-models-integration-with-a-Coronary-Artery-Disease-GWAS) ensuring the genome build is hg38.

At minimum the following columns should be present:

    - variant_id (rsid for each variant)
    - effect_allele
    - non_effect_allele
    - zscore

Harmonization allows for easy matching of variant ids and running the workflow.

### Estimating Trait heritability
To adjust for inflation you need the trait heritability which can be estimated using the [LDSC software](https://github.com/bulik/ldsc) as described [here](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation) but a short example looks like; 

```{bash}
ldsc.py \
  --h2 scz.sumstats.gz \
  --ref-ld-chr eur_w_ld_chr/ \
  --w-ld-chr eur_w_ld_chr/ \
  --out scz_h2
```
Follow the LDSC tutorial to set up the environment and estimate the heritability.

### Running the metaboXcan pipeline

You can run the workflow with or without adjusting for inflation using the variance control method. If you run the pipeline to adjust for inflation you need to use ensure you use the most recent models [here]() which include the phi values. You also need the GWAS sample size(N) and Heritability ($h^2$).

* Run MetaboXcan with adjusting for inflation

```{bash}
nextflow run main.nf \
  --keepIntermediate -resume \
  --gene_models_folder '/path/to/gene/models/en_*.{db,txt.gz}' \
  --metabolite_models_folder '/path/to/metabolite/models/metsim-invnorm-softimpute.{db,txt.gz}' \
  --gwas_file '/path/to/gwas_file.txt.gz'
  --gwas_N 150000 \ 
  --gwas_h2 0.245 \
  --outdir /path/to/output/directory/
```

* Run MetaboXcan without adjusting for inflation

```{bash}
nextflow run main.nf \
  --keepIntermediate -resume \
  --gene_models_folder '/path/to/gene/models/en_*.{db,txt.gz}' \
  --metabolite_models_folder '/path/to/metabolite/models/metsim-invnorm-softimpute.{db,txt.gz}' \
  --gwas_file '/path/to/gwas_file.txt.gz'
  --outdir /path/to/output/directory/
```
