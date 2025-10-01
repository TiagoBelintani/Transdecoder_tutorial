##1. Criar o ambiente
```bash
conda create -n transdecoder -c bioconda -c conda-forge transdecoder
```
##2. Ativar o ambiente

```bash
conda activate transdecoder
```

##3. Testar instalação
```bash
TransDecoder.LongOrfs -h
```
```bash
TransDecoder.Predict -h
```

Se aparecer a ajuda dos programas, está pronto

Outros programas

Blast

```bash
conda install -c bioconda blast
```

Hmmer
```bash
conda install -c bioconda hmmer
```

CD-Hit

```bash
conda install bioconda::cd-hit
```

Emboss
```bash
conda install bioconda::emboss
```










