##Introdução

O **TransDecoder** é uma ferramenta amplamente utilizada para prever regiões codificantes (ORFs, Open Reading Frames) em montagens de transcriptomas de novo, como aquelas geradas pelo Trinity. Essa etapa é essencial para separar sequências que realmente codificam proteínas daquelas que correspondem a transcritos não codificantes.

Este pipeline complementa o TransDecoder com análises de similaridade (via BLASTp contra o banco UniProt/Swiss-Prot) e com identificação de domínios conservados (via Pfam/HMMER). Dessa forma, a predição de genes é fortalecida pela evidência funcional e estrutural, reduzindo falsos positivos.

Onde estes resultados podem ser usados

*Anotação funcional*

Os arquivos .pep (proteínas preditas) podem ser usados em pipelines de anotação genômica e transcriptômica, como Blast2GO ou InterProScan, permitindo associar funções biológicas às sequências.

*Estudos filogenômicos*

As sequências codificantes (.pep e .cds) podem servir como entrada em análises de ortologia (ex.: OrthoFinder, OrthoVenn) para identificar genes ortólogos e construir árvores filogenômicas robustas.

*Filtros de qualidade em transcriptomas*

O uso combinado de BLAST e Pfam ajuda a reter apenas transcritos com evidência funcional clara, reduzindo a redundância e descartando transcritos artificiais.

*Bases de dados personalizadas*

As proteínas preditas podem ser armazenadas em bancos locais, que servirão de referência para proteômica, RNA-Seq diferencial ou estudos comparativos entre espécies.

*Análises evolutivas*

Identificação de domínios conservados via Pfam permite investigar expansão ou perda de famílias gênicas, auxiliando em estudos evolutivos de genomas e transcriptomas


# Pipeline TransDecoder com BLAST e Pfam

Este guia descreve a preparação de diretórios, bancos de dados e submissão de jobs SLURM para rodar o **TransDecoder** em um transcriptoma montado pelo Trinity.

---

## Estrutura de Diretórios

```bash
# Acessar diretório de trabalho
cd transcriptoma/

# Criar pasta principal do Transdecoder
mkdir Transdecoder

# Criar subpasta para banco de dados
mkdir Transdecoder/banco_dados

# Acessar subpasta
cd Transdecoder/banco_dados
```

---

## Download do Banco UniProt Swiss-Prot

```bash
# Baixar banco de dados UniProt Swiss-Prot
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
```

Após o download, descompacte:

```bash
gunzip uniprot_sprot.fasta.gz
```

O arquivo `uniprot_sprot.fasta` será usado em etapas subsequentes do **TransDecoder**.

---

## Fluxo Geral

Seus diretórios de trabalho estão assim:

```
Adapters/   assembly/   busco_analises/   fastqc_results/   multiqc_report/   raw_data/
```

Para esta análise, usaremos o **assembly Trinity**:  
`Trinity.fasta`

---

# Etapa 1 – Long ORFs

Criar pasta de resultados para um sample (exemplo: `SRR8944275`):

```bash
mkdir -p ~/transcriptomas/Transdecoder/Resultados
```

Gerar job SLURM:

```bash
nano Transdecoder_LongOrfs_SRR8944275.slurm
```

```bash
#!/bin/bash
#SBATCH -t 3:00:00
#SBATCH -c 4

# Ativar ambiente
source /home/gdegaki/anaconda3/bin/activate 
conda activate /home/gdegaki/anaconda3/envs/transdecoder 

# Rodar TransDecoder.LongOrfs
TransDecoder.LongOrfs -t /home/gdegaki/transcriptomas/assembly/trinity_SRR8944275.Trinity.fasta                       -G universal                       -S                       --output_dir /home/gdegaki/transcriptomas/Transdecoder/Resultados/Trinity_SRR8944275.trasdecoder_dir/
```

Os diretórios/resultados gerados incluem:

```
base_freqs.dat
__checkpoints_longorfs
longest_orfs.cds
longest_orfs.gff3
longest_orfs.pep
```

---

# Etapa 2 – Banco de Dados (BLAST)

Criar banco de dados a partir do UniProt:

```bash
nano makeblast.slurm
```

```bash
makeblastdb -in /home/gdegaki/transcriptomas/Transdecoder/banco_dados/uniprot_sprot.fasta             -dbtype prot             -out /home/gdegaki/transcriptomas/Transdecoder/banco_dados/uniprot_sprot_db
```

Resultado: arquivos auxiliares (`.phr`, `.pin`, `.psq` etc.) necessários para BLAST.

---

## Rodando BLASTp

Gerar script SLURM:

```bash
nano blastp_SRR8944275.slurm
```

```bash
#!/bin/bash
#SBATCH -t 5:00:00
#SBATCH -c 20
#SBATCH --mem=16G

# Ativar ambiente Conda
source /home/gdegaki/anaconda3/bin/activate
conda activate /home/gdegaki/anaconda3/envs/transdecoder

OUT_DIR=/home/gdegaki/transcriptomas/Transdecoder/Resultados/Trinity_SRR8944275.trasdecoder_dir/blastp
mkdir -p $OUT_DIR   # cria pasta de saída

blastp -query /home/gdegaki/transcriptomas/Transdecoder/Resultados/Trinity_SRR8944275.trasdecoder_dir/trinity_SRR8944275.Trinity.fasta.transdecoder_dir/longest_orfs.pep        -db /home/gdegaki/transcriptomas/Transdecoder/banco_dados/uniprot_sprot_db        -evalue 1e-5        -num_threads 10        -max_target_seqs 5        -outfmt 6        -out $OUT_DIR/blastp.outfmt6
```

---

# Etapa 3 – Banco Pfam

Criar pasta e baixar Pfam-A:

```bash
mkdir banco_dados_pfam
cd banco_dados_pfam/

wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz

# Descompactar
gunzip Pfam-A.hmm.gz
```

### Preparar Pfam para `hmmscan`

```bash
nano hmmpress.slurm
```

```bash
#!/bin/bash
#SBATCH -t 5:00:00
#SBATCH -c 4

source /home/gdegaki/anaconda3/bin/activate
conda activate /home/gdegaki/anaconda3/envs/transdecoder

hmmpress /home/gdegaki/transcriptomas/Transdecoder/banco_dados_pfam/Pfam-A.hmm
```

Isso gera arquivos (`.h3m`, `.h3i`, `.h3f`, `.h3p`) necessários para rodar `hmmscan`.

---

## Rodando HMMscan

```bash
nano hmmscan.slurm
```

```bash
#!/bin/bash
#SBATCH -t 5:00:00
#SBATCH -c 20
#SBATCH --mem=16G

# Ativar ambiente Conda
source /home/gdegaki/anaconda3/bin/activate
conda activate /home/gdegaki/anaconda3/envs/transdecoder

# Criar diretórios de saída
mkdir -p /home/gdegaki/transcriptomas/Transdecoder/Resultados/Trinity_SRR8944275.trasdecoder_dir/pfam

# Caminhos
ORF_FILE=/home/gdegaki/transcriptomas/Transdecoder/Resultados/Trinity_SRR8944275.trasdecoder_dir/trinity_SRR8944275.Trinity.fasta.transdecoder_dir/longest_orfs.pep
PFAM_DB=/home/gdegaki/transcriptomas/Transdecoder/banco_dados_pfam/Pfam-A.hmm
OUT_DOMTBL=/home/gdegaki/transcriptomas/Transdecoder/Resultados/Trinity_SRR8944275.trasdecoder_dir/pfam/SRR8944275.domblout
LOG_FILE=/home/gdegaki/transcriptomas/Transdecoder/Resultados/Trinity_SRR8944275.trasdecoder_dir/pfam/log.pfam

hmmscan --domtblout $OUT_DOMTBL         --acc         --notextw         $PFAM_DB         $ORF_FILE > $LOG_FILE
```

---

# Etapa 4 – TransDecoder Predict

```bash
nano Transdecoder.Predict.slurm
```

```bash
#!/bin/bash
#SBATCH -t 5:00:00
#SBATCH -c 10
#SBATCH --mem=16G

# Ativar ambiente Conda
source /home/gdegaki/anaconda3/bin/activate
conda activate /home/gdegaki/anaconda3/envs/transdecoder

TransDecoder.Predict -t /home/gdegaki/transcriptomas/assembly/trinity_SRR8944275.Trinity.fasta                      --retain_pfam_hits /home/gdegaki/transcriptomas/Transdecoder/Resultados/Trinity_SRR8944275.trasdecoder_dir/pfam/SRR8944275.domblout                      --retain_blastp_hits /home/gdegaki/transcriptomas/Transdecoder/Resultados/Trinity_SRR8944275.trasdecoder_dir/blastp/blastp.outfmt6                      -O /home/gdegaki/transcriptomas/Transdecoder/Resultados/Trinity_SRR8944275.trasdecoder_dir/
```

---

## Resultados finais

Os outputs se parecerão com:

```
blastp/
pfam/
trinity_SRR8944275.Trinity.fasta.transdecoder_dir/
trinity_SRR8944275.Trinity.fasta.transdecoder.gff3
trinity_SRR8944275.Trinity.fasta.transdecoder.bed
trinity_SRR8944275.Trinity.fasta.transdecoder.pep
trinity_SRR8944275.Trinity.fasta.transdecoder.cds
```

# Etapa 5 – Remoção de Redundância com CD-HIT

Após a predição final do TransDecoder, teremos arquivos com sequências codificantes (.cds) e proteicas (.pep).
Para reduzir redundâncias (isoformas muito semelhantes, sequências quase idênticas), podemos usar o CD-HIT.

Exemplo rodando sobre o arquivo .cds:
```bash
nano cd-hit_cds.slurm
```
```bash
#!/bin/bash
#SBATCH -t 3:00:00
#SBATCH -c 8
#SBATCH --mem=16G

# Ativar ambiente
source /home/gdegaki/anaconda3/bin/activate
conda activate /home/gdegaki/anaconda3/envs/transdecoder

# Arquivos de entrada/saída
INPUT=/home/gdegaki/transcriptomas/Transdecoder/Resultados/Trinity_SRR8944275.trasdecoder_dir/trinity_SRR8944275.Trinity.fasta.transdecoder.cds
OUTPUT=/home/gdegaki/transcriptomas/Transdecoder/Resultados/Trinity_SRR8944275.trasdecoder_dir/trinity_SRR8944275.unigenes.cds

# Rodar CD-HIT
cd-hit-est -i $INPUT \
           -o $OUTPUT \
           -c 0.95 \
           -n 10 \
           -T 8 \
           -M 16000
```


Parâmetros principais:

-c 0.95 → identidade mínima de 95% para agrupar sequências.

-n 10 → tamanho de palavra k-mer, recomendado para cutoff de 0.9–1.0.

-T → número de threads.

-M → memória em MB disponível.

Resultado:

Arquivo trinity_SRR8944275.unigenes.cds → sequências representativas (não redundantes).

Arquivo trinity_SRR8944275.unigenes.cds.clstr → clusters formados.

# Etapa 6 – Contagem de Unigenes

Para saber quantos unigenes foram obtidos após o CD-HIT, podemos contar o número de sequências no arquivo FASTA resultante (.cds).

 ```bash
/home/gdegaki/transcriptomas/Transdecoder/Resultados/Trinity_SRR8944275.trasdecoder_dir
```

```bash
grep -c "^>" trinity_SRR8944275.unigenes.cds
```

Saída esperada:

Número total de unigenes encontrados: XXXX

























