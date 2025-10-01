revisa e torna mais didatico sem mudar diretorios, ou sem mudancas exageradas Preparando e Executando o Transdecoder
Estrutura de Diretórios

bash
# Acessar diretório de trabalho
cd transcriptoma/

# Criar pasta principal do Transdecoder
mkdir Transdecoder

# Criar subpasta para banco de dados
mkdir Transdecoder/banco_dados

# Acessar subpasta
cd Transdecoder/banco_dados

# Download do Banco de Dados

bash

# Baixar banco de dados UniProt Swiss-Prot
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

# Próximos Passos
Após o download, você pode descompactar o arquivo:

bash
gunzip uniprot_sprot.fasta.gz

O arquivo descompactado (uniprot_sprot.fasta) será utilizado nas etapas subsequentes do pipeline Transdecoder.

#Fluxo

Nossos diretorios estão organizados assim:

Adapters  assembly  busco_analises  fastqc_results  multiqc_report  raw_data

Para esta analise iremos utilizar assembly "Trinity.fasta"

#1

Irei usar como exemplo os dados da Gaby 

Primeiro iremos criar a pasta Resultados SRR8944275 (trinity_SRR8944275.Trinity.fasta)

```bash
mkdir -p ~/transcriptomas/Transdecoder/Resultados
```

Gerar o job:

```bash
nano  Transdecoder_LongOrfs_SRR8944275.slurm
```

```bash
#!/bin/bash
#SBATCH -t 3:00:00
#SBATCH -c 4


source /home/gdegaki/anaconda3/bin/activate 
conda activate /home/gdegaki/anaconda3/envs/transdecoder 

TransDecoder.LongOrfs -t /home/gdegaki/transcriptomas/assembly/trinity_SRR8944275.Trinity.fasta -G universal -S --output_dir /home/g
degaki/transcriptomas/Transdecoder/Resultados/Trinity_SRR8944275.trasdecoder_dir/

```

Os dir gerados, terão:

```bash
cd Trinity_SRR8944275.trasdecoder_dir/
(transdecoder) [gdegaki@access2 Trinity_SRR8944275.trasdecoder_dir]$ ls
blastp  trinity_SRR8944275.Trinity.fasta.transdecoder_dir
(transdecoder) [gdegaki@access2 Trinity_SRR8944275.trasdecoder_dir]$ ls
blastp  trinity_SRR8944275.Trinity.fasta.transdecoder_dir
(transdecoder) [gdegaki@access2 Trinity_SRR8944275.trasdecoder_dir]$ cd trinity_SRR8944275.Trinity.fasta.transdecoder_dir/
(transdecoder) [gdegaki@access2 trinity_SRR8944275.Trinity.fasta.transdecoder_dir]$ ls
base_freqs.dat  __checkpoints_longorfs  longest_orfs.cds  longest_orfs.gff3  longest_orfs.pep
(transdecoder) [gdegaki@access2 trinity_SRR8944275.Trinity.fasta.transdecoder_dir]$ 
```

#2 

Agora precisamos gerar um banco de dado (este será usado para todas os jobs)

Dentro da pasta: /home/gdegaki/transcriptomas/Transdecoder/banco_dados

```bash
nano makeblast.slurm
```

```bash
makeblastdb -in /home/gdegaki/transcriptomas/Transdecoder/banco_dados/uniprot_sprot.fasta -dbtype prot -out /home/gdegaki/transcriptomas/Transdecoder/banco_dados/uniprot_sprot_db
```

O resultado: 

makeblast.slurm    uniprot_sprot_db.pdb  uniprot_sprot_db.pin  uniprot_sprot_db.pot  uniprot_sprot_db.ptf  uniprot_sprot.fasta
slurm-4247477.out  uniprot_sprot_db.phr  uniprot_sprot_db.pjs  uniprot_sprot_db.psq  uniprot_sprot_db.pto


Agora na pasta: /home/gdegaki/transcriptomas/Transdecoder

```bash
nano blastp_SRR8944275.slurm
```

```bash
#!/bin/bash
#SBATCH -t 5:00:00
#SBATCH -c 10
#SBATCH --mem=16G

# Ativar ambiente Conda
source /home/gdegaki/anaconda3/bin/activate
conda activate /home/gdegaki/anaconda3/envs/transdecoder

# Definir caminhos
PEP_FILE=/home/gdegaki/transcriptomas/Transdecoder/Resultados/Trinity_SRR8944275.trasdecoder_dir/trinity_SRR8944275.Trinity.fasta.tr
ansdecoder_dir/longest_orfs.cds
DB_FASTA=/home/gdegaki/transcriptomas/Transdecoder/banco_dados/uniprot_sprot.fasta
DB_NAME=/home/gdegaki/transcriptomas/Transdecoder/banco_dados/uniprot_sprot_db
OUT_DIR=/home/gdegaki/transcriptomas/Transdecoder/Resultados/Trinity_SRR8944275.trasdecoder_dir/blastp

# Criar diretório de saída se não existir
mkdir -p $OUT_DIR

# Rodar BLASTp
blastp -query $PEP_FILE \
       -db $DB_NAME \
       -evalue 1e-5 \
       -num_threads 10 \
       -max_target_seqs 5 \
       -outfmt 6 \
       -out $OUT_DIR/blastp.outfmt6
```


















