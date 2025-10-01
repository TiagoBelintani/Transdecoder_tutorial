Pipeline TransDecoder para ORFs em aranhas
#1 Estrutura de Diretórios
# Acessar diretório de trabalho
cd transcriptoma/

# Criar pasta principal do TransDecoder
mkdir Transdecoder

# Criar subpasta para banco de dados
mkdir Transdecoder/banco_dados

# Acessar subpasta do banco
cd Transdecoder/banco_dados

#2 Download do banco de dados UniProt Swiss-Prot
# Baixar o banco curado Swiss-Prot
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

# Descompactar o arquivo
gunzip uniprot_sprot.fasta.gz


O arquivo uniprot_sprot.fasta será usado nas etapas seguintes do pipeline TransDecoder.

#3 Organização do fluxo de trabalho

Diretórios do projeto:

Adapters  assembly  busco_analises  fastqc_results  multiqc_report  raw_data


Para esta análise, utilizaremos assembly: Trinity.fasta

#4 Etapa 1 – TransDecoder.LongOrfs

Criar a pasta de resultados:

mkdir -p ~/transcriptomas/Transdecoder/Resultados


Criar o job SLURM:

nano Transdecoder_LongOrfs_SRR8944275.slurm


Conteúdo do job:

#!/bin/bash
#SBATCH -t 3:00:00
#SBATCH -c 4

# Ativar ambiente Conda
source /home/gdegaki/anaconda3/bin/activate 
conda activate /home/gdegaki/anaconda3/envs/transdecoder 

# Executar TransDecoder.LongOrfs
TransDecoder.LongOrfs \
  -t /home/gdegaki/transcriptomas/assembly/trinity_SRR8944275.Trinity.fasta \
  -G universal -S \
  --output_dir /home/gdegaki/transcriptomas/Transdecoder/Resultados/Trinity_SRR8944275.trasdecoder_dir/

Resultado esperado

Dentro de:

cd /home/gdegaki/transcriptomas/Transdecoder/Resultados/Trinity_SRR8944275.trasdecoder_dir/trinity_SRR8944275.Trinity.fasta.transdecoder_dir/


Você encontrará:

base_freqs.dat
__checkpoints_longorfs
longest_orfs.cds
longest_orfs.gff3
longest_orfs.pep


O arquivo longest_orfs.pep será usado para BLASTp.

#5 Etapa 2 – Formatar banco BLAST

Dentro da pasta:

cd /home/gdegaki/transcriptomas/Transdecoder/banco_dados
nano makeblast.slurm


Conteúdo do job:

#!/bin/bash
#SBATCH -t 1:00:00

makeblastdb -in /home/gdegaki/transcriptomas/Transdecoder/banco_dados/uniprot_sprot.fasta \
            -dbtype prot \
            -out /home/gdegaki/transcriptomas/Transdecoder/banco_dados/uniprot_sprot_db

Resultado esperado

Arquivos criados:

makeblast.slurm
uniprot_sprot_db.pdb
uniprot_sprot_db.pin
uniprot_sprot_db.pot
uniprot_sprot_db.ptf
uniprot_sprot_db.phr
uniprot_sprot_db.pjs
uniprot_sprot_db.psq
uniprot_sprot_db.pto
uniprot_sprot.fasta
slurm-XXXXXX.out

#6 Etapa 3 – BLASTp

Criar o job SLURM:

cd /home/gdegaki/transcriptomas/Transdecoder
nano blastp_SRR8944275.slurm


Conteúdo do job:

#!/bin/bash
#SBATCH -t 5:00:00
#SBATCH -c 10
#SBATCH --mem=16G

# Ativar ambiente Conda
source /home/gdegaki/anaconda3/bin/activate
conda activate /home/gdegaki/anaconda3/envs/transdecoder

# Definir caminhos
PEP_FILE=/home/gdegaki/transcriptomas/Transdecoder/Resultados/Trinity_SRR8944275.trasdecoder_dir/trinity_SRR8944275.Trinity.fasta.transdecoder_dir/longest_orfs.pep
DB_FASTA=/home/gdegaki/transcriptomas/Transdecoder/banco_dados/uniprot_sprot.fasta
DB_NAME=/home/gdegaki/transcriptomas/Transdecoder/banco_dados/uniprot_sprot_db
OUT_DIR=/home/gdegaki/transcriptomas/Transdecoder/Resultados/Trinity_SRR8944275.trasdecoder_dir/blastp

# Criar diretório de saída se não existir
mkdir -p $OUT_DIR

# Formatar banco BLAST caso ainda não esteja pronto
if [ ! -f "${DB_NAME}.pin" ]; then
    makeblastdb -in $DB_FASTA -dbtype prot -out $DB_NAME
fi

# Executar BLASTp
blastp -query $PEP_FILE \
       -db $DB_NAME \
       -evalue 1e-5 \
       -num_threads 10 \
       -max_target_seqs 5 \
       -outfmt 6 \
       -out $OUT_DIR/blastp.outfmt6

✅
















