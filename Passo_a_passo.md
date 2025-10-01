Preparando e Executando o Transdecoder
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
Download do Banco de Dados
bash
# Baixar banco de dados UniProt Swiss-Prot
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
Próximos Passos
Após o download, você pode descompactar o arquivo:

bash
gunzip uniprot_sprot.fasta.gz
O arquivo descompactado (uniprot_sprot.fasta) será utilizado nas etapas subsequentes do pipeline Transdecoder.

