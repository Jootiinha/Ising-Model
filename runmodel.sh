#!/bin/bash
start=`date +%s`
echo "=================================="
echo "Iniciando simulação"
echo "=================================="
echo "=================================="
echo "=================================="
echo "Passo 1 ===> Instalação de dependências"
echo " * Necessário informar credências de super usuário - sudo"
#Instalar as dependências
sudo apt-get install -y git g++ libgsl-dev
echo "Dependências instaladas"
echo "=================================="
echo "=================================="

echo "Passo 2 ===> Buscando código da aplicação em repositório git (Github)"
echo " * Clonando repositório"
git clone https://github.com/Jootiinha/Ising-Model.git
cd Ising-Model/
echo "Repositório clonado"
echo "=================================="
echo "=================================="

echo "Passo 3 ==> Executando simulação"
g++ model.cpp -lgsl -lgslcblas -lm*


end=`date +%s`
runtime=$((end-start))

echo "=================================="
echo "=================================="
echo " * Simulação executada com sucesso *"
echo " * Tempo total de execução do script: $runtimes*"
echo " * Dados da simulação no arquivo DATA.1.dat *"
echo "=================================="
echo "=================================="



