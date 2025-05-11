<h1 align="center">Dinâmica Molecular da Crotoxina (PDB: 3R0L) em água</h1>

<div align="center">
  <strong>🚀 Objetivo 📚</strong>
</div>

<div align="center">
  <p>Simular a proteína Crotoxina em uma caixa cúbica de água com temperatura de 298 K e 1 bar de pressão. A Crotoxina (CTX) é uma proteína que compõe o veneno da cobra Cascavel.</p>
  <p>Explore, colabore e divirta-se! 😄</p>
</div>

## 📖 Índice

- [Visão Geral](#visão-geral)
- [Tecnologias](#tecnologias)
- [Configuração do Ambiente](#configuração-do-ambiente)
- [Como Contribuir](#como-contribuir)
- [Licença](#licença)

## Arquivos iniciais.

Inicialmente precisamos obter as coordenadas da nossa biomolécula, campos de forças e arquivos inputs para a dinâmica. Essa etapa faz parte do planejamento do projeto.

Vamos trabalhar com a biomolécula [Crotoxina](https://doi.org/10.1016/j.jmb.2011.07.027) que possui o codigo [3R0L](https://www.rcsb.org/structure/3R0L) no PDB. O PDB é um banco com várias biomoléculas depositadas e identificadas por códigos. Explore mais informações do PDB e da biomolécula.

<img src="./img/crotoxina.png" alt="CTX">

>[!TIP]
> Organize o diretório de trabalho criando as pastas `analysis` para os arquivos de analises e `inputs` para os arquivos .mdp da dinâmica molecular.
>

```
├── 3r0l.pdb
├── amber14sb_parmbsc1_cufix.ff
│   ├── aminoacids.arn
│   ├── aminoacids.c.tdb
│   ├── aminoacids.hdb
│   ├── aminoacids.n.tdb
│   ├── aminoacids.r2b
│   ├── aminoacids.rtp
│   ├── aminoacids.vsd
│   ├── atomtypes.atp
│   ├── ca-sol7.itp
│   ├── ca-sol7.pdb
│   ├── cufix.itp
│   ├── dna.arn
│   ├── dna.hdb
│   ├── dna.r2b
│   ├── dna.rtp
│   ├── ffbonded.itp
│   ├── ffnonbonded.itp
│   ├── ffnonbonded.itp~
│   ├── forcefield.doc
│   ├── forcefield.itp
│   ├── forcefield.itp~
│   ├── gbsa.itp
│   ├── ions.itp
│   ├── Makefile.am
│   ├── Makefile.in
│   ├── mg-sol6.itp
│   ├── mg-sol6.pdb
│   ├── README.md
│   ├── rna.arn
│   ├── rna.hdb
│   ├── rna.r2b
│   ├── rna.rtp
│   ├── spce.itp
│   ├── spc.itp
│   ├── tip3p.itp
│   ├── tip4pew.itp
│   ├── tip4p.itp
│   ├── tip5p.itp
│   ├── urea.itp
│   └── watermodels.dat
├── analysis
└── inputs
    ├── ions.mdp
    ├── md.mdp
    ├── minim.mdp
    ├── npt.mdp
    └── nvt.mdp
```

## Preparo da topologia da molécula: campos de forças.

Nessa etapa, é necessário escolher o modelo de água e o campo de força utilizado. O arquivo `3r0l.pdb` contém as coordenadas da biomolécula com moleculas de água e ligantes e será necessário remover as moléculas de água (`HOH`) e outros ligantes (`HETATM`) para evitar erros. Isso pode ser feito manualmente direto no arquivo ou pelo prompt de comando:

```
grep -v HETATM 3r0l.pdb > 3r0l_clean.pdb
```

Para escolher o campo de força e o modelo de água:

```
gmx pdb2gmx -v -f 3r0l_clean.pdb -o crotoxina.gro

# -v = verbose, para visualizar o processo.
# -f = file input, arquivo de coordenadas de entrada.
# -o = file output, arquivo de coordenadas de saída.
```
Quando solicitado, digite o número correspondente para selecionar o campo de força e o modelo de água.

O Gromacs assumirá valores canônicos para cada aminoácidos, levando em consideração valores de pH próximos da neutralidade.

>[!NOTE]
>Saiba mais sobre o comando [gmx2pdb](https://manual.gromacs.org/documentation/current/onlinehelp/gmx-pdb2gmx.html).
>Será gerado os seguintes arquivos:
> - crotoxina.gro = arquivo com as coordenadas de cada átomo da biomolécula compatível com o campo de força.
> - topol.top = arquivo com a topologia da biomolécula, ou seja, com os parâmetros necessários para o cálculo das forças.
> - posre.itp = arquivo de topologia auxiliar indicando os átomos com restrição por padrão.
>


## Definindo a caixa de simulação.

Forneça instruções claras e detalhadas sobre como configurar o ambiente de desenvolvimento localmente. Isso pode incluir:

- Pré-requisitos
- Instalação de dependências
- Configuração do banco de dados
- Configuração de variáveis de ambiente
- Execução de migrações ou scripts de inicialização
- ...

Certifique-se de fornecer exemplos de comandos ou scripts necessários para executar o projeto corretamente.

## Minimização do sistema

Se você deseja contribuir para o projeto, siga estas etapas:

1. Faça um fork do repositório e clone-o em sua máquina local.
2. Crie uma nova branch para suas modificações:
   ```
   git checkout -b minha-branch
   ```
3. Faça as modificações desejadas e adicione-as ao stage:
   ```
   git add .
   ```
4. Faça um commit das suas alterações:
   ```
   git commit -m "Minhas modificações"
   ```
5. Envie suas alterações para o repositório remoto:
   ```
   git push origin minha-branch
   ```
6. Abra um pull request para que suas modificações sejam revisadas e incorporadas ao projeto.

## Equilíbrio NVT e NPT: termostatos e barostatos.

Se você deseja contribuir para o projeto, siga estas etapas:

## Produção: integradores.

Se você deseja contribuir para o projeto, siga estas etapas:

---

## 📄 Licença

Este projeto está licenciado sob a [Nome da Licença]. Consulte o arquivo [LICENSE](LICENSE) para obter mais informações sobre os termos de licenciamento.
