<h1 align="center">Dinâmica Molecular da Crotamina (PDB: 1H5O) em água</h1>

<div align="center">
  <strong>🚀 Objetivo 📚</strong>
</div>

<div align="center">
  <p>Um repositório incrível com um projeto espetacular! 🎉</p>
  <p>Aqui você encontrará informações sobre o projeto, tecnologias utilizadas, instruções para configurar o ambiente de desenvolvimento e muito mais.</p>
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

Vamos trabalhar com a biomolécula [Crotamina](https://doi.org/10.1016/0003-9861(56)90444-1) que possui o codigo [1H5O](https://www.rcsb.org/structure/1H5O) no PDB. O PDB é um banco com várias biomoléculas depositadas e identificadas por códigos. Explore mais informações do PDB e da molécula.

<img src="./img/crotamina.jpg" alt="Crotamina">

>[!TIP]
> Organize o diretório de trabalho criando as pastas `analysis` para os arquivos de analises e `inputs` para os arquivos .mdp da dinâmica.
>

```
├── 1h5o.pdb
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

- Tecnologia 1
- Tecnologia 2
- Tecnologia 3
- ...

Descreva as tecnologias utilizadas no projeto e forneça links para as respectivas documentações ou sites oficiais.

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
