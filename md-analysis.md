<h1 align="center">Análises de Resultados</h1>

<div align="center">
  <strong>🚀 Objetivo 📚</strong>
</div>

<div align="center">
  <p>O objetivo deste tutorial é realizar algumas análises de uma dinâmica molecular básica, obtendo: RMSD, RMSF, DSSP, H-Bonds, SASA, Rg, Distribuição radial e propriedades termodinâmicas.</p>
  <p>Explore, colabore e estude! 😄 Dúvidas: <a href="mailto:patrick.faustino@unesp.br">patrick.faustino@unesp.br</a></p>
</div>

## 📖 Índice

- [Ajustes iniciais](#ajustes-iniciais)


## Ajustes iniciais

Para iniciar a simulação, obtenha os arquivos de topologia (campos de força), as coordenadas iniciais da biomolécula e os parâmetros de entrada para a dinâmica molecular.

Utilize a estrutura da [insulina humana](https://doi.org/10.1107/S1744309110000461) com o código [3I40](https://www.rcsb.org/structure/3I40) do PDB, que possui uma resolução de 1,85 Å. **Dê preferência a estruturas com resolução cristalográfica inferior a 2,5 Å**, pois isso garante uma geometria molecular mais confiável e detalhada, o que é fundamental para a qualidade da simulação. Uma resolução menor proporciona maior detalhamento cristalográfico.

Acesse a página da estrutura no [PDB (*Protein Data Bank*)](https://www.rcsb.org/) para uma análise aprofundada. Para garantir maior precisão e realismo, explore os detalhes complementares da estrutura. Verifique o método experimental usado para sua obtenção, a presença de ligantes, possíveis modificações estruturais e os estados de protonação dos resíduos.

<div align="center">
<img src="img/insulina.png" alt="insulina">
</div>

>PDB 3I40, insulina humana. O VMD (*Visual Molecular Dynamics*) possui esquema de cores para estruturas de biomoléculas: 🟣 violeta para alfa-hélices; 🟡 amarelo para beta-folhas; 🔵 ciano para voltas e ⚪ branco para superhélices ou cordas.

>[!TIP]
>Organize seu diretório de trabalho. Crie duas subpastas: `analysis`, destinada aos resultados das análises, e `inputs`, para armazenar os arquivos de parâmetros da dinâmica molecular (.mdp).
>



---

### 🧪⚗️ *Boas simulações moleculares!* 🦠🧬

---
## 📜 Citação

- FAUSTINO, Patrick Allan dos Santos. **Tutorials: Análise de Resultados de Dinâmica Molecular**. [*S. l.*]: Github, 18 jul. 2025. DOI 10.5281/zenodo.16062830. Disponível em: [https://github.com/patrickallanfaustino/tutorials-md/blob/main/md-analysis.md](https://github.com/patrickallanfaustino/tutorials-md/blob/main/md-analysis.md). Acesso em: 18 jul. 2025.
