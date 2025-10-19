<h1 align="center">Dinâmica Molecular de Biomolécula (PDB: 1S0Q) em água - modificado</h1>

<div align="center">
  <strong>🚀 Objetivo 📚</strong>
</div>

<div align="center">
  <p>O objetivo deste tutorial é simular a enzima digestiva tripsina pancreática bovina em uma caixa cúbica com água sob condições de 298 K e 1 bar, com modificações do fluxo de trabalho.</p>
  <p>Explore, colabore e estude! 😄 Dúvidas: <a href="mailto:patrick.faustino@unesp.br">patrick.faustino@unesp.br</a></p>
</div>

## Considerações

Após compreender perfeitamente como realizar uma dinâmica molecular simples da [tripsina pancreática bovina](https://github.com/patrickallanfaustino/tutorials-md/blob/main/md-easy.md), a evolução e destreza em realizar dinâmicas moleculares complexas exige etapas adicionais, que foram estudadas e são utilizadas por mim em artigos, dissertações, teses e quando necessário rigor acadêmico. Algumas observações importantes:

- Dependendo do sistema estudado, não torna os resultados mais precisos.
- Adiciona uma etapa de minimização antes da neutralização.
- Faz uso progressivo do termostato e barostato de Berendsen.

Assumindo que o estudante esta familializado com a dinâmica molecular anterior da [tripsina pancreática bovina](https://github.com/patrickallanfaustino/tutorials-md/blob/main/md-easy.md) e possui a prática necessária, segue abaixo o resumo das etapas:

**Preparo dos arquivos e escolha do campo de força**
```
grep -v HETATM 1S0Q.pdb > 1S0Q_clean.pdb
gmx pdb2gmx -v -f 1S0Q_clean.pdb -o tripsina.gro
vmd tripsina.gro
```

**Definição da caixa de simulação, solvatação e minimização**
```
gmx editconf -f tripsina.gro -o box.gro -c -d 2.0 -bt cubic
gmx solvate -cp box.gro -cs spc216.gro -o solvated.gro -p topol.top
gmx grompp -v -f inputs/minim.mdp -c solvated.gro -o solvated_em.tpr -p topol.top
gmx mdrun -v -deffnm solvated_em
```

**Neutralização**
```
gmx grompp -v -f inputs/ions.mdp -c solvated_em.gro -o ions.tpr -p topol.top
gmx genion -s ions.tpr -o solvated_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15
```

**Minimização**
```
gmx grompp -v -f inputs/minim.mdp -c solvated_ions.gro -o em.tpr -p topol.top
gmx mdrun -v -deffnm em
gmx energy -f em.edr -s em.tpr -o potential.xvg
xmgrace potential.xvg
```

**NVT Berendsen**
```
gmx grompp -v -f inputs/nvt_1.mdp -c em.gro -r em.gro -o nvt_1.tpr -p topol.top
gmx mdrun -v -deffnm nvt_1
```

**NVT V-rescale + NPT Berendsen**
```
gmx grompp -v -f inputs/npt_1.mdp -c nvt_1.gro -r nvt_1.gro -t nvt_1.cpt -o npt_1.tpr -p topol.top
gmx mdrun -v -deffnm npt_1
gmx energy -f npt_1.edr -s npt_1.tpr -o temperature.xvg
xmgrace temperature.xvg
```

**NVT V-rescale + NPT C-rescale**
```
gmx grompp -v -f inputs/npt_2.mdp -c npt_1.gro -r npt_1.gro -t npt_1.cpt -o npt_2.tpr -p topol.top
gmx mdrun -v -deffnm npt_2
gmx energy -f npt_2.edr -s npt_2.tpr -o pressure.xvg
xmgrace pressure.xvg
gmx energy -f npt_2.edr -s npt_2.tpr -o density.xvg
xmgrace density.xvg
```

**Produção**
```
gmx grompp -v -f inputs/md.mdp -c npt_2.gro -t npt_2.cpt -o md_5ns.tpr -p topol.top
gmx mdrun -v -deffnm md_5ns
```

>[!NOTE]
>Pode ocorrer alertas pelo Gromacs sobre o uso do termostato e barostato Berendsen. Para suprimir, use `-maxwarn [x]`, onde [x] é a quantidade de alertas emitidos.
>

Link para visualizar o video demonstrativo da dinâmica: [https://youtu.be/IQGiznRc0Xo](https://youtu.be/IQGiznRc0Xo).

- [Imagens e video](md-visual.md)
- [Análises de resultados](md-analysis.md)


---

### 🧪⚗️ *Boas simulações moleculares!* 🦠🧬

---

## 📜 Citação

- FAUSTINO, Patrick Allan dos Santos. **Tutorials: Dinâmica Molecular de Biomoléculas (PDB: 1S0Q) em água**. [*S. l.*]: Github, 18 jul. 2025. DOI 10.5281/zenodo.16062830. Disponível em: [https://github.com/patrickallanfaustino/tutorials-md/blob/main/md-easy.md](https://github.com/patrickallanfaustino/tutorials-md/blob/main/md-easy.md). Acesso em: 18 jul. 2025.
