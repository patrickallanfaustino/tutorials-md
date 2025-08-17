<h1 align="center">Dinâmica Molecular da Insulina Humana (PDB: 3I40) em água</h1>

<div align="center">
  <strong>🚀 Objetivo 📚</strong>
</div>

<div align="center">
  <p>O objetivo deste tutorial é simular a insulina humana, o hormônio que regula o metabolismo da glicose, em uma caixa com água cúbica sob condições de 298 K e 1 bar.</p>
  <p>Explore, colabore e estude! 😄</p>
</div>

## 📖 Índice

- [Arquivos iniciais](#arquivos-iniciais)
- [Preparo da topologia da molécula: campos de forças](#preparo-da-topologia-da-molécula-campos-de-forças)
- [Definindo a caixa de simulação: dimensões, solvatação e neutralização](#definindo-a-caixa-de-simulação-dimensões-solvatação-e-neutralização)
- [Minimização do sistema](#minimização-do-sistema)
- [Equilíbrio NVT e NPT: termostatos e barostatos](#equilíbrio-nvt-e-npt-termostatos-e-barostatos)
- [Produção: integradores](#produção-integradores)

## Arquivos iniciais

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

```
├── 3i40.pdb
├── amber14sb_parmbsc1_cufix.ff
├── analysis
└── inputs
    ├── ions.mdp
    ├── md.mdp
    ├── minim.mdp
    ├── npt.mdp
    └── nvt.mdp
```

## Preparo da topologia da molécula: campos de forças

O arquivo **3i40.pdb** contém, além das coordenadas da biomolécula, moléculas de água (`HOH`) e outros ligantes (`HETATM`). Remova esses componentes extras para evitar erros nas etapas subsequentes. Realize essa limpeza de duas maneiras: editando o arquivo manualmente ou utilizando os comandos de terminal apresentados a seguir.

```
grep -v HETATM 3i40.pdb > 3i40_clean.pdb

# grep -v HOH 3i40.pdb > 3i40_clean.pdb
```

Observe que algumas biomoléculas apresenta múltiplas cadeias, identificadas como `chain A`, `chain B`, e assim por diante. Remova as cadeias desnecessárias. No caso, a cadeia B foi removida com um editor de texto simples.

Em seguida, escolha o campo de força e o modelo de água que serão usados na simulação:

>[!NOTE]
>Caso utilize um campo de força externo, coloque a pasta correspondente no diretório de trabalho como <name>.ff.
>

```
gmx pdb2gmx -v -f 3i40_clean.pdb -o insulina.gro

# -v = verbose, para visualizar o processo.
# -f = file input, arquivo de entrada das coordenadas.
# -o = file output, arquivo de saída das coordenadas.
```

O programa solicitará duas escolhas em sequência. Responda a cada prompt da seguinte forma:
 - Para o campo de força, digite 1 para selecionar AMBER03.
 - Para o modelo de água, digite 1 novamente para selecionar TIP3P, o padrão recomendado para a família AMBER.

O GROMACS utiliza estados de protonação canônicos para cada aminoácido (assumindo pH neutro) e adiciona os hidrogênios correspondentes. Ao final do processo, o programa conserva a carga líquida total da biomolécula. Confirme este valor no terminal, procurando pela mensagem: `Total charge in system -2.000 e`.

Para visualizar no VMD, utilize:
```
vmd insulina.gro
```

>[!NOTE]
>Saiba mais sobre [gmx2pdb](https://manual.gromacs.org/documentation/current/onlinehelp/gmx-pdb2gmx.html).
>
>Serão criados os arquivos:
> - insulina.gro = arquivo com as posições iniciais de cada átomo da biomolécula compatível com o campo de força escolhido.
> - topol.top = arquivo com os parâmetros necessários para as integrações e derivações numéricas.
> - posre.itp = arquivo de topologia auxiliar indicando os átomos com restrições por padrão.
>

Campo de Força  |  Informações  |  Modelo de água  |  cut-off
------- | ---------- | -------- | -------- 
**OPLS**    | O campo de força OPLS-AA (Optimized Potentials for Liquid Simulations – All Atom) é amplamente usado para simulações de proteínas, pequenas moléculas, solventes, lipídios, dentre outros. | TIP4P recomendado, mas pode usar TIP3P. Não recomendado SPC. | 1.0~1.2 nm
**AMBER**   | A família de campos de força AMBER (como amber99sb, amber99sb-ildn, amber14, etc.) é amplamente usada para proteínas, DNA/RNA e simulações biomoleculares. | TIP3P, não recomendado TIP4P e SPC. | 1.0~1.2 nm
**CHARMM**  | O campo de força CHARMM (como charmm36-jul2022.ff) é extremamente detalhado, especialmente para lipídios, proteínas e açúcares, e foi parametrizado com *switching functions*, o que o diferencia das abordagens anteriores. | TIP3P modificado, não substituir por TIP3P comum. | 1.2 nm
**GROMOS**  | O campo de força GROMOS96 (como gromos54a7.ff) é uma escolha clássica para simulações de proteínas, sistemas aquosos e alguns tipos de estudos de bioenergia. Ele é o único desta lista a usar potencial truncado sem PME. | SPC | 1.4 nm

| Modelo | Tipo | Descrição |
|--------|---------|--------------------------------|
| **SPC** | 3 pontos | Modelo rígido, ângulo fixo de 109.47°, parametrizado para propriedades macroscópicas. |
| **SPC/E** | 3 pontos | Versão estendida do SPC, com correção de energia de polarização. Melhor densidade e constante dielétrica. |
| **TIP3P** | 3 pontos | Muito usado com AMBER e CHARMM. Simples e compatível com muitos campos de força. |
| **TIP4P** | 4 pontos | Inclui ponto virtual (M-site) para carga negativa fora do oxigênio, melhorando propriedades de fase. |
| **TIP5P** | 5 pontos | Dois pontos extra para os pares de elétrons do oxigênio. Mais preciso para estrutura tetraédrica, porém mais custoso. |

>[!IMPORTANT]
>Escolha o campo de força e o modelo de água com base na natureza do seu sistema e nas propriedades que você deseja investigar.
>
>É de extrema importância o conhecimento completo sobre os formatos de arquivos utilizados pelo GROMACS. Para estudos: [File formats topology](https://manual.gromacs.org/current/reference-manual/topologies/topology-file-formats.html) e [File formats](https://manual.gromacs.org/current/reference-manual/file-formats.html).
>

## Definindo a caixa de simulação: dimensões, solvatação e neutralização

Nesta etapa, procede-se à definição da caixa de simulação, ajustando-se suas dimensões, a distância da biomolécula em relação às bordas e outros parâmetros relevantes para a correta configuração do sistema.

```
gmx editconf -f insulina.gro -o box.gro -c -d 2.5 -bt cubic

# -c = center, para centralizar a biomolécula na caixa.
# -d = distance, distância em nm entre todas moléculas e a borda.
# -bt = box type, formato da caixa.
```

O formato da caixa pode ser definido em `cubic`, `triclinic`, `octahedron` ou `dodecahedron`. A escolha do formato da caixa de simulação é de responsabilidade do pesquisador e deve considerar o formato da biomolécula, com o objetivo de otimizar o número de moléculas no sistema. Essa estratégia permite economizar recursos computacionais, equilibrando tempo de simulação e demanda energética.

As dimensões da caixa selecionada podem ser verificadas na mensagem de saída. Recomenda-se que a distância entre a biomolécula e as bordas da caixa esteja entre 1,0 e 2,5 nm, pois esses valores são considerados ideais.

>[!NOTE]
>Saiba mais sobre [editconf](https://manual.gromacs.org/documentation/current/onlinehelp/gmx-editconf.html).
>
>Essa função é util para converter arquivos .pdb <--> .gro usando `gmx editconf -f <file>.gro -o <file>.pdb`.
>

>[!IMPORTANT]
>A tag `-box` pode ser utilizada para definir as dimensões da caixa de simulação. Por exemplo, ao executar `gmx editconf -f insulina.gro -o box.gro -c -d 2.5 -bt cubic -box 10 10 10`, obtém-se uma caixa cúbica com arestas de 10 nm. Nessa configuração, a distância da borda definida como 2,5 nm será considerada, resultando em um espaço útil de 7,5 nm para a acomodação das moléculas, garantindo o afastamento adequado entre a molécula e as bordas da caixa.
>
>**E quando não definimos `-box`?** Nessa configuração, o algoritmo do GROMACS definirá as dimensões da caixa com base no tamanho máximo da biomolécula, acrescido da distância especificada para a borda. Essa abordagem proporciona uma margem suficiente para garantir uma dinâmica molecular segura, ao mesmo tempo em que promove o uso eficiente dos recursos computacionais.
>

<div align="center">
<img src="img/box.png" alt="caixa de simulação">
</div>

>PDB 3I40, insulina humana em uma caixa de simulação cubica 7.8 x 7.8 x 7.8 nm.

### Solvatação
Na sequência, a caixa de simulação será preenchida com moléculas de água, uma vez que o objetivo é investigar a solvatação da insulina em meio aquoso. Esse procedimento assegura que a biomolécula esteja imersa em um ambiente que simule condições fisiológicas adequadas para a dinâmica molecular.

```
gmx solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top

# -cp = coordenates protein, coordenadas do nosso soluto (geralmente, proteina).
# -cs = coordenates solvent, coordenadas da molecula que será usada como solvente.
# -p = processing, para processar o arquivo de topologia do sistema.
```

O GROMACS preencherá toda a caixa de simulação com moléculas de água provenientes do arquivo `spc216.gro`, que é adequado para o modelo TIP3P. As moléculas de água serão identificadas pelo resíduo **SOL**. Na mensagem de saída, é possível observar a quantidade total de moléculas de solvente adicionadas por meio da linha `Number of solvent molecules`, sendo essa informação incorporada automaticamente ao `[ molecules ]` do arquivo de topologia do sistema.

>[!NOTE]
>Saiba mais sobre [solvate](https://manual.gromacs.org/documentation/current/onlinehelp/gmx-solvate.html).
>
>Adicionalmente, podemos definir **-box** para definir as dimensões de uma nova caixa de simulação e **-maxsol** para definir a quantidade máxima de moleculas adicionadas, útil nos cálculos de concentrações.
>

>[!IMPORTANT]
>Para modelos de água TIP4P, `-cs` utilize `tip4p.gro`.
>

<div align="center">
<img src="img/solvate.png" alt="proteina solvatada">
</div>

>PDB 3I40 solvatada com água modelo TIP3P

### Neutralização
Na última etapa do preparo da caixa de simulação, procede-se à neutralização do sistema por meio da adição de íons. Essa etapa é fundamental, pois os integradores utilizados nas simulações apresentam maior eficiência em sistemas eletricamente neutros. Considerando que a insulina apresenta carga total de -2,000 e, como evidenciado anteriormente na preparação da topologia, será necessária a adição de cátions de modo a compensar essa carga e garantir a neutralidade global do sistema.

Antes de neutralizar com a função `genion`, é necessário gerar um arquivo binário .tpr com as informações necessárias para o processamento:

```
gmx grompp -v -f inputs/ions.mdp -o ions.tpr -c solv.gro -p topol.top

# -c = coordenates, arquivo com as coordenadas do sistema.
```

Na tag -f está indicado o arquivo [ions.mdp](inputs-easy/ions.mdp) da pasta `/inputs`. Esse arquivo possui todos os parâmetros necessários para o processamento dessa etapa. Recomenda-se um [estudo intensivo](https://manual.gromacs.org/current/user-guide/mdp-options.html) sobre os parâmetros desse arquivo.

>[!NOTE]
>Saiba mais sobre [grompp](https://manual.gromacs.org/documentation/current/onlinehelp/gmx-grompp.html).
>
>Em algumas oportunidades, o GROMACS gera `warnings` que devem ser verificados e, se necessário, suprimidos com **-maxwarn [x]**, onde `x` é a quantidade de `warnings` a ser suprimidos.
>

Neste momento, procede-se à neutralização da caixa de simulação, adicionando-se os íons necessários para assegurar a neutralidade do sistema.

```
gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15

# -s = submit binary, arquivo binário criado anteriormente com todas informações do sistema.
# -pname = nome do cátion(+), nesse caso NA Sódio.
# -nname = nome do ânion(-), nesse caso CL Cloro.
# -neutral = para neutralizar completamente o sistema, às vezes desnecessário.
# -conc 0.15 = concentration, define a concentração em mol/L.
```

Por padrão, o GROMACS adiciona íons de sódio (NA) e cloreto (CL) em quantidade suficiente apenas para neutralizar a biomolécula e neste caso, considerando a carga líquida de -2,000 e, serão adicionados dois íons NA ao sistema. Entretanto, ao utilizar as opções `-conc 0.15` e, opcionalmente, `-neutral`, é possível garantir a adição de uma solução fisiológica a 0,9% m/m, simulando um ambiente semelhante ao sistema biológico humano, além de assegurar a neutralidade do sistema. Na mensagem de saída, pode-se observar a mensagem `Will try to add 45 NA ions and 43 CL ions`, indicando o número de íons a serem incorporados para atingir a concentração e a neutralidade desejados.

O `genion` solicitara selecionar qual o grupo de moléculas que serão substituídas pela adição dos ions. Por convenção, utiliza-se o grupo **SOL** para remover moléculas de água em troca dos ions. Selecione o número correspondente ao SOL.

>[!NOTE]
>Saiba mais sobre [genion](https://manual.gromacs.org/documentation/current/onlinehelp/gmx-genion.html).
>

<div align="center">
<img src="img/neutralization.png" alt="proteina solvatada e neutralizada">
</div>

>PDB 3I40 solvatada e neutralizada. Em 🔵 NA e 🟢 CL.

Pronto, agora nossa caixa de simulação está pronta!

## Minimização do sistema

Neste momento, realiza-se a minimização da energia potencial do sistema, processando as eventuais sobreposições entre as moléculas. Para isso, é gerado novamente o arquivo binário .tpr e, em seguida, procede-se à minimização energética do sistema, garantindo uma configuração estrutural inicial estável e adequada para as etapas subsequentes da simulação.

```
gmx grompp -v -f inputs/minim.mdp -c solv_ions.gro -o em.tpr -p topol.top
```
```
gmx mdrun -v -deffnm em

# -deffnm = define o nome padrão de todos arquivos de entrada e saida.
```

A função `mdrun` constitui o núcleo do processo de dinâmica molecular no GROMACS. Recomenda-se simplificar os nomes dos arquivos de entrada e saída utilizando a opção `-deffnm`. O nome utilizado em `grompp -o <name>.tpr` deve ser o mesmo especificado em `-deffnm`, garantindo consistência entre os arquivos utilizados. Para a etapa de minimização, é adotado o arquivo de parâmetros [minim.mdp](inputs-easy/minim.mdp), que contém as opções específicas para o procedimento de minimização energética.


>[!NOTE]
>Saiba mais sobre [mdrun](https://manual.gromacs.org/documentation/current/onlinehelp/gmx-mdrun.html).
>

Para o acompanhamento eficiente desta etapa, recomenda-se a análise do gráfico de energia potencial do sistema. Para isso, o arquivo .edr, que armazena as energias calculadas durante a simulação, deve ser lido e convertido em um arquivo .xvg. Esse procedimento permite a avaliação gráfica da convergência e da estabilidade energética do sistema ao longo do processo de minimização.

```
gmx energy -f em.edr -s em.tpr -o potential.xvg
```

Verifique na tabela o número correspondente ao 'Potential' e digite-o, seguindo por um espaço e pelo número 0 (zero). Exemplo: **10 0**.

>[!NOTE]
>Saiba mais sobre [energy](https://manual.gromacs.org/documentation/current/onlinehelp/gmx-energy.html).
>

Utiliza-se o `XMGrace` para visualizar o gráfico:

```
xmgrace potential.xvg
```

Observa-se a curva gerada no gráfico, a qual indica a minimização efetiva do sistema.

<div align="center">
<img src="img/minim.png" alt="gráfico da energia minimizada">
</div>

## Equilíbrio NVT e NPT: termostatos e barostatos

As próximas etapas consiste no ajuste da temperatura e da pressão do sistema, estabelecendo-se 298,15 K (25 ºC) para a temperatura e 1 bar (0,98 atm) para a pressão. Essas condições visam simular um ambiente termodinâmico semelhante ao meio biológico natural.

### NVT: ajustando a temperatura da caixa de simulação
Mantendo o número de moléculas (N), o volume (V) e a temperatura (T) constantes, gera-se o arquivo binário .tpr utilizando o arquivo de parâmetros [nvt.mdp](inputs-easy/nvt.mdp). No arquivo `nvt.mdp` define-se os parâmetros:

* Define-se a restrição da biomolécula, com `define = -DPOSRES`.
* Define-se o tempo para o ajuste da temperatura, em `nsteps = 50000` x 0,002 (dt) = 100 ps.
* Define-se o algoritmo para o ajuste da temperatura, em `tcoupl = V-rescale`.
* Define-se os grupos para o ajuste da temperatura, com `tc-grps = Protein   Non-Protein`.
* Define-se a constante de acoplamento da temperatura, com `tau-t = 1.0`.
* Define-se a temperatura de referência, em `ref-t = 298.15`.

Destaca-se algumas considerações específicas:

* A restrição de posição dos átomos não-hidrogênios da biomolécula nas etapas subsequentes é necessária para preservar a conformação da biomolécula enquanto se promove o ajuste do solvente ao seu redor. Caso algum átomo exceda o limite estabelecido no arquivo posre.itp (padrão 1000 kJ/mol/nm), será permitido o movimento apenas desse átomo, mantendo os demais restritos conforme os parâmetros definidos.
* A aplicação do termostato em grupos distintos, como definido em `tc-grps = Protein Non-Protein`, é mais eficiente e garante maior acurácia ao controle de temperatura. Essa abordagem permite que a proteína e o solvente sejam tratados separadamente, ajustando com precisão as variações térmicas de cada componente do sistema.
* A constante de acoplamento da temperatura, definida com `tau-t = 1.0`, assegura que o termostato seja aplicado nesse intervalo de tempo, medido em picossegundos. Esse valor pode variar entre **0,5 e 1,0 ps**, devendo garantir que permaneça sempre **menor que a constante de acoplamento da pressão**. Ressalta-se que valores demasiadamente pequenos para tau-t podem ocasionar instabilidade no sistema, levando à 'explosão' (colapso estrutural ou erros críticos durante a simulação).

```
gmx grompp -v -f inputs/nvt.mdp -c em.gro -r em.gro -o nvt.tpr -p topol.top

# -r = restrain file, arquivo de coordenadas com as restrinções iniciais (geralmente mesmo arquivo).
```
```
gmx mdrun -f -deffnm nvt
```

>[!NOTE]
>Nota-se a performance na mensagem de saída, pode ser útil para planejar o tempo da simulação baseado no seu computador. Exemplo: 210.65 ns/day ou 0.114 hour/ns.
>

Procede-se a visualização do gráfico para a verificação da temperatura do sistema. Essa análise permite confirmar se a temperatura média está de acordo com o valor estabelecido nos parâmetros de simulação, além de avaliar possíveis flutuações durante o processo.

```
gmx energy -f nvt.edr -s nvt.tpr -o temperature.xvg
```

Selecione o número correspondente a 'Temperature' seguida por espaço e 0 (zero).

```
xmgrace temperature.xvg
```

<div align="center">
<img src="img/temperature.png" alt="gráfico da temperatura">
</div>

Após 20 ps, observa-se que a temperatura do sistema estabilizou em 298,15 K. Caso a estabilização não seja alcançada, recomenda-se aumentar o valor de `nsteps` e repetir a etapa. Com a temperatura devidamente controlada, procede-se ao ajuste da pressão do sistema.

### NPT: ajustando a pressão da caixa de simulação
Mantendo o número de moléculas (N), o pressão (P) e a temperatura (T) constantes, gera-se o arquivo binário .tpr utilizando o arquivo de parâmetros [npt.mdp](inputs-easy/npt.mdp). Nesse arquivo `npt.mdp` define-se:

* O algoritmo responsável por ajustar a pressão, com `pcoul = C-rescale`.
* A constante de acoplamento da pressão, em `tau-p = 3.0`.
* A pressão de referência, em `ref-p = 1.0`.

Os demais parâmetros utilizados nesta etapa são idênticos ou semelhantes aos empregados na etapa NVT, contudo o tempo de equilíbrio costuma ser um pouco maior na etapa NPT. O arquivo binário .tpr será gerado a partir das coordenadas obtidas previamente na etapa NVT, garantindo a continuidade do processo de simulação sob o novo conjunto de condições.

```
gmx grompp -v -f inputs/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -o npt.tpr -p topol.top

# -t = time file, arquivo com checkpoint anterior (geralmente utilizado para indicar o ponto de partida com relação a dinâmica anterior)
```
```
gmx mdrun -v -deffnm npt
```

<div align="center">
<img src="img/pressure.png" alt="gráfico da pressão">
</div>
<div align="center">
<img src="img/density.png" alt="gráfico da densidade">
</div>

A análise do gráfico de pressão revela a presença de picos distintos, que não são representativos nem adequados para avaliar o desempenho do barostato. Para esse fim, o gráfico de densidade mostra-se mais apropriado, pois permite observar a estabilização da densidade do sistema, geralmente acompanhada de pequenas variações, indicando o equilíbrio adequado sob as condições simuladas.

A seguir, apresenta-se um breve resumo sobre os principais termostatos e barostatos utilizados em simulações de dinâmica molecular.

| Termostato | Características | Vantagens | Limitações
|--------|---------|-------------|---------------|
| **Berendsen** | Rápido para equilibrar temperatura | Simples e eficiente para equilíbrios | Não reproduz corretamente as flutuações canônicas |
| **V-rescale*** | Mantém temperatura média correta e flutuações realistas | Estável e mais preciso que Berendsen | Ligeiramente mais complexo |
| **Nose-Hoover** | Mantém distribuição canônica (NVT) | Correto estatisticamente | Pode ter acoplamento mais lento |

| Barostato | Características | Vantagens | Limitações
|--------|---------|-------------|---------------|
| **Berendsen** | Ajusta pressão rapidamente durante o equilíbrio | Simples, ideal para pré-produção | Não reproduz corretamente as flutuações canônicas |
| **Parrinello-Rahman** | Permite flutuações de volume e forma da caixa (NPT) | Correto para simulações de produção | Pode ser instável sem bom equilíbrio inicial |
| **C-rescale*** | Versão estocástica rigorosa de controle de pressão. Mantém flutuações canônicas corretas no ensemble NPT | Produz NPT canônico exato, mais robusto e estável que Parrinello-Rahman em algumas situações | Disponível a partir do GROMACS 2023, não testado quanto Parrinello-Rahman |

>[!IMPORTANT]
>A escolha do termostato e barostato deve considerar a natureza do sistema e as propriedades que se deseja investigar.
>
>O GROMACS recomenda: **V-rescale** e **C-rescale**.
>

Agora estamos prontos para nossa simulação!

## Produção: integradores

Working...
