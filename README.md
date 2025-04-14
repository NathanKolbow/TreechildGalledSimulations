# Simulation Study Design

This simulation study is split into two portions:

1. The true network is tree-child and galled and we search for networks under this restriction.
    - Here: true species network $\to$ simulated gene trees $\to$ simulated sequences $\to$ estimated gene trees $\to$ observed concordance factors $\to$ estimated networks
2. The true network is either (a) non tree-child but galled, (b) tree-child but not galled, or (c) neither tree-child nor galled and we search for networks still under the tree-child/galled restriction.
    - Here: true species network $\to$ expected concordance factors

## Varied Parameters

- Number of taxa (10/20/30) [only 20 for portion (2)]
- Reticulation density (low/high) [only low for portion (2)]
- Level of ILS (low/high) [only low for portion (2)]
- Sequence length (500/100) [only applicable for portion (1)]
- Number of genes (100/1,000/10,000) [only applicable for portion (1)]

