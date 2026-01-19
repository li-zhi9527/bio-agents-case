# Agent Performance Evaluation

Biological AI agents are increasingly used to automate single-cell RNA-seq analysis, yet their evaluation remains
challenging because open-ended workflows rely on many context-dependent scientific decisions. Here we present an end-
to-end usability case study of contemporary AI agents with distinct architectural designs, evaluated under a standardized
local setup using a hierarchical framework spanning scientific validity, technical reliability, and usability and efficiency.
Using single-cell datasets with well-established biological ground truths, we show that successful execution and well-formed
reports do not necessarily imply scientifically valid conclusions. Agents with large, unrestricted action spaces tend to
prioritize uninterrupted execution, making them prone to silent fallback behaviors and context misinterpretation that lead
to biologically incorrect results. In contrast, tool-oriented systems exhibit stronger scientific validity and reproducibility,
but require greater user involvement and reduced autonomy, while multi-agent systems improve robustness to low-
level technical failures yet often struggle to translate computational outputs into biologically grounded interpretations.
Together, our findings expose fundamental limitations of execution-centric benchmarks and motivate process-aware
evaluation strategies that assess the quality of intermediate decisions and reasoning, providing guidance for the principled
benchmarking and development of reliable AI agents for biological research.
This is a GitHub repository of several case studies that records the behavior of current Bio AI agents applied on classic datasets we choose.


## 1. Project Overview
This repository contains detailed experimental results and performance reports for our study. Currently we evaluate **three different methods** across **two classic datasets** .

## 2. Performance Matrix
The following table provides links to the detailed reports for each experimental configuration. 

| Method | EMTAB3929 | GSE94820 | GitHub Repository | Research Paper (bioRxiv/arXiv) |
| :--- | :--- | :--- | :--- | :--- |
| **Biomni** | [Download PDF](./reports/Biomni_embryo.pdf) | [Download PDF](./reports/Biomni_blood.pdf) |[snap-stanford/Biomni](https://github.com/snap-stanford/Biomni)|[Huang et al. (2025)](https://doi.org/10.1101/2025.05.30.656746)|
| **Pantheon-CLI** | [Download PDF](./reports/Pantheon_embryo.pdf) | [Download PDF](./reports/Pantheon_blood.pdf) |[aristoteleo/pantheon-cli](https://github.com/aristoteleo/pantheon-cli)|[PantheonOS Team](https://pantheonos.stanford.edu/cli/docs/intro/getting-started)|
| **STELLA** | [View Details (MD)](./reports/STELLA_embryo.md) | [View Details (MD)](./reports/STELLA_blood.md) |[zaixizhang/STELLA](https://github.com/zaixizhang/STELLA)|[Jin et al. (2025)](https://arxiv.org/abs/2507.02004)|


## 3. Key Findings


---
[Back to Main Paper](#) | [Contact Authors](mailto:liyuzhuo23@mails.tsinghua.edu.cn)
