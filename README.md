# Biological Agent Performance Evaluation

Biological AI agents are increasingly used to automate single-cell RNA-seq analysis, yet their evaluation remains challenging because open-ended workflows rely on many context-dependent scientific decisions. Here we present a usability case study of recent AI agents with distinct architectural designs, evaluated under a standardized local setup using a hierarchical framework spanning **scientific validity**, **technical reliability**, and **usability and efficiency**.


## 1. Project Overview
This repository contains detailed experimental results and performance reports for our study. Currently we evaluate **three different methods** across **two classic datasets**.

## 2. Methods & References
Below are the official resources and research papers for the agents used in this study.

| Method | GitHub Repository | Research Paper (bioRxiv/arXiv) |
| :--- | :--- | :--- |
| **Biomni** | [snap-stanford/Biomni](https://github.com/snap-stanford/Biomni) | [Huang et al. (2025)](https://doi.org/10.1101/2025.05.30.656746) |
| **Pantheon-CLI** | [aristoteleo/pantheon-cli](https://github.com/aristoteleo/pantheon-cli) | [PantheonOS Team](https://pantheonos.stanford.edu/cli/docs/intro/getting-started) |
| **STELLA** | [zaixizhang/STELLA](https://github.com/zaixizhang/STELLA) | [Jin et al. (2025)](https://arxiv.org/abs/2507.02004) |

## 3. Performance Matrix (Detailed Reports)
This table summarizes the performance of different agents across two benchmark datasets. Click on the links to view or download the detailed reports.

| Method | EMTAB3929 (Embryo) | GSE94820 (Blood) |
| :--- | :--- | :--- |
| **Biomni** | [Download PDF](./reports/Biomni_embryo.pdf) | [Download PDF](./reports/Biomni_blood.pdf) |
| **Pantheon-CLI** | [Download PDF](./reports/Pantheon_embryo.pdf) | [Download PDF](./reports/Pantheon_blood.pdf) |
| **STELLA** | [View Details (MD)](./reports/STELLA_embryo.md) | [View Details (MD)](./reports/STELLA_blood.md) |

---
[Contact Authors](mailto:liyuzhuo23@mails.tsinghua.edu.cn)
