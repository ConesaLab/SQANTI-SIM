<p align="center">
  <img src="https://github.com/ConesaLab/SQANTI-SIM/blob/main/docs/SQANTI-SIM_logo.png" alt="" width="300">
</p>

<p align="center">
  <a href="https://github.com/ConesaLab/SQANTI-SIM/wiki/Requirements-and-installation">Installation</a>
  ·
  <a href="https://github.com/ConesaLab/SQANTI-SIM/wiki">Documentation</a>
  ·
  <a href="https://github.com/ConesaLab/SQANTI-SIM/issues">Report bug</a>
  ·
  <a href="https://github.com/ConesaLab/SQANTI3">SQANTI3</a>
</p>

***

SQANTI-SIM is an open-source long-read RNA-seq (LRS) simulation tool designed to generate high-fidelity ONT and PacBio transcriptome long reads with precise control over transcript novelty based on [SQANTI3](https://github.com/ConesaLab/SQANTI3) structural categories, along with  orthogonal data supporting both known and novel transcripts.

By building upon the capabilities of [NanoSim](https://github.com/bcgsc/NanoSim), [PBSIM3](https://github.com/yukiteruono/pbsim3), and [IsoSeqSim](https://github.com/yunhaowang/IsoSeqSim) to simulate long reads using the provided GTF, SQANTI-SIM generates both simulated reads and a reduced GTF with excluded transcripts flagged as novel. Additionally, SQANTI-SIM goes a step further by generating matching Illumina and CAGE data for the simulated LRS, thus faithfully reproducing LRS transcriptome datasets.

![small_workflow](https://github.com/ConesaLab/SQANTI-SIM/blob/main/docs/small_workflow.png)

## Documentation

Please refer to [Wiki](https://github.com/ConesaLab/SQANTI-SIM/wiki) for how to install and use SQANTI-SIM:

- [Overview](https://github.com/ConesaLab/SQANTI-SIM/wiki/Overview)

- [Requirements and installation](https://github.com/ConesaLab/SQANTI-SIM/wiki/Requirements-and-installation)

- [Parameters summary](https://github.com/ConesaLab/SQANTI-SIM/wiki/Parameters-summary)

- [Run SQANTI-SIM](https://github.com/ConesaLab/SQANTI-SIM/wiki/Run-SQANTI-SIM)

- [SQANTI-SIM: classif](https://github.com/ConesaLab/SQANTI-SIM/wiki/SQANTI-SIM:-classif)

- [SQANTI-SIM: design](https://github.com/ConesaLab/SQANTI-SIM/wiki/SQANTI-SIM:-design)

- [SQANTI-SIM: sim](https://github.com/ConesaLab/SQANTI-SIM/wiki/SQANTI-SIM:-sim)

- [Tutorial](https://github.com/ConesaLab/SQANTI-SIM/wiki/Tutorial)

## Change Log

Current version (04/10/2023): **0.2.1**

Updates, patches and releases:

- [Change log](https://github.com/ConesaLab/SQANTI-SIM/wiki/Change-Log)

## How to cite SQANTI-SIM

The SQANTI-SIM paper is currently in progress. Meanwhile, you can access the preprint [here](https://www.biorxiv.org/content/10.1101/2023.08.23.554392v1).
