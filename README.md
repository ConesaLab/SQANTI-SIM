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

*SQANTI-SIM* is a simulator of controlled novelty and degradation of transcripts sequenced by long reads. It is a wrapper tool for RNA-Seq long-reads simulators such as [IsoSeqSim](https://github.com/yunhaowang/IsoSeqSim) and [NanoSim](https://github.com/bcgsc/NanoSim) (formerly Trans-NanoSim) to simulate transcripts based on the SQANTI3 structural categories ([publication](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5848618/) and [code repository](https://github.com/ConesaLab/SQANTI3)).

SQANTI-SIM uses NanoSim and IsoSeqSim to simulate PacBio cDNA reads and Nanopore cDNA and dRNA reads, and implements a strategy to simulate novel transcripts based on SQANTI structural categories. Unlike current lrRNA-seq simulators, SQANTI-SIM simulates well-grounded novel transcripts and permits the assessment of the capacity of any lrRNA-seq transcript reconstruction tool to detect different types of transcripts not yet present in reference annotations.

![small_workflow](https://github.com/ConesaLab/SQANTI-SIM/blob/main/docs/small_workflow.png)

## Documentation

Please refer to [Wiki](https://github.com/ConesaLab/SQANTI-SIM/wiki) for how to install and use SQANTI-SIM 

## Change Log

Current version (12/09/2022): 0.1.0

Updates, patches and releases:

**0.1.0**:
- :tada: Initial release!

## How to cite SQANTI-SIM

SQANTI-SIM will be included in the upcoming SQANTI3 paper
