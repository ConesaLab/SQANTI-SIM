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

Please refer to [Wiki](https://github.com/ConesaLab/SQANTI-SIM/wiki) for how to install and use SQANTI-SIM 

## Change Log

Current version (25/07/2023): **0.2.0**

Updates, patches and releases:

### Version 0.2.0 (25/07/2023)

#### :rocket: New Features and Major Changes

- Implemented the simulation of **CAGE-seq** orthogonal data (`--CAGE`) when running SQANTI-SIM in *sample* mode.
- Added PacBio long read simulation using **PBSIM3** (`--pbsim`).
- Added the `full-sim` mode, which executes the *classif*, *design* and *sim* modules, allowing users to run the complete simulation with just **one simple instruction**.
- Provided pre-trained models and characterized datasets of the WTC11 human cell line.
- Updated pre-trained NanoSim models with WTC11 human cell line dRNA and cDNA data. Additionally, users have the option to utilize their specific pre-trained NanoSim models using the `--nanosim_model` parameter.
- Enhanced the `--diff_exp` parameter by deprecating `--low_prob` and `--high_prob`. Now, `--diff_exp` is a numeric value that allows users to adjust with a single parameter the bias of known and novel expression levels.

#### :bug: Bug Fixes and Minor Changes

- Added the option `--expr_file` to provide an expression distribution file to `sqanti-sim.py design sample`  to avoid mapping the same reads in different executions.
- When using the `--long_reads` option in `sqanti-sim.py design sample` mode, SQANTI-SIM now generates an expression file named "sqanti-sim_expression.tsv." This file can be reused in future executions using `--expr_file`.
- Added a default number of novel transcripts to simulate when the user does not specify any.
- Added additional example files, providing users with practical examples.
- Updated the classification algorithm to the latest SQANTI3 release, version 5.1.2, ensuring up-to-date transcript classification results.
- Updated the `sqanti-sim.py design sample` random expression sampling by implementing inverse transform sampling from ECDFs.
- Fixed Partial True Positives (PTP) and Positive Detection Rate (PDR) metrics to consider only non-redundant transcript models.
- Fixed conda environment installation issues related to bcbio-gff and updated dependencies.
- Fixed SQANTI-SIM Report generating empty tables when not all structural categories were detected. 

***

### Version 0.1.0 (14/07/2022)

- :tada: Initial release!

## How to cite SQANTI-SIM

SQANTI-SIM paper in progress
