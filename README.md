# Matrix Three-Pass Regression Filter: A Nowcasting Framework for Mixed-Frequency Multidimensional Data

**Can preserving the natural matrix structure of macroeconomic data improve GDP nowcasting?**

This project compares **matrix** and **vector** factor models in a mixed-frequency nowcasting framework.  
The paper develops the **Matrix MF--TPRF** (Matrix Mixed-Frequency Three-Pass Regression Filter), a matrix-valued extension of the Mixed-Frequency Three-Pass Regression Filter, and evaluates its performance against standard vector benchmarks through both **Monte Carlo simulations** and **pseudo real-time euro-area GDP nowcasting**.

## Motivation

Macroeconomic datasets are naturally organized as **countries × indicators × time**.  
Standard vector approaches ignore this structure by stacking the panel into a long vector.  
The matrix approach instead preserves the bilinear organization of the data and directly exploits cross-country and cross-variable comovement.

## Main contribution

This project:

- develops a **matrix-valued mixed-frequency nowcasting model**;
- handles **ragged-edge missingness** and asynchronous releases;
- compares **Matrix MF--TPRF** with both pooled and country-specific vector benchmarks in the same empirical environment;
- studies pseudo real-time GDP nowcasting for the main euro-area countries in a **real-data nowcasting application**;
- provides both **Monte Carlo evidence** to analyse the finite simple properties.

## Results at a glance

The empirical evidence suggests that the **matrix formulation can improve nowcast accuracy**, especially when cross-country dependence is strong.  
The gains are most evident relative to pooled vectorized formulations and become especially visible in the **post-COVID period**, when cross-country linkages appear more informative.  
The matrix specification also tends to perform particularly well for **very early nowcasts**, especially at **M1**.  
At the same time, country-specific vector benchmarks remain competitive and are computationally simpler, especially when the information set becomes very large.

## Keywords

**Matrix Factor Models · GDP Nowcasting · Mixed Frequency · Euro Area · TPRF**

---

# Repository structure

The project folder is called:

`Matrix_vs_Vector_TS`

The most relevant subfolders are:

- `code/`
- `LaTeX_TPRF/`
- `literature/`

This README is organized in two main parts:

1. **How to navigate the code**
2. **How to navigate the LaTeX project**

---

# Part I — How to navigate the code

## Main code folders

Inside

`Matrix_vs_Vector_TS/code`

the project is organized as follows:

- `functions/`  
  contains all functions used across the models;

- `MonteCarlo_simulation/`  
  contains the code for the Monte Carlo exercises;

- `TPRF_Models_EA/`  
  contains the empirical euro-area nowcasting code for the matrix and vector models;

- `data/`  
  contains the datasets used by the scripts.

The folder

`Matrix_vs_Vector_TS/literature`

contains the reference material for both the matrix and vector approaches.

---

## General logic of the code

The empirical code is divided into three main model classes:

- **Matrix model (Matrix MF-TPRF)**
- **Pooled vectorized model (VEC-P)**
- **Country-specific vector model (VEC-C)**

For each model, the general workflow is:

1. run the **main** script to estimate the model and save outputs;
2. run the corresponding **results** script (mandatory for the Matrix MF-TPRF and optional for country dpecific vector models);
3. combine the outputs into final tables and graphs for the paper.

Generated outputs are typically saved inside the relevant `results/` folders and then collected into the final paper-ready material.

---

## Functions

All functions used in the project are stored in:

`code/functions`

This is the best place to start if the goal is to understand the implementation details.

---

## Monte Carlo

The Monte Carlo code is stored in:

`code/MonteCarlo_simulation`

This folder contains the scripts used to generate the simulation exercises reported in the paper.

---

## Empirical euro-area nowcasting

The main empirical code is stored in:

`code/TPRF_Models_EA`

Inside this folder, the relevant subfolders are:

- `Matrix_MF-TPRF/`
- `VecTensor_MF-TPRF/`
- `Vector_MF-TPRF/`
- `Final_Tab_Graph/`

---

## Recommended running order for the empirical results

### 1. Matrix MF--TPRF

Run first:

`code/TPRF_Models_EA/Matrix_MF-TPRF/matrix.mf.tprf.main.R`

Then run:

`code/TPRF_Models_EA/Matrix_MF-TPRF/matrix.mf.tprf.results.R`

This should be repeated for the desired parametrization:

- `small`, `medium`, `large`
- `corr` or `LASSO`

---

### 2. Pooled vectorized benchmark

For the pooled vectorized benchmark, run:

`code/TPRF_Models_EA/VecTensor_MF-TPRF/all_vectorized.mf.tprf.main.R`

---

### 3. Country-specific vector benchmark

For the country-specific vector benchmark, run:

`code/TPRF_Models_EA/Vector_MF-TPRF/all_mf.tprf.main.R`

---

### 4. Country-by-country versions

If country-specific runs are needed, use:

- `code/TPRF_Models_EA/VecTensor_MF-TPRF/cc_vectorized.mf.tprf.main.R`
- `code/TPRF_Models_EA/Vector_MF-TPRF/cc_mf.tprf.main.R`

and then the corresponding results scripts:

- `code/TPRF_Models_EA/VecTensor_MF-TPRF/cc_vectorized.mf.tprf.results.R`
- `code/TPRF_Models_EA/Vector_MF-TPRF/cc_mf.tprf.results.R`

---

### 5. Final paper tables and graphs

Once the relevant model results have been generated, run:

- `code/TPRF_Models_EA/paper.results.models.R`
- `code/TPRF_Models_EA/paper.results.variables.R`

These scripts produce:

- final tables and graphs for model comparison;
- tables on variable selection by country and model size.

The final empirical outputs are stored in:

`code/TPRF_Models_EA/Final_Tab_Graph`

---

## Practical note

Most scripts are designed to save outputs automatically into the relevant `results/` folders.  
The final paper-ready material is then assembled inside `Final_Tab_Graph/`.

---

# Part II — How to navigate the LaTeX project

The main paper folder is:

`Matrix_vs_Vector_TS/LaTeX_TPRF`

This is the most important folder if the goal is to read, modify, or compile the paper.

From the current structure, the main elements are:

- `chapter/`
- `figures/`
- `frontmatter/`
- `output/`
- `preamble/`
- `main.tex`

There are also several auxiliary LaTeX files (`.aux`, `.fls`, `.fdb_latexmk`, etc.), which can safely be ignored.

---

## Main LaTeX logic

The file `main.tex` is the entry point of the document.

It does **not** contain the whole paper text.  
Its role is simply to assemble the document by calling the different files corresponding to chapters, sections, subsections, and appendix material.

So:

- if you want to understand how the document is assembled, open `main.tex`;
- if you want to edit a specific section, go directly to the relevant file inside `chapter/`.

---

## `chapter/`

This folder contains the actual text of the paper.

The paper is written in modular form: chapters, sections, subsections, and appendix material are stored in separate files.  
Therefore, if you are interested in a specific part of the manuscript, the correct place to go is the corresponding file inside `chapter/`, not `main.tex`.

In short:

- `main.tex` only assembles the document;
- `chapter/` contains the substantive text.

---

## `figures/`

This folder contains the figures used in the paper.

It includes the plots generated for:

- Monte Carlo results;
- euro-area empirical results;
- appendix figures.

If a figure appears in the paper, the corresponding image file is usually stored here.

---

## `preamble/`

This folder contains all global formatting settings, such as:

- packages;
- custom commands;
- style settings;
- theorem environments;
- spacing and layout choices.

If the goal is to modify the appearance of the paper, this is the right place to intervene.

---

## `frontmatter/`

This folder contains the front material of the paper, such as title page elements and other introductory parts.

---

## `output/`

This folder contains the compiled PDF.

To read the current version of the paper, open:

`LaTeX_TPRF/output/main.pdf`

---

## Practical compilation note

The LaTeX project is easiest to manage in **VS Code**.  
Because the project is relatively large, Overleaf may run into timeout issues, especially with the free version.

For this reason, VS Code is the recommended environment for editing and compiling the manuscript.

---

# Suggested navigation guide

## If you want to understand the code

Start from:

1. `code/functions/`
2. `code/TPRF_Models_EA/`
3. `code/MonteCarlo_simulation/`

## If you want to reproduce the empirical results

Follow this order:

1. `Matrix_MF-TPRF/matrix.mf.tprf.main.R`
2. `Matrix_MF-TPRF/matrix.mf.tprf.results.R`
3. pooled vector scripts
4. country-specific vector scripts
5. `paper.results.models.R`
6. `paper.results.variables.R`

## If you want to read or edit the paper

Start from:

1. `LaTeX_TPRF/main.tex`
2. `LaTeX_TPRF/chapter/`
3. `LaTeX_TPRF/output/main.pdf`

---

# Final summary

In short:

- `code/functions/` contains all reusable functions;
- `code/MonteCarlo_simulation/` contains the Monte Carlo code;
- `code/TPRF_Models_EA/` contains the empirical euro-area nowcasting code;
- `code/TPRF_Models_EA/Final_Tab_Graph/` contains final tables and graphs;
- `literature/` contains the reference material;
- `LaTeX_TPRF/main.tex` assembles the paper;
- `LaTeX_TPRF/chapter/` contains the actual text;
- `LaTeX_TPRF/preamble/` controls formatting;
- `LaTeX_TPRF/figures/` contains the paper figures;
- `LaTeX_TPRF/output/main.pdf` is the compiled paper.