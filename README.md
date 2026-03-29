# Matrix vs Vector Time Series for GDP Nowcasting

**Can preserving the natural matrix structure of macroeconomic data improve GDP nowcasting?**

This project compares **Matrix** and **Vector** factor models in a mixed-frequency nowcasting framework.  
The paper develops the **Matrix MF–TPRF**, a matrix-valued extension of the Mixed-Frequency Three-Pass Regression Filter, and evaluates its performance against the standard vector benchmark through simulation exercises and pseudo real-time euro-area GDP nowcasting.

## Motivation

Macroeconomic datasets are naturally organized as **countries × indicators × time**.  
Standard vector approaches ignore this structure by stacking the panel into a long vector.  
The matrix approach instead preserves the bilinear organization of the data and exploits cross-country comovement directly.

## Main contribution

This project:

- develops a **matrix-valued mixed-frequency nowcasting model**;
- handles **ragged-edge missingness** and asynchronous releases;
- compares **Matrix MF–TPRF** and **Vector MF–TPRF** in the same empirical environment;
- studies euro-area GDP nowcasting for **Germany, France, Italy, and Spain**.

## Results at a glance

The empirical evidence suggests that the **matrix formulation can improve nowcast accuracy**, especially when cross-country dependence is strong.  
The gains are more visible for **France and Spain**, while **Germany and Italy** show stronger improvements in the post-COVID period.  
Vector methods remain useful benchmark models, but the matrix specification appears better suited for highly interconnected macroeconomic systems.

## Keywords

**Matrix Factor Models · GDP Nowcasting · Mixed Frequency · Euro Area · TPRF**

---

# Repository structure

The main project folder is:

`Matrix_vs_Vector_TS`

Inside this folder, the repository is organized into four main directories:

- `code`
- `LaTeX`
- `LaTeX_TPRF`
- `literature`

This README is meant to help the reader distinguish clearly between:

1. the **code side** of the project;
2. the **paper / LaTeX side** of the project.

---

# Part I — Code and research material

## 1. `code/`

This folder contains the computational core of the project.

Here you can find the scripts and functions used to:

- preprocess the data;
- build the predictor panels;
- estimate the vector and matrix models;
- run simulations;
- perform pseudo real-time nowcasting;
- select hyperparameters;
- generate empirical results used in the paper.

This is the right place to start if the goal is to understand:

- how the models are implemented;
- how the nowcasting exercises are run;
- how the empirical results are produced.

It is strongly recommended to read this folder carefully, since it contains the full logic of the implementation and is the natural bridge between the methodology in the paper and the actual estimation routines.

## 2. `literature/`

This folder collects the theoretical and documentary background of the project.

In particular, it includes:

- references related to **vector** methods;
- references related to **matrix** methods;
- notes, comments, and corrections received from the professors during the drafting phase of the paper.

This folder is useful both as a bibliography archive and as a record of how the project evolved during the writing and revision process.

---

# Part II — Navigating the LaTeX paper folder

## 3. `LaTeX_TPRF/`

This is the main folder for the **paper writing and compilation**.

It contains the full LaTeX structure of the paper in modular form.  
The most important subfolders and files inside `LaTeX_TPRF` are:

- `chapter/`
- `figures/`
- `frontmatter/`
- `output/`
- `preamble/`
- `references`
- `main.tex`
- LaTeX auxiliary compilation files (`.aux`, `.toc`, `.fls`, `.fdb_latexmk`, etc.)

---

## 4. `main.tex`

The file `main.tex` is the central file of the paper.

Its role is mainly organizational: it does **not** contain the whole text of the paper, but instead assembles the document by calling the different `.tex` files placed in the various subfolders.

In practice:

- `main.tex` is the entry point of the LaTeX project;
- it imports the different parts of the paper;
- it controls the order in which the sections appear in the final document.

So, if someone wants to understand how the document is built, `main.tex` is the first file to inspect.

---

## 5. `preamble/`

This folder contains everything related to the **global formatting** of the paper.

This is where one can manage:

- LaTeX packages;
- custom commands;
- theorem or proposition environments;
- layout settings;
- stylistic choices;
- numbering conventions;
- headers, spacing, and general formatting rules.

### Important note
If the goal is to modify the **appearance** of the document, the correct place to intervene is the `preamble/` folder, not the chapter files.

---

## 6. `output/`

This folder contains the compiled output of the paper.

Most importantly, this is where the final PDF is stored:

- `main.pdf`

This is the file to open if you simply want to read or share the final version of the paper.

So:

- to **read the final paper** → open `LaTeX_TPRF/output/main.pdf`
- to **change the content** → edit the `.tex` files called by `main.tex`
- to **change the formatting** → work in `LaTeX_TPRF/preamble/`

---

## 7. `chapter/`

This folder contains the **main body of the paper**, divided into thematic sections.

According to the current structure, it includes folders such as:

- `Introduction`
- `Methodology`
- `Empirics`
- `Appendix`

and files such as:

- `abstract`
- `conclusion`

### Internal organization

Each subfolder corresponds to one major part of the paper.  
Inside each of these folders, the various sections and subsections are stored in separate files and kept in the intended order.

This means that:

- the paper is not written in one long file;
- each chapter is broken down into smaller section files;
- `main.tex` calls these files in sequence.

This modular structure makes the project easier to read, edit, and maintain.

---

## 8. `figures/`

This folder contains the figures, plots, and graphical material used in the paper.

Whenever a figure is inserted into the LaTeX document, the corresponding file is typically stored here.

---

## 9. `frontmatter/`

This folder contains the material that belongs to the front part of the paper, before the main chapters.

Depending on the structure of the document, this may include:

- title page material;
- introductory formal elements;
- other preliminary parts of the manuscript.

---

# Quick navigation guide

Depending on the task, the relevant folder is:

- **To understand the implementation** → `code/`
- **To consult theory and comments from the writing process** → `literature/`
- **To understand how the paper is assembled** → `LaTeX_TPRF/main.tex`
- **To edit the text of the paper** → `LaTeX_TPRF/chapter/`
- **To edit global formatting** → `LaTeX_TPRF/preamble/`
- **To manage figures** → `LaTeX_TPRF/figures/`
- **To read the final compiled version** → `LaTeX_TPRF/output/main.pdf`

---

# Suggested reading order

For someone approaching the repository for the first time, the suggested order is:

1. read this `README`;
2. inspect the `code/` folder;
3. inspect `LaTeX_TPRF/main.tex`;
4. explore the `chapter/` folder;
5. open `LaTeX_TPRF/output/main.pdf`;
6. consult `literature/` for references and comments.

---

# Final summary

In short:

- `code/` contains the computational and empirical implementation;
- `literature/` contains references and revision material;
- `LaTeX_TPRF/` contains the full paper structure;
- `main.tex` assembles the paper;
- `preamble/` controls formatting;
- `chapter/` contains the text, split into chapters and sections;
- `output/main.pdf` is the final compiled paper.