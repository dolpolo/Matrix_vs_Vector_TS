# Matrix vs Vector Time Series for GDP Nowcasting

**Can preserving the natural matrix structure of macroeconomic data improve GDP nowcasting?**

This project compares **Matrix** and **Vector** factor models in a mixed-frequency nowcasting framework.  
The paper develops the **Matrix MF–TPRF**, a matrix-valued extension of the Mixed-Frequency Three-Pass Regression Filter, and evaluates its performance against the standard vector benchmark through simulation and pseudo real-time euro-area GDP nowcasting. :contentReference[oaicite:5]{index=5}

## Motivation

Macroeconomic datasets are often naturally organized as **countries × indicators × time**.  
Standard vector methods ignore this structure by stacking everything into a long vector.  
The matrix approach instead preserves the bilinear organization of the panel and exploits cross-country comovement directly. :contentReference[oaicite:6]{index=6}

## Contribution

This paper:

- develops a **matrix-valued mixed-frequency nowcasting model**;
- handles **ragged-edge missingness** and asynchronous releases;
- compares **Matrix MF–TPRF** and **Vector MF–TPRF** on the same data environment;
- studies euro-area GDP nowcasting for **Germany, France, Italy, and Spain**. :contentReference[oaicite:7]{index=7}

## Results at a glance

The evidence suggests that the **matrix formulation improves nowcast accuracy**, especially when cross-country linkages are strong.  
The gains are most evident for **France and Spain**, while **Germany and Italy** show stronger improvements after COVID. Vector methods remain useful and lighter benchmarks, but the matrix specification appears better suited for highly interconnected macroeconomic systems. :contentReference[oaicite:8]{index=8}

## Project contents

- `code/` — estimation, simulation, nowcasting, utilities  
- `LaTeX/` — paper source  
- `README.md` — overview of the project  

## Keywords

**Matrix Factor Models · GDP Nowcasting · Mixed Frequency · Euro Area · TPRF**