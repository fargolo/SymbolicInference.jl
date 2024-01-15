# SymbolicInference

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://fargolo.github.io/SymbolicInference.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://fargolo.github.io/SymbolicInference.jl/dev/)
[![Build Status](https://github.com/fargolo/SymbolicInference.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/fargolo/SymbolicInference.jl/actions/workflows/CI.yml?query=branch%3Amain)

Probability-based inferences based on the symbolic method.  

White-paper link: https://docs.google.com/document/d/1giF2vcG8dcEHf4SONktg21Kk9KU5SXEAyO0fpNnyINY/edit#heading=h.d898sqp17elx  

## Time-series  

### Chaotic, non-linear systems and recurrence analysis  

The recurrence of states, in the meaning that states are again arbitrarily close after some time of divergence, is a fundamental property of deterministic dynamical systems and is typical for nonlinear or chaotic systems.  

Poincaré discusse this property in 1890 and it was later proved by Constantin Carathéodory (see Poincaré recurrence theorem). Further on, several techniques address recurrences in dynamical systems for inference.  

Recurrence plots (RPs) were creature to capture such patterns and several parameters that caractherize the underlying time-series can be obtained with recurrence quantification analysis (RQA).  

Since each pair of states is mapped into a binary value ('close enough',1, or 'not close enough',0), making probabilistic inference with symbolic methods is straightforward.  

The procedure in `double_inference_weighted()` iterates over [recurrence matrices](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/recurrenceanalysis/stable/) diagonals and checks whether the size of the longest consecutive sequence is unexpectedly large. 

## Graphs  
