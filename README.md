Use the atoms_helper object (under atoms_helper.py), put in your input data. It'll generate a graph as an object attribute and assign each atom its position (adsorbate/slab, coordinated or not).

You can then make a subgraph of the adsorbate and however many extra slab atoms you want connected to the coordinated adsorbate atoms.

You can get all subgraphs from an input graph. This uses the RDkit functionality to do so. Think Benson's Group Additivity, but for all subgraphs instead of Benson's groups. 

You can put a bunch of adsorbed systems or whatever into the Subgraph Enumerator object and it'll enumerate those subgraphs for you and return the subgraph frequency matrix, or just the counter. 
Comparing full NetworkX graphs is slow, so I hash them first with the NetworkX Weisfeiler-Lehman has, and for whatever reason, that's faster. There's probably a trade-off with accuracy, but I leave that as an exercise for the reader.



Look at these papers for more on subgraph enumeration:

Thermochemistry of gas-phase and surface species via LASSO-assisted subgraph selection (Gu et al.)
React. Chem. Eng., 2018,3, 454-466

Thermochemical Data Fusion Using Graph Representation Learning (Bhattacharjee et al.)
J. Chem. Inf. Model. 2020, 60, 10, 4673–4683

Regularized machine learning on molecular graph model explains systematic error in DFT enthalpies(Bhattacharjee et al.)
Bhattacharjee, H. et al., Regularized machine learning on molecular graph model explains systematic error in DFT enthalpies. Sci Rep 11, 14372 (2021).
