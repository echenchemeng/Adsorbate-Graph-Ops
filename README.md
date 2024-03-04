Use the atoms_helper object (under atoms_helper.py), put in your input data. It'll generate a graph as an object attribute and assign each atom its position (adsorbate/slab, coordinated or not).

You can then make a subgraph of the adsorbate and however many extra slab atoms you want connected to the coordinated adsorbate atoms.

You can get all subgraphs from an input graph. This uses the RDkit functionality to do so. Think Benson's Group Additivity, but for all subgraphs instead of Benson's groups. 

Pending: with a set of molecules (or adsorbate systems), you can enumerate all the subgraphs. 
Comparing full NetworkX graphs is slow, so I hash them first with the NetworkX Weisfeiler-Lehman has, and for whatever reason, that's faster. There's probably a trade-off with accuracy, but I leave that as an exercise for the reader.



Look at these papers for more on subgraph enumeration:
https://pubs-rsc-org.udel.idm.oclc.org/en/content/articlelanding/2018/re/c7re00210f
https://pubs-acs-org.udel.idm.oclc.org/doi/full/10.1021/acs.jcim.0c00699
https://www-nature-com.udel.idm.oclc.org/articles/s41598-021-93854-w
