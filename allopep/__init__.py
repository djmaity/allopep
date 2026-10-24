"""AlloPep: allosteric peptide prediction pipeline.

Keep this module free of imports. It is loaded from three Conda
environments -- ``gaps`` for :mod:`allopep.gaps_predict`, ``PepGLAD`` for
:mod:`allopep.recover_relaxation`, and ``allopep`` for everything else --
and only the last one has Biopython, freesasa, psutil and OpenStructure.
Importing anything here breaks the other two environments.
"""
