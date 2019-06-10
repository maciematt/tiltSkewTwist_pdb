tiltSkewTwist_pdb
=================

This script takes in pdb's and spits out the skequential tilt, skew, and twist angles between component domains, as defined by the user.

Requires: numpy

Usage: `./tiltSkewTwist.py [infile] [dom1 residues] [dom1 atoms] [dom1 ref] [dom2 residues] [dom2 atoms] [dom2_ref]`

`[dom1/dom2 residues]` is the residue range of a given domain given in PyMOL format, i.e. `low1-high1+low2-high2` etc., same goes for `[dom1/dom2 atoms]`, eg. `n+ca+c` (write `all` to include all atom types), and `[dom1/dom2 ref]` is the seq. number of the reference residue.

Note: review the code to choose the coordinate set type (right- or left-handed; right-handed by default), the angle convention (0 - 360 degrees or -180 - 180 degrees; the latter by default) etc., as annotated within the code.


visualize_tilt_skew_twist
=========================

This script can be used to visualize the output of tiltSkewTwist_pdb.py.

Requires: SAGE

Usage: `sage ./visualize_tilt_skew_twist.sage`

Before running make sure, however, to swap the names in the infilename and outfilename variables to match filename of the output of tiltSkewTwist_pdb.py and the desired output filename.
