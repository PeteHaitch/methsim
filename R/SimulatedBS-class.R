### =========================================================================
### SimulatedBS: An S4 class to store data on a single simulated
### bisulfite-sequencing assay; BS = [WGBS | RRBS | eRRBS]
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Design
###
### list("chr:pos1:...:posm" = (z1, ..., zm))
###     "chr:pos1:...:posm": The co-ordinates of methylation loci assayed by
###                          the read.
###           (z1, ..., zm): The methylation state of each methylation loci
###                          assayed by the read.
###
### An object of this class is returned by the
### simulate,SimulateBSParam-method.
### This class will only retain those reads overlapping at least one
### methylation locus.
