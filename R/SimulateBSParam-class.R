### =========================================================================
### SimulateBSParam: An S4 class to store the parameters used to simulate
### bisulfite-sequencing data from a SimulatedMethylome object;
### BS = [WGBS | RRBS | eRRBS]
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Design
###
### list(SimulatedMethylome, aveCov, errorRate, target)
### SimulatedMethylome: A SimulatedMethylome object.
###             aveCov: The average sequencing coverage to simulate.
###          errorRate: The error rate (combined sequencing error and
###                     bisulfite-conversion error) to use in the simulation.
###             target: The co-ordinates of baits/enriched regions used in
###                     simulating RRBS/eRRBS data.
