# Copyright (C) 2015 Peter Hickey
#
# This file is part of methsim.
#
# methsim is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# methsim is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with methsim  If not, see <http://www.gnu.org/licenses/>.

### =========================================================================
### Prepare Lister data for use as internal data. This internal data is
### eventually made external via the call to utils::data()
### -------------------------------------------------------------------------
###

# WARNING: For internal use only and not to be run in directory containing
# package source code.
# WARNING: The use of this file requires that methsim is already installed on
# the system.
# WARNING: The creation of the example data is not yet fully reproducible. The
# raw data (output files of 'methtuple') are very large, and depend on their
# own set of parameters. More relevant, the PartitionedMethylome depends on
# random number generation, and beta_by_pm_region, lor_by_my_region,
# pattern_freqs_by_pm_region all depend on the original PartitionedMethylome.
#
# This script is was written to run in
# ~/hickey/methylation_m-tuples/analyses/methsim and assumes the following
# directory structure:
#    |-rds
#    |---Lister
#    |-----beta_by_pm_region
#    |-----lor_by_pm_region
#    |-----PartitionedMethylome
#    |-----pattern_freqs_by_pm_region

prepare <- function(sample_name) {

  # The dataset is currently hardcoded. If this is changed then this function
  # will require modification.
  dataset <- "Lister"
  BSgenomeName <- "BSgenome.Hsapiens.UCSC.hg18"
  l_pm <- readRDS(paste0("rds/", dataset, "/PartitionedMethylome/", dataset,
                         "_pm.rds"))
  meth_level <- readRDS(paste0("rds/", dataset, "/beta_by_pm_region/", dataset,
                               "_beta_by_pm_region.rds"))
  cometh <- readRDS(paste0("rds/", dataset, "/lor_by_pm_region/", dataset,
                           "_lor_by_pm_region.rds"))
  pattern_freqs <- readRDS(paste0("rds/", dataset,
                                  "/pattern_freqs_by_pm_region/", dataset,
                                  "_pattern_freqs_by_pm_region.rds"))
  pm <- l_pm[[sample_name]]
  meth_level <- meth_level[sample == sample_name, ][, sample := NULL]
  cometh <- cometh[sample == sample_name, ][, sample := NULL]
  pattern_freqs <- pattern_freqs[sample == sample_name][, sample := NULL]

  SimulateMethylomeParam(BSgenome = BSgenomeName,
                         PartitionedMethylome = pm,
                         MethLevelDT = meth_level,
                         ComethDT = cometh,
                         PatternFreqsDT = pattern_freqs,
                         SampleName = sample_name)
}

ADS <- prepare("ADS")
ADS_iPSC <- prepare("ADS-iPSC")
ADS_adipose <- prepare("ADS-adipose")

# These files should go in data/
save(ADS, file = "ADS.RData", compress = "xz")
save(ADS_iPSC, file = "ADS_iPSC.RData", compress = "xz")
save(ADS_adipose, file = "ADS_adipose.RData", compress = "xz")
