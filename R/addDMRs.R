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
### addDMRs: A method to add differentially methylated regions to a
### SimulatedMethylome.
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Design
###
### This should return a suitably modified SimulatedMethylome object and some
### information on the DMRs (e.g., location, effect size).
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### TODO
###
### This is a longer term goal that I haven't yet figured out how best to
### implement. The key question is how to control the effect size.

# addDMRs() will work on SimulatedMethylome2 objects by modifying the
# marginalProb assay. It should return a SimulatedMethylome2 object with either
# the original data as the first sample and the modified data as the second
# sample, or just the modified data.
