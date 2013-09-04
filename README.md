code and data in support of the fucci paper

=== Section 1: Processing microscopy data ===
Data processing workflow is as follows: Sample preparation --> Filming --> Automated selection of single cell wells --> Manual annotation --> Input tables --> Data analysis.

This section describes the last step. The input for this step is contained in files listed below. Each file corresponds to one experiment. Note that '20111118', '20120309', and '20120316' have a slightly different table format from '20120413', but this fact is taken care of during importing. In these tables and derived data structures, times are measured in hours.

data/
20111118_division_data.csv
20120309_division_data.csv
20120316_division_data.csv
20120413_division_data.csv


Script [runImportAnnotated.m] reads the *.csv files and yields Matlab data files for convenience of further processing. Note that this script imports one file at a time. To choose a different *.csv input, uncomment the corresponding line (lines 12--15, e.g., "par.sExpName = '20111118'; ..."). The script discards any inconsistent annotations by comparing event timings for individual cells (e.g., the cell is discarded if time to green off is earlier than time to green on). The output of the script is as follows.

data/
20111118-cycle-meas.mat
20120309-cycle-meas.mat
20120316-cycle-meas.mat
20120413-cycle-meas.mat


Script [runProcessFucciData.m] is then processes data, produces figures and prints some statistics. This script also produces intermediate data files.

data/
20111118-info.mat
20120309-info.mat
20120316-info.mat
20120413-info.mat


=== Section 2: Processing FACS data ===
File "data/SH1.93 G2M summary.xlsx" contains gated FACS results from the BrdU treatment experiment. "All in division" tab contains results for a dividing population (i.e., undivided cells are excluded). In this tab, the numbers are percentages out of total number of live dividing cells. For each gate, columns show triplicates.

Script [runStretchedG2M.m] reads this data file (and saves data structures in an intermediate file named "data/BrdU-data.mat"), fits the stretched G2M model to the data, and produces figures.


=== Legal stuff ===
Copyright (c) 2013 The Walter & Eliza Hall Institute of Medical Research

This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.

If a copy of the MPL was not distributed with this file, You can obtain one at http://mozilla.org/MPL/2.0/
