Scripts for the manuscript entitled "Rhythm in the premature neonate brain: Very early processing of auditory beat and meter". 

Note that the Matlab analysis script requires a running version of the EEGLAB toolbox (https://github.com/sccn/eeglab.git), Fieldtrip toolbox (https://github.com/fieldtrip/fieldtrip.git), and circStat toolbox (https://github.com/circstat/circstat-matlab.git).

There are three .mat files to run the Stat.m: 

DupleTripleRhythmMat.mat: subjects * frequencies (1,1.5,2,2.5,3 Hz, and averaged of un-related frequencies)

QuadrupleRhythmMat.mat: subjects * frequencies (0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3, and averaged of un-related frequencies)

Phase: two cells for Duple\Triple Rhythm and Quadruple Rhythm

Phase{1}: channels * subjects * frequencies (1,1.5,3 Hz)

Phase{2}: channels * subjects * frequencies (0.75,3 Hz)
