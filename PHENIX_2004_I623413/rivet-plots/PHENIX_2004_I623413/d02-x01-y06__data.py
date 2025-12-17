
import numpy as np
from numpy import nan, inf

add_legend_handle = [
  'Data',
  'Rivet.yoda [beam=AUAU130,cent=GEN]'
]

xpoints = {
                                  'Data' : [2.500000e+00, 1.000000e+01, 2.250000e+01, 4.500000e+01, 7.600000e+01],
    'Rivet.yoda [beam=AUAU130,cent=GEN]' : [2.500000e+00, 1.000000e+01, 2.250000e+01, 4.500000e+01, 7.600000e+01],
}

xedges = {
                                  'Data' : [0.000000e+00, 5.000000e+00, 1.500000e+01, 3.000000e+01, 6.000000e+01,
                                            9.200000e+01],
    'Rivet.yoda [beam=AUAU130,cent=GEN]' : [0.000000e+00, 5.000000e+00, 1.500000e+01, 3.000000e+01, 6.000000e+01,
                                            9.200000e+01],
}

ref_xerrs = [
  [abs(xpoints['Data'][i]   - xedges['Data'][i]) for i in range(len(xpoints['Data']))],
  [abs(xedges['Data'][i+1] - xpoints['Data'][i]) for i in range(len(xpoints['Data']))]
]

yvals = {
                                  'Data' : [4.050000e+01, 3.040000e+01, 1.550000e+01, 6.200000e+00, 1.000000e+00],
    'Rivet.yoda [beam=AUAU130,cent=GEN]' : [nan, nan, nan, nan, nan],
}

xerrs = {
                                  'Data' : [
                                              [2.500000e+00, 5.000000e+00, 7.500000e+00, 1.500000e+01, 1.600000e+01],
                                              [2.500000e+00, 5.000000e+00, 7.500000e+00, 1.500000e+01, 1.600000e+01],
                                           ],
    'Rivet.yoda [beam=AUAU130,cent=GEN]' : [
                                              [2.500000e+00, 5.000000e+00, 7.500000e+00, 1.500000e+01, 1.600000e+01],
                                              [2.500000e+00, 5.000000e+00, 7.500000e+00, 1.500000e+01, 1.600000e+01],
                                           ],
}

yerrs = {
                                  'Data' : [
                                              [6.519202e+00, 4.808326e+00, 2.404163e+00, 9.486833e-01, 2.236068e-01],
                                              [6.519202e+00, 4.808326e+00, 2.404163e+00, 9.486833e-01, 2.236068e-01],
                                           ],
    'Rivet.yoda [beam=AUAU130,cent=GEN]' : [
                                              [0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00],
                                              [0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00],
                                           ],
}

variation_yvals = {
}



# lists for ratio plot
ratio0_yvals = {
                                  'Data' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00],
    'Rivet.yoda [beam=AUAU130,cent=GEN]' : [nan, nan, nan, nan, nan],
}

ratio0_yerrs = {
                                  'Data' : [
                                              [1.609680e-01, 1.581686e-01, 1.551073e-01, 1.530134e-01, 2.236068e-01],
                                              [1.609680e-01, 1.581686e-01, 1.551073e-01, 1.530134e-01, 2.236068e-01],
                                           ],
    'Rivet.yoda [beam=AUAU130,cent=GEN]' : [
                                              [0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00],
                                              [0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00],
                                           ],
}

ratio0_variation_vals = {
}

ratio_band_edges = {
}
