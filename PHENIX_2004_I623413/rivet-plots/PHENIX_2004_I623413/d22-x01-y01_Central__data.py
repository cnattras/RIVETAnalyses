
import numpy as np
from numpy import nan, inf

add_legend_handle = [
  'Rivet.yoda [beam=AUAU130,cent=GEN]'
]

xpoints = {
    'Rivet.yoda [beam=AUAU130,cent=GEN]' : [4.500000e-01, 5.500000e-01, 7.000000e-01, 9.000000e-01, 1.200000e+00,
                                            1.600000e+00],
}

xedges = {
    'Rivet.yoda [beam=AUAU130,cent=GEN]' : [4.000000e-01, 5.000000e-01, 6.000000e-01, 8.000000e-01, 1.000000e+00,
                                            1.400000e+00, 1.800000e+00],
}

ref_xerrs = [
  [abs(xpoints['Rivet.yoda [beam=AUAU130,cent=GEN]'][i]   - xedges['Rivet.yoda [beam=AUAU130,cent=GEN]'][i]) for i in range(len(xpoints['Rivet.yoda [beam=AUAU130,cent=GEN]']))],
  [abs(xedges['Rivet.yoda [beam=AUAU130,cent=GEN]'][i+1] - xpoints['Rivet.yoda [beam=AUAU130,cent=GEN]'][i]) for i in range(len(xpoints['Rivet.yoda [beam=AUAU130,cent=GEN]']))]
]

yvals = {
    'Rivet.yoda [beam=AUAU130,cent=GEN]' : [0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                                            0.000000e+00],
}

xerrs = {
    'Rivet.yoda [beam=AUAU130,cent=GEN]' : [
                                              [5.000000e-02, 5.000000e-02, 1.000000e-01, 1.000000e-01, 2.000000e-01,
                                               2.000000e-01],
                                              [5.000000e-02, 5.000000e-02, 1.000000e-01, 1.000000e-01, 2.000000e-01,
                                               2.000000e-01],
                                           ],
}

yerrs = {
    'Rivet.yoda [beam=AUAU130,cent=GEN]' : [
                                              [0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                                               0.000000e+00],
                                              [0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
                                               0.000000e+00],
                                           ],
}

variation_yvals = {
}



# lists for ratio plot
ratio0_yvals = {
    'Rivet.yoda [beam=AUAU130,cent=GEN]' : [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                            1.000000e+00],
}

ratio0_yerrs = {
    'Rivet.yoda [beam=AUAU130,cent=GEN]' : [
                                              [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                               1.000000e+00],
                                              [1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00, 1.000000e+00,
                                               1.000000e+00],
                                           ],
}

ratio0_variation_vals = {
}

ratio_band_edges = {
}
