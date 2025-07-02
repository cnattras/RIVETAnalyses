#! /usr/bin/env python3
#This script is from Andy Buckley and is for converting YODA -> YODA2
#It should be used with some caution, as not all files will convert correctly.
#Preferred methods are to regenerate the file if it was generated from Rivet
#or to download the file again from HEPData if it was downloaded from HEPData

import yoda

aos = yoda.read("input.yoda.gz")
rtn = dict()
for aoname, ao in aos.items():
    newao = None
    if '1D' in ao.type():
        newao = yoda.Estimate0D(aoname)
        newao.setAnnotation('IsRef', 1)
        for p in ao.points():
            newao.set(p.x(), [ - p.xErrMinus(), p.xErrPlus() ])
    elif '2D' in ao.type():
        xedges = [ ao.edges(0)[0][0] ] + [ tup[1] for tup in ao.edges(0) ]
        newao = yoda.BinnedEstimate1D(xedges, aoname)
        newao.setAnnotation('IsRef', 1)
        for p in ao.points():
            newao.binAt(p.x()).set(p.y(), [ - p.yErrMinus(), p.yErrPlus() ])
    else:
        print('WARNING - file has %s' % ao.type())
    rtn[aoname] = newao
yoda.writeYODA(rtn, "output.yoda.gz")
