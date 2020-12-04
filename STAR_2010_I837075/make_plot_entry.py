#!/usr/bin/python3

class plot_printer:
    def update_dict(self, in_dict):
        for key in in_dict:
            self.dict[key] = in_dict[key]
    def __init__(self, in_dict=dict(), fout=False):
        self.dict = dict()
        for key in ('d','x','y','Title','XLabel','YLabel'):
            self.dict[key] = ''
        self.update_dict(in_dict)
        self.fout = fout
    def has_data(self):
        for x in self.dict.values():
            if x == '':
                return False
        return True
    def print(self, comments=[]):
        d = self.dict
        string = "\n"
        for C in comments:
            string += '# ' + C + '\n'
        if '#' in d:
            for x in d['#']:
                string += '# ' + x + '\n'
        string += f"""BEGIN PLOT /STAR_2010_I837075/d0{d['d']}-x0{d['x']}-y0{d['y']}
Title={d['Title']}
XLabel={d['XLabel']}
YLabel={d['YLabel']}
END PLOT
"""
        if self.fout:
            self.fout.write(string)
        else:
            print(string)

# fout = open('STAR_2010_I837075.plot','w')
plot = plot_printer({
    'XLabel' : r"Transverse Momentum $p_{\rm T}$ [GeV/c]",
    'YLabel' : r'$\frac{1}{2\pi p_{\rm T}}\,\frac{\rm d^2 N}{\mathrm{d}p_{\rm T}\mathrm{dy}} [(\mathrm{GeV}/c)^2]$',
    'x' : 1,
    'y' : 1,
    'd' : 1,
    'Title' : "Title" })

plot.fout = open('STAR_2010_I837075.plot','w')

# plot out the pi-
for (y,cent) in zip([1,3,5,7],['0-10','10-20','20-40','40-60']):
    print(y, cent)
    plot.update_dict({
        'y':y,
        'Title': r'FIG 1(a) $\pi^-$ spectra: $\sqrt{s_\mathrm{NN}}$=200 GeV  Cu+Cu '+cent+r'\% Centrality'
    })
    plot.print([f"pi- spectra {cent}"])

for (y,cent) in zip([2,4,6,8],['0-10','10-20','20-40','40-60']):
    print(y, cent)
    plot.update_dict({
        'y':y,
        'Title': r'FIG 1(a) $\pi^+$ spectra: $\sqrt{s_\mathrm{NN}}$=200 GeV  Cu+Cu '+cent+r'\% Centrality'
    })
    plot.print([f"pi+ spectra {cent}"])

plot.update_dict({'d':2})
for (y,cent) in zip([1,3,5,7],['0-10','10-20','20-40','40-60']):
    print(y, cent)
    plot.update_dict({
        'y':y,
        'Title': r'FIG 1(b) $\bar{p}$ spectra: $\sqrt{s_\mathrm{NN}}$=200 GeV  Cu+Cu '+cent+r'\% Centrality'
    })
    plot.print([r"\bar{p}"+f" spectra {cent}"])

for (y,cent) in zip([2,4,6,8],['0-10','10-20','20-40','40-60']):
    print(y, cent)
    plot.update_dict({
        'y':y,
        'Title': r'FIG 1(b) $p$ spectra: $\sqrt{s_\mathrm{NN}}$=200 GeV  Cu+Cu '+cent+r'\% Centrality'
    })
    plot.print([f"$p$ spectra {cent}"])

#------------------------------------------------------------
plot.update_dict({ 'd':3, 'YLabel':r'Ratio: $\frac{\pi^-}{\pi^+}$' })
for (y,cent) in zip([1,2,3,4],['0-10','10-20','20-40','40-60']):
    print(y, cent)
    plot.update_dict({
        'y':y,
        'Title': r'FIG 2(a) $\sqrt{s_\mathrm{NN}}$=200 GeV  Cu+Cu '+cent+r'\% Centrality'
    })
    plot.print([f"$pi+/pi-$ spectra {cent}"])

plot.update_dict({ 'd':4, 'YLabel':r'Ratio: $\frac{\bar{p}}{p}$' })
for (y,cent) in zip([1,2,3,4],['0-10','10-20','20-40','40-60']):
    print(y, cent)
    plot.update_dict({
        'y':y,
        'Title': r'FIG 2(b) $\sqrt{s_\mathrm{NN}}$=200 GeV  Cu+Cu '+cent+r'\% Centrality'
    })
    plot.print([f"pbar/p spectra {cent}"])


# FIG 3a
plot.update_dict({ 'd':5, 'YLabel':r'$\mathrm{R}^{\pi^{+} + \pi^{-}}_{\rm AA}$' })
for (y,cent) in zip([1,2,3,4],['0-10','10-20','20-40','40-60']):
    print(y, cent)
    plot.update_dict({
        'y':y,
        'Title': r'FIG 3(a) $\sqrt{s_\mathrm{NN}}$=200 GeV  Cu+Cu '+cent+r'\% Centrality'
    })
    plot.print([f"FIG 3(a) RAA pions: {cent}"])

# FIG 3b
plot.update_dict({ 
    'd':6, 
    'y':2,
    'YLabel':r'$\mathrm{R}^{\pi^{+} + \pi^{-}}_{\rm AA}\, (5<p_\mathrm{T}<8\,\mathrm{GeV}_{c})$',
    'XLabel':r'Centrality',
    'Title': r'FIG 3(b) $\sqrt{s_\mathrm{NN}}$=200 GeV  Cu+Cu'
    })
plot.print(['FIG 3b: pion RAA w.r.t centrality'])

# FIG 4a
plot.update_dict({ 'd':7, 'YLabel':r'$\mathrm{R}^{p + \bar{p}}_{\rm AA}$' })
for (y,cent) in zip([1,2,3,4],['0-10','10-20','20-40','40-60']):
    print(y, cent)
    plot.update_dict({
        'y':y,
        'Title': r'FIG 4(a) $\sqrt{s_\mathrm{NN}}$=200 GeV  Cu+Cu '+cent+r'\% Centrality'
    })
    plot.print([f"FIG 4(a) RAA pions: {cent}"])

# FIG 5a
plot.update_dict({ 
    'd':8, 
    'YLabel' : r'$(p + \bar{p})/(\pi^++\pi^-)$',
    'XLabel' : r"Transverse Momentum $p_{\rm T}$ [GeV/c]"
    })
for (y,cent) in zip([1,2,3,4],['0-10','10-20','20-40','40-60']):
    print(y, cent)
    plot.update_dict({
        'y':y,
        'Title': r'FIG 5(a) $\sqrt{s_\mathrm{NN}}$=200 GeV  Cu+Cu '+cent+r'\% Centrality'
    })
    plot.print([f"FIG 5(a) (pbar+p)/(pi+ + pi-): {cent}"])

# FIG 5b
plot.update_dict({ 
    'd':9, 
    'y':2,
    'YLabel' : r'$(p + \bar{p})/(\pi^++\pi^-)\;(3<p_\mathrm{T}<4\,\mathrm{GeV}/c)$',
    'XLabel' : r"Centrality",
    'Title': r'FIG 5(b) $\sqrt{s_\mathrm{NN}}$=200 GeV  Cu+Cu'
    })
plot.print([f"FIG 5(b) (pbar+p)/(pi+ + pi-): {cent}"])
