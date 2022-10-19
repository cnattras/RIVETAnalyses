BEGIN PLOT /PHENIX_2001_I562409/d01-x01-y01
Title=$\frac{1}{N_{Event}2\pi p_{T}}\frac{d^{2}N}{dp_{T}dy}$ of $\pi^0$ mesons in $PbSc$ in Au+Au at $\sqrt{s_{NN}}=$ 130GeV ($0-10\%$)
XLabel=$p_{T}$(GeV/c)
YLabel=$\frac{1}{N_{Event}2\pi p_{T}}\frac{d^{2}N}{dp_{T}dy}$
END PLOT

BEGIN PLOT /PHENIX_2001_I562409/d02-x01-y01
Title=$\frac{1}{N_{Event}2\pi p_{T}}\frac{d^{2}N}{dp_{T}dy}$ of $\pi^0$ mesons in $PbSc$ in Au+Au at $\sqrt{s_{NN}}=$ 130GeV ($60-80\%$)
XLabel=$p_{T}$(GeV/c)
YLabel=$\frac{1}{N_{Event}2\pi p_{T}}\frac{d^{2}N}{dp_{T}dy}$
END PLOT

BEGIN PLOT /PHENIX_2001_I562409/d03-x01-y01
Title=$\frac{1}{N_{Event}2\pi p_{T}}\frac{d^{2}N}{dp_{T}dy}$ of $\pi^0$ mesons in $PbGl$ in Au+Au at $\sqrt{s_{NN}}=$ 130GeV ($0-10\%$)
XLabel=$p_{T}$(GeV/c)
YLabel=$\frac{1}{N_{Event}2\pi p_{T}}\frac{d^{2}N}{dp_{T}dy}$
END PLOT

BEGIN PLOT /PHENIX_2001_I562409/d04-x01-y01
Title=$\frac{1}{N_{Event}2\pi p_{T}}\frac{d^{2}N}{dp_{T}dy}$ of $\pi^0$ mesons in $PbGl$ in Au+Au at $\sqrt{s_{NN}}=$ 130GeV ($60-80\%$)
XLabel=$p_{T}$(GeV/c)
YLabel=$\frac{1}{N_{Event}2\pi p_{T}}\frac{d^{2}N}{dp_{T}dy}$
END PLOT

BEGIN PLOT /PHENIX_2001_I562409/d05-x01-y01
Title=$\frac{1}{N_{Event}2\pi p_{T}}\frac{d^{2}N}{dp_{T}dy}$ of $h^{\pm}$ in in Au+Au at $\sqrt{s_{NN}}=$ 130GeV ($0-10\%$)
XLabel=$p_{T}$(GeV/c)
YLabel=$\frac{1}{N_{Event}2\pi p_{T}}\frac{d^{2}N}{dp_{T}dy}$
END PLOT

BEGIN PLOT /PHENIX_2001_I562409/d06-x01-y01
Title=$\frac{1}{N_{Event}2\pi p_{T}}\frac{d^{2}N}{dp_{T}dy}$ of $h^{\pm}$ in in Au+Au at $\sqrt{s_{NN}}=$ 130GeV ($60-80\%$)
XLabel=$p_{T}$(GeV/c)
YLabel=$\frac{1}{N_{Event}2\pi p_{T}}\frac{d^{2}N}{dp_{T}dy}$
END PLOT

BEGIN PLOT /PHENIX_2001_I562409/d07-x01-y01
Title=$PbSc$ versus $R_{AA}$  
XLabel=$PbSc$
YLabel=$R_{AA}$
END PLOT

BEGIN PLOT /PHENIX_2001_I562409/d08-x01-y01
Title=$PbGl$ versus $R_{AA}$ 
XLabel=$PbGl$
YLabel=$R_{AA}$
END PLOT

BEGIN PLOT /PHENIX_2001_I562409/d09-x01-y01
Title=$charged hadrons$ versus $R_{AA}$ 
XLabel=$charged hadrons$
YLabel=$R_{AA}$
END PLOT

BEGIN PLOT /PHENIX_2001_I562409/d10-x01-y01
Title=$PbSc$ ratio of $central$/$peripheral$
XLabel=$central$
YLabel=$peripheral$
END PLOT

BEGIN PLOT /PHENIX_2001_I562409/d11-x01-y01
Title=$PbGl$ ratio of $central$/$peripheral$
XLabel=$central$
YLabel=$peripheral$
END PLOT

BEGIN PLOT /PHENIX_2001_I562409/d12-x01-y01
Title=$charged hadrons$ ratio of $central$/$peripheral$
XLabel=$central$
YLabel=$peripheral$
END PLOT

list = [d03-x01-y01 d04-x01-y01 d05-x01-y01 d06-x01-y01 d07-x01-y01 d08-x01-y01 d09-x01-y01 d10-x01-y01 d11-x01-y01 d12-x01-y01]
DrawOnly=list

