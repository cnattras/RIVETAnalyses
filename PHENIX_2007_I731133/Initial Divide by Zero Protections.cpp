/*Replaced by lines 762 - 900. Change was done due to concerns that the 
scale factors kept being tacted onto the sow_pp historgram consectivelty*/
          //Figure 16: Rda
          //d07-x01-y01:
          //denominator: must do our process to our cross section as done in Figure 13 plus multiplying it by <Tda>
          binShift(*hCrossSec["ppEtadAuc0088"]);
          binShift(*hEtaPt["ptyieldsdAuc0088b"]);
          if(sow["sow_pp"]->sumW()>0){
          hCrossSec["ppEtadAuc0088"]->scaleW(1. / sow["sow_pp"]->sumW());
          hCrossSec["ppEtadAuc0088"]->scaleW(cross);
          hCrossSec["ppEtadAuc0088"]->scaleW(0.2); //scaling by <TdA> as shown in Table II
        }
         if(sow["sow_dAu"]->sumW()>0){
          //numerator: must do our process to our invariant yield like in Figure 14
          hEtaPt["ptyieldsdAuc0088b"]->scaleW(1. / sow["sow_dAu"]->sumW());
        }
        if(sow["sow_dAu"]->sumW()>0 && sow["sow_pp"]->sumW()>0){
          //Rda
          divide(hEtaPt["ptyieldsdAuc0088b"], hCrossSec["ppEtadAuc0088"], hRda["EtadAuc0088"]);
        }
          
          //d07-x01-y02
          //denominator
          binShift(*hCrossSec["ppEtadAuc0020"]);
        if(sow["sow_pp"]->sumW()>0){
          hCrossSec["ppEtadAuc0020"]->scaleW(1. / sow["sow_pp"]->sumW());
          hCrossSec["ppEtadAuc0020"]->scaleW(cross);
          hCrossSec["ppEtadAuc0020"]->scaleW(0.36); //scaling by <TdA> as shown in Table II
        }  
          //numerator
        if(sow["sow_dAuc0020"]->sumW()>0){
          binShift(*hEtaPt["ptyieldsdAuc0020b"]);
          hEtaPt["ptyieldsdAuc0020b"]->scaleW(1. / sow["sow_dAuc0020"]->sumW());
        }
          //Rda
        if(sow["sow_dAuc0020"]->sumW()>0 && sow["sow_pp"]->sumW()>0){
          divide(hEtaPt["ptyieldsdAuc0020b"], hCrossSec["ppEtadAuc0020"], hRda["EtadAuc0020"]);
        } 

          //d07-x01-y03
          //denominator
        if(sow["sow_pp"]->sumW()>0){
          binShift(*hCrossSec["ppEtadAuc2040"]);
          hCrossSec["ppEtadAuc2040"]->scaleW(1. / sow["sow_pp"]->sumW());
          hCrossSec["ppEtadAuc2040"]->scaleW(cross);
          hCrossSec["ppEtadAuc2040"]->scaleW(0.25); //scaling by <TdA> as shown in Table II
        }  
          //numerator
        if(sow["sow_dAuc2040"]->sumW()>0){
          binShift(*hEtaPt["ptyieldsdAuc2040b"]);
          hEtaPt["ptyieldsdAuc2040b"]->scaleW(1. / sow["sow_dAuc2040"]->sumW());
        }  
          //Rda
        if(sow["sow_dAuc2040"]->sumW()>0 && sow["sow_pp"]->sumW()>0){
          divide(hEtaPt["ptyieldsdAuc2040b"], hCrossSec["ppEtadAuc2040"], hRda["EtadAuc2040"]);
        }  

          //d07-x01-y04
          //denominator
        if(sow["sow_pp"]->sumW()>0){
          binShift(*hCrossSec["ppEtadAuc4060"]);
          hCrossSec["ppEtadAuc4060"]->scaleW(1. / sow["sow_pp"]->sumW());
          hCrossSec["ppEtadAuc4060"]->scaleW(cross);
          hCrossSec["ppEtadAuc4060"]->scaleW(0.17); //scaling by <TdA> as shown in Table II
        }
          //numerator
        if(sow["sow_dAuc4060"]->sumW()>0){
          binShift(*hEtaPt["ptyieldsdAuc4060b"]);
          hEtaPt["ptyieldsdAuc4060b"]->scaleW(1. / sow["sow_dAuc4060"]->sumW());
        }
          //Rda
        if(sow["sow_dAuc4060"]->sumW()>0 && sow["sow_pp"]->sumW()>0){
          divide(hEtaPt["ptyieldsdAuc4060b"], hCrossSec["ppEtadAuc4060"], hRda["EtadAuc4060"]);
        }

          //d07-x01-y05
          //denominator
        if(sow["sow_pp"]->sumW()>0){
          binShift(*hCrossSec["ppEtadAuc6088"]);
          hCrossSec["ppEtadAuc6088"]->scaleW(1. / sow["sow_pp"]->sumW());
          hCrossSec["ppEtadAuc6088"]->scaleW(cross);
          hCrossSec["ppEtadAuc6088"]->scaleW(0.073); //scaling by <TdA> as shown in Table II
        }
          //numerator
        if(sow["sow_dAuc6088"]->sumW()>0){
          binShift(*hEtaPt["ptyieldsdAuc6088b"]);
          hEtaPt["ptyieldsdAuc6088b"]->scaleW(1. / sow["sow_dAuc6088"]->sumW());
        }
          //Rda
        if(sow["sow_dAuc6088"]->sumW()>0 && sow["sow_pp"]->sumW()>0){
          divide(hEtaPt["ptyieldsdAuc6088b"], hCrossSec["ppEtadAuc6088"],hRda["EtadAuc6088"]);
        }

          //Figure 17: RAA  
          //d08-x01-y01:
          //denominator
        if (sow["sow_pp"]->sumW() > 0) {       
          binShift(*hCrossSec["ppEtaAuAuc0020"]);
          hCrossSec["ppEtaAuAuc0020"]->scaleW(1. / sow["sow_pp"]->sumW());
          hCrossSec["ppEtaAuAuc0020"]->scaleW(cross);
          hCrossSec["ppEtaAuAuc0020"]->scaleW(18.5); //scaling by <TdA> as shown in Table II
        }
          //numerator: must do our process to our invariant yield like in Figure 14
        if(sow["sow_AuAuc0020"]->sumW()>0){
          binShift(*hEtaPt["ptyieldsAuAuc0020b"]);
          hEtaPt["ptyieldsAuAuc0020b"]->scaleW(1. / sow["sow_AuAuc0020"]->sumW());
        }
          //RAA
        if(sow["sow_AuAuc0020"]->sumW()>0 && sow["sow_pp"]->sumW()>0){
          divide(hEtaPt["ptyieldsAuAuc0020b"], hCrossSec["ppEtaAuAuc0020"], hRaa["EtaAuAuc0020"]);
        }

          //d08-x01-y02
          //denominator
        if (sow["sow_pp"]->sumW() > 0) {
          binShift(*hCrossSec["ppEtaAuAuc2060"]);
          hCrossSec["ppEtaAuAuc2060"]->scaleW(1. / sow["sow_pp"]->sumW());
          hCrossSec["ppEtaAuAuc2060"]->scaleW(cross);
          hCrossSec["ppEtaAuAuc2060"]->scaleW(4.6); //scaling by <TdA> as shown in Table II
        }
          //numerator
        if(sow["sow_AuAuc2060"]->sumW()>0){
          binShift(*hEtaPt["ptyieldsAuAuc2060b"]);
          hEtaPt["ptyieldsAuAuc2060b"]->scaleW(1. / sow["sow_AuAuc2060"]->sumW());
        }
          //RAA
        if(sow["sow_AuAuc2060"]->sumW()>0 && sow["sow_pp"]->sumW()>0){
          divide(hEtaPt["ptyieldsAuAuc2060b"], hCrossSec["ppEtaAuAuc2060"], hRaa["EtaAuAuc2060"]);
        }

          //d08-x01-y03
          //denominator
        if (sow["sow_pp"]->sumW() > 0) {
          binShift(*hCrossSec["ppEtaAuAuc6092"]);
          hCrossSec["ppEtaAuAuc6092"]->scaleW(1. / sow["sow_pp"]->sumW());
          hCrossSec["ppEtaAuAuc6092"]->scaleW(cross);
          hCrossSec["ppEtaAuAuc6092"]->scaleW(0.3); //scaling by <TdA> as shown in Table II
        }
          //numerator
        if(sow["sow_AuAuc6092"]->sumW()>0){
          binShift(*hEtaPt["ptyieldsAuAuc6092b"]);
          hEtaPt["ptyieldsAuAuc6092b"]->scaleW(1. / sow["sow_AuAuc6092"]->sumW());
        }
          //RAA
        if(sow["sow_AuAuc6092"]->sumW()>0 && sow["sow_pp"]->sumW()>0){
          divide(hEtaPt["ptyieldsAuAuc6092b"], hCrossSec["ppEtaAuAuc6092"], hRaa["EtaAuAuc6092"]);
        }