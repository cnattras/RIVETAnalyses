//Include the following code in the scope that precedes your init() in your Rivet Analysis

//create binShift function
void binShift(YODA::Histo1D& histogram) {
    std::vector<YODA::HistoBin1D> binlist = histogram.bins();
    int n = 0;
    for (YODA::HistoBin1D bins : binlist) {
        double p_high = bins.xMax();
        double p_low = bins.xMin();
        //Now calculate f_corr
        if (bins.xMin() == binlist[0].xMin()) { //Check if we are working with first bin
            float b = 1 / (p_high - p_low) * log(binlist[0].sumW()/binlist[1].sumW());
            float f_corr = -b * (p_high - p_low) * pow(M_E, -b * (p_high+p_low) / 2) / (pow(M_E, -b * p_high) - pow(M_E, -b*p_low));
            histogram.bin(n).scaleW(f_corr);
            n += 1;
        } else if (bins.xMin() == binlist.back().xMin()){ //Check if we are working with last bin
            float b = 1 / (p_high - p_low) * log(binlist[binlist.size()-2].sumW() / binlist.back().sumW());
            float f_corr = -b * (p_high - p_low) * pow(M_E, -b * (p_high+p_low) / 2) / (pow(M_E, -b * p_high) - pow(M_E, -b*p_low));
            histogram.bin(n).scaleW(f_corr);
        } else { //Check if we are working with any middle bin
            float b = 1 / (p_high - p_low) * log(binlist[n-1].sumW() / binlist[n+1].sumW());
            float f_corr = -b * (p_high - p_low) * pow(M_E, -b * (p_high+p_low) / 2) / (pow(M_E, -b * p_high) - pow(M_E, -b*p_low));
            histogram.bin(n).scaleW(f_corr);
            n += 1;
        }
    }
}

/*
USAGE:
in finalize()

binShift(*HISTONAME);

Example pulled from PHENIX_2003_I619987 analysis:

binShift(*hProtonPt["AuAuc0010a"]);
binShift(*hPionPosPt["AuAuc0010"]);
divide(hProtonPt["AuAuc0010a"], hPionPosPt["AuAuc0010"], RatioPtoPiPos["AuAuc0010"]);

If you are doing a ratio paper and are using divide in the finalize(),
you will need to apply the binShift to both of the histograms before calling divide() (as done in the example above)
To understand the mathematics of the binShift correction, you can check out an early analysis note here:
https://drive.google.com/file/d/1auYEgcFhdk2gSoT3KrErnSDllQ-B_zE4/view?usp=sharing
*/
