# Current commands
    rivet-mkhtml --pwd -o plots_test  Rivet_CuCu.yoda:'legend label'

    rivet-merge  -O beam -o Rivet_final.yoda Cu.yoda p.yoda

    To plot it:
    rivet-mkhtml --pwd -o plots_test  Rivet_pp.yoda:'legend label'

    rivet-merge -O cent -O beam -O Rivet_final.yoda Rivet_CuCu.yoda Rivet_pp.yoda Rivet_hepmc_pp_200GeV_*.hepmc.yoda
    rivet-merge  -O beam -o Rivet_final.yoda Rivet_CuCu.yoda Rivet_pp.yoda 

    rivet-merge -o Rivet_final.yoda Rivet_CuCu.yoda Rivet_pp.yoda

# Questions to ask:

Commands: 


Note:
    When using `book(` function:
        - !! Any name with a "_" prepend, will not be plotted later in the rivet-mkhtml
        - Can us the `refData` function to select the data location without
          giving a default name
          i.e. (I think)
          `book(_h["name"], "_this_name", refData(7, 1, 2));`
          (The _ in _this_name above says -- "don't plot this")


 1.  Why `->fill()` appears to not take a user weighting?
    -- by design: will always use the weights in ->fill() from book,
       can't really get around it. Must fix *hepmc input file.
 1.  Why bin edges on manually binned histograms aren't matching up
    -- probably NA -- using "_" prefix, along with mkAxisCode() and refData()
     allows better flexibility in input binning
 1.  Am I having trouble from manually made histograms in the `finalize()` section
    -- This probably results from 


        * Q: Are there standard libraries for triggers in RIVET?
          A: No, these are corrected for in the experiment, prior to comparison.

 # Possible improvents:

 1.  Maybe modify my file with perl to change all weights to 1



Note:
	I can use .gzip input files, but it runs a little faster if I unzip the files locally.
