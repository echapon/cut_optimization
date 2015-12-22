How to use the code:

First create the histograms:
- root -l
- .L cutoptim.C+
- cutoptim toto(some_ttree)
- toto.Loop()
- .q

If needed, do this for both the sample containing signal (e.g. MC) and the sample containing background (e.g. data)

Then draw them:
- root -l plot.C'("signal.root","background.root")'
