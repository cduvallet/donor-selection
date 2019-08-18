# donor-selection

Code and data accompanying PLOS ONE publication "Framework for rational donor selection in fecal microbiota transplant trials." Duvallet et al. 2019.

## Code

All figures and analyses in the paper can be reproduced through the code in this repo. I've written a Makefile that outlines all of the steps I took, but it has unfortunately not been tested or debugged because I finished my PhD before I was able to finalize this work as fully reproducible. However, the Makefile and associated comments within it should hopefully provide enough guidance for anyone who wishes to reproduce the steps I took. It'll just be a bit hacky, that's all.

I did check that all the figure-generating notebooks run without errors. Hopefully all the files required to run the notebooks and generate the figures are clear.

I also didn't gather all the required packages for this project, but for the most part they are pretty standard. The only non-standard package I used is `feather format`, which makes reading and writing large files much faster. Also, I worked in Python 2. Again, I didn't do many fancy things here so I doubt there will be issues if you're using Python 3 but just an FYI.

A few notes on the specific parts of this work:

## Figure 2 (IBD case study)

The data that I used for the IBD case study was downloaded by another author on this paper, and I did not include that process in the Makefile. The instructions and relevant files are indicated in comments in the Makefile.

Similarly, I did not include the bash scripts with all the QIIME 2 steps to process the raw data in the Makefile, but they are in the `src/data/` folder and should be able to run from start to finish once you have activated your QIIME 2 environment.

Also note that the code sometimes changes from assuming you are in this main director vs. in the sub-directory that the data is in.

I didn't include these mostly because I didn't want to figure out how to encode the file dependencies in `make`, not because there are lots of intermediate steps that you have to take (hopefully).

## Figure 3 (liver cirrhosis / metabolomics case study)

The metabolomics data in the liver cirrhosis case study was acquired after personal communication with the authors. They have since made the data available to download, so you should grab it from their repo. It is quite possible that the data they've uploaded is not compatible with my code, so apologies in advance for any changes to the code you have to make there. Hopefully the steps I take in `src/exploration/2018-11-14.bn10_metabolomics.tidy_julians_new_data.ipynb` (and the copied over `src/data/bn10.tidy_data.py`) are clear enough for you to get your data into the same tidy format that I used to make my figures.

## Figure 4 (power simulation)

I think this one should actually run from start to finish without any problems. You can try running `make fig4`.

## Supplementary figures

These should also run using `make` once you've gotten the other figures sorted.

A note that in the two FMT studies that had multiple measured outcomes (i.e. remission and response; month 1 and month 6), I've also included the butyrate abundance vs. FMT response figures for all the other outcomes in this repo (but not in the paper). These are made by the same notebook that makes the patient matching Supp Fig 4, `src/exploration/2018-10-29.figures.ibd_butyrate_match_patient_donors.ipnyb` (which is also copied into the `src/figures/` folder).

## Old stuff

Anything you see that isn't in the Makefile should be related to old analyses that I did. Specifically, the Hsiao and David datasets were cholera datasets from our original intention with this project (help OB think about how to do donor selection for FMT studies in developing countries).
