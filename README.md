
# PAGES2k phase2 composite
Data and Matlab code for the second phase of the PAGES2k synthesis of temperature-sensitive proxies.
The main goal of this repository is to give project collaborators

## Organization
X pieces of code to know about:

- **pages2k_composite_prep.m** prepares the data for compositing, then we slice and dice the dataset in all kinds of ways:

- **pages2k_compositeByArchive**
- **pages2k_compositeByLatitudeBand
- **sensitivity to record length
- **YOUR_CONTRIBUTION_HERE

- **pages2k_composite_regression.m** produces the temperature-scaled composite with error estimates

Al the data are in "data". The figures are all exported to "figs"
Code dependencies needed to run the above codes are included in "utilities". If you are missing anything, we will be happy to add it. To comply with the license, some third-party packages may need to be installed separately:

- [export_fig](https://github.com/altmany/export_fig)


## To contribute
### Do you GitHub?
If you have gotten this far, it means that you successfully created an account. Fabulous.
The next step is to get familiar with [version control](https://backlogtool.com/git-guide/en/intro/intro1_1.html) in general, and [GitHub in particular](http://readwrite.com/2013/09/30/understanding-github-a-journey-for-beginners-part-1).

Note that the GitHub learning curve is much smoother with a [Desktop App](https://desktop.github.com/) (Mac and Windows only. If you are on Linux, you presumably can't get enough of the command line, so you'll be fine)

### Start a branch
[Start a new branch](https://help.github.com/articles/creating-and-deleting-branches-within-your-repository/), preferably indexed by your name so we can track what's going on. The project leader (Julien) will be notified that you have done so, but don't assume they have read the code unless you made a pull request (see below).

### Create a new file (or more)
Please use a very descriptive name so that it doesn't get lost with everyone else'e edits. Something like pages2k_composite_sensitivity_to_[feature of your choice].m.

### Hack away
Make figures, tables, and analyses. Explain what you've done using abundant in-line comments.

### Make a pull request
This will inform the project leader (Julien) that your changes are ready to be shared. He will review them and incorporate them to the master branch if appropriate.

### Repeat until you have assuaged all your inner doubts about the virtues of a composite
(this make take a while)
