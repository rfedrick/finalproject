# Phylogenetic Biology - Final Project

## Guidelines - you can delete this section before submission

This is a stub for your final project. Edit/ delete the text in this readme as needed.

There are two ways you can use this document:  
- You can download this file to a folder on your computer, edit this document and add other files (data, code, etc), and then zip up and submit the folder on canvas.
- You can for the [repository](finalproject) containing this document on gitub. Then commit and push your canges to the repository, and submit a link to the repository on canvas.

Github is a great way to work on projects, but also has a steep initial learning curve.


Some guidelines and tips:

- Use the stubs below to write up your final project. Alternatively, if you would like the writeup to be an executable document (with [knitr](http://yihui.name/knitr/), [jupytr](http://jupyter.org/), or other tools), you can create it as a separate file and put a link to it here in the readme.

- For information on formatting text files with markdown, see https://guides.github.com/features/mastering-markdown/ . You can use markdown to include images in this document by linking to files in the repository, eg `![GitHub Logo](/images/logo.png)`.

- The project must be entirely reproducible. In addition to the results, the repository must include all the data (or links to data) and code needed to reproduce the results.

- If you are working with unpublished data that you would prefer not to publicly share at this time, please contact me to discuss options. In most cases, the data can be anonymized in a way that putting them in a public repo does not compromise your other goals.

- Paste references (including urls) into the reference section, and cite them with the general format (Smith at al. 2003).

OK, here we go.

#  Genome Evolution and Karyotype in Amphibians

## Introduction and Goals

.  Genome size is a conspicuous, highly variable trait across the animal tree of life. Its functional significance isn't immediately obvious, as having more DNA does not necessarily mean the organism has a greater number of genes, nor does it mean that those genes serve a more complicated function. This is because the abundance of noncoding regions and transposons within the genome allows genome size to vary dramatically without ever changing the organism's coding genes. While not obvious, changing volumes of DNA do have significant impacts on the life histories of different organisms. Genome size is known to correlate positively with cell size, as well as with longer periods of DNA replication required before cell division. These consequences have far reaching implications on an organismal level, including increased development times, slower limb regeneration, and has even been correlated with decreasing brain complexity in frogs and salamanders (with higher cell sizes correlating with reduced complexity) (Sessions, 2008). Within vertebrates, amphibians have the greatest known variation in genome sizes, with genomes ranging from 0.95 pg in the Ornate Burrowing Frog to 120 pg in the Neuse Waterdog (Leidke et al., 2018). Within lissamphibians, Caudatan genomes are significantly larger than the other two clades, with the smallest Caudatan genome (13.87 pg) already larger than the largest frog or caecilian genome. 
.  Karyotype describes the number of chromosomes in a genome, describing elements of its 3-D structure and sub-dividing the genome into discrete chunks. While genome size variation is primarily driven by the dynamics of the DNA within chromosomes, I want to know if changing superstructural elements like chromosomes might also have some impact on genome evolution. One obvious connection between karyotype and genome size comes through polyploidy events. When a genome duplicates, it doubles both the size of the genome and its karyotype, and while I would be interested in testing whether selective pressure for smaller genomes intensifies immediately after a polyploidy event, I do not have the framework with which to test that hypothesis.
.  In my analysis, I used diploid number to sidestep the question of polyploidy entirely, focusing my investigation on processes that would result in more gradual karyotype change. Barring polyploidy events, karyotype has been proposed to change primarily through fusion and fission of chromosomes (with centric insertion not visible via change in diploid number) (Simakov et al., 2022; Yoshida and Kitano, 2021). One of the most interesting of these processes is the idea of fusion with mixing, documented by Simakov et al., 2022, in which the genes from the two fused chromosomes mix within the newly formed chromosome. A process such as this begs the question "in this kind of dynamic rearrangement, is it possible that the evolutionary rate of genome size evolution could also be impacted?". Though they weren't correlating chromosome number to genome size, a study in the genus Lilium has found a relationship between genome size and different karyological characters (Du et al., 2017). They found relationships between genome size and karyotype asymmetry, diversity in centromeric index, and total haploid length, and while it's possible this pattern only holds true within the genus, it indicates that there is the potential for karyotypic features to influence genome size. While I didn't find any papers that cited whole chromosomal deletions or duplication as significant drivers of karyotype evolution in animals, it is also worth considering that they would directly impact the size of the genome in a similar, if less pronounced, manner as whole genome duplication events.
.  In my project, I investigate the impacts of karyotype change on the evolution of genome size. I'm looking at changes in the diploid number of amphibians through time using ancestral state reconstruction, and correlating it with the rate of genome size evolution on specific branches of the phylogeny. To this end, MBASR was used to reconstruct the ancestral diploid number of amphibia, from which karyotype change was correlated with evolutionary rates produced though BAMM analysis. 
.  Understanding genome size variation within amphibia will allow for us to better understand the pressures driving genomic evolution. A potential interaction with karyotype would reinforce the significance of structural elements of the genome, and the idea that better understanding these structures is an integral part of genomics moving forward.


## Methods

  Data on C-value (measured in pg) was taken from Animal Genome Size Database (curated by Lietke at al., 2018), and Karyotype information was pulled from the Amphibian Karyotype Database (Perkins et al., n.d.). These values were then pruned into two datasets, one for each reference phylogeny (Pyron, 2014 and Hime et al., 2020), containing all species for which C-value and Karyotype were known.
  These datasets were then further split between salamanders and all other amphibians in order to analyze variability in evolutionary rates of genome size using BAMM (Bayesian Analysis of Macroevolutionary Mixtures, n.d.). This splitting was necessary to account for a rapid increase in genome size at the base of Caudata (Lietke at al., 2018). The package BAMMtools was used to estimate the priors and parameters needed to perform MCMC, and these analyses were then run on the Farnam Cluster. Once the analyses were completed, the results were visualized using functions from BAMMtools.
  The Ancestral state reconstruction from karyotype data was performed using the MBASR (MrBayes Ancestral States with R) script designed to streamline the process of ancestral state reconstruction (Heritage, 2021). This toolkit allows for karyotype to be reconstructed as an ordered variable, but is limiting in that it only allows for up to 6 discrete ordered states to be specified. This wasn't a catastrophic limitation, as the majority of Amphibian karyotypes had 2n numbers of 24, 26, or 28, but it did force me to flatten a significant amount of karyotype diversity above and below the mean into just three categories. If given more time, I would have preferred to use a more well-established program for karyotype ancestral state reconstruction such as Mesquite, but this MBASR reconstruction should be sufficient for identifying most branches likely to represent shifts in karyotype (Takagui et al., 2021). 
  To test whether karyotype shifts could be impacting genome size evolution, branches where MBASR predicted a change in karyotype were identified, then their corresponding branches from the BAMM analysis (minus the branch connecting Caudata to the rest of amphibia) were segregated from branches where no shift was predicted. One-way ANOVA and the non-parametric Kruskal–Wallis test were then performed on the two groups to test whether their means are significantly different from each other, which would indicate that karyotype has a relationship to genome size.
  

## Results
BAMM
  Once BAMM analyses were split into Caudata and the rest of amphibia, all runs resulted in zero predicted significant rate shifts being the most likely scenario. Within the Pyron phylogeny excluding salamanders (Figure 1A), the only group with noticeably accelerated evolution was the clade defined by Neobactrachus and Crinia, a group of Australian frogs. Salamanders, on the other hand, didn't have a particular clade with significantly accelerated rates of genome size evolution, but rather a general trend towards accelerated evolution at the root, which slows down as you get closer to the tips (Figure 1B). The Hime phylogeny analysis revealed a similar pattern in salamanders (but also showed an increased rate of genome size evolution in plethodontidae) (Figure 2B), while the group without salamanders showed a more varied spread of evolutionary rates across frogs and caecilians (Figure 2A).

MBASR
  Ancestral state reconstruction within the MBASR framework predicted an ancestral diploid number of 26 when run for both the Pyron (Figure 3) and Hime (Figure 4) phylogenies. 57 branches in the Pyron phylogeny were identified to have potential karyotype changes (12 in salamanders, 45 outside of salamanders), and 30 such branches (7 in salamanders, 23 outside of salamanders) were identified in the Pyron phylogeny.

ANOVA
  Interestingly, the ANOVA and Krukskal-Wallis tests between branches with and without predicted karyotype change both found that the Hime phylogeny (Figure 5A) showed no signs of karyotype impacting the rates of genome size evolution (p > .8), whereas the larger Pyron phylogeny (Figure 5B) returned a highly significant p-value (ANOVA: p = 0.000109; Krukskal-Wallis: p = 7.97e-06) from both tests, which would suggested the rate of genome size evolution was greater in branches that were predicted to have experienced karyotype change compared to those predicted not to. 
Alternate_universe ANOVA
After the noticing the distribution of rates returned by the BAMM analysis of the Hime phylogeny only including salamanders looked very different from the other three runs, I ran the exact same protocol to generate a second set of BAMM results for the Hime phylogeny only including salamanders. This time, the distribution of rates matched what was seen on the other three runs, and when added together with the rates from Hime phylogeny without salamanders. This resulted in a p value just under .05 (ANOVA: p = 0.0361; Krukskal-Wallis: p = 0.04332), now agreeing with the results of the ANOVA on the Pyron phylogeny that suggested the rate of genome size evolution was greater in branches that were predicted to have experienced karyotype change. 
*See full analysis in folder “Alternate_universe Hime Salamander results”
  
## Discussion

  Based on these results, it is very difficult to make an assertion about whether or not karyotype change is likely to have a relationship with genome size. It is striking that two different phylogenies (using C-value and karyotype data from the same two databases) could result in such dramatically different interpretations on karyotypic influence in genome size evolution. The p-values from ANOVAs performed on the Pyron and Hime phylogenies are over 5 orders of magnitude apart, a result hard to attribute just to the difference in data quality between the two phylogenies. 
  The most apparent difference between these two phylogenies is the number of taxa used. Pyron was building off the phylogeny he created in Pyron and Wiens, 2011, which already included 2,871 species, compared to the modest 286 taxa included in Hime et al., 2020. Even when all tips without full genome size and karyotype data were dropped, Pryon's phylogeny still retained 273 tips compared to Hime's 59. Another major difference between the approaches taken by these two phylogenies lies in the amount of data from each species used to construct the tree. Hime's phylogeny used a set of 291,282 aligned nucleotide positions across 22 loci to construct their tree, whereas Pyron's tree was only drawing from a maximum of 12,712 bp per species. Part of the reason for these contrasting approaches comes down to Pyron's main goal being to resolve family-level taxonomic classifications, whereas Hime's team set out to better resolve the deeper nodes of the Amphibian family tree.
  One follow-up analysis that I'd like to perform is to pare down both trees until only taxa shared by both trees remain. Then I would repeat all the above analyses to see if the results still show a significant difference between the two phylogenies. If the reason why one tree results in significance and the other doesn't is because of the greater amount of data in Pyron's tree leading to a more confident statistic, this should equalize the playing field and result in a similar p-value for both phylogenies. A similar result between these pared down phylogenies could also indicate the discrepancy was driven by taxon sampling differences. And if this discrepancy persists, it could point to a more fundamental difference about the way the phylogenies are being interpreted by BAMM, or a difference in topography that dramatically alters the reconstructed histories of these traits. 
  Personally, I expect that there was an error with how I set up my BAMM analyses. Particularly, when looking at the mean rates of all branches produced by BAMM, the shape of the data for the Hime phylogeny only including salamanders stood out as an outlier compared to the other 3 (see histogram "Hime rates with Salamanders" in Phylogenies_project.html), so it is possible that this run was not performed properly (whether that be because of the MCMC or the priors). Since, in my experience BAMM runs have been somewhat inconsistent, another point of future investigation might be to essentially bootstrap BAMM to give me a better idea of the variability I'm dealing with.

# Alternate Run
  After noticing the strange histogram formed by the rates in the Hime salamanders only BAMM run, I decided to run the exact same code in BAMM again, and sure enough I got a different result. In this new potential distribution of rates for Hime salamanders, the branches with predicted karyotype change were now also statistically significantly larger than those without predicted karyotype change. The data looks better visually, and the results agree more closely with the results from the Pyron phylogeny, but this markedly different result from same input data makes me skeptical about how much meaning that can be drawn from any single BAMM run. In my mind, this reinforces the need for something like a bootstrap test be developed to evaluate the results of any given BAMM run.


## References
Bayesian Analysis of Macroevolutionary Mixtures. Welcome - bamm 2.5.0 documentation. Available at: http://bamm-project.org/ (Accessed: December 15, 2022).

Du, Y. P., Bi, Y., Zhang, M. F., Yang, F. P., Jia, G. X., & Zhang, X. H. (2017). Genome size diversity in Lilium (Liliaceae) is correlated with karyotype and environmental traits. Frontiers in plant science, 8, 1303. https://www.frontiersin.org/articles/10.3389/fpls.2017.01303/full

Gregory, R.T. Animal Genome Size Database. Available at: http://www.genomesize.com/ (Accessed: December 15, 2022).

Heritage, S. (2021). MBASR: Workflow-simplified ancestral state reconstruction of discrete traits with MrBayes in the R environment. bioRxiv. https://www.biorxiv.org/content/10.1101/2021.01.10.426107v1.full

Liedtke, H.C., Gower, D.J., Wilkinson, M. et al. Macroevolutionary shift in the size of amphibian genomes and the role of life history and climate. Nat Ecol Evol 2, 1792–1799 (2018). https://doi.org/10.1038/s41559-018-0674-4 

Hime, P. M., Lemmon, A. R., Lemmon, E. C. M., Prendini, E., Brown, J. M., Thomson, R. C., ... & Weisrock, D. W. (2021). Phylogenomics reveals ancient gene tree discordance in the amphibian tree of life. Systematic biology, 70(1), 49-66. https://academic.oup.com/sysbio/article/70/1/49/5828236#226184871

Pyron, R. A., & Wiens, J. J. (2011). A large-scale phylogeny of Amphibia including over 2800 species, and a revised classification of extant frogs, salamanders, and caecilians. Molecular phylogenetics and evolution, 61(2), 543-583. https://www.sciencedirect.com/science/article/pii/S105579031100279X#f0015

Pyron, R. A. (2014). Biogeographic analysis reveals ancient continental vicariance and recent oceanic dispersal in amphibians. Systematic biology, 63(5), 779-797. https://academic.oup.com/sysbio/article/63/5/779/2847944#115269521

R.D. Perkins, J.R. Gamboa, M.M. Jonika, J. Lo, A. Shum, R.H. Adams, H. Blackmon. The Amphibian Karyotype Database. Chromosome Research. https://doi.org/10.1007/s10577-019-09613-1

Sessions, S. K. (2008). Evolutionary cytogenetics in salamanders. Chromosome Research, 16(1), 183-201. https://link.springer.com/article/10.1007/s10577-007-1205-3

Simakov, O., Bredeson, J., Berkoff, K., Marletaz, F., Mitros, T., Schultz, D. T., ... & Rokhsar, D. S. (2022). Deeply conserved synteny and the evolution of metazoan chromosomes. Science advances, 8(5), eabi5884. https://www.science.org/doi/10.1126/sciadv.abi5884

Takagui, F. H., Viana, P., Baumgärtner, L., Bitencourt, J. A., Margarido, V. P., Lui, R. L., ... & Giuliano-Caetano, L. (2021). Reconstruction of the Doradinae (Siluriformes-Doradidae) ancestral diploid number and NOR pattern reveals new insights about the karyotypic diversification of the Neotropical thorny catfishes. Genetics and Molecular Biology, 44. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8612126/

Yoshida, K., & Kitano, J. (2021). Tempo and mode in karyotype evolution revealed by a probabilistic model incorporating both chromosome number and morphology. PLoS genetics, 17(4), e1009502. https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009502
![image](https://user-images.githubusercontent.com/120404559/208016239-d7ad79ad-9809-40b1-aa95-6a6657a8bb2e.png)


