# Spatial Epidemics Dynamics and Synchronization
The final project for Dr Earn's Mathematical Biology course at McMaster Univresity.

By Nicole Dumont, Melody Fong, and Carolina Weishaar.

## Purpose
This model explores the spatial synchronization of epidemics. We use the metapopulation susceptible, infectious, and removed (SIR) model with sinusoidal seasonal forcing. Our aim is to determine which parameters lead to a synchronous state of the model so that diseases can be more easily eradicated. 

## Files
Main Report: Contains a summary, background information, description of methods and models used, results, plots, and conclusion.
  * `ModelStudentsProject.pdf`: The pdf of the main report for the project. 
  * `ModelStudentsProject.tex`: The tex file of the main report.
  * `4mbapreamble.tex`: A tex file written by Dr Earn for formtting of assignments and reports in the 4MB3 course.
  * `ModelStudents.bib`: The bibilography of the report.
  * `images`: A folder conatining images used in the report.
  * `texcount.pl`: A Perl script used for generating the word count in the main report. From http://app.uio.no/ifi/texcount/
  
Supplementary Material: The R code written for the project to produce plots and descriptions of the code. Results are reproducible. 
  * `ModelStudentsSupplementaryMaterial.pdf`: A pdf of the supplementary material.
  * `ModelStudentsSupplementaryMaterial.rnw`: The knitr file containing the code. 
  * `bifurcation.ode`: Code to generate the forced SIR model bifurcation diagram using XPPAUT. See ModelStudentsSupplementaryMaterial for description.
  * `SuppBib.bib`: The bibilography of the supplementary material.
  
Beamer Presentation: A presentation created using beamer for the project. Plots are from an older version of the code.
  * `Presentation.pdf`: A pdf of the presentation.
  * `Presentation.rnw`: A knitr file of the beamer presentation
  * `images`: A folder conatining images used in the presenation.
  * `Theme`: A folder conatining metropolis beamer theme files. From https://github.com/matze/mtheme
