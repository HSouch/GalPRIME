---
title: 'GalPRIME: A Python pipeline for massively parallel extraction of galaxy light profiles'
tags:
- Python
- astronomy
- galaxies
- galaxy morphology
- stellar haloes

authors:
-   name: Harrison Souchereau
    orcid: 0000-0001-5079-1865
-   name: Ivana Damjanov
    orcid: 0000-0003-4797-5246
-   name: Marcin Sawicki
    orcid: 0000-0002-7712-7857
-   name: Devin Williams
    orcid: 0009-0002-1526-0086

affiliations:
-   name: Saint Mary's University, Canada 
    index: 1
-   name: Yale University, United States
    index: 2
date: July 2026
bibliography: paper.bib

---

# Summary

Across cosmic time, galaxies assemble mass via a broad suite of physical processes, governed by both galaxy internal properties and their environment.
This rich assembly history is encoded within their observed radial light distributions  `[@Hopkins:2010; @Hirschmann:2014; @Pillepich:2014; @Sachdeva:2017; @Elias:2018; @Pasha:2025; @Williams:2025]`.
One-dimensional surface brightness profiles (i.e., radial light profiles, SBP) enhance the signal-to-noise for these distributions, enabling studies that drive galaxy growth out to extremely low surface brightness levels (i.e. galaxy haloes) `[@Buitrago:2017; @Huang:2018; @Wang:2019; @Merritt:2020; @Cheng:2021; @Li:2021]`.


# Statement of Need


GalPRIME is designed for large-scale galaxy surveys with millions of objects.


# State of the Field

A wealth of tools currently exist for the pre-processing and subsequent extraction of light profiles from astronomical images.


# Software Design

GalPRIME is designed to be used in a 2-step approach. 

The first step is to use the tools in GalPRIME to determine how best to apply
preprocessing methods, in particular bright interloper masking and background subtraction. 
This allows the user to configure certain parameters that light profile extraction can be sensitive to, and whose best-options may vary
depending on characteristics of the the imaging survey in question (survey depth, pixel scale, etc).



# Research Impact Statement




# AI Usage Disclosure

No AI tools were used in the development, documentation, or paper authoring of this work.


# Acknowledgements

The authors thank Chris Geroux and Ross Dickson at Compute Canada, for their helpful and timely advice to improve the reliability and speed of the code. The authors also thank Song Huang, Michael Keim, Isabel Medlock, and Viraj Pandya for suggestions and helpful discussions.

We acknowledge the support of the Canada Research Chair Program and the Natural Sciences and Engineering Research Council of Canada (NSERC) grants RGPIN-2018-05425, DDG-2024-00004, RGPIN-2020-06023, and RGPAS-2020-00065.