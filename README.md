## `colonyR`

The `colonyR` package provides an interface between R and the Colony software
for parentage analysis (see http://www.zsl.org/science/software/colony). It was
developed and has been used on a GNU/Linux platform.

I realised only after starting it that another similar project already existed
(see [here](https://r-forge.r-project.org/projects/rcolony/) for the project
page and [here](https://r-forge.r-project.org/scm/viewvc.php/?root=rcolony) for
the Subversion repository).

`colonyR` relies on having Colony installed on the system, independently from
R. Colony might be accessible either from `$PATH` or the executable can be
copied directly in the R script folder (see `?colRun` from R).

The functions documentation should be enough to get you up and running. Details
about the parameters are to be found in the pdf manual accompanying Colony
itself.

This R package has been working nicely for me until now, but be aware that I
also had troubles when trying to run it on simulated genotypes (resulting in
Colony crashes).
