## Lesson template

This lesson is derived from the
[Carpentries template](https://github.com/carpentries/lesson-example).
See [intructions](https://carpentries.github.io/lesson-example/setup.html)
for how to build and preview it on your local computer.

I found that for `make serve` to work on Windows, I had to create a user
environmental variable called `R_LIBS_USER` and point it to the path to
my package library in my Documents folder.  It also seemed to work better to
install all required packages before running `make serve` rather than letting
`make serve` attempt to do it.

## Genomics files

See the [setup document](setup.md) for the genome and annotation files that
are needed.  You should put these in `_episodes_rmd/data` in order to build the
site.
