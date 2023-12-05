<!-- This file is part of the COMANDO project which is released under the MIT
license. See file LICENSE for full license details.

AUTHOR: Marco Langiu -->
# Building the documentation

This directory contains everything you need to build the documentation of COMANDO, which is also hosted on [readthedocs](https://comando.readthedocs.io/en/latest) locally.

The tools required to generate the documentation can be installed by using the `dev` or `all` features upon installation (see [main README.md](../README.md)).
This documentation is built with sphinx and contains both hand-written and automatically generated content.

To build the documentation run
```shell
# in comando/docs...
python -m sphinx -T -b <target> -d build/doctrees -D language=en source/ build/<target>
```
where `<target>` is `html` or `singlehtml` for documentation you can view using a prowser or `latex` for the latex files required to build a pdf version of the documentation using `pdflatex`.

Leaving `<target>` empty will build the default target which is `html`.
The `html` documentation will be created in the directory `build\html` in which you will find the entry point `index.html`.