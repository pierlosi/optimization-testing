.. This file is part of the COMANDO project which is released under the MIT
.. license. See file LICENSE for full license details.
..
.. AUTHOR: Marco Langiu
Welcome to COMANDO's documentation!
===================================

COMANDO is a next generation modeling framework for **Component-Oriented
Modeling AND Optimization** of the design and operation of energy systems.
An energy system is considered to be a collection of different interconnected
components whose purpose is to satisfy demands of various commodities such as,
e.g., electric power, heat and cooling, in a variety of different operating
conditions.

When such a system is built (or extended), there are many design and
operational decisions that need to be made.
In this context, optimizing the design and operation means finding a set of
decisions that results in a minimal value for some generalized costs, taking
into account restrictions imposed by individual components, their connections
or by other safety-related, social, political, or economic considerations.

COMANDO provides means to...

- model existing energy systems and possible extensions in a flexible, component-oriented fashion.
- use the resulting system models to create mathematical optimization problems
- solve these problems directly or use tools to automatically reformulate them to a form, more amenable to solution

To determine whether COMANDO suits your needs and avoid confusion when
something doesn't work as expected, it is important to first understand some
basic concepts and terminology introduced by COMANDO.
We therefore recommend that you read the background section of this
documentation before attempting to start hacking away.
If you feel brave and know what you're doing, you can also skip ahead to the
quickstart section, instead.

When using COMANDO in an academic context please cite our preprint_ on arXiv.org:

.. code:: bibtex

   @Article{langiu2021comando,
     author        = {Marco Langiu and David Yang Shu and Florian Joseph Baader and Dominik Hering and Uwe Bau and Andr\'{e} Xhonneux and Dirk M\"uller and Andr\'e Bardow and Alexander Mitsos and Manuel Dahmen},
     title         = {COMANDO: A Next-Generation Open-Source Framework for Energy SystemsOptimization},
     howpublished  = {\url{https://arxiv.org/abs/2102.02057v1}},
     eprint        = {http://arxiv.org/abs/1707.02514v4},
     eprintclass   = {math.OC},
     eprinttype    = {arXiv},
     optauthor     = {M.~Langiu and D.~Y.~Shu and F.~J.~Baader and D.~Hering and U.~Bau and A.~Xhonneux and D.~M\"uller and A.~Bardow and A.~Mitsos and M.~Dahmen},
     optnote       = {submitted on 04Feb2021 to cace},
     primaryclass  = {math.OC},
     year          = {2021},
   }

.. todo::

  - Extend :ref:`quickstart` examples for all 3 modeling phases.
  - Find out how to change the name of the 'Parameters' and 'Variables' sections in the API to display 'Arguments' and 'Attributes' instead, to avoid confusion with :class:`comando.Parameter` and :class:`comando.Variable`!


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   background.rst
   install.rst
   quickstart.rst
   interfaces
   api/comando.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. _preprint: https://arxiv.org/abs/2102.02057
