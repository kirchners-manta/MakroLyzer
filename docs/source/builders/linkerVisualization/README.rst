:orphan:

Linker drawings
===============

Activate an environment with RDKit (MakroLyzer's ``builder`` extra). From
the repository root, generate the first drawing with:

.. code-block:: bash

   python docs/source/builders/linkerVisualization/draw_linker.py

The script reads ``beta-O-4.xyz`` beside itself and writes ``pics/beta-O-4.svg``.
It also accepts an XYZ path and an optional output path:

.. code-block:: bash

   python docs/source/builders/linkerVisualization/draw_linker.py input.xyz -o drawing.svg

The same script draws all supplied linker structures. For example, generate
the phenylcoumaran drawing from the repository root:

.. code-block:: bash

   python docs/source/builders/linkerVisualization/draw_linker.py docs/source/builders/linkerVisualization/phenylcoumaran-beta-5.xyz

This writes ``pics/phenylcoumaran-beta-5.svg`` relative to the input file's
directory, creating ``pics/`` if needed. ``--output`` overrides this location.

XYZ files must contain a complete, neutral, closed-shell molecule, including
hydrogens, with coordinates in angstroms. RDKit infers bonds and bond orders
from those coordinates; inspect the resulting chemistry before publishing
drawings from new inputs. See RDKit's `bond inference documentation
<https://www.rdkit.org/docs/source/rdkit.Chem.rdDetermineBonds.html>`_.

The supplied beta-O-4 example is a G chain unit joined to an incoming H unit.
The biphenyl ether example joins incoming G at C5 to chain S at O4.
Its yellow core covers C5--O--C4, with purple C5 and pink C4 labels.
Its panel excludes incoming S, which lacks C5-H, and requires a free chain
phenolic OH. The existing chain O4 becomes the ether oxygen.
Yellow highlights cover all three atoms and both bonds of its
C-beta--O--C4 path. Greek position labels identify the participating carbons.
For phenylcoumaran and resinol, the complete cyclic cores are highlighted,
including their ether closures. For spirodienone, the highlight covers the
spiro ether ring, dienone ring and carbonyl, and the existing beta-O-4 path.
For the acyclic diarylpropane diol, the highlight covers the incoming
alpha/beta/gamma backbone and the retained chain ring's C1. The supplied
example has an incoming S unit and a retained H ring. Its panel describes
direct product construction: only an isolated chain tail can be removed,
and coproducts are omitted.
Primes distinguish participating units: resinol uses primed/unprimed alpha,
beta and gamma; spirodienone uses primes for the dienone unit and double
primes for its beta-ether partner. These are local figure labels.

Position labels use purple (``#785EF0``) for the incoming unit and pink
(``#DC267F``) for the chain unit. In the spirodienone trimer, orange
(``#FE6100``) labels identify the third unit (the existing
beta-O-4 partner). Atom symbols remain black and highlights remain yellow.
For the topologically symmetric resinol core, the first matched half is
assigned to the chain, following the builder atom order of the supplied XYZ.
XYZ does not itself record which unit was added first.

Each recognized motif includes a text panel below the structure: H/G/S
availability for incoming and chain units, required free sites before adding
the linker, and special conditions. Header colours match the position labels.
The phenylcoumaran panel excludes chain-side S; the spirodienone panel records
its beta-O-4 precursor requirement and the third partner. Panel wording is
maintained in ``PANEL_CONTENT`` in the script and describes local eligibility;
an existing chain must still leave those functional groups available.

The drawings and panels request Arial, with Helvetica/sans-serif fallbacks.
Table rules use gently uneven pen-like strokes. Aligned H/G/S rows associate each availability tick with its unit type.
The left column has no heading, and there are no rules between these rows.
The lower cells repeat the role names as Incoming free sites and Chain free
sites, using the same colours as their column headings and structure labels.
All panel text is centred in its column. The top rule and separator before Special conditions use thicker strokes.
The thin divider below S stops at the Special conditions column.
Panel text is 20 px; availability ticks/crosses remain 25 px. Multiple
requirements or special conditions use one bullet per complete statement,
with long statements wrapped onto continuation lines. A single statement
is displayed without a bullet. Store complete statements in ``PANEL_CONTENT``;
the renderer handles line wrapping.

Motifs are recognized by connectivity rather than filenames or atom order.
Other structures can be drawn, but only the six supplied motifs have
highlight and label rules. The SVG has
a shared 920 x 690 pixel canvas with the same table position in each figure.
Structures are centred above the table at a consistent bond/label scale;
larger inputs are reduced uniformly to fit. It uses skeletal notation with OH groups and alternating aromatic double bonds,
without lone-pair dots or assigned stereochemistry. It preserves the
specific H/G substituents in the XYZ; generic X/Y labels are not yet added.

SVGs can be embedded directly in Sphinx with an ``image`` or ``figure``
directive. RDKit is needed to regenerate them, rather than to display them
or build the documentation. The gallery is embedded in ``../01Lignin.rst``;
its CSS/JavaScript live under ``docs/source/_static/`` and are registered in
``conf.py``. It supports arrow buttons, keyboard navigation and wraparound.
Without JavaScript and in print, all figures remain visible in sequence.
