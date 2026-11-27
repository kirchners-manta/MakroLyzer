"""Draw a linker XYZ structure with yellow motif highlights and position labels.

Run with no arguments to render beta-O-4.xyz beside this script. XYZ must
contain all hydrogens and coordinates in angstroms for a neutral molecule:
it stores no bonds, so connectivity and bond orders are inferred from geometry.
The drawing uses skeletal notation (implicit carbon-bound H, explicit OH),
without lone-pair dots or stereochemistry inferred from the 3D conformation.
"""

import argparse
from html import unescape
from pathlib import Path
import re
from textwrap import wrap
from xml.etree import ElementTree as ET

from rdkit import Chem
from rdkit.Chem import rdDepictor, rdDetermineBonds
from rdkit.Chem.Draw import rdMolDraw2D


# Mapped query atoms define the highlighted core, rather than relying on XYZ
# atom order. Unmapped atoms constrain its chemical context and stay uncoloured.
# Label tokens a/b/g are converted to Greek SVG text after RDKit places them.
MOTIFS = (
    (
        'beta-O-4',
        '[c:4]-[O:2]-[CH:1](-[CH2]-[OH])-[CH](-[O])-[#6]',
        {1: 'b', 2: None, 4: '4'},
        {1: 'incoming', 4: 'chain'},
    ),
    (
        'phenylcoumaran-beta-5',
        '[c:5]1[c:4][O:6][CH:2](-c)[CH:1]1-[CH2]-[OH]',
        {1: 'b', 2: 'a', 4: '4', 5: '5', 6: None},
        {1: 'incoming', 2: 'incoming', 4: 'chain', 5: 'chain'},
    ),
    (
        'resinol-beta-beta',
        '[CH:1]12[CH:2](-c)[O:3][CH2:4][CH:5]1[CH:6](-c)[O:7][CH2:8]2',
        # Resinol's two halves have equivalent topology. For the supplied XYZ,
        # the first matched half is the original chain unit (builder atom order).
        {1: 'b', 2: 'a', 3: None, 4: "g'", 5: "b'", 6: "a'", 7: None, 8: 'g'},
        {1: 'chain', 2: 'chain', 4: 'incoming', 5: 'incoming',
         6: 'incoming', 8: 'chain'},
    ),
    (
        'spirodienone-beta-1',
        '[C:1]12([CH:2]([CH:12](-[CH2]-[OH])-[O:13]-[c:14])'
        '[O:3][CH:4](-c)[CH:5]1-[CH2]-[OH])'
        '[CH:6]=[C:7]([O])-[C:8](=[O:9])-[C:10]([O])=[CH:11]2',
        # Include the spiro ring, dienone ring/carbonyl, and precursor beta ether.
        # Prime denotes the dienone unit; double prime denotes its ether partner.
        {1: "1'", 2: "a'", 3: None, 4: 'a', 5: 'b', 6: None,
         7: None, 8: None, 9: None, 10: None, 11: None,
         12: "b'", 13: None, 14: "4''"},
        {1: 'chain', 2: 'chain', 4: 'incoming', 5: 'incoming',
         12: 'chain', 14: 'partner'},
    ),
)

LABEL_COLOURS = {'incoming': '#785EF0', 'chain': '#DC267F', 'partner': '#FE6100'}
FONT_FAMILY = 'Arial, Helvetica, sans-serif'

# Documentation of local eligibility, before the depicted connection is made.
# H/G/S identify unit types; existing bonds must also leave the listed sites
# available. S lacks C5-H, so it cannot supply phenylcoumaran's chain-side C5.
PANEL_CONTENT = {
    'beta-O-4': {
        'chain_types': 'HGS',
        'incoming_sites': ('β (α–β double bond)',),
        'chain_sites': ('4 (free phenolic OH)',),
        'notes': (),
    },
    'phenylcoumaran-beta-5': {
        'chain_types': 'HG',
        'incoming_sites': ('β (α–β double bond)',),
        'chain_sites': ('5 and free 4-OH',),
        'notes': ('S has no free 5 position.', 'The 4-OH also participates in ring closure.'),
    },
    'resinol-beta-beta': {
        'chain_types': 'HGS',
        'incoming_sites': ('β and free γ-OH', 'α–β double bond'),
        'chain_sites': ('β and free γ-OH', 'α–β double bond'),
        'notes': ('Both γ-OH groups form the cyclic ether bridges.',),
    },
    'spirodienone-beta-1': {
        'chain_types': 'HGS',
        'incoming_sites': ('β (α–β double bond)',),
        'chain_sites': ('1, free α-OH and 4-OH on a β-O-4 precursor',),
        'notes': ('Requires preceding β-O-4.', 
                  'Third partner: H / G / S (orange labels).'),
    },
}


def annotate_linkers(molecule):
    """Annotate supported motifs in place; return highlighted atom/bond indices.

    Highlight beta-O-4's three-atom C-beta--O--C4 path and complete cyclic
    cores for phenylcoumaran, resinol, and spirodienone. Primed labels separate
    participating units; symmetric matches are labelled only once per core.
    """
    atoms = set()
    bonds = set()
    for name, smarts, labels, roles in MOTIFS:
        query = Chem.MolFromSmarts(smarts)
        seen = set()
        for match in molecule.GetSubstructMatches(query):
            mapped = {
                atom.GetAtomMapNum(): match[atom.GetIdx()]
                for atom in query.GetAtoms() if atom.GetAtomMapNum()
            }
            core = frozenset(mapped.values())
            if core in seen:
                continue
            seen.add(core)
            # Later, more specific motifs take precedence over their precursor
            # substructures when choosing the information panel (spirodienone).
            molecule.SetProp('_drawing_motif', name)
            atoms.update(core)
            # Transfer only bonds belonging to the mapped query core; this
            # avoids highlighting unrelated bonds between neighboring motifs.
            for bond in query.GetBonds():
                start, end = bond.GetBeginAtom(), bond.GetEndAtom()
                if start.GetAtomMapNum() and end.GetAtomMapNum():
                    bonds.add(molecule.GetBondBetweenAtoms(
                        match[start.GetIdx()], match[end.GetIdx()],
                    ).GetIdx())
            for number, label in labels.items():
                if label is not None:
                    atom = molecule.GetAtomWithIdx(mapped[number])
                    atom.SetProp('atomNote', label)
                    atom.SetProp('_linker_role', roles[number])
    return sorted(atoms), sorted(bonds)


def colour_position_labels(svg, molecule):
    """Colour SVG position notes by unit role and replace tokens with Greek.

    RDKit emits each note character separately, in atom order. Check that
    sequence before recolouring, so changes to its SVG format fail visibly
    instead of assigning a label to the wrong unit. Element symbols stay black.
    """
    notes = iter(
        (character, atom.GetProp('_linker_role'))
        for atom in molecule.GetAtoms() if atom.HasProp('atomNote')
        for character in atom.GetProp('atomNote')
    )
    greek = {'a': '&#945;', 'b': '&#946;', 'g': '&#947;', "'": '&#8242;'}

    def replace_note(match):
        """Replace one checked SVG note glyph with its coloured Greek glyph."""
        expected = next(notes, None)
        if expected is None or unescape(match[2]) != expected[0]:
            raise ValueError('RDKit SVG notes do not match the atom label sequence.')
        character, role = expected
        header = match[1].replace('fill:#000000', f'fill:{LABEL_COLOURS[role]}')
        header = header.replace("class='note'", f"class='note' data-role='{role}'")
        return header + greek.get(character, character) + match[3]

    svg = re.sub(r"(<text\b[^>]*class='note'[^>]*>)(.*?)(</text>)", replace_note, svg)
    if next(notes, None) is not None:
        raise ValueError('RDKit omitted position labels from the SVG.')
    return svg


def read_xyz(path):
    """Return a sanitized RDKit molecule inferred from a neutral XYZ file."""
    molecule = Chem.MolFromXYZBlock(path.read_text(encoding='utf-8'))
    if molecule is None or molecule.GetNumAtoms() == 0:
        raise ValueError(f'Cannot read a molecule from {path}.')

    # Recover single/double bonds before removing explicit hydrogens. Their
    # positions and valences are needed to distinguish OH, ether, and alkene.
    rdDetermineBonds.DetermineBonds(molecule, charge=0, embedChiral=False)
    Chem.SanitizeMol(molecule)
    if len(Chem.GetMolFrags(molecule)) != 1:
        raise ValueError('Expected one connected linker structure.')
    if any(atom.GetNumRadicalElectrons() for atom in molecule.GetAtoms()):
        raise ValueError('Expected a closed-shell structure with all hydrogens present.')

    # A relaxed model conformation is not an experimental stereo assignment.
    Chem.RemoveStereochemistry(molecule)
    return Chem.RemoveHs(molecule)


def add_information_panel(svg, motif):
    """Return SVG with a three-column eligibility panel below the structure.

    Every figure has the same canvas and panel position for stable gallery
    navigation. Structures retain their scale unless they exceed the shared
    drawing area; those are reduced uniformly to fit without clipping.
    """
    namespace = 'http://www.w3.org/2000/svg'
    ET.register_namespace('', namespace)
    root = ET.fromstring(svg)
    _, _, old_width, old_height = map(float, root.attrib['viewBox'].split())
    width, structure_height = 920, 400
    scale = min(1, (width - 24) / old_width, structure_height / old_height)
    top = structure_height + 12
    height = top + 278
    root.set('width', f'{width:g}px')
    root.set('height', f'{height:g}px')
    root.set('viewBox', f'0 0 {width:g} {height:g}')

    # Center each molecule inside the same area while keeping bond lengths and
    # label sizes consistent across the four supplied structures.
    drawing = ET.Element(f'{{{namespace}}}g', {
        'transform': (
            f'translate({(width - old_width * scale) / 2:g},'
            f'{(structure_height - old_height * scale) / 2:g}) scale({scale:g})'
        ),
    })
    drawing.extend(list(root))
    root[:] = []
    ET.SubElement(root, f'{{{namespace}}}rect', {
        'width': str(width), 'height': str(height), 'fill': 'white',
    })
    root.append(drawing)
    panel = ET.SubElement(root, f'{{{namespace}}}g', {
        'id': 'linker-information', 'font-family': FONT_FAMILY,
        'font-size': '20', 'fill': '#222222',
    })

    def text(x, y, value, colour='#222222', bold=False, centered=True, font_size=20):
        """Add one line of native SVG text at the given pixel position."""
        node = ET.SubElement(panel, f'{{{namespace}}}text', {
            'x': str(x), 'y': str(y), 'fill': colour,
            'font-weight': 'bold' if bold else 'normal',
            'text-anchor': 'middle' if centered else 'start',
            'font-size': str(font_size),
        })
        node.text = value

    def items(x, y, values, wrap_width):
        """Centre wrapped points, adding one bullet per distinct requirement.

        Content entries are complete statements, not individual display lines:
        wrapping one statement must never create additional bullet points.
        """
        multiple = len(values) > 1
        for value in values:
            colour = LABEL_COLOURS['partner'] if 'partner' in value.lower() else '#222222'
            for index, chunk in enumerate(wrap(value, width=wrap_width,
                                               break_long_words=False, break_on_hyphens=False)):
                prefix = '• ' if multiple and index == 0 else ''
                text(x, y, prefix + chunk, colour)
                y += 24
            y += 6

    def line(x1, y1, x2, y2, stroke_width=1.2):
        """Draw a gently uneven pencil-like rule, reproducible on every run."""
        dx, dy = x2 - x1, y2 - y1
        length = max(abs(dx), abs(dy))
        nx, ny = -dy / length, dx / length
        # A few small opposing bends resemble a pen stroke while keeping the
        # endpoints aligned. Fixed offsets avoid jitter between regenerations.
        path = f'M {x1},{y1} '
        for index in range(3):
            start, end = index / 3, (index + 1) / 3
            sign = 1 if index % 2 == 0 else -1
            first, second = start + 0.1, end - 0.1
            path += (f'C {x1 + dx * first + nx * 4 * sign},'
                     f'{y1 + dy * first + ny * 4 * sign} '
                     f'{x1 + dx * second - nx * 4 * sign},'
                     f'{y1 + dy * second - ny * 4 * sign} '
                     f'{x1 + dx * end},{y1 + dy * end} ')
        ET.SubElement(panel, f'{{{namespace}}}path', {
            'd': path, 'fill': 'none', 'stroke': '#A09B93',
            'stroke-width': str(stroke_width), 'stroke-linecap': 'round',
        })

    content = PANEL_CONTENT[motif]
    # Align H/G/S with both tick columns without rules between those rows.
    # Repeated role names below tie free-site requirements to their units.
    monomer_edge, chain_edge, special_edge = 112, width * 0.38, width * 0.67
    incoming_center = (monomer_edge + chain_edge) / 2
    chain_center = (chain_edge + special_edge) / 2
    special_center = (special_edge + width - 12) / 2
    line(12, top, width - 12, top, stroke_width=2.4)
    line(monomer_edge, top, monomer_edge, height - 12)
    line(chain_edge, top, chain_edge, height - 12)
    line(special_edge, top, special_edge, height - 12, stroke_width=2.4)
    text(incoming_center, top + 30, 'Incoming monomer', LABEL_COLOURS['incoming'], True)
    text(chain_center, top + 30, 'Chain monomer', LABEL_COLOURS['chain'], True)
    text(special_center, top + 30, 'Special conditions', bold=True)
    line(12, top + 46, special_edge, top + 46)
    for index, unit in enumerate('HGS'):
        y = top + 72 + index * 32
        text(62, y, unit, bold=True, centered=True)
        text(incoming_center, y + 2, '✓', centered=True, font_size=25)
        text(chain_center, y + 2, '✓' if unit in content['chain_types'] else '✗',
             centered=True, font_size=25)
    line(12, top + 150, special_edge, top + 150)
    text(incoming_center, top + 182, 'Incoming free sites', LABEL_COLOURS['incoming'], True)
    text(chain_center, top + 182, 'Chain free sites', LABEL_COLOURS['chain'], True)
    line(monomer_edge, top + 194, special_edge, top + 194)
    for x, field in ((incoming_center, 'incoming_sites'), (chain_center, 'chain_sites')):
        items(x, top + 218, content[field], wrap_width=22)
    items(special_center, top + 65, content['notes'], wrap_width=26)
    return ET.tostring(root, encoding='unicode')


def draw_linker(molecule):
    """Return an annotated SVG and eligibility panel without editing the input."""
    molecule = Chem.Mol(molecule)
    rdDepictor.Compute2DCoords(molecule)
    atoms, bonds = annotate_linkers(molecule)

    # Let RDKit fit the canvas to the molecule, including OH labels, instead
    # of leaving empty space around a fixed 900 x 600 rectangle.
    # Native SVG text supports Greek glyphs in the browser's sans-serif font.
    drawer = rdMolDraw2D.MolDraw2DSVG(-1, -1, noFreetype=True)
    options = drawer.drawOptions()
    options.useBWAtomPalette()
    options.bondLineWidth = 2
    options.scalingFactor = 40
    options.padding = 0.03
    options.annotationFontScale = 0.95
    yellow = (1.0, 0.9, 0.3)
    rdMolDraw2D.PrepareAndDrawMolecule(
        drawer, molecule, highlightAtoms=atoms, highlightBonds=bonds,
        highlightAtomColors={index: yellow for index in atoms},
        highlightBondColors={index: yellow for index in bonds},
    )
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText().replace('font-family:sans-serif', f'font-family:{FONT_FAMILY}')
    svg = colour_position_labels(svg, molecule)
    if molecule.HasProp('_drawing_motif'):
        svg = add_information_panel(svg, molecule.GetProp('_drawing_motif'))
    return svg


def main():
    """Read an XYZ file and write its SVG to pics/, or to --output."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        'xyz', nargs='?', type=Path,
        default=Path(__file__).resolve().with_name('beta-O-4.xyz'),
        help='Input XYZ file (default: beta-O-4.xyz beside this script).',
    )
    parser.add_argument('-o', '--output', type=Path, help='Output SVG filename.')
    args = parser.parse_args()
    output = args.output or args.xyz.parent / 'pics' / args.xyz.with_suffix('.svg').name
    if output.suffix.lower() != '.svg':
        parser.error('The output filename must end in .svg.')
    if output.resolve() == args.xyz.resolve():
        parser.error('The output must differ from the input file.')
    try:
        svg = draw_linker(read_xyz(args.xyz))
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(svg, encoding='utf-8')
    except (OSError, ValueError, RuntimeError) as exc:
        parser.exit(1, f'Cannot draw {args.xyz}: {exc}\n')
    print(f'Wrote {output}')


if __name__ == '__main__':
    main()
