"""
Central registry and factories for structure analyzers.
-> Factory: Function that creates and returns an analyzer instance or None.

Each entry in `ANALYZERS_REGISTRATION` maps a short key to a factory
callable that accepts `(args, **context)` and returns an analyzer
instance or `None` if it should not be created for the given args.
"""
from MakroLyzer.structure_modules.Hbonds import HBondsAnalyzer
from MakroLyzer.structure_modules.Anisotropy import AnisotropyAnalyzer
from MakroLyzer.structure_modules.Asphericity import AsphericityAnalyzer
from MakroLyzer.structure_modules.RadiusOfGyration import RadiusOfGyrationAnalyzer
from MakroLyzer.structure_modules.MoleculeCount import MoleculeCountAnalyzer
from MakroLyzer.structure_modules.EndToEndDistance import EndToEndDistanceAnalyzer
from MakroLyzer.structure_modules.OrderParameter import OrderParameterAnalyzer
from MakroLyzer.structure_modules.Ramachandran import RamachandranAnalyzer
from MakroLyzer.structure_modules.Dihedrals import DihedralsAnalyzer
from MakroLyzer.structure_modules.ChemicalFormula import ChemicalFormulaAnalyzer
from MakroLyzer.structure_modules.structureBase import OutputHandler


def create_hbonds(args, **context):
    if not args.get('hydrogenBonds'):
        return None
    cutoffs = args['hydrogenBonds']
    hbonds_file = args.get('hbonds_file') or 'hydrogenBonds.csv'
    output_handler = OutputHandler(hbonds_file, mode='streaming')
    return HBondsAnalyzer(cutoffs, output_handler)


def create_anisotropy(args, **context):
    val = args.get('anisotropyFactor')
    if val is None:
        return None
    out_file = args.get('anisotropy_file') if args.get('anisotropy_file') is not None else (val if isinstance(val, str) else 'anisotropyFactor.csv')
    output_handler = OutputHandler(out_file, mode='streaming')
    return AnisotropyAnalyzer(output_handler)


def create_asphericity(args, **context):
    val = args.get('asphericityParameter')
    if val is None:
        return None
    out_file = args.get('asphericity_file') if args.get('asphericity_file') is not None else (val if isinstance(val, str) else 'asphericityParameter.csv')
    output_handler = OutputHandler(out_file, mode='streaming')
    return AsphericityAnalyzer(output_handler)


def create_radius_of_gyration(args, **context):
    val = args.get('radiusOfGyration')
    if val is None:
        return None
    out_file = args.get('Rg_file') if args.get('Rg_file') is not None else (val if isinstance(val, str) else 'Rg.csv')
    output_handler = OutputHandler(out_file, mode='streaming')
    return RadiusOfGyrationAnalyzer(output_handler)


def create_molecule_count(args, **context):
    val = args.get('NoSubgraphs')
    if val is None:
        return None
    out_file = args.get('NoSub_file') if args.get('NoSub_file') is not None else (val if isinstance(val, str) else 'NoSubGraphs.csv')
    output_handler = OutputHandler(out_file, mode='streaming')
    return MoleculeCountAnalyzer(output_handler)


def create_end_to_end(args, **context):
    val = args.get('endToEndDistance')
    if val is None:
        return None
    out_file = args.get('e2e_file') if args.get('e2e_file') is not None else (val if isinstance(val, str) else 'e2eDistances.csv')
    output_handler = OutputHandler(out_file, mode='streaming')
    static_topology = context.get('static_topology', False)
    backbone_cache = context.get('backbone_cache')
    return EndToEndDistanceAnalyzer(output_handler, static_topology=static_topology, backbone_cache=backbone_cache)


def create_order_parameter(args, **context):
    if not args.get('orderParameter'):
        return None
    # Expect that main passes BoxSize, NoCellsPerDim, MolecularVectorLength in context
    BoxSize = context.get('BoxSize')
    NoCellsPerDim = context.get('NoCellsPerDim')
    MolecularVectorLength = context.get('MolecularVectorLength')
    static_topology = context.get('static_topology', False)
    backbone_cache = context.get('backbone_cache')
    order_file = args.get('order_file') or 'orderParameter.csv'
    output_handler = OutputHandler(order_file, mode='streaming')
    return OrderParameterAnalyzer(BoxSize, NoCellsPerDim, MolecularVectorLength, output_handler, static_topology=static_topology, backbone_cache=backbone_cache)


def create_ramachandran(args, **context):
    val = args.get('Ramachandran')
    if val is None:
        return None
    out_file = args.get('Ramachandran_file') if args.get('Ramachandran_file') is not None else (val if isinstance(val, str) else 'Ramachandran.csv')
    output_handler = OutputHandler(out_file, mode='streaming')
    return RamachandranAnalyzer(output_handler)


def create_dihedrals(args, **context):
    dihedral_val = args.get('dihedral')
    if dihedral_val is None:
        return None
    dihedral_handler = None
    dihedral_list_handler = None
    cistrans_handler = None
    if dihedral_val is not None:
        dihedral_file = 'dihedrals.csv'
        dihedral_handler = OutputHandler(dihedral_file, mode='collect')
        dihedral_list_handler = OutputHandler(dihedral_file.replace('.csv', '_list.csv'), mode='streaming')
        cistrans_file = 'CisTrans.csv'
        cistrans_handler = OutputHandler(cistrans_file, mode='streaming')
    static_topology = context.get('static_topology', False)
    backbone_cache = context.get('backbone_cache')
    return DihedralsAnalyzer(
        dihedral_output_handler=dihedral_handler,
        dihedral_list_output_handler=dihedral_list_handler,
        cistrans_output_handler=cistrans_handler,
        dihedral_range=args.get('dihedral_range'),
        special_dihedral=args.get('special_dihedral'),
        static_topology=static_topology,
        backbone_cache=backbone_cache
    )


def create_chemical_formula(args, **context):
    val = args.get('formula')
    if val is None:
        return None
    out_file = args.get('formula_file') if args.get('formula_file') is not None else (val if isinstance(val, str) else 'chemicalFormulas.csv')
    formula_handler = OutputHandler(out_file, mode='streaming')
    return ChemicalFormulaAnalyzer(formula_handler)


ANALYZERS_REGISTRATION = {
    'hbonds': create_hbonds,
    'anisotropy': create_anisotropy,
    'asphericity': create_asphericity,
    'radiusOfGyration': create_radius_of_gyration,
    'noSubgraphs': create_molecule_count,
    'endToEndDistance': create_end_to_end,
    'orderParameter': create_order_parameter,
    'AminoAcidRamachandran': create_ramachandran,
    'dihedrals': create_dihedrals,
    'chemicalFormula': create_chemical_formula,
}
