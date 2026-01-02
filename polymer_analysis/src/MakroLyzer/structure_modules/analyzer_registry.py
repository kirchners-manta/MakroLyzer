"""
Central registry and factories for structure analyzers.
-> Factory: Function that creates and returns an analyzer instance or None.

Each entry in `ANALYZERS_REGISTRATION` maps a short key to a factory
callable that accepts `(args, **context)` and returns an analyzer
instance or `None` if it should not be created for the given args.
"""
from src.MakroLyzer.structure_modules.Hbonds import HBondsAnalyzer
from src.MakroLyzer.structure_modules.Anisotropy import AnisotropyAnalyzer
from src.MakroLyzer.structure_modules.Asphericity import AsphericityAnalyzer
from src.MakroLyzer.structure_modules.RadiusOfGyration import RadiusOfGyrationAnalyzer
from src.MakroLyzer.structure_modules.MoleculeCount import MoleculeCountAnalyzer
from src.MakroLyzer.structure_modules.EndToEndDistance import EndToEndDistanceAnalyzer
from src.MakroLyzer.structure_modules.OrderParameter import OrderParameterAnalyzer
from src.MakroLyzer.structure_modules.Ramachandran import RamachandranAnalyzer
from src.MakroLyzer.structure_modules.Dihedrals import DihedralsAnalyzer
from src.MakroLyzer.structure_modules.ChemicalFormula import ChemicalFormulaAnalyzer
from src.MakroLyzer.structure_modules.structureBase import OutputHandler


def create_hbonds(args, **context):
    if not args.get('hydrogenBonds'):
        return None
    cutoffs = args['hydrogenBonds']
    output_handler = OutputHandler(args['hbonds_file'], mode='streaming')
    return HBondsAnalyzer(cutoffs, output_handler)


def create_anisotropy(args, **context):
    if not args.get('anisotropyFactor'):
        return None
    output_handler = OutputHandler(args['anisotropy_file'], mode='streaming')
    return AnisotropyAnalyzer(output_handler)


def create_asphericity(args, **context):
    if not args.get('asphericityParameter'):
        return None
    output_handler = OutputHandler(args['asphericity_file'], mode='streaming')
    return AsphericityAnalyzer(output_handler)


def create_radius_of_gyration(args, **context):
    if not args.get('radiusOfGyration'):
        return None
    output_handler = OutputHandler(args['Rg_file'], mode='streaming')
    return RadiusOfGyrationAnalyzer(output_handler)


def create_molecule_count(args, **context):
    if not args.get('noSubgraphs'):
        return None
    output_handler = OutputHandler(args['noSub_file'], mode='streaming')
    return MoleculeCountAnalyzer(output_handler)


def create_end_to_end(args, **context):
    if not args.get('endToEndDistance'):
        return None
    output_handler = OutputHandler(args['e2e_file'], mode='streaming')
    return EndToEndDistanceAnalyzer(output_handler)


def create_order_parameter(args, **context):
    if not args.get('orderParameter'):
        return None
    # Expect that main passes BoxSize, NoCellsPerDim, MolecularVectorLength in context
    BoxSize = context.get('BoxSize')
    NoCellsPerDim = context.get('NoCellsPerDim')
    MolecularVectorLength = context.get('MolecularVectorLength')
    output_handler = OutputHandler(args['order_file'], mode='streaming')
    return OrderParameterAnalyzer(BoxSize, NoCellsPerDim, MolecularVectorLength, output_handler)


def create_ramachandran(args, **context):
    if not args.get('AminoAcidRamachandran'):
        return None
    output_handler = OutputHandler(args['AARamachandran_file'], mode='streaming')
    return RamachandranAnalyzer(output_handler)


def create_dihedrals(args, **context):
    if not (args.get('dihedral') or args.get('cisTrans')):
        return None
    dihedral_handler = OutputHandler(args['dihedral_file'], mode='collect') if args.get('dihedral') else None
    dihedral_list_handler = OutputHandler(args['dihedral_file'].replace('.csv', '_list.csv'), mode='streaming') if args.get('dihedral') else None
    cistrans_handler = OutputHandler(args['CisTrans_file'], mode='streaming') if args.get('dihedral') else None
    return DihedralsAnalyzer(
        dihedral_output_handler=dihedral_handler,
        dihedral_list_output_handler=dihedral_list_handler,
        cistrans_output_handler=cistrans_handler,
        dihedral_range=args.get('dihedral_range')
    )


def create_chemical_formula(args, **context):
    if not args.get('formula'):
        return None
    formula_handler = OutputHandler(args['formula_file'], mode='streaming')
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
