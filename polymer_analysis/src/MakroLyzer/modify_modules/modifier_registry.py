"""
Central registry and factories for structure modifiers.
-> Factory: Function that creates and returns an modifier instance or None.

Each entry in `MODIFIERS_REGISTRATION` maps a short key to a factory
callable that accepts `(args, **context)` and returns an modifier
instance or `None` if it should not be created for the given args.
"""
from MakroLyzer.modify_modules.Patterns import PatternModifier
from MakroLyzer.modify_modules.structureModifierBase import ModifyOutputHandler
from MakroLyzer.modify_modules.AminoSaturation import AminoSaturationModifier
from MakroLyzer.modify_modules.SubgraphPrint import SubgraphPrintModifier
from MakroLyzer.modify_modules.WrapGraph import WrapGraphModifier
from MakroLyzer.modify_modules.FunctionalizePE import FunctionalizePEModifier
from MakroLyzer.modify_modules.FunctionalizePEsurface import FunctionalizePEsurfaceModifier


def create_patterns(args, **context):
    pattern_path = args.get('patternFile')
    if not pattern_path:
        return None
    out_file = args.get('patternID_file') if args.get('patternID_file') is not None else 'patternIDs.csv'
    output_handler = ModifyOutputHandler(out_file)
    return PatternModifier(pattern_path, output_handler)

def create_AminoSaturation(args, **context):
    val = args.get('saturation')
    if val is None:
        return None
    out_file = val if isinstance(val, str) else 'saturatedPolymers.xyz'
    output_handler = ModifyOutputHandler(out_file)
    return AminoSaturationModifier(output_handler)

def create_subgraph_print(args, **context):
    val = args.get('subgraph_coords')
    if val is None:
        return None
    out_file = val if isinstance(val, str) else 'subgraphCoordinates.xyz'
    output_handler = ModifyOutputHandler(out_file)
    return SubgraphPrintModifier(output_handler)

def create_wrap_graph_print(args, **context):
    val = args.get('wrapGraphPrint')
    if val is None:
        return None
    out_file = val if isinstance(val, str) else 'wrappedGraph.xyz'
    output_handler = ModifyOutputHandler(out_file)
    return WrapGraphModifier(output_handler)

def create_functionalize_PE(args, **context):
    val = args.get('functionalizePE')
    if val is None:
        return None
    mode = val['mode']
    func_type = val['func_type']
    percentage = val.get('percentage')
    neighbor_exclusion = val.get('neighbor_exclusion')
    distance = val.get('distance')
    out_file = args.get('functionalizePE_file') if args.get('functionalizePE_file') is not None else 'functionalizedPE.xyz'
    output_handler = ModifyOutputHandler(out_file)
    return FunctionalizePEModifier(mode, func_type, percentage, neighbor_exclusion, distance, output_handler)

def create_functionalize_PEsurf(args, **context):
    val = args.get('functionalizePEsurface')
    if val is None:
        return None
    percentage = val['percentage']
    func_type = val['func_type']
    neighbor_exclusion = val.get('neighbor_exclusion')
    out_file = args.get('functionalizePE_file') if args.get('functionalizePE_file') is not None else 'functionalizedPE.xyz'
    output_handler = ModifyOutputHandler(out_file)
    box_size = context.get('boxSize')
    return FunctionalizePEsurfaceModifier(
        percentage,
        func_type,
        output_handler,
        box_size=box_size,
        neighbor_exclusion=neighbor_exclusion,
    )


MODIFIERS_REGISTRATION = {
    'functionalizePE': create_functionalize_PE,
    'functionalizePEsurface': create_functionalize_PEsurf,    
    'patternFile': create_patterns,
    'AminoSaturation': create_AminoSaturation,
    'subgraphPrint': create_subgraph_print,
    'wrapGraphPrint': create_wrap_graph_print,
}
