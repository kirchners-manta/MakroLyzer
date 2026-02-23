import argparse
import re
import sys

ASCII_ART = (
    "\x1b[1m"
    r"""
                             
                                @    @@ @            
                              @@ @@ @    @           
                              @   @@    @           
                             @     @      @@         
                             @         @@    @        
                            @@         @@     @      
                            @                  @     
                            @                 @      
                            @                @@      
                            @              @@        
                           @               @         
                         @                 @         
                       @                   @         
                     @                      @    +---------------------------|-----------------------------+    
                    @                       @    |   __  __       _          |   _                         |  
                   @                        @    |  |  \/  | __ _| | ___ __ _|_ | |   _   _ _______ _ __   |   
                  @                         @    |  | |\/| |/ _` | |/ / '__/ _ \| |  | | | |_  / _ \ '__|  |
                 @           @              @    |  | |  | | (_| |   <| | | (_) | |__| |_| |/ /  __/ |     |  
     @@ @        @           @              @    |  |_|  |_|\__,_|_|\_\_|  \___/|_____\__, /___\___|_|     |
   @@    @       @           @       @      @    |                          /|        |___/                |   
  @@    @@      @            @       @      @    |                                                         |
  @    @@       @            @       @     @     |    MakroLyzer - Makromolecule Analyzer                  |   
  @    @@       @             @      @     @     +---------------------------------------------------------+   
  @     @@      @             @      @    @         
   @        @@                 @     @    @          
    @@                          @    @    @                                         
      @@@                       @    @    @          
          @ @ @ @ @ @ @ @ @     @ @  @  @ @      
    """
    "\x1b[0m"
    r"""                                                                
--------------------------------------------------------------------------------------------------
   Group of Prof. Dr. Barbara Kirchner
   Maintainer: Katrin Drysch
   
   Please cite:
   - J. Phys. Chem. B 2025, 129, 50, 12997-13008 (DOI: 10.1021/acs.jpcb.5c06175)
   
   For the polyethylene surface functionalization module as well as the surface atoms
   analysis module, please additionally cite:
   - J. Comput. Chem. 2018, 39, 2118-2125 (DOI: 10.1002/jcc.25384)
   - J. Comput. Chem. 2011, 32, 2319-2327 (DOI: 10.1002/jcc.21787)
--------------------------------------------------------------------------------------------------
"""
)

class ArtParser(argparse.ArgumentParser):
    def format_help(self):
        help_text = super().format_help()
        help_text = re.sub(r'\n(?= {2}-)', '\n\n', help_text)
        return ASCII_ART + "\n" + help_text

def format_header(title, width=90):
    core = f" {title} "
    if len(core) >= width:
        return f"\x1b[1m{core}\x1b[0m"
    dash_total = width - len(core)
    left = dash_total // 2
    right = dash_total - left
    return f"\x1b[1m{'-' * left}{core}{'-' * right}\x1b[0m"

def readCommandLine() -> dict:
    
    """
    Read input parameters from command line arguments.
    
    Returns:
        dict: Dictionary containing the input parameters.
    """
    
    # Create the argument parser
    parser = ArtParser(prog='MakroLyzer',
                       formatter_class=argparse.RawTextHelpFormatter,
                       description='MakroLyzer - Makromolecule Analysis Tool')
    
    # Input Group # ----------------------------------------------------------------------------------------
    input_group = parser.add_argument_group(
        format_header("Arguments for specifying input files and settings")
    )
    file_group = input_group.add_mutually_exclusive_group(required=True)
    
    file_group.add_argument('-xyz', '--xyzFile', 
                        help='Path to the XYZ-trajectory file', 
                        default=None) 
    file_group.add_argument('-lmp', '--lmpFile', 
                        help='Path to the LAMMPS-trajectory file',
                        default=None)
    
    input_group.add_argument('-nth', '--nthStep',
                        help='Analyze/Modify every nth step from the trajectory (default: 1)',
                        type=int,
                        default=1)
    input_group.add_argument('-bs', '--BoxSize', 
                        help='Box size for periodic boundary conditions. (default: None)',
                        type=float)
    input_group.add_argument('-vf', '--vibFactor',
                        dest='vibFactor',
                        help='Vibration factor to scale covalent bond cutoff. (default: 1.15)',
                        type=float,
                        default=None)
    input_group.add_argument(
                        '--static-topology',
                        dest='staticTopology',
                        action='store_true',
                        help='Build graph topology only on the first frame and only update coordinates afterwards.\n'
                              'NOTE: This is not recommended for usage with modules that modify the structure of the macromolecule.')
    input_group.add_argument(
                        '-sel', '--subgraph-selection',
                        dest='subgraphSelection',
                        nargs='+',
                        default=None,
                        help='Select subgraphs by exact formula or element set.\n'
                        'Exact composition: include digits (e.g. C2H6, CH4, C1H4).\n'
                        'Element set only: no digits (e.g. CHO) -> any counts, only these elements.\n'
                        'Provide multiple selections separated by space or comma.\n'
                        'Example: -sel C2H6 C4H10 C1H4 CHO')
   
    
    # Structure Analyzers Group # -------------------------------------------------------------------------------
    analyzer_group = parser.add_argument_group(title=format_header("Arguments for structure analyzers"),
                                                description='This module includes various analyzers to compute structural properties of macromolecules.')
    
    analyzer_group.add_argument(
                        '-d', '--dihedral',
                        nargs='?', const='dihedrals.csv', default=None,
                        help='Calculate dihedral angles.\n'
                        '1. Dihedral angle counts (file: dihedrals.csv)\n'
                        '2. Dihedral angles per along each polymer strand (file: dihedrals_list.csv)\n'
                        '3. Cis-trans counts (file: cisTrans.csv)\n'
                        'Optionally provide a base-output filename (default: dihedrals.csv)')
    
    analyzer_group.add_argument(
                        '--dihedral-range',
                        help='Range of dihedral angles.\n'
                        '"abs": (0-180), "nonabs": (-180-180), (default: 0-180)',
                        choices=['abs', 'nonabs'],
                        default='abs')
    
    analyzer_group.add_argument(
                        '--special-dihedral',
                        dest='special_dihedral',
                        nargs='+',
                        default=None,
                        help='Take only certain dihedrals into account.\n'
                        'Provide the atom types of the four consecutive elements along the backbone you\n'
                        'want to calculate the dihedral angles for.\n'
                        'Usage example: -d --special-dihedral OCCO')
    
    analyzer_group.add_argument(
                        '-Rg', '--radiusOfGyration',
                        nargs='?', const='Rg.csv', default=None,
                        help='Calculate radius of gyration for the whole polymer.\n'
                             'Optionally provide an output filename (default: Rg.csv)')
    
    analyzer_group.add_argument(
                        '-af', '--anisotropyFactor',
                        nargs='?', const='anisotropyFactor.csv', default=None,
                        help='Calculate anisotropy factor.\n'
                        'Optionally provide an output filename (default: anisotropyFactor.csv)')
    
    analyzer_group.add_argument(
                        '-as', '--asphericityParameter',
                        nargs='?', const='asphericityParameter.csv', default=None,
                        help='Calculate asphericity parameter.\n'
                        'Optionally provide an output filename (default: asphericityParameter.csv)')
    
    analyzer_group.add_argument(
                    '-hb', '--hydrogenBonds',
                    nargs='+',
                    help='Calculate the number of hydrogen bonds in the given polymer with specific cutoffs.\n'
                    'Data to be specified:   Acceptor Element Type : b : c : alpha\n'
                    '                        b : maximum hydrogen atom (H) - acceptor atom (A) distance (float)\n'
                    '                        c : maximum donor atom (D) - acceptor atom distance (float)\n'
                    '                        alpha : maximum angle between the D-H and H-A vectors (float)\n'
                    'Usage example: -hb O:2.4:3.4:30 N:2.8:3.9:25',
                    type=element_distance_tuple)
    
    analyzer_group.add_argument(
                        '-hb-file', '--hbonds-file',
                        nargs='?', const='hydrogenBonds.csv', default='hydrogenBonds.csv',
                        help='Output filename for hydrogen bonds (default: hydrogenBonds.csv)')
    
    analyzer_group.add_argument(
                        '-op', '--orderParameter',
                        help='Calculate the nematic order parameter S* for the polymer structure\n'
                        'The order parameter S* can either be calculated for the entire graph, or the space can be divided into\n'
                        'smaller cells, and S* can be calculated for each cell individually.\n'
                        'The overall order parameter is then the average of the order parameters of all cells.\n'
                        'Parameters to be specified:  <BoxSize>:<NoCellsPerDim>[:<MolecularVectorLength>]\n'
                        '                    BoxSize : Size of the simulation box in x, y, z directions (float,float,float) or one value if symmetrical\n'
                        '                    NoCellsPerDim : Number of cells in each direction (int, int, int) or only one value if symmetrical\n'
                        '                    MolecularVectorLength : Number of atoms used to calculate a molecular vector (int)\n'
                        '                                            A value of 2 means that the vector is calculated between two consecutive atoms / ring ceters along the backbone.\n'
                        '                                            Required for --op-vector-source atom/ring-center.\n'
                        '                                            Ignored for --op-vector-source ring-normal.\n'
                        'Usage examples: -op 100:1:3   or   -op 100:1',
                        type=OrderParam)
    
    analyzer_group.add_argument(
                        '-op-file', '--order-file',
                        nargs='?', const='orderParameter.csv', default='orderParameter.csv',
                        help='Output filename for the nematic order parameter S* (default: orderParameter.csv)')
    
    analyzer_group.add_argument(
                        '--op-vector-source',
                        choices=['atom', 'ring-center', 'ring', 'ring-normal'],
                        default='atom',
                        help='Vector source for order parameter calculation: atom (default), ring-center, or ring-normal. "ring" is an alias for ring-center.\n'
                        'If "atom" is chosen, molecular vectors are assigned via consecutive atoms along the polymer backbone.\n'
                        'If "ring-center" is chosen, molecular vectors are assigned via the consecutive centers of rings along the polymer backbone.\n'
                        'If "ring-normal" is chosen, molecular vectors are assigned via the normal vectors of all rings with matching size.')
    
    analyzer_group.add_argument(
                        '--op-ring-size',
                        type=RingCycleSize,
                        default=None,
                        help='Optional ring size filter for ring-based order-parameter vectors.\n'
                        'Accepts a single size (e.g. 6) or a range [min,max] (e.g. [4,7]).\n'
                        'Minimum allowed ring size is 3.')
    
    analyzer_group.add_argument(
                        '-e2e', '--endToEndDistance',
                        nargs='?', const='endToEndDistances.csv', default=None,
                        help='Calculate end-to-end distance.\n'
                        'Optionally provide an output filename (default: e2eDistances.csv)')
    
    analyzer_group.add_argument(
                        '-Ramachandran', '--Ramachandran',
                        nargs='?', const='AARamachandran.csv', default=None,
                        help='Calculate Ramachandran matrix.\n'
                        'Optionally provide an output filename (default: Ramachandran.csv)')

    analyzer_group.add_argument(
                        '-NoSub', '--NoSubgraphs',
                        nargs='?', const='noSubGraphs.csv', default=None,
                        help=('Calculate the number of polymer strands, as well as the number of chains and rings.\n'
                              'Optionally provide an output filename (default: noSubGraphs.csv)'))
    
    analyzer_group.add_argument(
                        '-f', '--formula',
                        nargs='?', const='chemicalFormulas.csv', default=None,
                        help='Get chemical formulas of the polymer.\n'
                        'Optionally provide an output filename (default: chemicalFormulas.csv)')
    
    analyzer_group.add_argument(
                        '-sA', '--surfaceAtoms',
                        nargs='?', const='surfaceAtoms.csv', default=None,
                        help='Identify surface atoms of the polymer and classify them as nonpolar, polar, and ring atoms.\n'
                        'Ring atoms are carbon atoms that are part of rings matching --sA-ring-size.\n'
                        'Output columns: Frame, Nonpolar Surface Atoms, Polar Surface Atoms, Ring Surface Atoms.\n'
                        'Optionally provide an output filename (default: surfaceAtoms.csv).\n'
                        'Usage example: -sA surfaceAtoms.csv'
    )
    
    analyzer_group.add_argument(
                        '--sA-ring-size',
                        dest='sA_ring_size',
                        type=RingCycleSize,
                        default=6,
                        help='Ring size filter for ring atom assignment in surface-atom analysis.\n'
                        'Accepts a single size (e.g. 6) or a range [min,max] (e.g. [4,7]).\n'
                        'Minimum allowed ring size is 3. (default: 6).\n'
                        'Usage examples: --sA-ring-size 6   or   --sA-ring-size [4,7]')
    
    # Structure Modifiers Group # ------------------------------------------------------------------------------
    modifier_group = parser.add_argument_group(
        format_header("Arguments for structure modifiers")
    )
    
    modifier_group.add_argument(
                        '-funcPE', '--functionalizePE',
                        default=None,
                        help='Functionalize Polyethylene (PE) with CO/COH/ester/amide groups.\n'
                        'Modes:\n'
                        '  random:<percentage>:<func_type>:<neighbor_exclusion>\n'
                        '  periodic:<distance>:<func_type>\n'
                        'The random mode ramdomly distibutes the functional groups <func_type> over the polymerstrands,\n'
                        'with a minimum distance of <neighbor_exclusion>. The <percentage> value gives the percentage of\n'
                        'adjusted carbon atoms in the backbone, excluding end group carbons (CH3).\n'
                        'The periodic mode adds functional groups of type <func_type> in distances of <distance> along the\n'
                        'backbones of the PE strands.\n'
                        'Usage examples:\n'
                        '  -funcPE random:15.8:CO:2\n'
                        '  -funcPE periodic:6:ester',
                        type=functionalize_pe)
    
    modifier_group.add_argument(
                        '-funcPEsurf', '--functionalizePEsurface',
                        default=None,
                        help='Functionalize Polyethylene (PE) surface with CO/COH/ester/amide groups (random only).\n'
                        'Modes:\n'
                            'random:<percentage>:<func_type>:<neighbor_exclusion>\n'
                        'Usage example:\n'
                            '-funcPEsurf random:15:CO:2',
                        type=functionalize_pe_surface)
    
    modifier_group.add_argument(
                        '-funcPE-file', '--functionalizePE-file',
                        nargs='?', const='functionalizedPE.xyz', default=None,
                        help='Output filename for functionalized polyethylene (default: functionalizedPE.xyz)')
    
    modifier_group.add_argument(
                        '-p', '--patternFile',
                        dest='patternFile',
                        help='Assigns specific pattern IDs to fragments based on predefined patterns.\n'
                        'Required argument: Path to the TXT file containing the patterns in list of lists format.\n'
                        'Example pattern file content: [[‘O_CC’, ‘C_CO’, ‘C_CO’], [‘C_CC, ‘C_CC’]]')

    modifier_group.add_argument(
                        '-pID-file', '--patternID-file',
                        nargs='?', const='patternIDs.csv', default=None,
                        help='Output file name for assigned patterns (default: patternIDs.csv)')
    
    modifier_group.add_argument(
                        '-s', '--saturation',
                        nargs='?', const='saturatedPolymers.xyz', default=None,
                        help='Saturate the ends polymers.\n'
                        'Optionally provide an output filename (default: saturatedPolymers.xyz)')
    
    modifier_group.add_argument(
                    '-sub', '--subgraph-coords',
                    nargs='?', const='subgraphCoordinates.xyz', default=None,
                    help='Get subgraph-coordinates.\n'
                    'Optionally provide an output filename (default: subgraphCoordinates.xyz)')

    modifier_group.add_argument(
                    '-wrap', '--wrapGraphPrint',
                    nargs='?', const='wrappedGraph.xyz', default=None,
                    help='Write the wrapped graph coordinates. Requires a box size.\n'
                    'Optionally provide an output filename (default: wrappedGraph.xyz)')

    
    args = vars(parser.parse_args())
    
    # Check if the required arguments are provided
    if len(sys.argv) == 1:
        parser.print_help()
        print("\nNo arguments provided. Try again ...")
        sys.exit(1)
    
    return args

def element_distance_tuple(value):
    parts = value.split(':')
    if len(parts) != 4:
        raise argparse.ArgumentTypeError(
            f"Invalid format: '{value}'. Expected format: <element>:H-Acceptor-Distance:Donor-Acceptor-Distance:Angle-Cutoff"
        )
    element, HAcceptor_dist, DonorAcceptor_dist, Angle_cut = parts
    try:
        HAcceptor_dist = float(HAcceptor_dist)
        DonorAcceptor_dist = float(DonorAcceptor_dist)
        Angle_cut = float(Angle_cut)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"Invalid numeric values in '{value}'. Expected format: <element>:H-Acceptor-Distance:Donor-Acceptor-Distance:Angle-Cutoff"
        )
    return (element, HAcceptor_dist, DonorAcceptor_dist, Angle_cut)


def OrderParam(value):
    parts = value.split(':')
    if len(parts) not in (2, 3):
        raise argparse.ArgumentTypeError(
            f"Invalid format: '{value}'. Expected format: BoxSize(x, y, z):n(nx, ny, nz)[:unitSize]"
        )
    boxSize = tuple(map(float, parts[0].split(',')))
    if len(boxSize) != 3:
        boxSize = (boxSize[0], boxSize[0], boxSize[0]) 
    n = tuple(map(int, parts[1].split(',')))
    if len(n) != 3:
        n = (n[0], n[0], n[0])
    unitSize = int(parts[2]) if len(parts) == 3 else None
    return (boxSize, n, unitSize)

def RingCycleSize(value):
    value = value.strip()

    if value.startswith('[') and value.endswith(']'):
        inner = value[1:-1].strip()
        bounds = [item.strip() for item in inner.split(',')]
        if len(bounds) != 2:
            raise argparse.ArgumentTypeError(
                f"Invalid range '{value}'. Expected format: [min,max] with integers and min ring size >= 3."
            )
        try:
            min_size = int(bounds[0])
            max_size = int(bounds[1])
        except ValueError:
            raise argparse.ArgumentTypeError(
                f"Invalid range '{value}'. Expected integers in format [min,max]."
            )
        if min_size < 3:
            raise argparse.ArgumentTypeError("Minimum ring size is 3.")
        if max_size < min_size:
            raise argparse.ArgumentTypeError(
                f"Invalid range '{value}'. Expected max >= min."
            )
        return (min_size, max_size)

    try:
        size = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"Invalid ring size '{value}'. Expected an integer (e.g. 6) or a range [min,max] (e.g. [4,7])."
        )
    if size < 3:
        raise argparse.ArgumentTypeError("Minimum ring size is 3.")
    return size

def _parse_func_type(func_type, value):
    if func_type not in ['CO', 'COH', 'ester', 'amide']:
        raise argparse.ArgumentTypeError(
            f"Invalid func_type in '{value}'. Expected 'CO', 'COH', 'ester' or 'amide'."
        )
    return func_type

def functionalize_pe(value):
    parts = value.split(':')
    if len(parts) < 3:
        raise argparse.ArgumentTypeError(
            "Invalid format. Expected random:<percentage>:<func_type>:<neighbor_exclusion> "
            "or periodic:<distance>:<func_type>."
        )
    mode = parts[0].lower()
    if mode == 'random':
        if len(parts) != 4:
            raise argparse.ArgumentTypeError(
                "Invalid random format. Expected random:<percentage>:<func_type>:<neighbor_exclusion>."
            )
        try:
            percentage = float(parts[1])
        except ValueError:
            raise argparse.ArgumentTypeError(
                f"Invalid percentage value in '{value}'. Expected a float between 0 and 100."
            )
        if not (0 <= percentage <= 100):
            raise argparse.ArgumentTypeError(
                f"Percentage must be between 0 and 100 in '{value}'."
            )
        func_type = _parse_func_type(parts[2], value)
        try:
            neighbor_exclusion = int(parts[3])
        except ValueError:
            raise argparse.ArgumentTypeError(
                f"Invalid neighbor_exclusion in '{value}'. Expected an integer >= 0."
            )
        if neighbor_exclusion < 0:
            raise argparse.ArgumentTypeError(
                f"neighbor_exclusion must be >= 0 in '{value}'."
            )
        return {
            'mode': 'random',
            'percentage': percentage,
            'func_type': func_type,
            'neighbor_exclusion': neighbor_exclusion,
        }
    if mode == 'periodic':
        if len(parts) != 3:
            raise argparse.ArgumentTypeError(
                "Invalid periodic format. Expected periodic:<distance>:<func_type>."
            )
        try:
            distance = int(parts[1])
        except ValueError:
            raise argparse.ArgumentTypeError(
                f"Invalid distance in '{value}'. Expected an integer >= 1."
            )
        if distance < 1:
            raise argparse.ArgumentTypeError(
                f"distance must be >= 1 in '{value}'."
            )
        func_type = _parse_func_type(parts[2], value)
        return {
            'mode': 'periodic',
            'distance': distance,
            'func_type': func_type,
        }
    raise argparse.ArgumentTypeError(
        f"Invalid mode in '{value}'. Expected 'random' or 'periodic'."
    )

def functionalize_pe_surface(value):
    data = functionalize_pe(value)
    if data['mode'] != 'random':
        raise argparse.ArgumentTypeError(
            "Surface functionalization only supports random mode."
        )
    return data
