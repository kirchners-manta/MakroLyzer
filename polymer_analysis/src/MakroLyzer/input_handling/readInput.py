import argparse
import sys

ASCII_ART = r"""
                             
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

class ArtParser(argparse.ArgumentParser):
    def format_help(self):
        return ASCII_ART + "\n" + super().format_help()

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
    input_group = parser.add_argument_group('----------------| Arguments for specifying input files and settings |----------------')
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
    
    
    # Structure Analyzers Group # -------------------------------------------------------------------------------
    analyzer_group = parser.add_argument_group(title='----------------| Arguments for structure analyzers |----------------',
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
                        'Parameters that need to be specified:   BoxSize : NoCellsPerDim : MolecularVectorLength\n'
                        '                    BoxSize : Size of the simulation box in x, y, z directions (float,float,float) or one value if symmetrical\n'
                        '                    NoCellsPerDim : Number of cells in each direction (int, int, int) or only one value if symmetrical\n'
                        '                    MolecularVectorLength : Number of atoms that are used to calculate a molecular vector from (int)\n'
                        'Usage example: -op 100:1:3',
                        type=OrderParam)
    
    analyzer_group.add_argument(
                        '-op-file', '--order-file',
                        nargs='?', const='orderParameter.csv', default='orderParameter.csv',
                        help='Output filename for the nematic order parameter S* (default: orderParameter.csv)')
    
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
    
    
    # Structure Modifiers Group # ------------------------------------------------------------------------------
    modifier_group = parser.add_argument_group('----------------| Arguments for structure modifiers |----------------')
    
    modifier_group.add_argument(
                        '-p', '--patternFile',
                        dest='patternFile',
                        help='Assigns specific pattern IDs to fragments based on predefined patterns.\n'
                        'Required argument: Path to the TXT file containing the patterns in list of lists format.\n'
                        'Example pattern file content: [[C,C,C],[C,C,O]]')

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
                    help='Get subgraph-coordinates; optionally provide output filename (default: subgraphCoordinates.xyz)')

    
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
    if len(parts) != 3:
        raise argparse.ArgumentTypeError(
            f"Invalid format: '{value}'. Expected format: BoxSize(x, y, z):n(nx, ny, nz):unitSize"
        )
    boxSize = tuple(map(float, parts[0].split(',')))
    if len(boxSize) != 3:
        boxSize = (boxSize[0], boxSize[0], boxSize[0]) 
    n = tuple(map(int, parts[1].split(',')))
    if len(n) != 3:
        n = (n[0], n[0], n[0])
    unitSize = int(parts[2])
    return (boxSize, n, unitSize)
