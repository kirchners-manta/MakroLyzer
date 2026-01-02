from MakroLyzer.input_handling import readXYZ
from MakroLyzer.input_handling import readLMP
from MakroLyzer.input_handling import estimateFrames
from MakroLyzer.structure_modules import graphs

from MakroLyzer.structure_modules.analyzer_registry import ANALYZERS_REGISTRATION

from tqdm import tqdm

def main(args):
    """
    Perfoms the main analysis of the polymer structure.
    """
    # Get the trajectory file path 
    if args['xyzFile']:
        trajectoryFilePath = args['xyzFile']
        n_frames = estimateFrames.EstimateFrames.estimateFramesXYZ(trajectoryFilePath)
        read = readXYZ.readXYZ
    elif args['lmpFile']:
        trajectoryFilePath = args['lmpFile']
        n_frames = estimateFrames.EstimateFrames.estimateFramesLMP(trajectoryFilePath)
        read = readLMP.readLMP
        
    # Get the modulo for reading frames
    nthStep = args.get('nthStep', 1)
    
    # Get the box size
    boxSize = args.get('BoxSize', None)
        # check if the order parameter is provided
    if args['orderParameter']:
        BoxSize, NoCellsPerDim, MolecularVectorLength = args['orderParameter']
    if not boxSize:
        boxSize = BoxSize[0]
            
    
    # Build analyzers from registry ---------------------------------------------------------------
    analyzers = {}
    # prepare context for factories
    context = {'boxSize': boxSize}
    if args.get('orderParameter'):
        context.update({'BoxSize': BoxSize, 'NoCellsPerDim': NoCellsPerDim, 'MolecularVectorLength': MolecularVectorLength})

    for key, factory in ANALYZERS_REGISTRATION.items():
        try:
            instance = factory(args, **context)
        except Exception:
            instance = None
        if instance is not None:
            analyzers[key] = instance
            if hasattr(instance, 'initialize_output'):
                instance.initialize_output()

    # Loop ----------------------------------------------------
    for i, xyz_frame in enumerate(tqdm(read(trajectoryFilePath), total=n_frames, desc="Creating something magical", unit="frame", ncols=100)):
        if i % nthStep != 0:
            continue
        
        # Get Graph object of the polymer box
        boxGraph = graphs.GraphManager(xyz_frame, boxSize=boxSize)
        
        # Run all registered analyzers for this frame
        for analyzer in analyzers.values():
            analyzer.run(boxGraph, i)
            
    # Finalize output for class-based analyzers ---------------------------------
    for analyzer_name, analyzer in analyzers.items():
        analyzer.finalize_output()
