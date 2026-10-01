"""
Central registry and factories for structure analyzers.
-> Factory: Function that creates and returns an analyzer instance or None.

Each entry in `ANALYZERS_REGISTRATION` maps a short key to a factory
callable that accepts `(args, **context)` and returns an analyzer
instance or `None` if it should not be created for the given args.
"""
from MakroLyzer.dynamic_modules.MSD import MSDAnalyzer
from MakroLyzer.dynamic_modules.dynamicBase import OutputHandler


def create_MSD(args, **context):
    if not args.get('MSD'):
        return None
    correlation_depth = args['MSD']
    if correlation_depth is None:
        return None
    MSD_file = args.get('MSD_file') or 'MSD.csv'
    output_handler = OutputHandler(MSD_file, mode='collect')
    return MSDAnalyzer(
        correlation_depth=correlation_depth,
        timestep=context['timestep'],
        output_handler=output_handler,
    )

ANALYZERS_REGISTRATION = {
    'MSD': create_MSD,
}
