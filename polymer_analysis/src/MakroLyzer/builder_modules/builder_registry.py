"""
Factories for builders selected through the input argument dictionary.

Each factory accepts ``(args, **context)`` and returns a builder instance,
or ``None`` when that builder was not requested.
"""

from MakroLyzer.builder_modules.BuilderBase import BuilderOutputHandler
from MakroLyzer.builder_modules.ligninBuilder import LigninBuilder


def create_lignin(args, **context):
    """Create a lignin builder from ``lignin`` (recipe) and ``lignin_file``.

    The output path defaults to ``lignin.xyz`` when not supplied.
    """
    recipe_file = args.get('lignin')
    if recipe_file is None:
        return None
    output_file = args.get('lignin_file') or 'lignin.xyz'
    output_handler = BuilderOutputHandler(output_file)
    return LigninBuilder(recipe_file, output_handler,
                         max_iterations=args.get('builder_max_iterations', 1500))


BUILDERS_REGISTRATION = {
    'lignin': create_lignin,
}
