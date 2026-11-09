"""
Dispatch recipe-based structure construction to registered builders.
"""

from MakroLyzer.builder_modules.builder_registry import BUILDERS_REGISTRATION


def main(args):
    """
    Run requested builders and return their graphs keyed by registry name.

    Factories return ``None`` for builders not selected in ``args``.
    Each selected builder validates its recipe, constructs its graph, and
    writes any configured output through its ``run()`` method.
    Initialization and construction errors propagate to the caller.
    """
    builders = {}
    for key, factory in BUILDERS_REGISTRATION.items():
        builder = factory(args)
        if builder is not None:
            builders[key] = builder

    results = {}
    for key, builder in builders.items():
        results[key] = builder.run()
    return results
