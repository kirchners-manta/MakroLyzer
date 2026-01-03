"""
Central registry and factories for structure modifiers.
-> Factory: Function that creates and returns an modifier instance or None.

Each entry in `MODIFIERS_REGISTRATION` maps a short key to a factory
callable that accepts `(args, **context)` and returns an modifier
instance or `None` if it should not be created for the given args.
"""
from MakroLyzer.modify_modules.Patterns import PatternModifier
from MakroLyzer.modify_modules.structureModifierBase import ModifyOutputHandler


def create_patterns(args, **context):
    pattern_path = args.get('patternFile')
    if not pattern_path:
        return None
    output_handler = ModifyOutputHandler('Patterns.csv')
    return PatternModifier(pattern_path, output_handler)


MODIFIERS_REGISTRATION = {
    'patternFile': create_patterns,
}
