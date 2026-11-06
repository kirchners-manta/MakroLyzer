from MakroLyzer.builder_modules.BuilderBase import Builder

class LigninBuilder(Builder):
    """
    Builds a lignin strand corresponding to the users given sequence of monomers and likages.
    Inherits from Builder.
    """
    
    def __init__(self, recipe_file, output_handler=None):
        super().__init__(output_handler)
        self.recipe = self._read_recipe(recipe_file)