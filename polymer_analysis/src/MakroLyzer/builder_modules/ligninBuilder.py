from MakroLyzer.builder_modules.BuilderBase import Builder
from MakroLyzer.builder_modules import globalHelpers, liginConnectionLogic
from MakroLyzer.builder_modules.template_io import load_template_graph
from pathlib import Path

class LigninBuilder(Builder):
    """
    Builds a lignin strand corresponding to the users given sequence of monomers and likages.
    Inherits from Builder.
    """
    
    def __init__(self, recipe_file, output_handler=None, max_iterations=1500):
        super().__init__(output_handler)
        if isinstance(max_iterations, bool) or not isinstance(max_iterations, int) or max_iterations < 0:
            raise ValueError("max_iterations must be a positive integer.")
        self.max_iterations = max_iterations
        self.recipe = globalHelpers.read_recipe(recipe_file) # globally for all builders
        
    SUPPORTED_MONOMERS = frozenset({'H', 'G', 'S'})
    SUPPORTED_LINKERS = frozenset({
        'beta-O-4', 'beta-5', 'beta-beta', '5-5', '4-O-beta',
        '5-O-4', 'beta-1',
    })

    def load_monomer(self, monomer, monomer_id=0, start_index=0):
        """Load an H, G, or S template with explicit bonds and fresh node data."""
        if not isinstance(monomer, str) or monomer not in self.SUPPORTED_MONOMERS:
            raise ValueError(f"Unsupported monomer: {monomer!r}.")
        directory = Path(__file__).resolve().parent / 'buildingBlocks' / 'lignin'
        return load_template_graph(
            directory / f'{monomer}-topo.xyz', directory / 'bonding-info.txt',
            monomer_type=monomer, monomer_id=monomer_id, start_index=start_index,
            allowed_elements={'C', 'H', 'O'}, skip_missing_bonds=True,
        )

    def _connect_monomers(self, chain, linker, monomer, *, chain_monomer_id):
        """Dispatch a linkage, returning a new graph without editing inputs.

        ``chain_monomer_id`` identifies the preceding recipe monomer, rather
        than selecting an arbitrary matching atom anywhere in the chain.
        """
        if linker == '4-O-beta':
            return liginConnectionLogic.connect_4_O_beta(chain, monomer, chain_monomer_id)
        elif linker == 'beta-O-4':
            return liginConnectionLogic.connect_beta_O_4(chain, monomer, chain_monomer_id)
        elif linker == 'beta-5':
            return liginConnectionLogic.connect_beta_5(chain, monomer, chain_monomer_id)
        elif linker == 'beta-beta':
            return liginConnectionLogic.connect_beta_beta(
                chain, monomer, chain_monomer_id, max_iterations=self.max_iterations,
            )
        else:
            raise NotImplementedError(f"Linkage {linker!r} is not implemented yet.")

    @staticmethod
    def _forbidden_linker_pairs():
        """Return forbidden ordered pairs of consecutive linkers.

        The first linker uses a site on the incoming middle monomer; the
        next linker must not reuse that site when extending the chain.
        These rules apply to all monomer types and are not a complete
        chemical compatibility check.
        """
        return {
            ('beta-O-4', '4-O-beta'): 'beta position of the shared monomer is already occupied',
            ('4-O-beta', 'beta-O-4'): 'O-4 position of the shared monomer is already occupied',
            ('4-O-beta', 'beta-5'): 'O-4 position of the shared monomer is already occupied',
            ('beta-5', '4-O-beta'): 'beta position of the shared monomer is already occupied',
            ('beta-5', 'beta-beta'): 'beta position of the shared monomer is already occupied',
            ('beta-O-4', 'beta-beta'): 'beta position of the shared monomer is already occupied',
            ('beta-beta', '4-O-beta'): 'beta position of the shared monomer is already occupied',
        }

    @staticmethod
    def _forbidden_monomer_linkages():
        """Map ordered (existing monomer, linker) triples to reasons.

        Rules are directional: add the reverse triple separately if needed.
        Each monomer value is an explicit H, G, or S; no wildcards are used.
        """
        return {
            ('S', 'beta-5'): 'S does not have a free 5 position',
        }
        
    def validate_recipe(self):
        """Validate a linear recipe, raising ValueError for invalid input."""
        if not isinstance(self.recipe, dict):
            raise ValueError("Recipe must be a dictionary.")

        for key in ('monomer', 'linker'):
            if key not in self.recipe:
                raise ValueError(f"Recipe is missing required field '{key}'.")
            if not isinstance(self.recipe[key], list):
                raise ValueError(f"Recipe field '{key}' must be a list.")

        monomers = self.recipe['monomer']
        linkers = self.recipe['linker']
        if not monomers:
            raise ValueError("Recipe must contain at least one monomer.")
        if len(monomers) != len(linkers) + 1:
            raise ValueError(
                "A linear recipe requires exactly one fewer linker than monomers "
                f"(got {len(monomers)} monomers and {len(linkers)} linkers)."
            )

        for key, values, supported in (
            ('monomer', monomers, self.SUPPORTED_MONOMERS),
            ('linker', linkers, self.SUPPORTED_LINKERS),
        ):
            for i, value in enumerate(values):
                if not isinstance(value, str) or value not in supported:
                    raise ValueError(
                        f"Unsupported {key} at index {i}: {value!r}. "
                        f"Supported values: {', '.join(sorted(supported))}."
                    )

        forbidden = self._forbidden_monomer_linkages()
        for i, linker in enumerate(linkers):
            combination = (monomers[i], linker)
            if combination in forbidden:
                raise ValueError(
                    f"Forbidden monomer linkage at linker index {i}: "
                    f"{combination!r}. {forbidden[combination]}"
                )

        forbidden = self._forbidden_linker_pairs()
        for i, combination in enumerate(zip(linkers, linkers[1:])):
            if combination in forbidden:
                raise ValueError(
                    f"Forbidden combination at linker indices {i} and {i + 1} "
                    f"(shared monomer index {i + 1}): "
                    f"{combination!r}. {forbidden[combination]}"
                )
        
    def build_structure(self):
        """
        1) Load first monomer and call it chain
        2) For each -linker-monomer tuple:
            - load monomer
            - connect monober to chain via the linker
            - energy optimization        
        """
        
        chain = self.load_monomer(self.recipe["monomer"][0])
        for i, linker in enumerate(self.recipe["linker"]):
            
            monomer = self.load_monomer(self.recipe["monomer"][i + 1], monomer_id=i + 1, start_index=max(chain.nodes) + 1)
            chain = self._connect_monomers(chain, linker, monomer, chain_monomer_id=i)
            chain = globalHelpers.optimize_structure(chain, max_iterations=self.max_iterations)

        # If only one monomer is given, optimize it
        if not self.recipe["linker"]:
            return globalHelpers.optimize_structure(chain, max_iterations=self.max_iterations)

        # Retry relaxation of the completed structure, retaining the previous
        # geometry if the new run (which adds small coordinate noise) is worse.
        previous_energy = chain.graph['uff_optimization']['final_energy']
        relaxed = globalHelpers.optimize_structure(chain, max_iterations=self.max_iterations)
        if relaxed.graph['uff_optimization']['final_energy'] < previous_energy:
            return relaxed
        
        return chain
    
    def render_output(self, data):
        """Write the built lignin graph as XYZ and return the output path."""
        if self.output_handler is None:
            return None
        header = f"{data.number_of_nodes()}\nLignin structure"
        return self.output_handler.write_output(header, data)
    
