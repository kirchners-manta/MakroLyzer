from numbers import Integral
from pathlib import Path

from MakroLyzer.modify_modules.structureModifierBase import StructureModifier

def write_functionalization_output(output_handler, graph, report, original_input):
    # Use the path assigned by the writer so both files share any collision suffix.
    path = output_handler.write_output(f"{graph.number_of_nodes()}\nelement  x         y        z", graph)
    if not isinstance(report, dict) or path is None:
        return

    lines = [f"Mode: {report['mode']}"]
    if original_input is not None:
        lines.insert(0, f"Input: {original_input}")
    if report['mode'] == 'random':
        lines.append(f"Neighbor exclusion (C-C hops): {report['neighbor_exclusion']}")
    else:
        lines.append(f"Periodic distance: {report['distance']}")
    if report['surface_atom_count'] is None:
        lines.append(f"Original carbon atoms: {report['carbon_count']}")
    else:
        lines.append(f"Detected surface atoms: {report['surface_atom_count']}")
        lines.append(f"Original surface carbon atoms: {report['carbon_count']}")
    for group in report['groups']:
        lines.extend(('', group['type']))
        if group['requested_percentage'] is not None:
            lines.append(f"Requested percentage: {group['requested_percentage']:g}%")
            lines.append(f"Target groups: {group['target_groups']}")
            lines.append(f"Group shortfall: {group['target_groups'] - group['actual_groups']}")
        lines.extend((
            f"Actual groups: {group['actual_groups']}",
            f"Actual percentage: {group['actual_percentage']:.2f}%",
            f"Maximum groups (group used alone): {group['maximum_groups']}",
            f"Maximum percentage (group used alone): {group['maximum_percentage']:.2f}%",
        ))
    Path(path).with_suffix('.info.txt').write_text('\n'.join(lines) + '\n', encoding='utf-8')


class FunctionalizePEModifier(StructureModifier):
    """Functionalize polyethylene by random or periodic placement.

    Random targets and report percentages use the original carbon count.
    Chain ends contribute to that count and are excluded from placement.
    """

    def __init__(self, mode, func_type, percentage=None, neighbor_exclusion=None, distance=None, output_handler=None,
                 *, requests=None, original_input=None):
        """Configure functionalization and optional output.

        Args:
            mode (str): Random percentage-based placement or periodic spacing.
                Use `"random"` or `"periodic"`.
            func_type (str): One of CO, COH, ester, or amide.
            percentage (float, optional): Percentage of original carbons to
                functionalize (0-100). Required in random mode.
            neighbor_exclusion (int, optional): Inclusive exclusion distance
                in C-C bonds for random placement. Defaults to one.
            distance (int, optional): Number of backbone carbons to leave
                between consecutive groups in periodic mode. Must be positive.
            output_handler (ModifyOutputHandler, optional): Writer for the XYZ output.
            requests (list[tuple[str, float]], optional): Ordered
                (group type, percentage) pairs for mixed random placement.
                Each group reserves sites in input order using the same
                exclusion distance. When supplied, these pairs set all group
                targets. The single-group arguments still undergo validation.
            original_input (str, optional): Raw user input to record in the
                companion .info.txt report.
        """
        super().__init__(output_handler)
        if func_type not in ['CO', 'COH', 'ester', 'amide']:
            raise ValueError("func_type must be 'CO', 'COH', 'ester' or 'amide'")
        if mode not in ['random', 'periodic']:
            raise ValueError("mode must be 'random' or 'periodic'.")
        if mode == 'random':
            if percentage is None:
                raise ValueError("percentage is required for random mode.")
            if not (0 <= percentage <= 100):
                raise ValueError("percentage must be between 0 and 100")
        else:
            if not isinstance(distance, Integral) or distance < 1:
                raise ValueError("distance must be a positive integer.")
        
        self.mode = mode
        self.func_type = func_type
        self.percentage = percentage
        self.neighbor_exclusion = neighbor_exclusion
        self.distance = distance
        self.requests = requests
        self.original_input = original_input
        self.report = None
        
    def compute(self, graph):
        """Functionalize eligible backbone carbons in place.

        Store the counts and percentages in `self.report` and return the
        input graph with the modifications applied.
        """
        options = {} if self.requests is None else {'requests': self.requests}
        self.report = graph.functionalizePE(
            self.mode,
            self.func_type,
            percentage=self.percentage,
            neighbor_exclusion=self.neighbor_exclusion,
            distance=self.distance,
            **options,
        )
        return graph

    def render_output(self, graph, frame_idx):
        """Write the modified XYZ and its accompanying .info.txt report.

        Args:
            graph (GraphManager): Graph returned by compute().
            frame_idx (int): Frame index supplied by StructureModifier.run().
                The output writer chooses the filename.
        """
        write_functionalization_output(self.output_handler, graph, self.report, self.original_input)
