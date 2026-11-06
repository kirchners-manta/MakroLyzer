"""Base class and coordinate output handler for builder modules."""

from abc import ABC, abstractmethod
import os


class BuilderOutputHandler:
    """Write a built structure's coordinates and optional atom attributes."""

    def __init__(self, file_path):
        self.file_path = file_path

    def write_output(self, header, graph, attributes=None, ID=None):
        """Write XYZ-style output without overwriting an existing file.

        Args:
            header (str): Atom-count line and comment line, separated by a newline.
            graph: Molecular graph with element, x, y, and z node attributes.
            attributes: Optional node attribute names to append as extra columns.
            ID: Optional fragment identifier, including zero.

        Returns:
            str: Path of the written file.
        """
        base, ext = os.path.splitext(self.file_path)
        new_path = f"{base}_frag{ID}{ext}" if ID is not None else os.fspath(self.file_path)
        base, ext = os.path.splitext(new_path)
        counter = 1
        while os.path.exists(new_path):
            new_path = f"{base}_{counter}{ext}"
            counter += 1

        attributes = tuple(attributes) if attributes is not None else ()
        with open(new_path, 'w') as file:
            file.write(header + '\n')
            for _, data in graph.nodes(data=True):
                element = data['element']
                x, y, z = data['x'], data['y'], data['z']
                row = f"{element:6} {x:10.4f} {y:10.4f} {z:10.4f}"
                for attribute in attributes:
                    row += f"    {data.get(attribute, '')}"
                file.write(row + '\n')
        return new_path


class Builder(ABC):
    """Validate a recipe, build a structure, and optionally write its output.

    Concrete builders store their recipe in their own constructor and implement
    the validation, construction, and output formatting methods below.
    """

    def __init__(self, output_handler=None):
        self.output_handler = output_handler

    @abstractmethod
    def validate_recipe(self):
        """Validate the recipe, raising an exception for invalid input."""
        pass

    @abstractmethod
    def build_structure(self):
        """Construct and return the molecular graph from the stored recipe."""
        pass

    @abstractmethod
    def render_output(self, data):
        """Format and write the graph returned by build_structure()."""
        pass

    def run(self):
        """Validate, construct, optionally write, and return the structure."""
        self.validate_recipe()
        result = self.build_structure()
        if self.output_handler is not None:
            self.render_output(result)
        return result
