"""Check alignment edge cases shared by lignin placement functions."""

import numpy as np
import pytest

from MakroLyzer.builder_modules import ligninHelpers as helpers


@pytest.mark.parametrize('source,target', [
    ([2, 0, 0], [4, 0, 0]),
    ([2, 0, 0], [-4, 0, 0]),
    ([0, 3, 0], [0, 0, -2]),
    ([1, -2, 3], [-4, 5, 2]),
])
def test_alignment_maps_directions_with_proper_rotation(source, target):
    rotation = helpers.alignment_rotation(source, target)
    np.testing.assert_allclose(
        rotation.apply(helpers.unit_vector(source)), helpers.unit_vector(target), atol=1e-12,
    )
    assert np.linalg.det(rotation.as_matrix()) == pytest.approx(1)


@pytest.mark.parametrize('vector', [[0, 0, 0], [np.nan, 1, 0], [np.inf, 0, 0]])
def test_alignment_rejects_invalid_directions(vector):
    with pytest.raises(ValueError):
        helpers.alignment_rotation(vector, [1, 0, 0])
