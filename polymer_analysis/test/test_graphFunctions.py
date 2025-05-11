import pytest

from src.PolyLyzer.input_handling import readInput
from src.PolyLyzer.structure_modules import graphs
from src.PolyLyzer.structure_modules.endToEndDistance import end_to_end_dist

@pytest.fixture
def sample_data1():
    return 'test_structures/01.xyz'

@pytest.fixture
def sample_data2():
    return 'test_structures/02.xyz'

@pytest.fixture
def sample_data3():
    return 'test_structures/03.xyz'

@pytest.fixture
def sample_data4():
    return 'test_structures/04.xyz'

@pytest.fixture
def pattern_data():
    return 'test_structures/04pattern.txt'
    
def test_end_to_end(sample_data1):
    xyz = readInput.readXYZ(sample_data1)
    testGraph = graphs.GraphManager(xyz)
    distance = end_to_end_dist(testGraph)
    correctDistance = 5.595 
    assert distance[0] == pytest.approx(correctDistance, abs=1e-3)
    
def test_remove_1order(sample_data1):
    xyz = readInput.readXYZ(sample_data1)
    testGraph = graphs.GraphManager(xyz)
    
    # number of nodes
    assert testGraph.number_of_nodes() == 26
    newGraph = testGraph.remove_1order()
    assert newGraph.number_of_nodes() == 8
    
    # type of nodes
    for node in newGraph.nodes():
        assert newGraph.nodes[node]['element'] == 'C'
        
def test_update_degree(sample_data1):
    xyz = readInput.readXYZ(sample_data1)
    testGraph = graphs.GraphManager(xyz)
    testGraph.surrounding()
    
    # degree of C nodes
    for node in testGraph.nodes():
        if testGraph.nodes[node]['element'] == 'C':
            assert testGraph.nodes[node]['degree'] == 4
            
    newGraph = testGraph.remove_1order()
    newGraph.update_degree()
    
    for node in newGraph.nodes():
        if newGraph.nodes[node]['element'] == 'C' and newGraph.nodes[node]['surroundingAtoms'] == 'C_CCHH':
            assert newGraph.nodes[node]['degree'] == 2
        elif newGraph.nodes[node]['element'] == 'C' and newGraph.nodes[node]['surroundingAtoms'] == 'C_CHHH':
            assert newGraph.nodes[node]['degree'] == 1
        
def test_longest_path(sample_data1):
    xyz = readInput.readXYZ(sample_data1)
    testGraph = graphs.GraphManager(xyz)
    testGraph.surrounding()
    
    # get longest path 
    longestPath = testGraph.find_longest_path()
    assert len(longestPath) == 10
    path = (['H_C', 'C_CHHH', 'C_CCHH', 'C_CCHH', 'C_CCHH', 'C_CCHH', 'C_CCHH', 'C_CCHH', 'C_CHHH', 'H_C'])
    for i in range(len(longestPath)):
        assert testGraph.nodes[longestPath[i]]['surroundingAtoms'] == path[i]
    
    newGraph = testGraph.remove_1order()
    newGraph.update_degree()
    longestPath = newGraph.find_longest_path()
    assert len(longestPath) == 8
    path = (['C_CHHH', 'C_CCHH', 'C_CCHH', 'C_CCHH', 'C_CCHH', 'C_CCHH', 'C_CCHH', 'C_CHHH'])
    for i in range(len(longestPath)):
        assert testGraph.nodes[longestPath[i]]['surroundingAtoms'] == path[i]
    
def test_surrounding(sample_data2):
    xyz = readInput.readXYZ(sample_data2)
    testGraph = graphs.GraphManager(xyz)
    testGraph.surrounding()
    
    # check surrounding atoms
    for node in testGraph.nodes():
        if testGraph.nodes[node]['element'] == 'C':
            assert testGraph.nodes[node]['surroundingAtoms'] == 'C_MgNOP'
        elif testGraph.nodes[node]['element'] == 'O':
            assert testGraph.nodes[node]['surroundingAtoms'] == 'O_CH'
        elif testGraph.nodes[node]['element'] == 'N':
            assert testGraph.nodes[node]['surroundingAtoms'] == 'N_CHH'
        elif testGraph.nodes[node]['element'] == 'P':
            assert testGraph.nodes[node]['surroundingAtoms'] == 'P_CHH'
        elif testGraph.nodes[node]['element'] == 'Mg':
            assert testGraph.nodes[node]['surroundingAtoms'] == 'Mg_CH'
            
def test_chemicalFormula(sample_data3):
    xyz = readInput.readXYZ(sample_data3)
    testGraph = graphs.GraphManager(xyz)
    testGraph.surrounding()
    refformula = [('C1H6Mg1N1O1P1', 2), ('C1H5Mg1N1O1P1', 1)]
    
    # check chemical formula
    formula = testGraph.get_chemicalFormulas()
    for i in range(len(refformula)):
        assert formula[i][0] == refformula[i][0]
        assert formula[i][1] == refformula[i][1]
        

def test_find_patterns(sample_data4):
    print("ANANAS")