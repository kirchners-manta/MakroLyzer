import pytest

from src.PolyLyzer.input_handling import readInput
from src.PolyLyzer.structure_modules import graphs
from src.PolyLyzer.structure_modules.endToEndDistance import end_to_end_dist
from src.PolyLyzer.structure_modules.dihedrals import get_all_dihedrals, get_CisTrans

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
def pattern_data04():
    return 'test_structures/04pattern.txt'

@pytest.fixture
def sample_data5():
    return 'test_structures/05.xyz'

@pytest.fixture
def sample_data5cis():
    return 'test_structures/05cis.xyz'

@pytest.fixture
def sample_data5cis2():
    return 'test_structures/05cis2.xyz'

@pytest.fixture
def sample_data5trans():
    return 'test_structures/05trans.xyz'

@pytest.fixture
def sample_data6():
    return 'test_structures/06.xyz'

    
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
    xyz = readInput.readXYZ(sample_data4)
    testGraph = graphs.GraphManager(xyz)
    testGraph.find_and_tag_patterns([['C_CCC', 'C_CCC', 'C_CC', 'C_CC', 'C_CC', 'C_CC', 'C_CC', 'C_CC'], ['C_CC', 'C_CCC', 'C_CC', 'C_CC', 'C_CC', 'C_CC', 'C_CC'], ['C_C']])
    
    # iterate over nodes and check fragment ID
    i = 0
    for node in testGraph.nodes():
        if i in range(0,4):
            assert testGraph.nodes[node]['fragmentID'] == 4
        elif i in range(4, 20):
            assert testGraph.nodes[node]['fragmentID'] == 0
        elif i in range(20, 36):
            assert testGraph.nodes[node]['fragmentID'] == 1
        elif i in range(36, 52):
            assert testGraph.nodes[node]['fragmentID'] == 2
        elif i in range(52, 66):
            assert testGraph.nodes[node]['fragmentID'] == 3
        i += 1
        
def test_dihedrals1(sample_data5):
    xyz = readInput.readXYZ(sample_data5)
    testGraph = graphs.GraphManager(xyz)
    
    # absolute dihedrals
    dihedrals = get_all_dihedrals(testGraph, file="dihedrals.csv")
    for angle, count in dihedrals:
        if angle == 65:
            assert count == 1
        elif angle == 177:
            assert count == 2
        elif angle == 180:
            assert count == 1
        else:
            assert count == 0
            
    # cis trans counts
    cisTrans = get_CisTrans(testGraph, file="CisTrans.csv")
    assert cisTrans[0][1] == 1
    assert cisTrans[1][1] == 3
            
    # relative dihedrals
    dihedrals = get_all_dihedrals(testGraph, file="dihedrals.csv", sign=True)
    for angle, count in dihedrals:
        if angle == -65:
            assert count == 1
        elif angle == -177:
            assert count == 2
        elif angle == 180:
            assert count == 1
        else:
            assert count == 0
            
def test_dihedrals2(sample_data5cis):
    xyz = readInput.readXYZ(sample_data5cis)
    testGraph = graphs.GraphManager(xyz)
    
    # absolute dihedrals
    dihedrals = get_all_dihedrals(testGraph, file="dihedrals.csv")
    for angle, count in dihedrals:
        if angle == 64:
            assert count == 1
        else:
            assert count == 0
            
    # cis trans counts
    cisTrans = get_CisTrans(testGraph, file="CisTrans.csv")
    assert cisTrans[0][1] == 1
    assert cisTrans[1][1] == 0
            
    # relative dihedrals
    dihedrals = get_all_dihedrals(testGraph, file="dihedrals.csv", sign=True)
    for angle, count in dihedrals:
        if angle == -64:
            assert count == 1
        else:
            assert count == 0
            
def test_dihedrals3(sample_data5cis2):
    xyz = readInput.readXYZ(sample_data5cis2)
    testGraph = graphs.GraphManager(xyz)
    
    # absolute dihedrals
    dihedrals = get_all_dihedrals(testGraph, file="dihedrals.csv")
    for angle, count in dihedrals:
        if angle == 4:
            assert count == 1
        else:
            assert count == 0
            
    # cis trans counts
    cisTrans = get_CisTrans(testGraph, file="CisTrans.csv")
    assert cisTrans[0][1] == 1
    assert cisTrans[1][1] == 0
            
    # relative dihedrals
    dihedrals = get_all_dihedrals(testGraph, file="dihedrals.csv", sign=True)
    for angle, count in dihedrals:
        if angle == 4:
            assert count == 1
        else:
            assert count == 0
    
def test_dihedrals4(sample_data5trans):
    xyz = readInput.readXYZ(sample_data5trans)
    testGraph = graphs.GraphManager(xyz)
    
    # absolute dihedrals
    dihedrals = get_all_dihedrals(testGraph, file="dihedrals.csv")
    for angle, count in dihedrals:
        if angle == 180:
            assert count == 1
        else:
            assert count == 0
            
    # cis trans counts
    cisTrans = get_CisTrans(testGraph, file="CisTrans.csv")
    assert cisTrans[0][1] == 0
    assert cisTrans[1][1] == 1
            
    # relative dihedrals
    dihedrals = get_all_dihedrals(testGraph, file="dihedrals.csv", sign=True)
    for angle, count in dihedrals:
        if angle == 180:
            assert count == 1
        else:
            assert count == 0

def test_dihedrals5(sample_data6):
    xyz = readInput.readXYZ(sample_data6)
    testGraph = graphs.GraphManager(xyz)
    
    # absolute dihedrals
    dihedrals = get_all_dihedrals(testGraph, file="dihedrals.csv")
    for angle, count in dihedrals:
        if angle == 2:
            assert count == 1
        elif angle == 64:
            assert count == 1
        elif angle == 66:
            assert count == 1
        elif angle == 81:
            assert count == 1
        elif angle == 172:
            assert count == 1
        elif angle == 175:
            assert count == 1
        elif angle == 177:
            assert count == 2
        elif angle == 179:
            assert count == 1
        elif angle == 180:
            assert count == 4
        else:
            assert count == 0
            
    # cis trans counts
    cisTrans = get_CisTrans(testGraph, file="CisTrans.csv")
    assert cisTrans[0][1] == 4
    assert cisTrans[1][1] == 9
    
    # relative dihedrals
    dihedrals = get_all_dihedrals(testGraph, file="dihedrals.csv", sign=True)
    for angle, count in dihedrals:
        if angle == -180:
            assert count == 2
        elif angle == -175:
            assert count == 1
        elif angle == -81:
            assert count == 1
        elif angle == -2:
            assert count == 1
        elif angle == 172:
            assert count == 1
        elif angle == 64:
            assert count == 1
        elif angle == 66:
            assert count == 1
        elif angle == 177:
            assert count == 2
        elif angle == 179:
            assert count == 1
        elif angle == 180:
            assert count == 2
        else:
            assert count == 0
        
    
    
    
    
        
    