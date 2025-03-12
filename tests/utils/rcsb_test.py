import os
import pytest
from unittest.mock import patch, MagicMock

from cchtools.utils.rcsb import RcsbPdbClusters, pdb_check_obsolete, get_pdb_entities, match_pdb_chain_to_entity


class TestRcsbPdbClusters:
    @pytest.fixture
    def mock_cluster_file(self, tmp_path):
        """Create a mock cluster file for testing"""
        cluster_dir = tmp_path / "clusters"
        cluster_dir.mkdir()
        cluster_file = cluster_dir / "pdb_clusters_30.out"
        
        # Create sample cluster data
        with open(cluster_file, "w") as f:
            f.write("1ABC_A 1XYZ_B\n")
            f.write("2DEF_C 2UVW_D\n")
        
        return str(cluster_dir)
    
    @pytest.fixture
    def clusters(self, mock_cluster_file):
        with patch("cchtools.utils.rcsb.RcsbPdbClusters._download_cluster_sets"):
            return RcsbPdbClusters(identity=30, cluster_dir=mock_cluster_file)
    
    def test_init(self, clusters):
        """Test initialization of RcsbPdbClusters"""
        assert clusters.identity == 30
        assert len(clusters.clusters) == 4
        assert clusters.clusters["1ABC_A"] == 0
        assert clusters.clusters["1XYZ_B"] == 0
        assert clusters.clusters["2DEF_C"] == 1
        assert clusters.clusters["2UVW_D"] == 1
    
    def test_download_cluster_sets(self, tmp_path):
        """Test downloading cluster sets"""

        clusters = RcsbPdbClusters(identity=30, cluster_dir=tmp_path)

        assert clusters.clusters["12E8_1"] == 0

        # Verify the cluster file was created
        cluster_file = os.path.join(tmp_path, "pdb_clusters_30.out")
        assert os.path.exists(cluster_file)

        # Verify the content of the cluster file
        with open(cluster_file, "r") as f:
            content = f.read()
        
        assert "12E8_1" in content
        


    @patch("cchtools.utils.rcsb.match_pdb_chain_to_entity")
    def test_get_seqclust_with_chain(self, mock_match, clusters):
        """Test getting sequence cluster with chain ID"""
        mock_match.return_value = "A"
        
        # Test existing cluster
        clusters.clusters["1ABC_A"] = 42
        assert clusters.get_seqclust("1ABC", chain_id="X") == 42
        mock_match.assert_called_with("1ABC", "X")
        
        # Test non-existing cluster
        # Mock the API response for non-existent PDB
        with patch("cchtools.utils.rcsb.pdb_check_obsolete") as mock_obsolete:
            mock_obsolete.return_value = None
            assert clusters.get_seqclust("9999", chain_id="Z") == "None"
    
    # @patch("cchtools.utils.rcsb.pdb_check_obsolete")
    # @patch("cchtools.utils.rcsb.match_pdb_chain_to_entity")
    # def test_get_seqclust_obsolete(self, mock_match, mock_obsolete, clusters):
    #     """Test getting sequence cluster for obsolete PDB"""
    #     mock_match.return_value = "A"
    #     mock_obsolete.return_value = "1ABC"
        
    #     # Set up clusters
    #     clusters.clusters["1ABC_A"] = 42
        
    #     # Test obsolete PDB redirects to existing one
    #     result = clusters.get_seqclust("9999", chain_id="X", check_obsolete=True)
    #     assert result == 42
    #     mock_obsolete.assert_called_with("9999")
    
    def test_get_pdbs_in_cluster(self, clusters):
        """Test getting all PDBs in a cluster"""
        # Set up clusters
        clusters.clusters = {
            "1ABC_A": 42,
            "2DEF_B": 42,
            "AF_XYZ_C": 42,
            "3GHI_D": 43
        }
        
        # Test without AlphaFold
        pdbs = clusters.get_pdbs_in_cluster(42, include_alphafold=False)
        assert len(pdbs) == 2
        assert "1ABC_A" in pdbs
        assert "2DEF_B" in pdbs
        assert "AF_XYZ_C" not in pdbs
        
        # Test with AlphaFold
        pdbs = clusters.get_pdbs_in_cluster(42, include_alphafold=True)
        assert len(pdbs) == 3
        assert "AF_XYZ_C" in pdbs


@patch("requests.get")
def test_pdb_check_obsolete(mock_get):
    """Test checking if a PDB is obsolete"""
    # Mock response for obsolete PDB
    mock_response = MagicMock()
    mock_response.json.return_value = {
        "1abc": [{
            "status_code": "OBS",
            "superceded_by": ["2xyz"]
        }]
    }
    mock_get.return_value = mock_response
    
    assert pdb_check_obsolete("1abc") == "2xyz"
    mock_get.assert_called_with("https://www.ebi.ac.uk/pdbe/api/pdb/entry/status/1abc")
    
    # Mock response for active PDB
    mock_response.json.return_value = {
        "2xyz": [{
            "status_code": "REL"
        }]
    }
    
    assert pdb_check_obsolete("2xyz") is None


@patch("requests.get")
def test_get_pdb_entities(mock_get):
    """Test getting PDB entities"""
    # Mock response
    mock_response = MagicMock()
    mock_response.json.return_value = {
        "1abc": [
            {"entity_id": "1", "molecule_type": "polypeptide(L)", "in_chains": ["A", "B"]},
            {"entity_id": "2", "molecule_type": "polypeptide(L)", "in_chains": ["C"]}
        ]
    }
    mock_get.return_value = mock_response
    
    entities = get_pdb_entities("1abc")
    assert len(entities) == 2
    assert entities[0]["entity_id"] == "1"
    assert entities[1]["entity_id"] == "2"
    
    # Test error handling
    mock_get.side_effect = Exception("API error")
    assert get_pdb_entities("error") is None


@patch("cchtools.utils.rcsb.get_pdb_entities")
def test_match_pdb_chain_to_entity(mock_get_entities):
    """Test matching PDB chain to entity"""
    # Mock entities
    mock_get_entities.return_value = [
        {"entity_id": "1", "molecule_type": "polypeptide(L)", "in_chains": ["A", "B"]},
        {"entity_id": "2", "molecule_type": "polypeptide(L)", "in_chains": ["C"]},
        {"entity_id": "3", "molecule_type": "DNA", "in_chains": ["D"]}
    ]
    
    # Test matching chains
    assert match_pdb_chain_to_entity("1abc", "A") == "1"
    assert match_pdb_chain_to_entity("1abc", "C") == "2"
    
    # Test non-protein chain
    assert match_pdb_chain_to_entity("1abc", "D") is None
    
    # Test non-existent chain
    assert match_pdb_chain_to_entity("1abc", "Z") is None
