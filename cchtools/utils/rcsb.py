import logging
import os
import urllib.request as request
from functools import lru_cache
from typing import List, Optional

import requests

logger = logging.getLogger(__name__)


class RcsbPdbClusters:
    def __init__(self, identity: int = 30, cluster_dir: str = os.path.expanduser("~/.cache/cchtools")):
        """
        Class for fetching sequence cluster IDs for a given PDB code and chain ID
        using RCSB mmseq2/blastclust predefined clusters.

        Clusters info at https://www.rcsb.org/docs/programmatic-access/file-download-services

        Args:
            identity (int): Identity threshold for clustering
            cluster_dir (str): Directory to store cluster files (default: ~/.cache/cchtools)
        """
        self.cluster_dir = cluster_dir
        self.identity = identity
        self.clusters: dict[str, int] = {}
        self._fetch_cluster_file()

    def _download_cluster_sets(self, cluster_file_path: str):
        """Download cluster file from RCSB"""
        os.makedirs(os.path.dirname(cluster_file_path), exist_ok=True)

        # Note that the files changes frequently as do the ordering of cluster within
        request.urlretrieve(
            f"https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-{self.identity}.txt",
            cluster_file_path,
        )

    def _fetch_cluster_file(self):
        """Load cluster file if found else download and load"""

        cluster_file_path = os.path.join(self.cluster_dir, f"pdb_clusters_{self.identity}.out")
        logging.info(f"cluster file path: {cluster_file_path}")

        # Fetch cluster file
        if not os.path.exists(cluster_file_path):
            logging.warning("Cluster definition not found, will download a fresh one.")
            logging.warning(
                "However, this will very likely lead to silent incompatibilities with any old 'pdbcode_mappings.pickle' files! Please better remove those manually."
            )
            self._download_cluster_sets(cluster_file_path)

        # Extract cluster IDs
        for n, line in enumerate(open(cluster_file_path, "rb")):
            for id in line.decode("ascii").split():
                self.clusters[id] = n

    def get_seqclust(
        self,
        pdb_code: str,
        entity_id: Optional[str] = None,
        chain_id: Optional[str] = None,
        check_obsolete: bool = True,
    ) -> Optional[str]:
        """Get sequence cluster ID for a pdb_code chain using RCSB mmseq2/blastclust predefined clusters

        When `check_obsolete` is True, the function will check if the PDB code is obsolete and if so, will return the cluster ID for the superceding PDB code.

        Args:
            pdb_code (str): PDB code
            entity_id (str): Entity ID
            chain_id (str): Chain ID
            check_obsolete (bool): Check if PDB code is obsolete

        Returns:
            str: Cluster ID as a string
            None: If unable to match to entity_id or chain_id
        """

        if entity_id and chain_id:
            raise Exception("Only define one of either `chain_id` or `entity_id`")

        # Fetch entity_id if only given chain id
        if chain_id:
            entity_id = match_pdb_chain_to_entity(pdb_code, chain_id)
        if entity_id is None:
            logging.info(f"unable to match to entity_id for {pdb_code}_{chain_id}")
            return None

        query_str = f"{pdb_code.upper()}_{str(entity_id).upper()}"  # e.g. 1ATP_I
        seqclust = self.clusters.get(query_str, "None")

        if check_obsolete and seqclust == "None":
            superceded = pdb_check_obsolete(pdb_code)
            if superceded is not None:
                logging.info(
                    f"Assigning cluster for obsolete entry via superceding: {pdb_code}->{superceded} {chain_id}"
                )
                return self.get_seqclust(
                    superceded, chain_id=chain_id, check_obsolete=False
                )  # assumes chain is same in superceding entry
        if seqclust == "None":
            logging.info(f"unable to assign cluster to {pdb_code}{chain_id}")
            return None
        return str(seqclust)

    def get_pdbs_in_cluster(self, cluster_id: str | int, include_alphafold: bool = False) -> List[str]:
        """Get all PDBs in a given cluster"""
        if isinstance(cluster_id, int):
            cluster_id = str(cluster_id)
        pdb_ids = [pdb for pdb, id in self.clusters.items() if str(id) == cluster_id]
        if not include_alphafold:
            pdb_ids = [pdb for pdb in pdb_ids if not pdb.startswith("AF_")]
        return pdb_ids


@lru_cache()
def pdb_check_obsolete(pdb_code: str) -> Optional[str]:
    """Check the status of a pdb, if it is obsolete return the superceding PDB ID else return None"""
    pdb_code = pdb_code.lower()

    try:
        r = requests.get(f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/status/{pdb_code}").json()
    except Exception as e:
        logger.info(f"Could not check obsolete status of {pdb_code}: {e}")
        return None
    if r[pdb_code][0]["status_code"] == "OBS":
        pdb_code = r[pdb_code][0]["superceded_by"][0]
        return pdb_code
    else:
        return None


@lru_cache()
def get_pdb_entities(pdb_code: str):
    """
    Tries to fetch the macromolecular entities of a PDB code from the EBI
    else returns None. See https://www.rcsb.org/docs/general-help/identifiers-in-pdb
    """
    try:
        pdb_code = pdb_code.lower()
        response = requests.get(f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{pdb_code}")
        entities = response.json()[pdb_code]
        return entities
    except Exception as e:
        logging.error(f"Error fetching PDB entities for {pdb_code}: {e}")
        return None


@lru_cache()
def match_pdb_chain_to_entity(pdb_code: str, chain_id: str) -> Optional[str]:
    """
    Tries to match a PDB file chain ID to an macromolecular entitiy, else returns None.
    See https://www.rcsb.org/docs/general-help/identifiers-in-pdb
    """

    pdb_code = pdb_code.lower()
    chain_id = chain_id.upper()

    # Get PDB entities
    entities = get_pdb_entities(pdb_code)

    for entity in entities:
        # Check if chain present in entity and protein then return 'entity_id'
        if chain_id in entity["in_chains"] and entity["molecule_type"] == "polypeptide(L)":
            return entity["entity_id"]
    else:
        return None


if __name__ == "__main__":
    clusters = RcsbPdbClusters(identity=100)
    cluster_id = clusters.get_seqclust("6QR8", chain_id="A")
    print(cluster_id)
    pdbs = clusters.get_pdbs_in_cluster(str(cluster_id))
    print(pdbs)
    print(len(pdbs))
