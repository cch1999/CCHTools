

#!NOTE! WARNING: These are not correct yet. NOTE!!!
RING_SKELETONS = {
        # Simple rings
        "5": "c1cccc1",
        "6": "c1ccccc1",
        
        # Bicyclic systems
        "55": "c1cc2cccc2c1",
        "65": "c1ccc2cccc2c1", 
        "66": "c1ccc2ccccc2c1",
        
        # Tricyclic systems
        "555": "c1cc2cc3cccc3c2c1",
        "565a": "c1cc2c3cccccc3c2c1",
        "565b": "c1cc2ccc3cccc3c2c1",
        "655": "c1ccc2cc3cccc3cc2c1",
        "656": "c1ccc2ccc3cccc3c2c1",
        "665a": "c1ccc2cccc3ccccc3c2c1",
        "665b": "c1ccc2cccc3ccccc3c2c1",
        "666a": "c1ccc2cccc3cccccc3c2c1",
        "666b": "c1ccc2ccccc2cc3ccccc3c1",
    }

#!NOTE! WARNING: These are not correct yet. NOTE!!!

RING_ATOMS = {
    'C': 'c',  # aromatic carbon
    'N': 'N',  # pyridine nitrogen
    'NH': 'N',  # pyrrole nitrogen
    'O': 'O',  # aromatic oxygen
    'S': 'S',  # aromatic sulfur
    'C=O': 'C(=O)',  # carbonyl carbon
}