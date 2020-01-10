import os
import json
import pandas as pd

import redis
from neo4j import GraphDatabase

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.rdmolops import RemoveStereochemistry
from rdkit.Chem.SaltRemover import SaltRemover


class ChemNormalization:
    # Storage for the configuration params
    config = None

    def __init__(self):
        self.config = self.get_config()
        pass

    def load_redis(self):
        """ Entry point for this processing. """
        config = self.get_config()

    def get_config(self):
        """ Loads the configuration settings. """

        # get the location of the configuration file
        cname = os.path.join(os.path.dirname(__file__), '..', 'config.json')

        # open the config file
        with open(cname, 'r') as json_file:
            # parse the contents of the configuration settingx
            data = json.load(json_file)

            # return to the caller
        return data

    def get_redis(config, dbid: str):
        """ Return a connection to a redis instance """
        return redis.StrictRedis(host=config['redis_host'],
                                 port=int(config['redis_port']),
                                 db=dbid,
                                 password=config['redis_password'])

    def get_simplified_smiles_for_chemicals(self) -> pd.DataFrame:
        """ This method gets a simplified smiles for every chemical substance in the robokop neo4j graph database.

            It sets up a query to the graph database to get all the chemical substances and their SMILES value.
            It then spools through each record and gets an associated simplified SMILES for the molecule using rdkit.
            The simplified SMILES values will be used as a grouping mechanism in the data saved in the redis database.
        """
        # Declare the columns in the output data frame
        cols:set = ['orig_smiles', 'simplified_smiles']

        # Create a target dataframe for the processed data
        df = pd.DataFrame(columns = cols)

        try:
            # Create the query. This is of course robokop specific
            cquery = f'''match (c:chemical_substance) where c.smiles is not NULL and c.smiles <> "" and c.smiles <> "**" and c.smiles <> "*" RETURN c.id, c.smiles order by c.smiles limit 10000'''

            # execute the query
            records = self.run_neo4j_query(url, cquery)

            # did we get some records
            if len(records) > 0:
                # loop through the records
                for r in records:
                    # do not process records with no SMILES value
                    if r.smiles == '[empty]':
                        print(f"Error - Got an empty SMILES for id: {r.id}. Proceeding.")
                        continue

                    try:
                        # Construct a molecule from a SMILES string
                        molecule = Chem.MolFromSmiles(r.smiles)
                    except Exception as e:
                        # alert the user there was an issue and continue
                        print(f"Error - Exception trying to get a molecule for id: {r.id} with smiles: {r.smiles}. Proceeding.")
                        continue

                    # did we get the molecule
                    if molecule is None:
                        # Couldn't parse the molecule
                        print(f"Error - Got an empty molecule for id: {r.id} with smiles: {r.smiles}. Proceeding.")
                        continue
                    try:
                        # get the uncharged version of the largest fragment
                        molecule_uncharged = rdMolStandardize.ChargeParent(molecule)

                        # Remove all stereo-chemistry info from the molecule
                        RemoveStereochemistry(molecule_uncharged)

                        # get the simplified SMILES value
                        simplified_smile = Chem.MolToSmiles(molecule_uncharged)

                        # append the record to the data frame
                        df = df.append({'orig_smiles': r.smiles, 'simplified_smile': simplified_smile}, ignore_index=True)
                    except Exception as e:
                        # alert the user that something was discovered in the original graph record
                        print(f'Error - Could not get a simplified SMILES for id: {r.id}, SMILES: {r.smiles}')

                # get the simplified SMILES in groups
                df = df.groupby('simplified_smile')

                # loop through each record and create the proper dict for redis insertion

        except Exception as e:
            print(f'Error - Exception thrown getting the data from the neo4j database.')

        # return to the caller
        return df

    def get_driver():
        # get a connection to the graph database
        driver = GraphDatabase.driver(config.neo4j_url, auth=("neo4j", "ncatsgamma"))

        # return to the caller
        return driver

    def run_neo4j_query(self, url: str, cypherquery: str):
        """ runs a cypher query """

        # get a connection to the graph db
        driver = self.get_driver(url)

        # start a session
        with driver.session() as session:
            # get the results of the query
            results = session.run(cypherquery)

        # return the data to the caller
        return list(results)
