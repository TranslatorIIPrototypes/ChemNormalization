import os
import json
import pandas as pd

import redis
from neo4j import GraphDatabase

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.rdmolops import RemoveStereochemistry


##############
# Class: ChemNormalization
#
# By: Phil Owen
# Date: 1/10/2020
# Desc: Class that gets all chemical substances from a neo4j graph database and:
#   gets a simplified SMILEs value for each
#   groups the simplified SMILES and inserts them into a redis database
##############
class ChemNormalization:
    # Storage for the configuration params
    config: json = None

    def __init__(self):
        self.config = self.get_config()
        pass

    def load_redis(self):
        try:
            # get the grouped simplified SMILES
            smiles_list: list = self.get_simplified_smiles_for_chemicals()

            # get a connection to the redis instance for the (SMILE -> Simple SMILES) and (Simple SMILES -> [Similar SMILES, ..]) data
            # smiles_to_simple_smiles_redis = self.get_redis(self.config, '0')
            # simple_smiles_to_similar_smiles_redis = self.get_redis(self.config, '1')

            # get the pipelines for redis loading
            # smiles_to_simple_smiles_redis_pipeline = smiles_to_simple_smiles_redis.pipeline()
            # simple_smiles_to_similar_smiles_pipeline = simple_smiles_to_similar_smiles_redis.pipeline()

            # start loading redis with the data
            for el in smiles_list:
                print(f'test {el}')

            # print(f'Dumping to simple smiles to similar smiles db ...')
            # simple_smiles_to_similar_smiles_pipeline.execute()

            # print(f'Dumping to simple smiles to similar smiles db ...')
            # simple_smiles_to_similar_smiles_pipeline.execute()

        except Exception as e:
            print(f'Bad news {e}')

        pass

    @staticmethod
    def get_config() -> json:
        """ Loads the configuration settings. """

        # get the location of the configuration file
        cname = os.path.join(os.path.dirname(__file__), '..', 'config.json')

        # open the config file
        with open(cname, 'r') as json_file:
            # parse the contents of the configuration settingx
            data = json.load(json_file)

        # return to the caller
        return data

    @staticmethod
    def get_redis(config: json, dbid: str) -> object:
        """ This method returns a connection to a redis instance """
        return redis.StrictRedis(host=config['redis_host'],
                                 port=int(config['redis_port']),
                                 db=dbid,
                                 password=config['redis_password'])

    def get_simplified_smiles_for_chemicals(self) -> list:
        """ This method gets smiles for every chemical substance in the robokop neo4j graph database and creates a simplified SMILES from each.
            The simplified SMILES values will be used as a grouping mechanism and saved in the redis database.
        """

        # init the dicts that will be used to load redis
        terms: list = []

        # Declare the columns in the output data frame
        cols: list = ['similar_smiles', 'simplified_smiles']

        # Create a target dataframe for the processed data
        df = pd.DataFrame(columns=cols)

        try:
            # Create the query. This is of course robokop specific
            # to return 2 that have the same simplified SMILES use this in the where clause ->  and (c.id="CHEBI:140593" or c.id="CHEBI:140451")
            # to return 2 that have the exact same SMILES use this in the where clause ->   and (c.id="CHEBI:85764" or c.id="CHEBI:140773")
            c_query: str = f'match (c:chemical_substance) where c.smiles is not NULL and c.smiles <> "" and c.smiles <> "**" and c.smiles <> "*" RETURN c.id, c.smiles order by c.smiles'  # limit 50000'''

            # execute the query
            records = self.run_neo4j_query(c_query)

            # did we get some records
            if len(records) > 0:
                # loop through the records
                for r in records:
                    try:
                        # Construct a molecule from a SMILES string
                        molecule = Chem.MolFromSmiles(r['c.smiles'])
                    except Exception as e:
                        # alert the user there was an issue and continue
                        print(f"Error - Exception trying to get a molecule for id: {r['c.id']} with smiles: {r['c.smiles']}, Execption {e}. Proceeding.")
                        continue

                    # did we get the molecule
                    if molecule is None:
                        # Couldn't parse the molecule
                        print(f"Error - Got an empty molecule for id: {r['c.id']} with smiles: {r['c.smiles']}. Proceeding.")
                        continue
                    try:
                        # get the uncharged version of the largest fragment
                        molecule_uncharged = rdMolStandardize.ChargeParent(molecule)

                        # Remove all stereo-chemistry info from the molecule
                        RemoveStereochemistry(molecule_uncharged)

                        # get the simplified SMILES value
                        simplified_smiles = Chem.MolToSmiles(molecule_uncharged)

                        # append the record to the data frame
                        df = df.append({'similar_smiles': r['c.smiles'], 'simplified_smiles': simplified_smiles}, ignore_index=True)
                    except Exception as e:
                        # alert the user that something was discovered in the original graph record
                        print(f"Error - Could not get a simplified SMILES for id: {r['c.id']}, SMILES: {r['c.smiles']}, Exception: {e}")

                # remove all duplicates
                df = df.drop_duplicates()

                # get the simplified SMILES in groups
                grps = df.groupby('simplified_smiles')

                # loop through each record and create the proper dict for redis insertion
                for _simplified_smiles, similar_smiles_group in grps:
                    # build up the record to insert
                    term: dict = {'simplified_smiles': _simplified_smiles, 'similar_smiles': [_similar_smiles for _similar_smiles in similar_smiles_group.similar_smiles]}

                    # add it to the dict
                    terms.append(term)
        except Exception as e:
            raise e

        # return to the caller
        return terms

    def get_driver(self) -> GraphDatabase.driver:
        # Gets a connection to the graph database
        driver = GraphDatabase.driver(self.config['neo4j_uri'], auth=(self.config['neo4j_user'], self.config['neo4j_password']))

        # return to the caller
        return driver

    def run_neo4j_query(self, cypherquery: str) -> list:
        """ Runs a cypher query """

        # get a connection to the graph db
        driver = self.get_driver()

        # start a session
        with driver.session() as session:
            # get the results of the query
            results = session.run(cypherquery)

        # return the data to the caller
        return list(results)
