# test BERN NER
import requests
from bs4 import BeautifulSoup
import re
from bioregistry import parse_curie

class Paper:
    def __init__(self, title=None, abstract=None, pmids=None):
        self.title = title
        self.abstract = abstract
        self.pmids = pmids
        self.bern_res = []
        self.adj_list = []

        self.gen_funcs = {
            'synonyms': lambda a: 'https://id.nlm.nih.gov/mesh/sparql?query=PREFIX%20rdf%3A%20%3Chttp%3A%2F%2Fwww.w3.org%2F1999%2F02%2F22-rdf-syntax-ns%23%3E%0D%0APREFIX%20rdfs%3A%20%3Chttp%3A%2F%2Fwww.w3.org%2F2000%2F01%2Frdf-schema%23%3E%0D%0APREFIX%20xsd%3A%20%3Chttp%3A%2F%2Fwww.w3.org%2F2001%2FXMLSchema%23%3E%0D%0APREFIX%20owl%3A%20%3Chttp%3A%2F%2Fwww.w3.org%2F2002%2F07%2Fowl%23%3E%0D%0APREFIX%20meshv%3A%20%3Chttp%3A%2F%2Fid.nlm.nih.gov%2Fmesh%2Fvocab%23%3E%0D%0APREFIX%20mesh%3A%20%3Chttp%3A%2F%2Fid.nlm.nih.gov%2Fmesh%2F%3E%0D%0APREFIX%20mesh2015%3A%20%3Chttp%3A%2F%2Fid.nlm.nih.gov%2Fmesh%2F2015%2F%3E%0D%0APREFIX%20mesh2016%3A%20%3Chttp%3A%2F%2Fid.nlm.nih.gov%2Fmesh%2F2016%2F%3E%0D%0APREFIX%20mesh2017%3A%20%3Chttp%3A%2F%2Fid.nlm.nih.gov%2Fmesh%2F2017%2F%3E%0D%0A%0D%0ASELECT%20%3FEntryTerm%0D%0AFROM%20%3Chttp%3A%2F%2Fid.nlm.nih.gov%2Fmesh%3E%0D%0AWHERE%20%7B%0D%0A%20%20mesh%3A{}%20meshv%3Aconcept%20%3Fconcept%20.%0D%0A%20%20%3Fconcept%20meshv%3Aterm%20%3Fterm%20.%0D%0A%20%20%3Fterm%20rdfs%3Alabel%20%3FEntryTerm%0D%0A%7D%20%0D%0A%0D%0A%23Entry%20Terms%20for%20Ofloxacin%20(D015242)%0D%0A&format=JSON&limit=50&offset=0&inference=true'.format(a),
            'descriptor': lambda a: 'https://id.nlm.nih.gov/mesh/sparql?query=PREFIX%20rdf%3A%20%3Chttp%3A%2F%2Fwww.w3.org%2F1999%2F02%2F22-rdf-syntax-ns%23%3E%0D%0APREFIX%20rdfs%3A%20%3Chttp%3A%2F%2Fwww.w3.org%2F2000%2F01%2Frdf-schema%23%3E%0D%0APREFIX%20xsd%3A%20%3Chttp%3A%2F%2Fwww.w3.org%2F2001%2FXMLSchema%23%3E%0D%0APREFIX%20owl%3A%20%3Chttp%3A%2F%2Fwww.w3.org%2F2002%2F07%2Fowl%23%3E%0D%0APREFIX%20meshv%3A%20%3Chttp%3A%2F%2Fid.nlm.nih.gov%2Fmesh%2Fvocab%23%3E%0D%0APREFIX%20mesh%3A%20%3Chttp%3A%2F%2Fid.nlm.nih.gov%2Fmesh%2F%3E%0D%0APREFIX%20mesh2015%3A%20%3Chttp%3A%2F%2Fid.nlm.nih.gov%2Fmesh%2F2015%2F%3E%0D%0APREFIX%20mesh2016%3A%20%3Chttp%3A%2F%2Fid.nlm.nih.gov%2Fmesh%2F2016%2F%3E%0D%0APREFIX%20mesh2017%3A%20%3Chttp%3A%2F%2Fid.nlm.nih.gov%2Fmesh%2F2017%2F%3E%0D%0A%0D%0A%20SELECT%20%3Fd%20%3FdName%0D%0A%20FROM%20%3Chttp%3A%2F%2Fid.nlm.nih.gov%2Fmesh%3E%0D%0A%20%0D%0A%20WHERE%20%7B%0D%0A%20%20%20%0D%0A%20%3Fd%20a%20meshv%3ADescriptor%20.%0D%0A%20%3Fd%20meshv%3Aactive%201%20.%0D%0A%20%3Fd%20rdfs%3Alabel%20%3FdName%20.%0D%0A%20FILTER(REGEX(%3FdName%2C%22{}%22%2C%22i%22))%0D%0A%20%0D%0A%20%7D%0D%0A%20ORDER%20BY%20%3Fd%0D%0A%0D%0A%23Find%20a%20set%20of%20descriptors%20and%20terms%20that%20contain%20%22infection.%22&format=JSON&limit=50&offset=0&inference=true'.format(a),
            'NCBIGene': lambda a: "https://pubchem.ncbi.nlm.nih.gov/rest/pug/gene/geneid/{}/summary/JSON".format(a),
            'NCBITaxon': lambda a: "https://pubchem.ncbi.nlm.nih.gov/rest/pug/taxonomy/taxid/{}/summary/JSON".format(a), 
            'OLS': lambda a: 'https://www.ebi.ac.uk/ols/api/ontologies/cl/terms?obo_id={}'.format(a),
            'cellosaurus': lambda a: 'https://api.cellosaurus.org/cell-line/{}?format=json'.format(a),
            'OMIM': lambda a: 'https://omim.org/clinicalSynopsis/{}'.format(a),
        }

        self.mesh_terms = {}
        self.ncbigene_terms = {}
        self.ncbitaxon_terms = {}
        self.omim_terms = {}
        self.cell_types = {}
        self.cell_lines = {}

    def __node__(self):
        return {
            self.pmids if self.pmids else 'id': {
                'title': self.title,
                'label': ':Paper',
                'ref_list': [], #doi
                'NER_terms':  {
                    'gene': self.__get_terms__('gene'), # (id, position in text)
                    'disease': self.__get_terms__('disease'), 
                    'drug': self.__get_terms__('drug'), 
                    'species': self.__get_terms__('species'), 
                    'mutation': self.__get_terms__('mutation'), 
                    'cell_type': self.__get_terms__('cell_type'), 
                    'cell_line': self.__get_terms__('cell_line'), 
                    'DNA': self.__get_terms__('DNA'), 
                    'RNA': self.__get_terms__('RNA'), 
                },
                'authors': [],
                'keywords': [],
                'publication_type': '',
            }
        }
    
    def __term_obj__(self):
        return {
            'mesh': self.mesh_terms,
            'ncbigene': self.ncbigene_terms,
            'ncbitaxon': self.ncbitaxon_terms,
            'omim': self.omim_terms,
            'cell_types': self.cell_types,
            'cell_lines': self.cell_lines,
        }
    
    def __get_terms__(self, term_type):
        term_list = []
        for term in self.bern_res:
            if term['obj'] == term_type:
                term_list.append({'id': term['id'][0], 'position': term['span']})
        
        return term_list

    def __adjacency__(self):
        # loop over mesh terms and add edge to adj_list
        for idx in self.mesh_terms.keys():
            self.adj_list.append({
                'from': self.pmids if self.pmids else 'id',
                'to': idx,
                'type': 'IN_MESH',
            })
        # do this for the rest of the term dictionaries 
        for idx in self.ncbigene_terms.keys():
            self.adj_list.append({
                'from': self.pmids if self.pmids else 'id',
                'to': idx,
                'type': 'IN_GENE',
            })
        for idx in self.ncbitaxon_terms.keys():
            self.adj_list.append({
                'from': self.pmids if self.pmids else 'id',
                'to': idx,
                'type': 'IN_SPECIES',
            })
        for idx in self.omim_terms.keys():
            self.adj_list.append({
                'from': self.pmids if self.pmids else 'id',
                'to': idx,
                'type': 'IN_DISEASE',
            })
        for idx in self.cell_types.keys():
            self.adj_list.append({
                'from': self.pmids if self.pmids else 'id',
                'to': idx,
                'type': 'IN_CELL_TYPE',
            })
        for idx in self.cell_lines.keys():
            self.adj_list.append({
                'from': self.pmids if self.pmids else 'id',
                'to': idx,
                'type': 'IN_CELL_LINE',
            })

        return self.adj_list

    def init_NER(self):
        if self.abstract:
            self.bern_res = [i for i in self.query_plain(self.abstract)['annotations'] if i['prob'] > 0.5 and i['id'][0] != 'CUI-less'] 
        elif self.pmids:
            self.bern_res = [i for i in self.query_pmid(self.pmids)['annotations'] if i['prob'] > 0.5 and i['id'][0] != 'CUI-less'] 

    def query_plain(self, abstract):
        url="http://bern2.korea.ac.kr/plain"
        return requests.post(url, json={'text': abstract}).json()

    def query_pmid(self, pmids):
        url="http://bern2.korea.ac.kr/pubmed"
        return requests.get(url + "/" + ",".join(pmids)).json()

    def populate_terms(self):
        for term in self.bern_res:
            id_raw = term['id'][0]
            parsed_curie = id_raw.split(':')
            db, idx = parsed_curie[0].lower(), parsed_curie[1]
            print(parsed_curie)
            if 'mesh' == db and idx not in self.mesh_terms.keys(): 
                # print(term['id'])
                mesh_res = requests.get(f'https://id.nlm.nih.gov/mesh/{idx}.json').json()
                synonyms = [i['EntryTerm']['value'] for i in requests.get(self.gen_funcs['synonyms'](idx)).json()['results']['bindings']]
                # print(mesh_res)
                # break
                self.mesh_terms[idx] = {
                    'db': db,
                    'heading': mesh_res['label']['@value'],
                    'qualifiers': [i[-6:] for i in mesh_res['allowableQualifier']] if 'allowableQualifier' in mesh_res.keys() else None,
                    'concepts': [i[-6:] for i in mesh_res['concept']] if 'concept' in mesh_res.keys() else None,
                    'pharmacological_actions': mesh_res['pharmacologicalAction'] if 'pharmacologicalAction' in mesh_res.keys() else None,
                    'synonyms': synonyms if synonyms else None,
                    'pred_type': term['obj'],
                }

            elif 'ncbigene' == db and idx not in self.ncbigene_terms.keys():
                ncbi_res = requests.get(self.gen_funcs['NCBIGene'](idx)).json()
                # print(ncbi_res)
                summary = ncbi_res['GeneSummaries']['GeneSummary'][0]
                self.ncbigene_terms[idx] = {
                    **summary,
                    'db': db,
                    'pred_type': term['obj'],
                }

            elif 'ncbitaxon' == db and idx not in self.ncbitaxon_terms.keys():
                ncbi_res = requests.get(self.gen_funcs['NCBITaxon'](idx)).json()
                # print(ncbi_res)
                summary = ncbi_res['TaxonomySummaries']['TaxonomySummary'][0]
                self.ncbitaxon_terms[idx] = {
                    **summary,
                    'db': db,
                    'pred_type': term['obj'],
                }

            # elif 'mim' == parsed_curie[0]: 
            #     omim_id = parsed_curie[1]
            #     page = requests.get(gen_funcs['OMIM'](omim_id))

            #     soup = BeautifulSoup(page.content, "html.parser")

            elif 'cl' == db and idx not in self.cell_types.keys():
                cl_res = requests.get(self.gen_funcs['OLS'](idx)).json()
                # print(cl_res)
                self.cell_types[idx] = {
                    'db': db,
                    'pred_type': term['obj'],
                }
                
            elif 'cellosaurus' == db and idx not in self.cell_lines.keys():
                cellosaurus_res = requests.get(self.gen_funcs['cellosaurus']('CVCL_J'+idx)).json()
                # print(cellosaurus_res)
                self.cell_lines[idx] = {
                    'db': db,
                    'pred_type': term['obj'],
                }
