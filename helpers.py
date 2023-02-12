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
            'title': self.title,
            'ref_list': [], #doi
            'NER_terms':  {
                'Gene': [], # (id, position in text)
                'Disease': [], 
                'Chemical': [], 
                'Species': [], 
                'Mutation': [], 
                'CellType': [], 
                'CellLine': [], 
                'DNA': [], 
                'RNA': [], 
            },
            'authors': [],
            'keywords': [],
            'publication_type': '',
        }

    def init_NER(self):
        if self.abstract:
            self.bern_res = [i for i in self.query_plain(self.abstract)['annotations'] if i['prob'] > 0.5] 
        elif self.pmids:
            self.bern_res = [i for i in self.query_pmid(self.pmids)['annotations'] if i['prob'] > 0.5] 

    def query_plain(abstract):
        url="http://bern2.korea.ac.kr/plain"
        return requests.post(url, json={'text': abstract}).json()

    def query_pmid(pmids):
        url="http://bern2.korea.ac.kr/pubmed"
        return requests.get(url + "/" + ",".join(pmids)).json()

    def populate_terms(self):
        for term in self.bern_res:
            id_raw = term['id'][0]
            parsed_curie = parse_curie(id_raw)
            print(parsed_curie)
            if 'mesh' == parsed_curie[0]: 
                # print(term['id'])
                idx = parsed_curie[1]
                mesh_res = requests.get(f'https://id.nlm.nih.gov/mesh/{idx}.json').json()
                synonyms = [i['EntryTerm']['value'] for i in requests.get(self.gen_funcs['synonyms'](idx)).json()['results']['bindings']]
                # print(mesh_res)
                # break
                self.mesh_terms[idx] = {
                    'heading': mesh_res['label']['@value'],
                    'qualifiers': [i[-6:] for i in mesh_res['allowableQualifier']] if 'allowableQualifier' in mesh_res.keys() else None,
                    'concepts': [i[-6:] for i in mesh_res['concept']] if 'concept' in mesh_res.keys() else None,
                    'pharmacological_actions': mesh_res['pharmacologicalAction'] if 'pharmacologicalAction' in mesh_res.keys() else None,
                    'synonyms': synonyms if synonyms else None,
                    'pred_type': term['obj'],
                }

            elif 'ncbigene' == parsed_curie[0]:
                gene_id = parsed_curie[1]
                ncbi_res = requests.get(self.gen_funcs['NCBIGene'](gene_id)).json()
                # print(ncbi_res)
                self.ncbigene_terms[gene_id] = {
                    'pred_type': term['obj'],
                }

            elif 'ncbitaxon' == parsed_curie[0]:
                tax_id = parsed_curie[1]
                ncbi_res = requests.get(self.gen_funcs['NCBITaxon'](tax_id)).json()
                # print(ncbi_res)
                self.ncbitaxon_terms[tax_id] = {
                    'pred_type': term['obj'],
                }

            # elif 'omim' == parsed_curie[0]: 
            #     omim_id = parsed_curie[1]
            #     page = requests.get(gen_funcs['OMIM'](omim_id))

            #     soup = BeautifulSoup(page.content, "html.parser")

            elif 'cl' == parsed_curie[0]:
                cl_id = parsed_curie[1]
                cl_res = requests.get(self.gen_funcs['OLS'](cl_id)).json()
                # print(cl_res)
                self.cell_types[cl_id] = {
                    'pred_type': term['obj'],
                }
                
            elif 'cellosaurus' == parsed_curie[0]:
                cellosaurus_id = parsed_curie[1]
                cellosaurus_res = requests.get(self.gen_funcs['cellosaurus'](cellosaurus_id)).json()
                # print(cellosaurus_res)
                self.cell_lines[cellosaurus_id] = {
                    'pred_type': term['obj'],
                }
