from helpers import Paper

from flask import Flask, render_template, request, redirect, url_for

app = Flask(__name__)

@app.route('/')
def main():
    return '<h1>connnexsci NER api</h1>'

@app.route('/compile_node', methods=['POST', 'GET'])
def compile_node():
    if request.method == 'POST':
        paper = Paper(title=request.form['title'], abstract=request.form['abstract'])
        paper.init_NER()
        paper.populate_terms()
        
        paper_node = paper.__node__()
        term_nodes = paper.__term_obj__()
        adj_list = paper.__adjacency__()
        return {'paper_node': paper_node, 'term_nodes': term_nodes, 'adj_list': adj_list}
    return render_template('form.html')