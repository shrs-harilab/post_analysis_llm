from Bio import Entrez
import pandas as pd
import numpy as np
import requests


# Get searched results
def search(query):
    Entrez.email = 'hah90@pitt.edu'
    handle = Entrez.esearch(db='pubmed',
                            sort='relevance',
                            retmax='5',
                            retmode='xml',
                            term=query)
    results = Entrez.read(handle)
    return results


# Get papers based on the results ids
def fetch_details(id_list):
    ids = ','.join(id_list)
    Entrez.email = 'hah90@pitt.edu'
    handle = Entrez.efetch(db='pubmed',
                           retmode='xml',
                           id=ids)
    results = Entrez.read(handle)
    return results


def get_pmc_id(pm_id):
    try:
        Entrez.email = 'hah90@pitt.edu'
        handle = Entrez.elink(dbfrom="pubmed",
                              db="pmc",
                              linkname="pubmed_pmc",
                              id=pm_id,
                              retmode="text")
        result = Entrez.read(handle)
        return result[0]['LinkSetDb'][0]['Link'][0]['Id']
    except:
        return '0'


def get_papers_summary(id_list):
    papers = fetch_details(id_list)
    paper_arr = []
    for paper in papers['PubmedArticle']:
        pm_id = ''.join(paper['MedlineCitation']['PMID'])
        pmc_id = get_pmc_id(pm_id)
        article = {
            'title': paper['MedlineCitation']['Article']['ArticleTitle'],
            'abstract': ''.join(paper['MedlineCitation']['Article']['Abstract']['AbstractText']),
            'doi': ''.join(paper['MedlineCitation']['Article']['ELocationID']),
            'full_text': f'http://www.ncbi.nlm.nih.gov/pmc/articles/{pmc_id}/pdf/' if pmc_id != '0' else ''
        }
        paper_arr.append(article)
    paper_df = pd.DataFrame(paper_arr)
    paper_df.index += 1
    paper_df.to_excel('./pubmed_search.xlsx')
    return paper_arr


def load_pdfs_to_knowledge_base(df_summary):
    headers = {
        "User-Agent": "Chrome/111.0.0.0"
    }

    for index, df_item in enumerate(df_summary):
        if df_item['full_text'] != '':
            response = requests.get(df_item['full_text'], headers=headers)
            # Save the PDF
            if response.status_code == 200:
                with open(f"knowledge_base/pdfs/{index + 1}.pdf", "wb") as f:
                    f.write(response.content)
            else:
                print(response.status_code)


results = search('delirium')
id_list = results['IdList']
summary = get_papers_summary(id_list)
load_pdfs_to_knowledge_base(summary)
