from Bio import Entrez
from dotenv import load_dotenv
import os
from mpire import WorkerPool, cpu_count
from mpire.utils import make_single_arguments
import requests
import shutil
from tqdm import tqdm
import logging
import pandas as pd
import numpy as np
from math import ceil
from argparse import ArgumentParser

class PubMedCrawler:
    def __init__(self, term: str):
        Entrez.email = os.environ.get("ENTREZ_EMAIL")
        Entrez.api_key = os.environ.get("ENTREZ_API_KEY")
        self._term = term
        self._tasks = list()
        self._details = list()
        self._DIR = os.path.dirname(__file__)
        self._PDFDIR = os.path.join(self._DIR, "pdfs")
        self._create_dir()
        self.logger = logging.Logger(__name__)
        self.logger.addHandler(logging.FileHandler(os.path.join(self._DIR, f"{self.__class__.__name__}.log")))

    def crawl(self, crawl_num: int):
        lack = crawl_num
        step = 200
        with tqdm(total=lack, desc="Fetching Records") as pbar:
            i = 0
            while lack > 0:
                result = self._search(i, step)
                lack = crawl_num - len(self._tasks)
                i += step
                pbar.update(n=result if lack >= 0 else crawl_num - pbar.last_print_n)
        with WorkerPool(cpu_count()+1) as pool:
            pool.map(
                self._consume,
                make_single_arguments(self._tasks[:crawl_num], generator=False),
                progress_bar=True,
                progress_bar_options={"desc": "Downloading Files"},
            )
        for batch in tqdm(
            np.array_split(self._tasks, ceil(len(self._tasks) / step)),
            desc="Recording Details",
        ):
            self._record_details(batch)
        
        self._to_file()

    def _consume(self, task: dict):
        try:
            response = requests.get(
                task["url"], headers={"User-Agent": "Chrome/111.0.0.0"}
            )
            if response.status_code == requests.codes.ok:
                with open(os.path.join(self._PDFDIR, f"{task['pmid']}.pdf"), "wb") as file:
                    file.write(response.content)
        except Exception as e:
            self.logger.error(f"URL: {task['url']}, EXCEPTION: {e}")

    def _search(self, retstart: int, retmax: int) -> int:
        handle = Entrez.esearch(
            db="pubmed",
            term=self._term,
            sort="relevance",
            retstart=retstart,
            retmax=retmax,
            retmode="xml",
        )
        records = Entrez.read(handle)
        handle.close()
        return self._link(records["IdList"])

    def _record_details(self, tasks: list[dict]):
        handle = Entrez.efetch(
            db="pubmed",
            retmode="xml",
            id=",".join([task["pmid"] for task in tasks]),
        )
        papers = Entrez.read(handle)
        handle.close()
        for task, paper in zip(tasks, papers["PubmedArticle"]):
            doi_list = list(
                filter(
                    lambda item: item.attributes["EIdType"] == "doi",
                    paper["MedlineCitation"]["Article"]["ELocationID"],
                )
            )
            task.update(
                {
                    "title": paper["MedlineCitation"]["Article"]["ArticleTitle"],
                    "doi": str(doi_list[0]) if len(doi_list) > 0 else None,
                }
            )

    def _link(self, pmid: list) -> int:
        ids = ",".join(pmid)
        handle = Entrez.elink(
            dbfrom="pubmed", db="pmc", linkname="pubmed_pmc", id=ids, cmd="prlinks"
        )
        result = Entrez.read(handle)
        handle.close()
        get = 0
        for idurl in result[0]["IdUrlList"]["IdUrlSet"]:
            for obj in idurl["ObjUrl"]:
                if obj["Provider"]["NameAbbr"] == "PMC":
                    self._tasks.append(
                        {
                            "pmid": idurl["Id"],
                            "url": obj["Url"] + "pdf",
                        }
                    )
                    get += 1
        return get

    def _to_file(self):
        df = pd.DataFrame(self._tasks)
        df.to_excel(os.path.join(self._DIR, "pubmed_search.xlsx"))

    def _create_dir(self):
        if os.path.exists(self._PDFDIR):
            shutil.rmtree(self._PDFDIR)
        os.mkdir(self._PDFDIR)


if __name__ == "__main__":
    load_dotenv()
    parser = ArgumentParser()
    parser.add_argument("-t", "--term", default="delirium")
    parser.add_argument("-n", "--num", type=int, default=30)
    args = parser.parse_args()
    crawler = PubMedCrawler(args.term)
    crawler.crawl(args.num)
