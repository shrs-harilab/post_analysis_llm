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
import math


class PubMedCrawler:
    def __init__(self, dir: str, max_crawl_num: int = 10000):
        Entrez.email = os.environ.get("ENTREZ_EMAIL")
        Entrez.api_key = os.environ.get("ENTREZ_API_KEY")
        self.max_crawl_num = max_crawl_num
        self.tasks = list()
        self.details = list()
        self.DIR = os.path.dirname(__file__)
        self.PDFDIR = os.path.join(self.DIR, "pdfs", dir)
        self._create_dir()
        self.logger = logging.Logger(__name__)
        self.logger.addHandler(
            logging.FileHandler(
                os.path.join(self.DIR, f"{self.__class__.__name__}.log"), mode="w"
            )
        )

    def crawl(self, term: str):
        self.term = term
        step = 200
        with tqdm(total=self.max_crawl_num, desc="Searching Records") as pbar:
            i = 0
            for _ in range(math.floor(self.max_crawl_num / step)):
                self._search(i, step)
                i += step
                pbar.update(step)
            mod = self.max_crawl_num % step
            if mod > 0:
                self._search(i, mod)
                pbar.update(mod)
        with WorkerPool(math.floor(cpu_count() / 2) + 1) as pool:
            pool.map(
                self._consume,
                make_single_arguments(self.tasks, generator=False),
                progress_bar=True,
                progress_bar_options={"desc": "Downloading Files"},
            )
        for batch in tqdm(
            np.array_split(self.tasks, ceil(len(self.tasks) / step)),
            desc="Recording Details",
        ):
            self._record_details(batch)

        self.save_records_to_file()

    def _consume(self, task: dict):
        try:
            response = requests.get(
                task["url"], headers={"User-Agent": "Chrome/111.0.0.0"}
            )
            if response.status_code == requests.codes.ok:
                with open(
                    os.path.join(self.PDFDIR, f"{task['pmid']}.pdf"), "wb"
                ) as file:
                    file.write(response.content)
        except Exception as e:
            self.logger.error(f"URL: {task['url']}, EXCEPTION: {e}")

    def _search(self, retstart: int, retmax: int) -> int:
        handle = Entrez.esearch(
            db="pubmed",
            term=self.term,
            sort="relevance",
            retstart=retstart,
            retmax=retmax,
            retmode="xml",
        )
        records = Entrez.read(handle)
        handle.close()
        self._link(records["IdList"])

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
        handle = Entrez.elink(
            dbfrom="pubmed",
            db="pmc",
            linkname="pubmed_pmc",
            id=",".join(pmid),
            cmd="prlinks",
        )
        result = Entrez.read(handle)
        handle.close()
        for idurl in result[0]["IdUrlList"]["IdUrlSet"]:
            for obj in idurl["ObjUrl"]:
                if obj["Provider"]["NameAbbr"] == "PMC":
                    self.tasks.append(
                        {
                            "pmid": idurl["Id"],
                            "url": obj["Url"] + "pdf",
                        }
                    )

    def save_records_to_file(self):
        df = pd.DataFrame(self.tasks)
        df.to_excel(os.path.join(self.DIR, "pubmed_search.xlsx"))

    def _create_dir(self):
        if os.path.exists(self.PDFDIR):
            shutil.rmtree(self.PDFDIR)
        os.makedirs(self.PDFDIR)


if __name__ == "__main__":
    load_dotenv()
    parser = ArgumentParser()
    parser.add_argument("-t", "--term", default="delirium")
    parser.add_argument("-d", "--dir")
    parser.add_argument("-n", "--num", type=int, default=30)
    args = parser.parse_args()
    crawler = PubMedCrawler(
        dir=args.dir if args.dir else "_".join(args.term.split()),
        max_crawl_num=args.num,
    )
    crawler.crawl(args.term)
