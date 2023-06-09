import os
import numpy as np
from tqdm import tqdm
from parse_paragraph_from_pdf import parse_paragraphs_from_pdf
from pymilvus import (
    connections,
    FieldSchema,
    DataType,
    CollectionSchema,
    Collection,
    MilvusException,
    utility,
)
from sentence_transformers import SentenceTransformer
from mpire import WorkerPool
from multiprocessing import cpu_count
from pympler.asizeof import asizeof


def parse_one_pdf_file(filename):
    pdf_filepath = os.path.join(PDF_ROOT_PATH, filename)
    paragraphs = parse_paragraphs_from_pdf(pdf_filepath)
    sources = []
    for i in range(len(paragraphs)):
        paragraphs[i] = paragraphs[i][:65534]
        sources.append(f"{filename}-{i}")
    return sources, paragraphs


def load_sources(pdf_rootpath):
    for _, _, filenames in os.walk(pdf_rootpath):
        n_jobs = round(cpu_count() * 0.8) + 1
        with WorkerPool(n_jobs) as pool:
            data = pool.map(
                parse_one_pdf_file,
                [f for f in filenames if f.endswith(".pdf")],
                progress_bar=True,
            )

    all_sources, all_paragraphs = [], []
    for sources, paragraphs in data:
        all_sources.extend(sources)
        all_paragraphs.extend(paragraphs)

    return all_sources, all_paragraphs


if __name__ == "__main__":
    dir_path = os.path.dirname(os.path.realpath(__file__))
    PDF_ROOT_PATH = os.path.join(dir_path, "pdfs")
    BATCH_SIZE = 64
    COLLECTION_NAME = "pdf"

    connections.connect(alias="default", host="localhost", port="19530")

    # create schema if not exists
    try:
        Collection(COLLECTION_NAME).drop()
    except:
        pass
    # creat schema
    collection = Collection(
        COLLECTION_NAME,
        CollectionSchema(
            [
                FieldSchema(
                    name="pk", dtype=DataType.INT64, is_primary=True, auto_id=True
                ),
                FieldSchema(name="source", dtype=DataType.VARCHAR, max_length=1024),
                FieldSchema(name="text", dtype=DataType.VARCHAR, max_length=65535),
                FieldSchema(name="vector", dtype=DataType.FLOAT_VECTOR, dim=384),
            ]
        ),
    )
    index = {
        "index_type": "IVF_FLAT",
        "metric_type": "L2",
        "params": {"nlist": 100},
    }

    sources, paragraphs = load_sources(PDF_ROOT_PATH)
    model = SentenceTransformer("all-MiniLM-L6-v2")
    embeddings = model.encode(
        paragraphs,
        normalize_embeddings=True,
        batch_size=BATCH_SIZE,
        show_progress_bar=True,
    )
    data = [sources, paragraphs, embeddings]
    # split data to insert into milvus, can't exceed 64MB per insert
    partition_num = round(asizeof(data) / (pow(1024, 2) * 64)) + 1
    print(f"Total data to insert: {asizeof(data)/pow(1024, 3)} GB")
    for i in range(len(data)):
        data[i] = np.array_split(data[i], partition_num)
    try:
        for d in tqdm(zip(*data)):
            collection.insert(list(d))
            collection.flush()
        collection.create_index("vector", index)
    except MilvusException as e:
        print(e)
