import os

import numpy as np
from tqdm import tqdm

from knowledge_base.parse_paragraph_from_pdf import parse_paragraphs_from_pdf
from pymilvus import connections, FieldSchema, DataType, CollectionSchema, Collection
from sentence_transformers import SentenceTransformer




def load_sources(pdf_rootpath):
    sources = []
    for (_, _, filenames) in os.walk(pdf_rootpath):
        for filename in filenames:
            if not filename.endswith(".pdf"):
                continue
            pdf_filepath = os.path.join(PDF_ROOT_PATH, filename)
            paragraphs = parse_paragraphs_from_pdf(pdf_filepath)
            sources.append(
                (filename, paragraphs)
            )
    return sources


if __name__ == '__main__':
    PDF_ROOT_PATH = "./pdfs"
    BATCH_SIZE = 64


    connect = connections.connect(
        alias="default",
        host='localhost',
        port='19530'
    )

    COLLECTION_NAME = "pdf"

    # create schema if not exists
    try:
        collection = Collection(COLLECTION_NAME)
        collection.drop()
    except:
        pass
    ## creat schema
    fields = [
        FieldSchema(name="id", dtype=DataType.VARCHAR, max_length=1024, is_primary=True, auto_id=False),
        FieldSchema(name="text", dtype=DataType.VARCHAR, max_length=65535),
        # FieldSchema(name="position", dtype=DataType.INT16),
        FieldSchema(name="embeddings", dtype=DataType.FLOAT_VECTOR, dim=384)
    ]
    schema = CollectionSchema(fields, COLLECTION_NAME)
    collection = Collection(COLLECTION_NAME, schema)

    sources = load_sources(PDF_ROOT_PATH)
    source_paragraph_pair_list = [ ("{}-{}".format(filename,i),paragraph) for (filename,paragraphs) in sources for i,paragraph in enumerate(paragraphs)]

    model = SentenceTransformer('all-MiniLM-L6-v2')

    pool = model.start_multi_process_pool()

    for batch in tqdm(np.array_split(source_paragraph_pair_list, len(source_paragraph_pair_list) / BATCH_SIZE)):
        ids = [x[0] for x in batch]
        texts = [x[1][:65534] for x in batch]
        model.mul
        collection.insert([
            ids,
            texts,
            model.encode_multi_process(texts,normalize_embeddings=True, batch_size=BATCH_SIZE)
        ])
        collection.flush()

index = {
    "index_type": "IVF_FLAT",
    "metric_type": "L2",
    "params": {"nlist": 100},
}
collection.create_index("embeddings", index)
