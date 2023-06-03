import numpy as np

from pymilvus import (
    connections,
    FieldSchema,
    DataType,
    CollectionSchema,
    Collection,
    MilvusException,
)
from sentence_transformers import SentenceTransformer
import pickle
from tqdm import tqdm
import os

model = SentenceTransformer("all-MiniLM-L6-v2")

connect = connections.connect(alias="default", host="localhost", port="19530")

COLLECTION_NAME = "alz_v1"

# create schema if not exists
try:
    Collection(COLLECTION_NAME).drop()
    Collection("als").drop()
except MilvusException as e:
    print(e)

## creat schema

collection = Collection(
    COLLECTION_NAME,
    CollectionSchema(
        [
            FieldSchema(name="pk", dtype=DataType.INT64, is_primary=True, auto_id=True),
            FieldSchema(name="title", dtype=DataType.VARCHAR, max_length=65535),
            FieldSchema(name="source", dtype=DataType.VARCHAR, max_length=65535),
            FieldSchema(name="text", dtype=DataType.VARCHAR, max_length=65535),
            FieldSchema(name="vector", dtype=DataType.FLOAT_VECTOR, dim=384),
        ]
    ),
)

data_file = os.path.join(os.path.dirname(__file__), "AlzConnected.pkl")
with open(data_file, "rb") as input_file:
    data = np.array([(k, v) for k, v in pickle.load(input_file).items()], dtype=object)


for batch in tqdm(np.array_split(data, 128)):
    titles = [x[1][0]["title"][:65534] for x in batch]
    sources = [x[0] for x in batch]
    texts = [
        "\n\n=======\n\n".join([item["body"] for item in x[1]])[:65534] for x in batch
    ]
    embedding_texts = [
        (x[1][0]["title"] + "\n\n" + x[1][0]["body"])[:65534] for x in batch
    ]
    collection.insert(
        [
            titles,
            sources,
            texts,
            model.encode(embedding_texts, normalize_embeddings=True, batch_size=32),
        ]
    )
    collection.flush()


index = {
    "index_type": "IVF_FLAT",
    "metric_type": "L2",
    "params": {"nlist": 100},
}
collection.create_index("vector", index)
