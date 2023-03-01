import json
import numpy as np

from pymilvus import connections, FieldSchema, DataType, CollectionSchema, Collection
from sentence_transformers import SentenceTransformer
import pickle
from tqdm import tqdm

model = SentenceTransformer('all-MiniLM-L6-v2')

connect = connections.connect(
  alias="default",
  host='localhost',
  port='19530'
)

ALS_COLLECTION_NAME = "als"

# create schema if not exists
try:
    collection = Collection(ALS_COLLECTION_NAME)
    collection.drop()
except:

    pass
## creat schema
fields = [
        FieldSchema(name="id", dtype=DataType.VARCHAR,max_length=128, is_primary=True, auto_id=False),
        FieldSchema(name="title", dtype=DataType.VARCHAR,max_length=65535),
        FieldSchema(name="url", dtype=DataType.VARCHAR,max_length=65535),
        FieldSchema(name="text", dtype=DataType.VARCHAR,max_length=65535),
        FieldSchema(name="embeddings", dtype=DataType.FLOAT_VECTOR, dim=384)
    ]
schema = CollectionSchema(fields, "als")
collection = Collection("als", schema)

data_file = "./AlzConnected.pkl"
with open(data_file,"rb") as input_file:
    data = [(k,v) for k,v in pickle.load(input_file).items()]


for batch in tqdm(np.array_split(data,128)):
    ids = [x[0] for x in batch]
    titles = [x[1][0]["title"][:65534] for x in batch]
    embedding_texts = [(x[1][0]["title"] + "\n\n" + x[1][0]["body"])[:65534] for x in batch]
    texts = ["\n\n=======\n\n".join([item["body"] for item in x[1]])[:65534] for x in batch]

    collection.insert([
        ids,
        titles,
        ids,
        texts,
        model.encode(embedding_texts,normalize_embeddings=True,batch_size=32)
    ])
    collection.flush()







index = {
    "index_type": "IVF_FLAT",
    "metric_type": "L2",
    "params": {"nlist": 100},
}
collection.create_index("embeddings", index)



# query
# collection = Collection("als")
# collection.load()
#
#
# def search(text:str,k:int = 5):
#     vectors_to_search = model.encode([text])
#     search_params = {
#         "metric_type": "L2",
#         "params": {"nprobe": 10},
#     }
#     result = collection.search(vectors_to_search, "embeddings", search_params, limit=k, output_fields=["id","title", "text"])
#     hits = result[0]
#     print(f"hits: {len(hits)}")
#     for hit in hits:
#         print(f"{hit.id} ({hit.score}): {hit.entity.title}")
#
#
# search("Smart device")




