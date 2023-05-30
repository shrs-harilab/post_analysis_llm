from langchain.embeddings import SentenceTransformerEmbeddings
from langchain.vectorstores import Milvus
import streamlit as st
@st.cache_resource
def load_collections() -> dict[str, Milvus]:
    embeddings = SentenceTransformerEmbeddings(model_name='all-MiniLM-L6-v2')
    embeddings.encode_kwargs = dict(normalize_embeddings = True)
    connection_args = {"host": "localhost", "port": "19530"}
    search_params = {"metric_type": "L2", "params": {"nprobe": 10}}
    vector_stores = dict.fromkeys(["alz", "pdf"])
    for collection in vector_stores.keys():
        vector_store = Milvus(embedding_function=embeddings, collection_name=collection, connection_args = connection_args, search_params=search_params)
        vector_stores[collection] = vector_store
    return vector_stores