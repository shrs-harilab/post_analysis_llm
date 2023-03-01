import typing

import streamlit as st
from pymilvus import connections, Collection
from sentence_transformers import SentenceTransformer



def use_state(key: str, default=None):
    if key not in st.session_state:
        st.session_state[key] = default

    def set_value(value):
        st.session_state[key] = value
    return st.session_state[key], set_value




@st.cache_resource
def load_collection():
    connections.connect(
        alias="default",
        host='localhost',
        port='19530'
    )

    collection = Collection("als")
    collection.load()
    return collection

@st.cache_resource
def load_model():
    model = SentenceTransformer('all-MiniLM-L6-v2')
    return model

collection = load_collection()
model = load_model()


def search(text:str,k:int = 5):
    if not text:
        return []
    vectors_to_search = model.encode([text],normalize_embeddings=True)
    search_params = {
        "metric_type": "L2",
        "params": {"nprobe": 10},
    }
    result = collection.search(vectors_to_search, "embeddings", search_params, limit=k, output_fields=["id","title", "text","url"])
    hits = result[0]
    return hits

page_size, set_page_size = use_state("page_size",10)
is_expand_set, set_is_expand_set = use_state("is_expand",set())

# UI


search_query = st.text_input("search from AlzConnected forum data")
hits = search(search_query,page_size)


def on_load_more():
    set_page_size(page_size + 5)

def on_expand_text(id:str):
    print(is_expand_set)
    is_expand_set.add(id)
    set_is_expand_set(is_expand_set)


TEXT_LIMIT = 500
for i,hit in enumerate(hits):
    st.markdown(f"{i+1}. **{hit.entity.title}**")
    fulltext = hit.entity.text
    if(len(fulltext) > TEXT_LIMIT and not(hit.entity.url in is_expand_set)):
        st.markdown(f"{fulltext[:TEXT_LIMIT]}...")
        st.button("READ MORE", key=hit.entity.url, on_click=lambda id=hit.entity.url:on_expand_text(id))
    else:
        st.markdown(f"{fulltext}")
    st.markdown("-------")

if hits:
    st.button("load more",on_click=on_load_more)

